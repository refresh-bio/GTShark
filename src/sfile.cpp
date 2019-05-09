// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

#include "defs.h"
#include "sfile.h"
#include "utils.h"
#include "lzma_wrapper.h"

// ************************************************************************************
CSampleFile::CSampleFile()
{
	rc = nullptr;
	rce = nullptr;
	rcd = nullptr;

	vios = new CVectorIOStream(v_uint8);

	mode = mode_t::none;
}

// ************************************************************************************
CSampleFile::~CSampleFile()
{
	Close();

	if (rc)
		delete rc;

	delete vios;
}

// ************************************************************************************
bool CSampleFile::OpenForReading(string file_name, bool& _extra_variants)
{
	if (!fi_sample.Open(file_name))
	{
		cerr << "Cannot open " << file_name << " file\n";
		exit(1);
	}

	rc = new CRangeDecoder<CVectorIOStream>(*vios);
	rcd = (CRangeDecoder<CVectorIOStream>*) rc;

	size_t file_size = fi_sample.FileSize();
	_extra_variants = (bool) fi_sample.GetByte();
	extra_variants = _extra_variants;

	file_size -= read_header_data();
	file_size--;						// extra variants info

	if (extra_variants)
		file_size -= read_extra_variants();

	v_uint8.clear();
	v_uint8.reserve(file_size);

	while (!fi_sample.Eof())
		v_uint8.push_back(fi_sample.GetByte());

	rcd->Start();

	mode = mode_t::decompress;
	ctx_flag = 0;

	return true;
}

// ************************************************************************************
bool CSampleFile::OpenForWriting(string file_name, bool _extra_variants)
{
	if (!fo_sample.Open(file_name))
	{
		cerr << "Cannot open " << file_name << " file\n";
		exit(1);
	}

	mode = mode_t::compress;
	extra_variants = _extra_variants;
	ctx_flag = 0;

	fo_sample.PutByte((uint8_t)extra_variants);

	rc = new CRangeEncoder<CVectorIOStream>(*vios);
	rce = (CRangeEncoder<CVectorIOStream>*) rc;

	rce->Start();
	
	return true;
}

// ************************************************************************************
bool CSampleFile::Close()
{
	if (mode == mode_t::compress)
	{
		rce->End();
		fo_sample.Write(vios->Data(), vios->Size());
		fo_sample.Close();
	}
	else if (mode == mode_t::decompress)
	{
		rcd->End();
		fi_sample.Close();
	}

	return true;
}

// ************************************************************************************
uint32_t CSampleFile::read_header_data()
{
	uint32_t readed_bytes = 0;

	uint8_t header_present = fi_sample.GetByte();
	++readed_bytes;

	if (header_present)
	{
		uint32_t header_compressed_size = fi_sample.ReadUInt(4);
		readed_bytes += 4;

		input_vec_rel_header.resize(header_compressed_size);
		
		fi_sample.Read(input_vec_rel_header.data(), header_compressed_size);
		readed_bytes += header_compressed_size;
	}

	uint32_t sample_name_size = fi_sample.ReadUInt(2);
	readed_bytes += 2;
	
	input_sample_name.reserve(sample_name_size + 1);
	
	for(uint32_t i = 0; i < sample_name_size; ++i)
		input_sample_name.push_back(fi_sample.Get());

	readed_bytes += sample_name_size;

	return readed_bytes;
}

// ************************************************************************************
uint32_t CSampleFile::read_extra_variants()
{
	uint32_t ev_present = fi_sample.GetByte();
	uint32_t no_bytes = 0;

	if (!ev_present)
		return 1;

	vector<uint8_t> vr_chrom, vc_chrom;
	vector<uint8_t> vr_pos, vc_pos;
	vector<uint8_t> vr_id, vc_id;
	vector<uint8_t> vr_ref, vc_ref;
	vector<uint8_t> vr_alt, vc_alt;
	vector<uint8_t> vr_qual, vc_qual;
	vector<uint8_t> vr_filter, vc_filter;
	vector<uint8_t> vr_info, vc_info;
	vector<uint8_t> vr_gt, vc_gt;

	for (auto d : {
		make_tuple(ref(vr_chrom), ref(vc_chrom), 9, "chrom"),
		make_tuple(ref(vr_pos), ref(vc_pos), 9, "pos"),
		make_tuple(ref(vr_id), ref(vc_id), 9, "id"),
		make_tuple(ref(vr_ref), ref(vc_ref), 9, "ref"),
		make_tuple(ref(vr_alt), ref(vc_alt), 9, "alt"),
		make_tuple(ref(vr_qual), ref(vc_qual), 9, "qual"),
		make_tuple(ref(vr_filter), ref(vc_filter), 9, "filter"),
		make_tuple(ref(vr_info), ref(vc_info), 9, "info"),
		make_tuple(ref(vr_gt), ref(vc_gt), 9, "gt")
		})
	{
		uint32_t comp_size = fi_sample.ReadUInt(4);
		get<1>(d).resize(comp_size);
		fi_sample.Read(get<1>(d).data(), get<1>(d).size());

		no_bytes += 4 + comp_size;

		CLZMAWrapper::Decompress(get<1>(d), get<0>(d));
//		cout << get<3>(d) << " size: " << get<1>(d).size() << endl;
	}

	auto p_chrom = vr_chrom.begin();
	auto p_id = vr_id.begin();
	auto p_ref = vr_ref.begin();
	auto p_alt = vr_alt.begin();
	auto p_qual = vr_qual.begin();
	auto p_filter = vr_filter.begin();
	auto p_info = vr_info.begin();
	auto p_pos = vr_pos.begin();
	auto p_gt = vr_gt.begin();

	int prev_pos = 0;
	while(true)
	{
		loc_v_desc.push_back(make_pair(variant_desc_t(), vector<uint8_t>()));

		loc_v_desc.back().first.chrom = get_string_from_vector(p_chrom);
		loc_v_desc.back().first.id = get_string_from_vector(p_id);
		loc_v_desc.back().first.ref = get_string_from_vector(p_ref);
		loc_v_desc.back().first.alt = get_string_from_vector(p_alt);
		loc_v_desc.back().first.qual = get_string_from_vector(p_qual);
		loc_v_desc.back().first.filter = get_string_from_vector(p_filter);
		loc_v_desc.back().first.info = get_string_from_vector(p_info);
		
		uint32_t dif_pos = 0;
		for (int i = 0; i < 4; ++i)
			dif_pos += ((uint32_t) * (p_pos + i)) << (8 * i);
		p_pos += 4;

		loc_v_desc.back().first.pos = (int32_t)(dif_pos + (uint32_t) prev_pos);
		prev_pos = (int32_t)(dif_pos + (uint32_t)prev_pos);
		
		loc_v_desc.back().second.push_back(*p_gt++);

		if (p_chrom == vr_chrom.end())
			break;
	}

	return no_bytes;
}

// ************************************************************************************
string CSampleFile::get_string_from_vector(vector<uint8_t>::iterator& p)
{
	string r;

	for (; *p; ++p)
		r.push_back((char)* p);

	++p;			// skip terminator

	return r;
}

// ************************************************************************************
void CSampleFile::ReadExtraVariants(vector<pair<variant_desc_t, vector<uint8_t>>>& v_desc)
{
	v_desc = move(loc_v_desc);
}

// ************************************************************************************
uint32_t CSampleFile::WriteExtraVariants(const vector<pair<variant_desc_t, vector<uint8_t>>>& v_desc)
{
	if (v_desc.empty())
	{
		fo_sample.PutByte(0);			// no extra data
		return 1;
	}

	uint32_t no_ev = (uint32_t)v_desc.size();

	vector<uint8_t> vr_chrom, vc_chrom;
	vector<uint8_t> vr_pos, vc_pos;
	vector<uint8_t> vr_id, vc_id;
	vector<uint8_t> vr_ref, vc_ref;
	vector<uint8_t> vr_alt, vc_alt;
	vector<uint8_t> vr_qual, vc_qual;
	vector<uint8_t> vr_filter, vc_filter;
	vector<uint8_t> vr_info, vc_info;
	vector<uint8_t> vr_gt, vc_gt;
	string empty_info = ".";

	int prev_pos = 0;
	for (auto& x : v_desc)
	{
		vr_chrom.insert(vr_chrom.end(), x.first.chrom.begin(), x.first.chrom.end());		vr_chrom.push_back(0);
		vr_id.insert(vr_id.end(), x.first.id.begin(), x.first.id.end());					vr_id.push_back(0);
		vr_ref.insert(vr_ref.end(), x.first.ref.begin(), x.first.ref.end());				vr_ref.push_back(0);
		vr_alt.insert(vr_alt.end(), x.first.alt.begin(), x.first.alt.end());				vr_alt.push_back(0);
		vr_qual.insert(vr_qual.end(), x.first.qual.begin(), x.first.qual.end());			vr_qual.push_back(0);
		vr_filter.insert(vr_filter.end(), x.first.filter.begin(), x.first.filter.end());	vr_filter.push_back(0);
//		vr_info.insert(vr_info.end(), x.first.info.begin(), x.first.info.end());			vr_info.push_back(0);
		vr_info.insert(vr_info.end(), empty_info.begin(), empty_info.end());				vr_info.push_back(0);

		uint32_t dif_pos = (x.first.pos - prev_pos);
		prev_pos = x.first.pos;

		for (int i = 0; i < 4; ++i)
			vr_pos.push_back((uint8_t)((dif_pos >> (8 * i)) & 0xff));

		for (auto y : x.second)
			vr_gt.push_back(y);
	}

	// Save variant descriptions
	uint32_t no_bytes = 0;
	
	fo_sample.PutByte(1);
	++no_bytes;

	for (auto d : {
		make_tuple(ref(vr_chrom), ref(vc_chrom), 9, "chrom"),
		make_tuple(ref(vr_pos), ref(vc_pos), 9, "pos"),
		make_tuple(ref(vr_id), ref(vc_id), 9, "id"),
		make_tuple(ref(vr_ref), ref(vc_ref), 9, "ref"),
		make_tuple(ref(vr_alt), ref(vc_alt), 9, "alt"),
		make_tuple(ref(vr_qual), ref(vc_qual), 9, "qual"),
		make_tuple(ref(vr_filter), ref(vc_filter), 9, "filter"),
		make_tuple(ref(vr_info), ref(vc_info), 9, "info"),
		make_tuple(ref(vr_gt), ref(vc_gt), 9, "gt")
		})
	{
		CLZMAWrapper::Compress(get<0>(d), get<1>(d), get<2>(d));
//		cout << get<3>(d) << " size: " << get<1>(d).size() << endl;
		fo_sample.WriteUInt(get<1>(d).size(), 4);
		fo_sample.Write(get<1>(d).data(), get<1>(d).size());

		no_bytes += 4 + get<1>(d).size();
	}

	return no_bytes;
}

// ************************************************************************************
bool CSampleFile::WriteHeaderAndSample(const string &db_header, const string &v_header, const string &sample_name)
{
	if (input_vec_rel_header.empty())
		fo_sample.PutByte(0);		// header absent
	else
	{
		fo_sample.PutByte(1);		// header present

		vector<uint8_t> vec_header(db_header.begin(), db_header.end());
		vector<uint8_t> vec_v_header(v_header.begin(), v_header.end());
		vector<uint8_t> vec_rel_header;

		CLZMAWrapper::CompressWithHistory(vec_header, vec_v_header, vec_rel_header, 9);
		fo_sample.WriteUInt(vec_rel_header.size(), 4);
		fo_sample.Write(vec_rel_header.data(), vec_rel_header.size());
	}

	fo_sample.WriteUInt(sample_name.size(), 2);
	fo_sample.Write((uint8_t *) sample_name.data(), sample_name.size());

	return true;
}

// ************************************************************************************
bool CSampleFile::ReadHeaderAndSample(const string &db_header, string &v_header, string &sample_name)
{
	if (input_vec_rel_header.empty())
		v_header = db_header;
	else
	{
		vector<uint8_t> vec_header(db_header.begin(), db_header.end());
		vector<uint8_t> vec_v_header;
		vector<uint8_t> vec_rel_header;

		CLZMAWrapper::DecompressWithHistory(vec_header, vec_rel_header, vec_v_header);

		v_header.assign(vec_v_header.begin(), vec_v_header.end());
	}

	sample_name = input_sample_name;

	return true;
}

// ************************************************************************************
uint64_t CSampleFile::determine_context(const run_t &run, uint32_t no_pred_same, uint32_t no_succ_same)
{
	context_t ctx = 0;

	for (int i = 0; i < 2; ++i)
	{
		ctx += ((context_t)run[i].first) << (16 * i + 8);
		ctx += ((ilog2(run[i].second) + 1) / 4) << (16 * i);
	}

	if (no_pred_same > no_succ_same)
		ctx += 1ull << 62;
	else if (no_pred_same < no_succ_same)
		ctx += 2ull << 62;

	if (no_pred_same > no_succ_same)
		ctx += (((context_t)ilog2(no_pred_same) + 3) / 4) << 32;
	else
		ctx += (((context_t)ilog2(no_succ_same) + 3) / 4) << 40;

	return ctx;
}

// ************************************************************************************
inline CSampleFile::ctx_map_e_t::value_type CSampleFile::find_rc_coder(const run_t &run, uint32_t no_pred_same, uint32_t no_succ_same)
{
	context_t ctx = determine_context(run, no_pred_same, no_succ_same);
	auto p = rc_coders.find(ctx);

	if (p == nullptr)
	{
		int init_stat[SIGMA] = { 1, 1, 1, 1 };
		rc_coders.insert(ctx, p = new CRangeCoderModel<CVectorIOStream>(rc, 4, 13, 1 << 13, init_stat, 4, mode == mode_t::compress));
	}

	return p;
}

// ************************************************************************************
inline CSampleFile::ctx_map_e_t::value_type CSampleFile::find_rc_coder(context_t ctx)
{
	ctx += 1ull << 62;

	auto p = rc_coders.find(ctx);

	if (p == nullptr)
	{
		int init_stat[] = { 1, 1, 1, 1, 1 };
		rc_coders.insert(ctx, p = new CRangeCoderModel<CVectorIOStream>(rc, 5, 15, 1 << 15, init_stat, 4, mode == mode_t::compress));
	}

	return p;
}

// ************************************************************************************
bool CSampleFile::Put(uint8_t value, run_t &run, uint32_t no_pred_same, uint32_t no_succ_same)
{
	auto p = find_rc_coder(run, no_pred_same, no_succ_same);

	p->Encode(value);

	return true;
}

// ************************************************************************************
bool CSampleFile::Get(uint8_t &value, run_t &run, uint32_t no_pred_same, uint32_t no_succ_same)
{
	auto p = find_rc_coder(run, no_pred_same, no_succ_same);

	value = (uint8_t) p->Decode();

	return true;
}

// ************************************************************************************
bool CSampleFile::PutFlag(uint8_t flag) 
{
	auto p = find_rc_coder(ctx_flag);

	p->Encode(flag);

	ctx_flag = (ctx_flag << 2) + flag;
	ctx_flag &= ctx_flag_mask;

	return true;
}

// ************************************************************************************
bool CSampleFile::GetFlag(uint8_t &flag)
{
	auto p = find_rc_coder(ctx_flag);

	flag = p->Decode();

	ctx_flag = (ctx_flag << 2) + flag;
	ctx_flag &= ctx_flag_mask;

	return true;
}

// EOF
