// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include <memory>
#include <iostream>

using namespace std;

#include "cfile.h"
#include "lzma_wrapper.h"
#include "utils.h"

// ************************************************************************************
void CCompressedFile::append(vector<uint8_t> &v_comp, string x)
{
	v_comp.insert(v_comp.end(), x.begin(), x.end());
	v_comp.push_back(0);
}

// ************************************************************************************
void CCompressedFile::append(vector<uint8_t> &v_comp, int64_t x)
{
	append(v_comp, to_string(x));
}

// ************************************************************************************
void CCompressedFile::read(vector<uint8_t> &v_comp, size_t &pos, string &x)
{
	x.clear();

	for (; pos < v_comp.size() && v_comp[pos] != 0; ++pos)
		x.push_back(v_comp[pos]);
	++pos;
}

// ************************************************************************************
void CCompressedFile::read(vector<uint8_t> &v_comp, size_t &pos, int64_t &x)
{
	string t;

	read(v_comp, pos, t);
	x = stoll(t);
}

// ************************************************************************************
bool CCompressedFile::load_descriptions()
{
	string s;

	// Load file header and characteristics
	no_variants = (uint32_t) fi_db.ReadUInt(4);
	i_variant = 0;
	no_samples = (uint32_t) fi_db.ReadUInt(4);
	ploidy = (uint8_t) fi_db.ReadUInt(1);
	neglect_limit = (uint32_t) fi_db.ReadUInt(4);

	// Load variant descriptions
	for (auto d : {
		make_tuple(ref(v_rd_meta), ref(v_cd_meta), ref(p_meta), "meta"),
		make_tuple(ref(v_rd_header), ref(v_cd_header), ref(p_header), "header"),
		make_tuple(ref(v_rd_samples), ref(v_cd_samples), ref(p_samples), "samples"),
		make_tuple(ref(v_rd_chrom), ref(v_cd_chrom), ref(p_chrom), "chrom"),
		make_tuple(ref(v_rd_pos), ref(v_cd_pos), ref(p_pos), "pos"),
		make_tuple(ref(v_rd_id), ref(v_cd_id), ref(p_id), "id"),
		make_tuple(ref(v_rd_ref), ref(v_cd_ref), ref(p_ref), "ref"),
		make_tuple(ref(v_rd_alt), ref(v_cd_alt), ref(p_alt), "alt"),
		make_tuple(ref(v_rd_qual), ref(v_cd_qual), ref(p_qual), "qual"),
		make_tuple(ref(v_rd_filter), ref(v_cd_filter), ref(p_filter), "filter"),
		make_tuple(ref(v_rd_info), ref(v_cd_info), ref(p_info), "info")
		})
	{
		size_t field_len = fi_db.ReadUInt(4);

		get<1>(d).resize(field_len);
		fi_db.Read(get<1>(d).data(), field_len);
		CLZMAWrapper::Decompress(get<1>(d), get<0>(d));

		get<2>(d) = 0;
//		cout << get<3>(d) << " size: " << get<1>(d).size() << endl;
	}
	
	v_meta.clear();
	read(v_rd_meta, p_meta, v_meta);

	v_header.clear();
	read(v_rd_header, p_header, v_header);

	v_samples.clear();
	string sample;
	for (uint32_t i = 0; i < no_samples; ++i)
	{
		read(v_rd_samples, p_samples, sample);
		v_samples.push_back(sample);
	}

	return true;
}

// ************************************************************************************
bool CCompressedFile::save_descriptions()
{
	// Save file header and characteristics
	fo_db.WriteUInt(no_variants, 4);
	fo_db.WriteUInt(no_samples, 4);
	fo_db.WriteUInt(ploidy, 1);
	fo_db.WriteUInt(neglect_limit, 4);

	append(v_rd_meta, v_meta);
	append(v_rd_header, v_header);

	for (auto &x : v_samples)
		append(v_rd_samples, x);

	// Save variant descriptions
	for (auto d : {
		make_tuple(ref(v_rd_meta), ref(v_cd_meta), 9, "meta"),
		make_tuple(ref(v_rd_header), ref(v_cd_header), 9, "header"),
		make_tuple(ref(v_rd_samples), ref(v_cd_samples), 9, "samples"),
		make_tuple(ref(v_rd_chrom), ref(v_cd_chrom), 9, "chrom"),
		make_tuple(ref(v_rd_pos), ref(v_cd_pos), 9, "pos"),
		make_tuple(ref(v_rd_id), ref(v_cd_id), 9, "id"),
		make_tuple(ref(v_rd_ref), ref(v_cd_ref), 9, "ref"),
		make_tuple(ref(v_rd_alt), ref(v_cd_alt), 9, "alt"),
		make_tuple(ref(v_rd_qual), ref(v_cd_qual), 9, "qual"),
		make_tuple(ref(v_rd_filter), ref(v_cd_filter), 9, "filter"),
		make_tuple(ref(v_rd_info), ref(v_cd_info), 9, "info")
		})
	{
		CLZMAWrapper::Compress(get<0>(d), get<1>(d), get<2>(d));
		cout << get<3>(d) << " size: " << get<1>(d).size() << endl;
		fo_db.WriteUInt(get<1>(d).size(), 4);
		fo_db.Write(get<1>(d).data(), get<1>(d).size());
	}

	return true;
}

// ************************************************************************************
CCompressedFile::CCompressedFile()
{
	open_mode = open_mode_t::none;

	rce = nullptr;
	rcd = nullptr;
}

// ************************************************************************************
CCompressedFile::~CCompressedFile()
{
	if (rce)
		delete rce;
	if (rcd)
		delete rcd;

	fo_gt.Close();
	fo_db.Close();
}

// ************************************************************************************
bool CCompressedFile::OpenForReading(string file_name)
{
	prev_pos = 0;

	if (!fi_db.Open(file_name + "_db"))
	{
		cerr << "Cannot open " << file_name << "_db file\n";
		return false;
	}

	if (!fi_gt.Open(file_name + "_gt"))
	{
		cerr << "Cannot open " << file_name << "_gt file\n";
		return false;
	}

	open_mode = open_mode_t::reading;

	load_descriptions();
	pbwt_initialised = false;

	rcd = new CRangeDecoder<CInFile>(fi_gt);
	rcd->Start();

	return true;
}

// ************************************************************************************
bool CCompressedFile::OpenForWriting(string file_name)
{
	prev_pos = 0;

	if(!fo_db.Open(file_name + "_db"))
	{
		cerr << "Cannot open " << file_name << "_db file\n";
		return false;
	}

	if (!fo_gt.Open(file_name + "_gt"))
	{
		cerr << "Cannot open " << file_name << "_gt file\n";
		return false;
	}

	open_mode = open_mode_t::writing;
	pbwt_initialised = false;

	rce = new CRangeEncoder<COutFile>(fo_gt);
	rce->Start();

	no_variants = 0;

	return true;
}

// ************************************************************************************
bool CCompressedFile::Close()
{
	if (open_mode == open_mode_t::writing)
	{
		save_descriptions();
		rce->End();
		delete rce;
		rce = nullptr;

		fo_db.Close();
		fo_gt.Close();
	}
	else if (open_mode == open_mode_t::reading)
	{
		fi_db.Close();
		fi_gt.Close();
	}

	open_mode = open_mode_t::none;
	
	return true;
}

// ************************************************************************************
int CCompressedFile::GetNoSamples()
{
	return no_samples;
}

// ************************************************************************************
void CCompressedFile::SetNoSamples(uint32_t _no_samples)
{
	no_samples = _no_samples;
}

// ************************************************************************************
int CCompressedFile::GetNoVariants()
{
	return no_variants;
}

// ************************************************************************************
bool CCompressedFile::GetMeta(string &_v_meta)
{
	_v_meta = v_meta;

	return true;
}

// ************************************************************************************
bool CCompressedFile::SetMeta(string &_v_meta)
{
	v_meta = _v_meta;

	return true;
}

// ************************************************************************************
bool CCompressedFile::GetHeader(string &_v_header)
{
	_v_header = v_header;

	return true;
}

// ************************************************************************************
bool CCompressedFile::SetHeader(string &_v_header)
{
	v_header = _v_header;

	return true;
}

// ************************************************************************************
bool CCompressedFile::AddSamples(vector<string> &_v_samples)
{
	v_samples = _v_samples;

	return true;
}

// ************************************************************************************
bool CCompressedFile::GetSamples(vector<string> &_v_samples)
{
	_v_samples = v_samples;

	return true;
}

// ************************************************************************************
// Zrwaca ploidy
int CCompressedFile::GetPloidy()
{
	return ploidy;
}

// ************************************************************************************
// Ustawia ploidy
void CCompressedFile::SetPloidy(int _ploidy)
{
	ploidy = (uint8_t) _ploidy;
}

// ************************************************************************************
int CCompressedFile::GetNeglectLimit()
{
	return neglect_limit;
}

// ************************************************************************************
void CCompressedFile::SetNeglectLimit(uint32_t _neglect_limit)
{
	neglect_limit = _neglect_limit;
}

// ************************************************************************************
bool CCompressedFile::Eof()
{
	if (open_mode == open_mode_t::reading)
	{
		return i_variant >= no_variants;
	}

	return false;
}

// ************************************************************************************
bool CCompressedFile::GetVariant(variant_desc_t &desc, vector<uint8_t> &data)
{
	if (i_variant >= no_variants)
		return false;

	int64_t pos;

	// Load variant description
	read(v_rd_chrom, p_chrom, desc.chrom);
	read(v_rd_pos, p_pos, pos);
	pos += prev_pos;
	prev_pos = pos;
	desc.pos = pos;
	
	read(v_rd_id, p_id, desc.id);
	read(v_rd_ref, p_ref, desc.ref);
	read(v_rd_alt, p_alt, desc.alt);
	read(v_rd_qual, p_qual, desc.qual);
	read(v_rd_filter, p_filter, desc.filter);
	read(v_rd_info, p_info, desc.info);

	// Load genotypes
//	v_rd_gt.clear();
	v_rle_gt.clear();
	v_rle_gt_large.clear();
//	data.clear();

	uint32_t total_len = 0;
	ctx_prefix = context_prefix_mask;
	ctx_symbol = context_symbol_mask;

	while (total_len < no_samples * ploidy)
	{
		uint8_t symbol;
		uint32_t len;
		decode_run_len(symbol, len);

		if (len == 0u)
			len = no_samples * ploidy - total_len;

		v_rle_gt_large.push_back(make_pair(symbol, len));

		total_len += len;
	}

	pbwt.Decode(v_rle_gt_large, v_rd_gt);

	data.resize(no_samples);

	if (ploidy == 1)
	{ 
		for (uint32_t i = 0; i < no_samples; ++i)
			data[i] = v_rd_gt[i];
	}
	else if (ploidy == 2)
	{
		for(uint32_t i = 0; i < no_samples * ploidy; ++i)

		if (i % 2 == 0)
			data[i / 2] = 0b00010000 + v_rd_gt[i];
		else
			data[i / 2] += (uint8_t)(v_rd_gt[i] << 2u);
	}

	++i_variant;

	return true;
}

// ************************************************************************************
bool CCompressedFile::SetVariant(variant_desc_t &desc, vector<uint8_t> &data)
{
	// Store variant description
	append(v_rd_chrom, desc.chrom);
	append(v_rd_pos, desc.pos - prev_pos);
	prev_pos = desc.pos;
	append(v_rd_id, desc.id);
	append(v_rd_ref, desc.ref);
	append(v_rd_alt, desc.alt);
	append(v_rd_qual, desc.qual);
	append(v_rd_filter, desc.filter);
	append(v_rd_info, desc.info);

	// Store genotypes
	v_rd_gt.resize(no_samples * ploidy);

	if (ploidy == 1)
		for (uint32_t i = 0; i < no_samples; ++i)
			v_rd_gt[i] = (data[i] & 0b00000011);
	else if(ploidy == 2)
		for (uint32_t i = 0; i < no_samples; ++i)
		{
			v_rd_gt[2 * i + 0] = (data[i] & 0b00000011);
			v_rd_gt[2 * i + 1] = ((data[i] >> 2) & 0b00000011);
		}

	pbwt.Encode(v_rd_gt, v_rle_gt_large);
	ctx_prefix = context_prefix_mask;
	ctx_symbol = context_symbol_mask;

	v_rle_gt_large.back().second = 0u;

	for (auto x : v_rle_gt_large)
		encode_run_len(x.first, x.second);

	++no_variants;

	return true;
}

// ************************************************************************************
bool CCompressedFile::GetVariantGenotypesRaw(vector<pair<uint8_t, uint32_t>> &rle_genotypes)
{
	if (i_variant >= no_variants)
		return false;

	rle_genotypes.clear();

	uint32_t total_len = 0;
	ctx_prefix = context_prefix_mask;
	ctx_symbol = context_symbol_mask;

	while (total_len < no_samples * ploidy)
	{
		uint8_t symbol;
		uint32_t len;

		decode_run_len(symbol, len);
		if (len == 0)
			len = no_samples * ploidy - total_len;
		
		rle_genotypes.push_back(make_pair(symbol, len));

		total_len += len;
	}

	++i_variant;

	return true;
}

// ************************************************************************************
bool CCompressedFile::GetVariantGenotypesRawAndDesc(variant_desc_t &desc, vector<pair<uint8_t, uint32_t>> &rle_genotypes)
{
	if (i_variant >= no_variants)
		return false;

	rle_genotypes.clear();

	int64_t pos;

	// Load variant description
	read(v_rd_chrom, p_chrom, desc.chrom);
	read(v_rd_pos, p_pos, pos);
	pos += prev_pos;
	prev_pos = pos;
	desc.pos = pos;

	read(v_rd_id, p_id, desc.id);
	read(v_rd_ref, p_ref, desc.ref);
	read(v_rd_alt, p_alt, desc.alt);
	read(v_rd_qual, p_qual, desc.qual);
	read(v_rd_filter, p_filter, desc.filter);
	read(v_rd_info, p_info, desc.info);

	uint32_t total_len = 0;
	ctx_prefix = context_prefix_mask;
	ctx_symbol = context_symbol_mask;

	while (total_len < no_samples * ploidy)
	{
		uint8_t symbol;
		uint32_t len;

		decode_run_len(symbol, len);
		if (len == 0)
			len = no_samples * ploidy - total_len;

		rle_genotypes.push_back(make_pair(symbol, len));

		total_len += len;
	}

	++i_variant;

	return true;
}

// ************************************************************************************
bool CCompressedFile::InitPBWT()
{
	if (open_mode == open_mode_t::reading)
	{
		pbwt.StartReverse(no_samples * ploidy, neglect_limit);
		pbwt_initialised = true;
	}
	else if(open_mode == open_mode_t::writing)
	{
		pbwt.StartForward(no_samples * ploidy, neglect_limit);
		pbwt_initialised = true;

		v_rd_gt.reserve(no_samples * ploidy);
	}

	return pbwt_initialised;
}

// ************************************************************************************
bool CCompressedFile::TrackItem(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, 
	uint8_t &value, uint32_t &item_new_pos)
{
	return pbwt.TrackItem(v_rle, item_prev_pos, value, item_new_pos);
}

// ************************************************************************************
bool CCompressedFile::TrackItems(const vector<pair<uint8_t, uint32_t>> &v_rle, array<uint32_t, 2> item_prev_pos, 
	array<uint8_t, 2> &value, array<uint32_t, 2> &item_new_pos)
{
	return pbwt.TrackItems(v_rle, item_prev_pos, value, item_new_pos);
}

// ************************************************************************************
bool CCompressedFile::EstimateValue(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, 
	uint8_t value, array<pair<uint8_t, uint32_t>, 2> &runs, uint32_t &item_new_pos)
{
	return pbwt.EstimateValue(v_rle, item_prev_pos, value, runs, item_new_pos);
}

// ************************************************************************************
bool CCompressedFile::RevertDecode(uint32_t &pos_sample_to_trace, const vector<pair<uint8_t, uint32_t>> &hist_rle_genotypes, const uint8_t reference_value)
{
	return pbwt.RevertDecode(pos_sample_to_trace, hist_rle_genotypes, reference_value);
}

// ************************************************************************************
CCompressedFile::ctx_map_e_t::value_type CCompressedFile::find_rce_coder(context_t ctx, uint32_t no_symbols, uint32_t max_log_counter)
{
	auto p = rce_coders.find(ctx);

	if (p == nullptr)
		rce_coders.insert(ctx, p = new CRangeCoderModel<COutFile>(rce, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, 1, true));

	return p;
}

// ************************************************************************************
CCompressedFile::ctx_map_d_t::value_type CCompressedFile::find_rcd_coder(context_t ctx, uint32_t no_symbols, uint32_t max_log_counter)
{
	auto p = rcd_coders.find(ctx);

	if (p == nullptr)
		rcd_coders.insert(ctx, p = new CRangeCoderModel<CInFile>(rcd, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, 1, false));

	return p;
}

// ************************************************************************************
void CCompressedFile::encode_run_len(uint8_t symbol, uint32_t len)
{
	// Encode symbol
	auto rc_sym = find_rce_coder(ctx_symbol + context_symbol_flag, 4, 15);
	rc_sym->Encode(symbol);
	ctx_symbol <<= 4;
	ctx_symbol += symbol;
	ctx_symbol &= context_symbol_mask;

	rce_coders.prefetch(ctx_symbol + context_symbol_flag);

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)symbol;
	ctx_prefix &= context_prefix_mask;

	// Encode run length
	auto rc_p = find_rce_coder(ctx_prefix + context_prefix_flag, 11, 10);

	uint32_t prefix = ilog2(len);

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)prefix;
	ctx_prefix &= context_prefix_mask;

	rce_coders.prefetch(ctx_prefix + context_prefix_flag);

	if (prefix < 2)
		rc_p->Encode(prefix);
	else if (prefix < 10)
	{
		rc_p->Encode(prefix);
		uint64_t ctx_suf = context_suffix_flag;
		ctx_suf += ((context_t)symbol) << 8;
		ctx_suf += (context_t) prefix;
		uint32_t max_value_for_this_prefix = 1u << (prefix - 1);

		auto rc_s = find_rce_coder(ctx_suf, max_value_for_this_prefix, 15);
		rc_s->Encode(len - max_value_for_this_prefix);
	}
	else
	{
		rc_p->Encode(10);		// flag for large value

		context_t ctx_large1 = context_large_value1_flag;
		ctx_large1 += ((context_t)symbol) << 16;
		auto rc_l1 = find_rce_coder(ctx_large1, 256, 15);
		uint32_t lv1 = (len >> 16) & 0xff;
		rc_l1->Encode(lv1);

		context_t ctx_large2 = context_large_value2_flag;
		ctx_large2 += ((context_t)symbol) << 16;
		ctx_large2 += (context_t) lv1;
		auto rc_l2 = find_rce_coder(ctx_large2, 256, 15);
		uint32_t lv2 = (len >> 8) & 0xff;
		rc_l2->Encode(lv2);

		context_t ctx_large3 = context_large_value3_flag;
		ctx_large3 += ((context_t)symbol) << 16;
		ctx_large3 += ((context_t) lv1) << 8;
		ctx_large3 += (context_t) lv2;
		auto rc_l3 = find_rce_coder(ctx_large3, 256, 15);
		uint32_t lv3 = len & 0xff;
		rc_l3->Encode(lv3);
	}
}

// ************************************************************************************
void CCompressedFile::decode_run_len(uint8_t &symbol, uint32_t &len)
{
	// Decode symbol
	auto rc_sym = find_rcd_coder(ctx_symbol + context_symbol_flag, 4, 15);
	symbol = (uint8_t) rc_sym->Decode();
	ctx_symbol <<= 4;
	ctx_symbol += (context_t) symbol;
	ctx_symbol &= context_symbol_mask;

	rcd_coders.prefetch(ctx_symbol + context_symbol_flag);

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)symbol;
	ctx_prefix &= context_prefix_mask;

	// Decode run length
	auto rc_p = find_rcd_coder(ctx_prefix + context_prefix_flag, 11, 10);

	uint32_t prefix = rc_p->Decode();

	if (prefix < 2)
		len = prefix;
	else if (prefix < 10)
	{
		uint64_t ctx_suf = context_suffix_flag;
		ctx_suf += ((context_t)symbol) << 8;
		ctx_suf += (context_t)prefix;
		uint32_t max_value_for_this_prefix = 1u << (prefix - 1);

		auto rc_s = find_rcd_coder(ctx_suf, max_value_for_this_prefix, 15);
		len = max_value_for_this_prefix + rc_s->Decode();
	}
	else
	{
		context_t ctx_large1 = context_large_value1_flag;
		ctx_large1 += ((context_t)symbol) << 16;
		auto rc_l1 = find_rcd_coder(ctx_large1, 256, 15);
		uint32_t lv1 = rc_l1->Decode();

		context_t ctx_large2 = context_large_value2_flag;
		ctx_large2 += ((context_t)symbol) << 16;
		ctx_large2 += (context_t) lv1;
		auto rc_l2 = find_rcd_coder(ctx_large2, 256, 15);
		uint32_t lv2 = rc_l2->Decode();

		context_t ctx_large3 = context_large_value3_flag;
		ctx_large3 += ((context_t)symbol) << 16;
		ctx_large3 += ((context_t) lv1) << 8;
		ctx_large3 += (context_t) lv2;
		auto rc_l3 = find_rcd_coder(ctx_large3, 256, 15);
		uint32_t lv3 = rc_l3->Decode();

		len = (lv1 << 16) + (lv2 << 8) + lv3;

		prefix = ilog2(len);
	}

	ctx_prefix <<= 4;
	ctx_prefix += (context_t)prefix;
	ctx_prefix &= context_prefix_mask;

	rcd_coders.prefetch(ctx_prefix + context_prefix_flag);
}

// EOF
