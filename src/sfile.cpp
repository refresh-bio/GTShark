// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
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

	cout << "No. RC coders: " << rc_coders.get_size() << endl;

	if (rc)
		delete rc;

	delete vios;
}

// ************************************************************************************
bool CSampleFile::OpenForReading(string file_name)
{
	if (!fi_sample.Open(file_name))
	{
		cerr << "Cannot open " << file_name << " file\n";
		exit(1);
	}

	rc = new CRangeDecoder<CVectorIOStream>(*vios);
	rcd = (CRangeDecoder<CVectorIOStream>*) rc;

	size_t file_size = fi_sample.FileSize();

	file_size -= read_header_data();

	v_uint8.clear();
	v_uint8.reserve(file_size);

	while (!fi_sample.Eof())
		v_uint8.push_back(fi_sample.GetByte());

	rcd->Start();

	mode = mode_t::decompress;

	return true;
}

// ************************************************************************************
bool CSampleFile::OpenForWriting(string file_name)
{
	if (!fo_sample.Open(file_name))
	{
		cerr << "Cannot open " << file_name << " file\n";
		exit(1);
	}

	rc = new CRangeEncoder<CVectorIOStream>(*vios);
	rce = (CRangeEncoder<CVectorIOStream>*) rc;

	rce->Start();

	mode = mode_t::compress;

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

// EOF
