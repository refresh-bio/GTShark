#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

#include <unordered_map> 

#include "io.h"
#include "vios.h"
#include "rc.h"
#include "vcf.h"
#include "sub_rc.h"
#include "context_hm.h"

class CSampleFile
{
	CInFile fi_sample;
	COutFile fo_sample;

	string input_sample_name;
	vector<uint8_t> input_vec_rel_header;
	bool extra_variants;

	enum class mode_t {none, compress, decompress} mode;

	CBasicRangeCoder<CVectorIOStream> *rc;
	CRangeEncoder<CVectorIOStream> *rce;
	CRangeDecoder<CVectorIOStream> *rcd;

	vector<uint8_t> v_uint8;
	CVectorIOStream *vios;

	typedef CContextHM<CRangeCoderModel<CVectorIOStream>> ctx_map_e_t;

	ctx_map_e_t rc_coders;
	uint64_t ctx_flag;
	const uint64_t ctx_flag_mask = 0xfffull;

	uint64_t determine_context(const run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);

	inline ctx_map_e_t::value_type find_rc_coder(const run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);
	inline ctx_map_e_t::value_type find_rc_coder(context_t ctx);

	uint32_t read_header_data();
	uint32_t read_extra_variants();

	vector<pair<variant_desc_t, vector<uint8_t>>> loc_v_desc;
	string get_string_from_vector(vector<uint8_t>::iterator& p);

public:
	CSampleFile();
	~CSampleFile();

	bool OpenForReading(string file_name, bool &_extra_variants);
	bool OpenForWriting(string file_name, bool _extra_variants);
	bool Close();

	bool ReadHeaderAndSample(const string &db_header, string &v_header, string &sample_name);
	bool WriteHeaderAndSample(const string &db_header, const string &v_header, const string &sample_name);

	bool Put(uint8_t value, run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);
	bool Get(uint8_t &value, run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);

	bool PutFlag(uint8_t flag);
	bool GetFlag(uint8_t &flag);

	void ReadExtraVariants(vector<pair<variant_desc_t, vector<uint8_t>>>& v_desc);
	uint32_t WriteExtraVariants(const vector<pair<variant_desc_t, vector<uint8_t>>>& v_desc);
};

// EOF
