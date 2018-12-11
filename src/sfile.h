#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include <unordered_map> 

#include "io.h"
#include "vios.h"
#include "rc.h"
#include "sub_rc.h"
#include "context_hm.h"

class CSampleFile
{
	CInFile fi_sample;
	COutFile fo_sample;

	string input_sample_name;
	vector<uint8_t> input_vec_rel_header;

	enum class mode_t {none, compress, decompress} mode;

	CBasicRangeCoder<CVectorIOStream> *rc;
	CRangeEncoder<CVectorIOStream> *rce;
	CRangeDecoder<CVectorIOStream> *rcd;

	vector<uint8_t> v_uint8;
	CVectorIOStream *vios;

	typedef CContextHM<CRangeCoderModel<CVectorIOStream>> ctx_map_e_t;

	ctx_map_e_t rc_coders;

	uint64_t determine_context(const run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);

	inline ctx_map_e_t::value_type find_rc_coder(const run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);

	uint32_t read_header_data();

public:
	CSampleFile();
	~CSampleFile();

	bool OpenForReading(string file_name);
	bool OpenForWriting(string file_name);
	bool Close();

	bool ReadHeaderAndSample(const string &db_header, string &v_header, string &sample_name);
	bool WriteHeaderAndSample(const string &db_header, const string &v_header, const string &sample_name);

	bool Put(uint8_t value, run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);
	bool Get(uint8_t &value, run_t &run, uint32_t no_pred_same, uint32_t no_succ_same);
};

// EOF
