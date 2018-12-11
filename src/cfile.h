#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include <string>
#include <vector>
#include "vcf.h"
#include "io.h"
#include "pbwt.h"
#include "rc.h"
#include "sub_rc.h"
#include <unordered_map>
#include "context_hm.h"

using namespace std;

class CCompressedFile
{
	CInFile fi_db;
	CInFile fi_gt;
	COutFile fo_db;
	COutFile fo_gt;

	CRangeEncoder<COutFile> *rce;
	CRangeDecoder<CInFile> *rcd;

	CPBWT pbwt;
	bool pbwt_initialised;

	enum class open_mode_t {none, reading, writing} open_mode;

	vector<uint8_t> v_rd_header, v_cd_header;
	vector<uint8_t> v_rd_meta, v_cd_meta;
	vector<uint8_t> v_rd_samples, v_cd_samples;

	vector<uint8_t> v_rd_chrom, v_cd_chrom;
	vector<uint8_t> v_rd_pos, v_cd_pos;
	vector<uint8_t> v_rd_id, v_cd_id;
	vector<uint8_t> v_rd_ref, v_cd_ref;
	vector<uint8_t> v_rd_alt, v_cd_alt;
	vector<uint8_t> v_rd_qual, v_cd_qual;
	vector<uint8_t> v_rd_filter, v_cd_filter;
	vector<uint8_t> v_rd_info, v_cd_info;

	vector<uint8_t> v_rd_gt;
	vector<uint32_t> v_rle_gt;
	vector<pair<uint8_t, uint32_t>> v_rle_gt_large;

	size_t p_meta;
	size_t p_header;
	size_t p_samples;

	size_t p_chrom;
	size_t p_pos;
	size_t p_id;
	size_t p_ref;
	size_t p_alt;
	size_t p_qual;
	size_t p_filter;
	size_t p_info;

	uint32_t no_variants;
	uint32_t i_variant;
	uint32_t no_samples;
	uint8_t ploidy;
	uint32_t neglect_limit;
	string v_meta;
	string v_header;
	vector<string> v_samples;

	int64_t prev_pos;

	const context_t context_symbol_flag = 1ull << 60;
	const context_t context_symbol_mask = 0xffff;

	const context_t context_prefix_mask = 0xfffff;
	const context_t context_prefix_flag = 2ull << 60;
	const context_t context_suffix_flag = 3ull << 60;
	const context_t context_large_value1_flag = 4ull << 60;
	const context_t context_large_value2_flag = 5ull << 60;
	const context_t context_large_value3_flag = 6ull << 60;
	
	context_t ctx_prefix;
	context_t ctx_symbol;
	
	typedef CContextHM<CRangeCoderModel<COutFile>> ctx_map_e_t;
	typedef CContextHM<CRangeCoderModel<CInFile>> ctx_map_d_t;

	ctx_map_e_t rce_coders;
	ctx_map_d_t rcd_coders;

	inline ctx_map_e_t::value_type find_rce_coder(context_t ctx, uint32_t no_symbols, uint32_t max_log_counter);
	inline ctx_map_d_t::value_type find_rcd_coder(context_t ctx, uint32_t no_symbols, uint32_t max_log_counter);

	inline void encode_run_len(uint8_t symbol, uint32_t len);
	inline void decode_run_len(uint8_t &symbol, uint32_t &len);

	void append(vector<uint8_t> &v_comp, string x);
	void append(vector<uint8_t> &v_comp, int64_t x);

	void read(vector<uint8_t> &v_comp, size_t &pos, string &x);
	void read(vector<uint8_t> &v_comp, size_t &pos, int64_t &x);

	bool load_descriptions();
	bool save_descriptions();

public:
	CCompressedFile();
	~CCompressedFile();

	bool OpenForReading(string file_name);
	bool OpenForWriting(string file_name);

	bool Close();

	int GetNoSamples();
	void SetNoSamples(uint32_t _no_samples);
	int GetNoVariants();
	bool GetMeta(string &_v_meta);
	bool SetMeta(string &_v_meta);
	bool GetHeader(string &_v_header);
	bool SetHeader(string &_v_header);

	bool AddSamples(vector<string> &_v_samples);

	bool GetSamples(vector<string> &_v_samples);

	int GetPloidy();
	void SetPloidy(int _ploidy);

	int GetNeglectLimit();
	void SetNeglectLimit(uint32_t _neglect_limit);

	bool Eof();

	bool GetVariant(variant_desc_t &desc, vector<uint8_t> &data);
	bool SetVariant(variant_desc_t &desc, vector<uint8_t> &data);
	bool GetVariantGenotypesRaw(vector<pair<uint8_t, uint32_t>> &rle_genotypes);
	bool GetVariantGenotypesRawAndDesc(variant_desc_t &desc, vector<pair<uint8_t, uint32_t>> &rle_genotypes);

	bool InitPBWT();

	bool TrackItem(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, uint8_t &value, uint32_t &item_new_pos);
	bool TrackItems(const vector<pair<uint8_t, uint32_t>> &v_rle, array<uint32_t, 2> item_prev_pos, array<uint8_t, 2> &value, array<uint32_t, 2> &item_new_pos);

	bool EstimateValue(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, uint8_t value, array<pair<uint8_t, uint32_t>, 2> &runs, uint32_t &item_new_pos);

	bool RevertDecode(uint32_t &pos_sample_to_trace, const vector<pair<uint8_t, uint32_t>> &hist_rle_genotypes, const uint8_t reference_value);
};

// EOF
