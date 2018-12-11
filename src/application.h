#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <list>
#include <deque>
#include <memory>

#include "params.h"
#include "vcf.h"
#include "cfile.h"
#include "sfile.h"

using namespace std;

// ******************************************************************************
class CApplication
{
	const size_t no_variants_in_buf = 8192u;
	typedef pair<uint8_t, uint32_t> run_desc_t;

	list<vector<run_desc_t>> l_hist_rle_genotypes;
	list<array<uint8_t, 2>> l_hist_tracked_samples;

	const size_t max_tracked_dist = 2048;

	deque<vector<run_desc_t>> d_hist_rle_genotypes;
	deque<array<uint8_t, 2>> d_hist_tracked_samples;

	CParams params;

	vector<pair<variant_desc_t, vector<uint8_t>>> v_vcf_data_compress, v_vcf_data_io;

	vector<tuple<uint8_t, run_t, uint32_t, uint32_t>> v_sample_data_compress, v_sample_data_io;
	vector<pair<variant_desc_t, uint8_t>> v_sample_d_data_compress, v_sample_d_data_io;

	mutex mtx;
	condition_variable cv;

	bool find_prev_value(const vector<run_desc_t> &v_rle_genotypes, const uint32_t max_pos, const uint8_t value, uint32_t &found_pos);
	bool find_next_value(const vector<run_desc_t> &v_rle_genotypes, const uint32_t min_pos, const uint8_t value, uint32_t &found_pos);

public:
	CApplication(const CParams &_params);
	~CApplication();

	bool CompressDB();
	bool DecompressDB();
	bool CompressSample();
	bool DecompressSample();
	bool ExtractSample();
};

// EOF
