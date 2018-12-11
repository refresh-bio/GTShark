#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include <vector>
#include "defs.h"

using namespace std;

class CPBWT
{
	size_t no_items;
	size_t neglect_limit;

	vector<int> v_perm_cur;
	vector<int> v_perm_prev;
	vector<uint8_t> v_tmp;

	vector<uint32_t> v_hist_complete;

public:
	CPBWT();
	~CPBWT();

	bool StartForward(const size_t _no_items, const size_t _neglect_limit);
	bool StartReverse(const size_t _no_items, const size_t _neglect_limit);

	bool Encode(vector<uint8_t> &v_input, vector<pair<uint8_t, uint32_t>> &v_rle);
	bool Decode(const vector<pair<uint8_t, uint32_t>> &v_rle, vector<uint8_t> &v_output);

	bool TrackItem(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, uint8_t &value, uint32_t &item_new_pos);
	bool TrackItems(const vector<pair<uint8_t, uint32_t>> &v_rle, array<uint32_t, 2> item_prev_pos, array<uint8_t, 2> &value, array<uint32_t, 2> &item_new_pos);

	bool RevertDecode(uint32_t &pos_sample_to_trace, const vector<pair<uint8_t, uint32_t>> &hist_rle_genotypes, const uint8_t reference_value);

	bool EstimateValue(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, uint8_t value, array<pair<uint8_t, uint32_t>, 2> &runs, uint32_t &item_new_pos);
};

// EOF
