// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include "pbwt.h"
#include <numeric>
#include <algorithm>
#include <iterator>
#include "utils.h"

// ************************************************************************************
CPBWT::CPBWT()
{

}

// ************************************************************************************
CPBWT::~CPBWT()
{

}

// ************************************************************************************
bool CPBWT::StartForward(const size_t _no_items, const size_t _neglect_limit)
{
	no_items = _no_items;
	neglect_limit = _neglect_limit;

	v_perm_cur.resize(no_items);
	v_tmp.resize(no_items);

	iota(v_perm_cur.begin(), v_perm_cur.end(), 0);
	v_perm_prev = v_perm_cur;
	
	return true;
}

// ************************************************************************************
bool CPBWT::StartReverse(const size_t _no_items, const size_t _neglect_limit)
{
	no_items = _no_items;
	neglect_limit = _neglect_limit;

	v_perm_cur.resize(no_items);

	iota(v_perm_cur.begin(), v_perm_cur.end(), 0);
	v_perm_prev = v_perm_cur;
	
	v_tmp.clear();
	v_tmp.resize(no_items, 0u);

	return true;
}

// ************************************************************************************
// Forward PBWT for non-binary alphabet
bool CPBWT::Encode(vector<uint8_t> &v_input, vector<pair<uint8_t, uint32_t>> &v_rle)
{
	vector<uint32_t> v_hist(SIGMA);
	uint32_t max_count;

	v_rle.clear();

	// Determine histogram of symbols
	calc_cumulate_histogram(v_input, v_hist, max_count);

	uint8_t prev_symbol = v_input[v_perm_prev[0]];
	uint32_t run_len = 0;

	// Make PBWT
	for (size_t i = 0; i < no_items; ++i)
	{
		uint8_t cur_symbol = v_input[v_perm_prev[i]];

		if (cur_symbol == prev_symbol)
			++run_len;
		else
		{

			v_rle.push_back(make_pair(prev_symbol, run_len));
			prev_symbol = cur_symbol;
			run_len = 1;
		}

		v_perm_cur[v_hist[cur_symbol]] = v_perm_prev[i];
		++v_hist[cur_symbol];
	}

	v_rle.push_back(make_pair(prev_symbol, run_len));

	// Swap only if no. of non-zeros is larger than neglect_limit
	if (no_items - max_count >= neglect_limit)
		swap(v_perm_prev, v_perm_cur);

	return true;
}

// ************************************************************************************
// Reverse PBWT for non-binary alphabet
bool CPBWT::Decode(const vector<pair<uint8_t, uint32_t>> &v_rle, vector<uint8_t> &v_output)
{
	vector<uint32_t> v_hist(SIGMA);
	uint32_t max_count;

	v_output.resize(no_items);

	calc_cumulate_histogram(v_rle, v_hist, max_count);

	auto p_rle = v_rle.begin();
	uint8_t cur_symbol = p_rle->first;
	uint32_t cur_cnt = p_rle->second;

	// Make PBWT
	for (size_t i = 0; i < no_items; ++i)
	{
		v_output[v_perm_prev[i]] = cur_symbol;

		v_perm_cur[v_hist[cur_symbol]] = v_perm_prev[i];
		++v_hist[cur_symbol];

		if (--cur_cnt == 0)
		{
			++p_rle;
			cur_symbol = p_rle->first;
			cur_cnt = p_rle->second;
		}
	}

	// Swap only if no. of non-zeros is larger than neglect_limit
	if (no_items - max_count >= neglect_limit)
		swap(v_perm_prev, v_perm_cur);

	return true;
}

// ************************************************************************************
bool CPBWT::TrackItem(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, uint8_t &value, uint32_t &item_new_pos)
{
	// Determine histogram of symbols in run-encoded data
	vector<uint32_t> v_hist_complete(SIGMA, 0u);
	vector<uint32_t> v_hist_partial(SIGMA, 0u);
	uint32_t max_count;

	for (auto x : v_rle)
		v_hist_complete[x.first] += x.second;

	cumulate_sums(v_hist_complete, max_count);

	int cur_pos = 0;
	for (auto x : v_rle)
	{
		if (item_prev_pos < cur_pos + x.second)
		{
			value = x.first;
			v_hist_partial[x.first] += item_prev_pos - cur_pos;
			break;
		}
		else
			v_hist_partial[x.first] += x.second;

		cur_pos += x.second;
	}

	// Swap only if no. of non-zeros is larger than neglect_limit
	if (no_items - max_count >= neglect_limit)
		item_new_pos = v_hist_complete[value] + v_hist_partial[value];
	else
		item_new_pos = item_prev_pos;

	return true;
}

// ************************************************************************************
bool CPBWT::TrackItems(const vector<pair<uint8_t, uint32_t>> &v_rle, array<uint32_t, 2> item_prev_pos, array<uint8_t, 2> &value, array<uint32_t, 2> &item_new_pos)
{
	// Determine histogram of symbols in run-encoded data
	vector<uint32_t> v_hist_complete(SIGMA, 0u);
	vector<uint32_t> v_hist_partial;
	uint32_t max_count;

	for (auto x : v_rle)
		v_hist_complete[x.first] += x.second;

	cumulate_sums(v_hist_complete, max_count);

	for (uint32_t i = 0; i < 2; ++i)
	{
		int cur_pos = 0;
		v_hist_partial.clear();
		v_hist_partial.resize(SIGMA, 0u);

		for (auto x : v_rle)
		{
			if (item_prev_pos[i] < cur_pos + x.second)
			{
				value[i] = x.first;
				v_hist_partial[x.first] += item_prev_pos[i] - cur_pos;
				break;
			}
			else
				v_hist_partial[x.first] += x.second;

			cur_pos += x.second;
		}

		// Swap only if no. of non-zeros is larger than neglect_limit
		if (no_items - max_count >= neglect_limit)
			item_new_pos[i] = v_hist_complete[value[i]] + v_hist_partial[value[i]];
		else
			item_new_pos[i] = item_prev_pos[i];
	}

	return true;
}

// ************************************************************************************
bool CPBWT::EstimateValue(const vector<pair<uint8_t, uint32_t>> &v_rle, uint32_t item_prev_pos, uint8_t value, array<pair<uint8_t, uint32_t>, 2> &runs, uint32_t &item_new_pos)
{
	uint32_t max_count;

	if (v_hist_complete.size() != SIGMA)
	{
		v_hist_complete.clear();
		v_hist_complete.resize(SIGMA, 0u);
	}
	else
		fill(v_hist_complete.begin(), v_hist_complete.end(), 0u);

	for (auto x : v_rle)
		v_hist_complete[x.first] += x.second;

	cumulate_sums(v_hist_complete, max_count);

	fill_n(runs.begin(), 2, make_pair(0u, 0u));
	int cur_pos = 0;
	uint32_t counter_for_value = 0;

	if (item_prev_pos == 0)
	{
		runs[0] = make_pair(0u, 0u);
		runs[1] = v_rle[0];
	}
	else
		for (uint32_t i = 0; i < v_rle.size(); ++i)
		{
			if (item_prev_pos == cur_pos + v_rle[i].second)
			{
				runs[0] = v_rle[i];
				if (i + 1 < v_rle.size())
					runs[1] = v_rle[i+1];

				if(value == v_rle[i].first)
					counter_for_value += v_rle[i].second;

				break;
			}
			else if (item_prev_pos < cur_pos + v_rle[i].second)
			{
				runs[0].first = runs[1].first = v_rle[i].first;
				runs[0].second = item_prev_pos - cur_pos;
				runs[1].second = v_rle[i].second - runs[0].second;
			
				if (value == v_rle[i].first)
					counter_for_value += runs[0].second;
			
				break;
			}

			cur_pos += v_rle[i].second;
			if (value == v_rle[i].first)
				counter_for_value += v_rle[i].second;
		}

	// Swap only if no. of non-zeros is larger than neglect_limit
	if (no_items - max_count >= neglect_limit)
		item_new_pos = v_hist_complete[value] + counter_for_value;
	else
		item_new_pos = item_prev_pos;

	return true;
}

// ************************************************************************************
bool CPBWT::RevertDecode(uint32_t &pos_sample_to_trace, const vector<pair<uint8_t, uint32_t>> &hist_rle_genotypes, const uint8_t reference_value)
{
	array<uint32_t, SIGMA> a_hist;
	uint32_t max_count;

	fill_n(a_hist.begin(), SIGMA, 0u);
	for (auto &x : hist_rle_genotypes)
		a_hist[x.first] += x.second;
	cumulate_sums(a_hist, max_count);

	uint8_t value = SIGMA - 1;

	for(uint32_t i = 1; i < SIGMA; ++i)
		if (pos_sample_to_trace < a_hist[i])
		{
			value = (uint8_t) (i - 1);
			break;
		}

	if (value != reference_value)
		return false;

	uint32_t new_pos_sample = a_hist[value];

	uint32_t cur_pos = 0;

	for (auto &x : hist_rle_genotypes)
	{
		if (x.first != value)
			cur_pos += x.second;
		else if (cur_pos + x.second < pos_sample_to_trace)
		{
			cur_pos += x.second;
			new_pos_sample += x.second;
		}
		else
		{
			new_pos_sample += pos_sample_to_trace - cur_pos;
			cur_pos += x.second;
		}
		
		if (cur_pos >= pos_sample_to_trace)
			break;
	}

	pos_sample_to_trace = new_pos_sample;

	return true;
}

// EOF
