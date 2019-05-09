// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

// This code is based on the example from XZ library

#include "lzma_wrapper.h"
#include <algorithm>
#include <iostream>

// *******************************************************************************************
void CLZMAWrapper::Compress(const vector<uint8_t> &v_text, vector<uint8_t> &v_text_compressed, int compression_mode)
{
	if (v_text.empty())
		return;

	lzma_stream strm = LZMA_STREAM_INIT;

	bool success = init_encoder(&strm, compression_mode);
	if (success)
		success = compress_impl(&strm, v_text, v_text_compressed);

	lzma_end(&strm);
}

// *******************************************************************************************
void CLZMAWrapper::Decompress(const vector<uint8_t> &v_text_compressed, vector<uint8_t> &v_text)
{
	lzma_stream strm = LZMA_STREAM_INIT;

	bool success = init_decoder(&strm);
	if (success)
		success = decompress_impl(&strm, v_text_compressed, v_text);

	lzma_end(&strm);
}

// *******************************************************************************************
void CLZMAWrapper::CompressWithHistory(const vector<uint8_t> &v_history, const vector<uint8_t> &v_text, vector<uint8_t> &v_text_compressed, int compression_mode)
{
	vector<uint8_t> v_history_compressed;
	vector<uint8_t> v_combined_input;
	vector<uint8_t> v_combined_compressed;

	Compress(v_history, v_history_compressed, compression_mode);

	v_combined_input = v_history;
	v_combined_input.push_back(0u);
	v_combined_input.insert(v_combined_input.end(), v_text.begin(), v_text.end());

	Compress(v_combined_input, v_combined_compressed, compression_mode);

	size_t same_symbols;

	for (same_symbols = 60; same_symbols < v_history_compressed.size(); ++same_symbols)
		if (v_history_compressed[same_symbols] != v_combined_compressed[same_symbols])
			break;

	v_text_compressed.clear();
	v_text_compressed.push_back((same_symbols >> 24) & 0xff);
	v_text_compressed.push_back((same_symbols >> 16) & 0xff);
	v_text_compressed.push_back((same_symbols >>  8) & 0xff);
	v_text_compressed.push_back((same_symbols      ) & 0xff);

	v_text_compressed.insert(v_text_compressed.end(), v_combined_compressed.begin(), v_combined_compressed.begin() + 60);
	v_text_compressed.insert(v_text_compressed.end(), v_combined_compressed.begin() + same_symbols, v_combined_compressed.end());
}

// *******************************************************************************************
void CLZMAWrapper::DecompressWithHistory(const vector<uint8_t> &v_history, const vector<uint8_t> &v_text_compressed, vector<uint8_t> &v_text, int compression_mode)
{
	vector<uint8_t> v_history_compressed;

	vector<uint8_t> v_combined_input;
	vector<uint8_t> v_combined_compressed;

	Compress(v_history, v_history_compressed, compression_mode);

	size_t same_symbols = 0u;

	same_symbols += ((size_t) v_text_compressed[0]) << 24;
	same_symbols += ((size_t) v_text_compressed[1]) << 16;
	same_symbols += ((size_t) v_text_compressed[2]) <<  8;
	same_symbols += ((size_t) v_text_compressed[3]);

	v_combined_compressed.assign(v_text_compressed.begin() + 4, v_text_compressed.begin() + 64);
	v_combined_compressed.insert(v_combined_compressed.end(), v_history_compressed.begin() + 60, v_history_compressed.begin() + same_symbols);
	v_combined_compressed.insert(v_combined_compressed.end(), v_text_compressed.begin() + 64, v_text_compressed.end());

	Decompress(v_combined_compressed, v_combined_input);

	v_text.assign(v_combined_input.begin() + v_history.size() + 1u, v_combined_input.end());
}

// *******************************************************************************************
// Initialize LZMA coder
bool CLZMAWrapper::init_encoder(lzma_stream *strm, uint32_t preset)
{
	lzma_ret ret = lzma_easy_encoder(strm, preset, LZMA_CHECK_CRC64);

	if (ret == LZMA_OK)
		return true;

	cerr << "Some bug in LZMA stage\n";

	return false;
}

// *******************************************************************************************
// Do actual compression
bool CLZMAWrapper::compress_impl(lzma_stream *strm, const vector<uint8_t> &v_in, vector<uint8_t> &v_out)
{
	lzma_action action = LZMA_RUN;

	uint8_t inbuf[BUFSIZ];
	uint8_t outbuf[BUFSIZ];

	strm->next_in = NULL;
	strm->avail_in = 0;
	strm->next_out = outbuf;
	strm->avail_out = sizeof(outbuf);

	size_t in_size = v_in.size();
	size_t in_pos = 0;

	while (true) {
		if (strm->avail_in == 0 && in_pos < in_size) {
			strm->next_in = inbuf;
			size_t to_read = std::min((size_t) sizeof(inbuf), in_size - in_pos);
			copy_n(v_in.begin() + in_pos, to_read, inbuf);
			in_pos += to_read;
			strm->avail_in = to_read;

			if(in_pos == in_size)
				action = LZMA_FINISH;
		}

		lzma_ret ret = lzma_code(strm, action);

		if (strm->avail_out == 0 || ret == LZMA_STREAM_END) {
			size_t write_size = sizeof(outbuf) - strm->avail_out;

			for (size_t i = 0; i < write_size; ++i)
				v_out.push_back(outbuf[i]);

			// Reset next_out and avail_out.
			strm->next_out = outbuf;
			strm->avail_out = sizeof(outbuf);
		}

		if (ret != LZMA_OK) {
			if (ret == LZMA_STREAM_END)
				return true;
	
			cerr << "Some bug in LZMA compression\n";

			return false;
		}
	}
}

// *******************************************************************************************
// Initialize LZMA decoder
bool CLZMAWrapper::init_decoder(lzma_stream *strm)
{
	lzma_ret ret = lzma_stream_decoder(strm, UINT64_MAX, LZMA_CONCATENATED);

	if (ret == LZMA_OK)
		return true;

	cerr << "Some bug in LZMA compression\n";

	return false;
}

// *******************************************************************************************
// Do actual decompression
bool CLZMAWrapper::decompress_impl(lzma_stream *strm, const vector<uint8_t> &v_in, vector<uint8_t> &v_out)
{
	lzma_action action = LZMA_RUN;

	uint8_t inbuf[BUFSIZ];
	uint8_t outbuf[BUFSIZ];

	strm->next_in = NULL;
	strm->avail_in = 0;
	strm->next_out = outbuf;
	strm->avail_out = sizeof(outbuf);

	size_t in_pos = 0;
	size_t in_size = v_in.size();

	while (true) {
		if (strm->avail_in == 0 && in_pos < in_size) {
			strm->next_in = inbuf;

			size_t to_read = std::min((size_t) sizeof(inbuf), in_size - in_pos);
			copy_n(v_in.begin() + in_pos, to_read, inbuf);
			in_pos += to_read;
			strm->avail_in = to_read;
			
			if(in_pos == in_size)
				action = LZMA_FINISH;
		}

		lzma_ret ret = lzma_code(strm, action);

		if (strm->avail_out == 0 || ret == LZMA_STREAM_END) {
			size_t write_size = sizeof(outbuf) - strm->avail_out;

			for (size_t i = 0; i < write_size; ++i)
				v_out.push_back(outbuf[i]);

			strm->next_out = outbuf;
			strm->avail_out = sizeof(outbuf);
		}

		if (ret != LZMA_OK) {
			if (ret == LZMA_STREAM_END)
				return true;

			return false;
		}
	}

	return true;
}

// EOF
