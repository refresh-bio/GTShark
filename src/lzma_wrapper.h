#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

#include <vector>
#include <string>

#define LZMA_API_STATIC
#include "lzma.h"

using namespace std;

// *******************************************************************************************
//
// *******************************************************************************************
class CLZMAWrapper 
{
	void forward();
	void reverse();

	static bool init_encoder(lzma_stream *strm, uint32_t preset);
	static bool init_decoder(lzma_stream *strm);

	static bool compress_impl(lzma_stream *strm, const vector<uint8_t> &v_in, vector<uint8_t> &v_out);
	static bool decompress_impl(lzma_stream *strm, const vector<uint8_t> &v_in, vector<uint8_t> &v_out);

public:
	CLZMAWrapper() {};
	~CLZMAWrapper() {};

	static void Compress(const vector<uint8_t> &v_text, vector<uint8_t> &v_text_compressed, int compression_mode = 0);
	static void Decompress(const vector<uint8_t> &v_text_compressed, vector<uint8_t> &v_text);

	static void CompressWithHistory(const vector<uint8_t> &v_history, const vector<uint8_t> &v_text, vector<uint8_t> &v_text_compressed, int compression_mode = 0);
	static void DecompressWithHistory(const vector<uint8_t> &v_history, const vector<uint8_t> &v_text_compressed, vector<uint8_t> &v_text, int compression_mode = 0);
};

// EOF
