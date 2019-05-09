#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

#include "params.h"
#include <vector>
#include <string>
#include <htslib/hts.h>
#include <htslib/vcf.h>

using namespace std;

typedef struct variant_desc_tag {
	string chrom;
	int64_t pos;
	string id;
	string ref;
	string alt;
	string qual;
	string filter;
	string info;

	bool operator==(const struct variant_desc_tag &x)
	{
		if (chrom == "" && x.chrom == "")
			return true;

		return chrom == x.chrom &&
			pos == x.pos;
//		&&
//			id == x.id;
//		&&
//			ref == x.ref &&
//			alt == x.alt &&
//			qual == x.qual &&
//			filter == x.filter;
//			info == x.info;
	}

	bool operator !=(const struct variant_desc_tag &x)
	{
		return !operator==(x);
	}

	bool operator<(const struct variant_desc_tag &x)
	{
		if (chrom != x.chrom)
		{
			if (chrom.empty())
				return false;
			if (x.chrom.empty())
				return true;

			return chrom < x.chrom;
		}
		if (pos != x.pos)
			return pos < x.pos;
		return alt < x.alt;
	}
} variant_desc_t;

class CVCF
{
    htsFile * vcf_file;
    bcf_hdr_t * vcf_hdr;
    bcf1_t * rec;
    int ploidy;
    
    bool first_variant;
    int curr_alt_number; //allela from ALT field (from 1)
    int32_t *tmpia; //to set genotypes in new variant    

public:
	CVCF();
	~CVCF();

	// Open VCF file for reading
	bool OpenForReading(string & file_name);

	// Open VCF file for writing
	bool OpenForWriting(string & file_name, file_type type, char bcf_compression_level);

	// Close VCF file
	bool Close();

	// If open, return no. of samples
	int GetNoSamples();

	// Jeśli plik jest otwarty do odczytu, to zwraca wiersze z metadanymi - każdy wiersz w osobnym stringu
	//bool GetMeta(vector<string> &v_meta);

	// Jeśli plik jest otwarty do zapisu, to zapisuje wiersze metadanych
	// Ta funkcja może być wykonana tylko zaraz po otwarciu plik a przed zapisem jakichkolwiek wariantów
	//bool SetMeta(vector<string> &v_meta);

	// Jeśli plik jest otwarty do odczytu, to zwraca wiersz nagłówka - każde pole w jednym stringu
	//bool GetHeader(vector<string> &v_header);

	// Jeśli plik jest otwarty do zapisu, to zapisuje wiersz nagłówka
	// Ta funkcja może być wykonana tylko zaraz po otwarciu plik a przed zapisem jakichkolwiek wariantów
	//bool SetHeader(vector<string> &v_header);

	// Jeśli plik jest otwarty do odczytu, to zwraca informację czy doszliśmy do końca
//	bool Eof();

	// Get complete header as a string
	bool GetHeader(string &v_header);
	
	// Parse complete header given as a string
	bool SetHeader(string & v_header);
	
	// If open for writing store the header
	bool WriteHeader();
	
	// Get ploidy
	int GetPloidy();
    
	// Set ploidy
	void SetPloidy(int _ploidy);
	
	// If file open give the next variant:
	// desc - variant description
	// data - genotypes encoded (1B for each sample) as:
	// * bits 0-1 information about 1st haplotype:
	//     00 - absent
	//     01 - present
	//     10 - multi-allele
	//     11 - unknown
	// * bits 2-3 - information about 2nd haplotype, encoded as above
	// * bit 4:
	//     0 - no phasing info
	//     1 - data phase
	bool GetVariant(variant_desc_t &desc, vector<uint8_t> &data);

	// Store info about variant - parameters the same as for GetVariant
	bool SetVariant(variant_desc_t &desc, vector<uint8_t> &data);
	
	// Get vector with sample names
	bool GetSamplesList(vector<string> &s_list);
	
	// Add sample names to file
	bool AddSamples(vector<string> &s_list);
    
	// Add a single sample name
	bool AddSample(string &s_name);
};

// EOF
