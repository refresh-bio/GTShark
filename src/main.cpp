// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <nmmintrin.h>
#include <chrono>

#include "params.h"
#include "application.h"

#include "vios.h"
#include "rc.h"
#include "sub_rc.h"
#include "io.h"
#include "utils.h"

using namespace std;
using namespace std::chrono;

CParams params;
CApplication *app;

int old_main(int argc, char **argv);

bool parse_params(int argc, char **argv);
void usage_main();
void usage_compress_db();
void usage_decompress_db();
void usage_compress_sample();
void usage_decompress_sample();
void usage_extract_sample();

// ******************************************************************************
void usage_main()
{
	cerr << "gtshark <mode>\n";
	cerr << "Parameters:\n";
	cerr << "  mode - one of:\n";
	cerr << "    compress-db       - compress VCF file with collection of samples\n";
	cerr << "    decompress-db     - decompress VCF file with collection of samples\n";
	cerr << "    compress-sample   - compress VCF file containing a single sample\n";
	cerr << "    decompress-sample - decompress VCF file containing a single sample\n";
	cerr << "    extract-sample    - extract a single sample from database\n";
}

// ******************************************************************************
void usage_compress_db()
{
	cerr << "gtshark compress-db [options] <input_vcf> <output_db>\n";
	cerr << "Parameters:\n";
	cerr << "  input_vcf - path to input VCF (or VCF.GZ or BCF) file\n";
	cerr << "  output_db - path to output database file\n";
	cerr << "Options:\n";
    cerr << "  -nl <value> - ignore rare variants; value is a limit of alternative alleles (default: " << params.neglect_limit << ")\n";
}

// ******************************************************************************
void usage_decompress_db()
{
	cerr << "gtshark decompress-db [options] <input_db> <output_vcf>\n";
	cerr << "Parameters:\n";
	cerr << "  input_db   - path to input database file\n";
	cerr << "  output_vcf - path to output VCF file\n";
    cerr << "Options:\n";
    cerr << "  -b - output BCF file (VCF file by default)\n";
    cerr << "  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)\t"<< endl;
}

// ******************************************************************************
void usage_compress_sample()
{
	cerr << "gtshark compress-sample [options] <database> <input_sample> <compressed_sample>\n";
	cerr << "Parameters:\n";
	cerr << "  database          - path to database file obtained using `compress-db' command\n";
	cerr << "  input_sample      - path to input VCF (or VCF.GZ or BCF) file containing a single sample\n";
	cerr << "  compressed_sample - path to output compressed file containing a single sample\n";
	cerr << "Options:\n";
	cerr << "  -sh               - store header of compressed_sample file\n";
	cerr << "  -ev               - allow differnt variant sets in sample file and database\n";
}

// ******************************************************************************
void usage_decompress_sample()
{
	cerr << "gtshark decompress-sample [options] <database> <compressed_sample> <output_sample>\n";
	cerr << "Parameters:\n";
	cerr << "  database          - path to database file obtained using `compress-db' command\n";
	cerr << "  compressed_sample - path to compressed file containing a single sample\n";
	cerr << "  output_sample     - path to output VCF file containing a single sample\n";
    cerr << "Options:\n";
    cerr << "  -b - output BCF file (VCF file by default)\n";
    cerr << "  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)\t"<< endl;
}

// ******************************************************************************
void usage_extract_sample()
{
	cerr << "gtshark extract-sample [options] <database> <sample_id> <output_sample>\n";
	cerr << "Parameters:\n";
	cerr << "  database      - path to database file obtained using `compress-db' command\n";
	cerr << "  sample id     - id of sample to decompress\n";
	cerr << "  output_sample - path to output VCF file containing a single sample\n";
    cerr << "Options:\n";
    cerr << "  -b - output BCF file (VCF file by default)\n";
    cerr << "  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)\t"<< endl;
}

// ******************************************************************************
bool parse_params(int argc, char **argv)
{
	if (argc < 2)
	{
		usage_main();

		return false;
	}

	if (string(argv[1]) == "compress-db")
		params.work_mode = work_mode_t::compress_db;
	else if (string(argv[1]) == "decompress-db")
		params.work_mode = work_mode_t::decompress_db;
	else if (string(argv[1]) == "compress-sample")
		params.work_mode = work_mode_t::compress_sample;
	else if (string(argv[1]) == "decompress-sample")
		params.work_mode = work_mode_t::decompress_sample;
	else if (string(argv[1]) == "extract-sample")
		params.work_mode = work_mode_t::extract_sample;

	// Compress-db
	if (params.work_mode == work_mode_t::compress_db)
	{
		if (argc < 4)
		{
			usage_compress_db();
			return false;
		}

		int i = 2;
		while (i < argc - 2)
		{
			if (string(argv[i]) == "-nl" && i + 1 < argc - 2)
			{
				params.neglect_limit = atoi(argv[i + 1]);
				i += 2;
			}
        }

		params.vcf_file_name = string(argv[i]);
		params.db_file_name = string(argv[i+1]);
	}
	else if (params.work_mode == work_mode_t::decompress_db)
	{
		if (argc < 4)
		{
			usage_decompress_db();
			return false;
		}

        int i = 2;
        while (i < argc - 2)
        {
            if (string(argv[i]) == "-b")
            {
                params.out_type = file_type::BCF;
                i ++;
            }
            else if (string(argv[i]) == "-c")
            {
                i++;
                if(i >= argc - 2)
                {
                    usage_compress_db();
                    return false;
                }
                int tmp = atoi(argv[i]);
                if(tmp < 0 || tmp > 9)
                {
                    usage_compress_db();
                    return false;
                }
                else
                {
                    if(tmp)
                        params.bcf_compression_level = argv[i][0];
                    else
                        params.bcf_compression_level = 'u';
                }
                i++;
            }
            else
            {
                cerr << "Unknown option : " << argv[i] << endl;
                usage_compress_db();
                return false;
            }
            

        }
		params.db_file_name = string(argv[i]);
		params.vcf_file_name = string(argv[i+1]);
	}
	else if (params.work_mode == work_mode_t::compress_sample)
	{
		if (argc < 5)
		{
			usage_compress_sample();
			return false;
		}

		int i = 2;
		while (i < argc - 3)
		{
			if (string(argv[i]) == "-sh")
			{
				params.store_sample_header = true;
				++i;
			}
			else if (string(argv[i]) == "-ev")
			{
				params.extra_variants = true;
				++i;
			}
			else
            {
                cerr << "Unknown option : " << argv[i] << endl;
                usage_compress_sample();
                return false;
            }
        }

		params.db_file_name = string(argv[i]);
		params.vcf_file_name = string(argv[i+1]);
		params.sample_file_name = string(argv[i+2]);
	}
	else if (params.work_mode == work_mode_t::decompress_sample)
	{
		if (argc < 5)
		{
			usage_decompress_sample();
			return false;
		}

        int i = 2;
        while (i < argc - 3)
        {
            if (string(argv[i]) == "-b")
            {
                params.out_type = file_type::BCF;
                i ++;
            }
            else if (string(argv[i]) == "-c")
            {
                i++;
                if(i >= argc - 2)
                {
                    usage_compress_db();
                    return false;
                }
                int tmp = atoi(argv[i]);
                if(tmp < 0 || tmp > 9)
                {
                    usage_compress_db();
                    return false;
                }
                else
                {
                    if(tmp)
                        params.bcf_compression_level = argv[i][0];
                    else
                        params.bcf_compression_level = 'u';
                }
                i++;
            }
            else
            {
                cerr << "Unknown option : " << argv[i] << endl;
                usage_decompress_sample();
                return false;
            }

        }
		params.db_file_name = string(argv[i]);
		params.sample_file_name = string(argv[i+1]);
		params.vcf_file_name = string(argv[i+2]);
	}
	else if (params.work_mode == work_mode_t::extract_sample)
	{
		if (argc < 5)
		{
			usage_extract_sample();
			return false;
		}
        int i = 2;
        while (i < argc - 3)
        {
            if (string(argv[i]) == "-b")
            {
                params.out_type = file_type::BCF;
                i ++;
            }
            else if (string(argv[i]) == "-c")
            {
                i++;
                if(i >= argc - 2)
                {
                    usage_compress_db();
                    return false;
                }
                int tmp = atoi(argv[i]);
                if(tmp < 0 || tmp > 9)
                {
                    usage_compress_db();
                    return false;
                }
                else
                {
                    if(tmp)
                        params.bcf_compression_level = argv[i][0];
                    else
                        params.bcf_compression_level = 'u';
                }
                i++;
            }
            else
            {
                cerr << "Unknown option : " << argv[i] << endl;
                usage_extract_sample();
                return false;
            }
        }

		params.db_file_name = string(argv[i]);
		params.id_sample = string(argv[i+1]);
		params.vcf_file_name = string(argv[i+2]);
	}
	else
	{
		cerr << "Unknown mode : " << argv[2] << endl;
		usage_main();

		return false;
	}

	return true;
}

// ******************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	app = new CApplication(params);

	bool result = true;

	if (params.work_mode == work_mode_t::compress_db)
		result = app->CompressDB();
	else if (params.work_mode == work_mode_t::decompress_db)
		result = app->DecompressDB();
	else if (params.work_mode == work_mode_t::extract_sample)
		result = app->ExtractSample();
	else if (params.work_mode == work_mode_t::compress_sample)
		result = app->CompressSample();
	else if (params.work_mode == work_mode_t::decompress_sample)
		result = app->DecompressSample();

	delete app;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	if (!result)
		std::cout << "Critical error!\n";

	std::cout << "Processing time: " << time_span.count() << " seconds.\n";

	fflush(stdout);

	return 0;
}

// EOF
