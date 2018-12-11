// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include "application.h"
#include "utils.h"
#include "lzma_wrapper.h"

#include <iostream>

using namespace std;

// ******************************************************************************
CApplication::CApplication(const CParams &_params)
{
	params = _params;

	// To avoid unnecessary reallocation of buffers
	v_vcf_data_compress.reserve(no_variants_in_buf);
	v_vcf_data_io.reserve(no_variants_in_buf);
}

// ******************************************************************************
CApplication::~CApplication()
{
}

// ******************************************************************************
bool CApplication::CompressDB()
{
	CBarrier barrier(3);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	bool end_of_processing = false;

	if (!vcf->OpenForReading(params.vcf_file_name))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	if (!cfile->OpenForWriting(params.db_file_name))
		return false;

	bool ploidy_initialised = false;

	cfile->SetNeglectLimit(params.neglect_limit);
	cfile->SetNoSamples(vcf->GetNoSamples());

	string header;
	vector<string> v_samples;

	vcf->GetHeader(header);
	vcf->GetSamplesList(v_samples);
	cfile->SetHeader(header);
	cfile->AddSamples(v_samples);

	// Thread reading VCF files in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			v_vcf_data_io.clear();

			for (size_t i = 0; i < no_variants_in_buf; ++i)
			{
  				v_vcf_data_io.push_back(make_pair(variant_desc_t(), vector<uint8_t>()));
				if (!vcf->GetVariant(v_vcf_data_io.back().first, v_vcf_data_io.back().second))
				{
					v_vcf_data_io.pop_back();
					break;
				}
			}

			if (!ploidy_initialised)
			{
				cfile->SetPloidy(vcf->GetPloidy());
				cfile->InitPBWT();
				ploidy_initialised = true;
			}

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread making PBWT and compressing data
	unique_ptr<thread> t_vcf(new thread([&] {
		while (!end_of_processing)
		{
			for (size_t i = 0; i < v_vcf_data_compress.size(); ++i)
				cfile->SetVariant(v_vcf_data_compress[i].first, v_vcf_data_compress[i].second);

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	size_t no_variants = 0;
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_compress.empty())
			end_of_processing = true;

		no_variants += v_vcf_data_compress.size();
		cout << no_variants << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_vcf->join();
	t_io->join();

	cfile->Close();
	vcf->Close();
	cout << endl;

	return true;
}

// ******************************************************************************
bool CApplication::DecompressDB()
{
	CBarrier barrier(3);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	bool end_of_processing = false;

	if (!vcf->OpenForWriting(params.vcf_file_name, params.out_type, params.bcf_compression_level))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	if (!cfile->OpenForReading(params.db_file_name))
		return false;

	cfile->InitPBWT();
	params.neglect_limit = cfile->GetNeglectLimit();

	uint32_t no_variants = cfile->GetNoVariants();
	uint32_t i_variant = 0;

	string header;
	vector<string> v_samples;

	cfile->GetHeader(header);
	cfile->GetSamples(v_samples);
	vcf->SetHeader(header);
	vcf->AddSamples(v_samples);
	vcf->WriteHeader();

	vcf->SetPloidy(cfile->GetPloidy());

	// Thread making rev-PBWT and decompressing data
	unique_ptr<thread> t_vcf(new thread([&] {
		while (!end_of_processing)
		{
			v_vcf_data_compress.clear();

			for (size_t i = 0; i < no_variants_in_buf && i_variant < no_variants; ++i, ++i_variant)
			{
				v_vcf_data_compress.push_back(make_pair(variant_desc_t(), vector<uint8_t>()));
				cfile->GetVariant(v_vcf_data_compress.back().first, v_vcf_data_compress.back().second);
			}
			
			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread writing VCF files in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			for (size_t i = 0; i < v_vcf_data_io.size(); ++i)
				vcf->SetVariant(v_vcf_data_io[i].first, v_vcf_data_io[i].second);
			v_vcf_data_io.clear();

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_io.empty())
			end_of_processing = true;

		cout << i_variant << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_vcf->join();
	t_io->join();

	cfile->Close();
	vcf->Close();
	cout << endl;

	return true;
}

// ******************************************************************************
bool CApplication::ExtractSample()
{
	CBarrier barrier(3);
	unique_ptr<CVCF> vcf(new CVCF());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	bool end_of_processing = false;
	vector<pair<uint8_t, uint32_t>> rle_genotypes;

	if (!vcf->OpenForWriting(params.vcf_file_name, params.out_type, params.bcf_compression_level))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	if (!cfile->OpenForReading(params.db_file_name))
		return false;

	cfile->InitPBWT();
	params.neglect_limit = cfile->GetNeglectLimit();

	uint32_t no_variants = cfile->GetNoVariants();
	uint32_t i_variant = 0;

	string header;
	vector<string> v_samples;

	cfile->GetHeader(header);
	cfile->GetSamples(v_samples);
	vcf->SetHeader(header);
	vcf->AddSample(params.id_sample);
	vcf->WriteHeader();

	auto p = find(v_samples.begin(), v_samples.end(), params.id_sample);
	if (p == v_samples.end())
	{
		cerr << "Sample: " << params.id_sample << " does not exist\n";
		return false;
	}

	uint32_t ploidy = cfile->GetPloidy();
	array<uint32_t, 2> sample_pos;
	array<uint32_t, 2> sample_pos_perm;
	array<uint8_t, 2> val{ 0,0 };

	if (ploidy == 1)
		sample_pos[0] = (uint32_t) (p - v_samples.begin());
	else if (ploidy == 2)
	{
		sample_pos[0] = 2u * (uint32_t) (p - v_samples.begin());
		sample_pos[1] = sample_pos[0] + 1;
	}
	sample_pos_perm = sample_pos;

	vcf->SetPloidy(ploidy);

	// Thread making rev-PBWT and decompressing data
	unique_ptr<thread> t_vcf(new thread([&] {
		while (!end_of_processing)
		{
			v_vcf_data_compress.clear();

			for (size_t i = 0; i < no_variants_in_buf && i_variant < no_variants; ++i, ++i_variant)
			{
				v_vcf_data_compress.push_back(make_pair(variant_desc_t(), vector<uint8_t>()));
				cfile->GetVariantGenotypesRawAndDesc(v_vcf_data_compress.back().first, rle_genotypes);

				uint8_t variant_data = 0u;

				if (ploidy == 1)
					cfile->TrackItem(rle_genotypes, sample_pos_perm[0], val[0], sample_pos_perm[0]);
				else if (ploidy == 2)
				{
					cfile->TrackItems(rle_genotypes, sample_pos_perm, val, sample_pos_perm);
					variant_data = 0b00010000u;		// Data phased
				}

//				variant_data += (uint8_t) (val[0] + (val[1] << 2));
				variant_data += val[0];
				variant_data += (uint8_t) (val[1] << 2);

				v_vcf_data_compress.back().second.push_back(variant_data);
			}

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread writing VCF files in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			for (size_t i = 0; i < v_vcf_data_io.size(); ++i)
				vcf->SetVariant(v_vcf_data_io[i].first, v_vcf_data_io[i].second);
			v_vcf_data_io.clear();

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		swap(v_vcf_data_compress, v_vcf_data_io);
		if (v_vcf_data_io.empty())
			end_of_processing = true;

		cout << i_variant << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_vcf->join();
	t_io->join();

	cfile->Close();
	vcf->Close();
	cout << endl;
	
	return true;
}

// ******************************************************************************
bool CApplication::find_prev_value(const vector<run_desc_t> &v_rle_genotypes, const uint32_t max_pos, const uint8_t value, uint32_t &found_pos)
{
	bool found = false;
	found_pos = 0;
	uint32_t cur_pos = 0;

	for (auto run : v_rle_genotypes)
	{
		if (run.first == value)
		{
			if (cur_pos + run.second >= max_pos)
			{
				found_pos = max_pos - 1;
				cur_pos = max_pos;
				found = true;
			}
			else
			{
				found_pos = cur_pos + run.second - 1;
				cur_pos += run.second;
				found = true;
			}
		}
		else
			cur_pos += run.second;

		if (cur_pos >= max_pos)
			break;
	}

	return found;
}

// ******************************************************************************
bool CApplication::find_next_value(const vector<run_desc_t> &v_rle_genotypes, const uint32_t min_pos, const uint8_t value, uint32_t &found_pos)
{
	bool found = false;
	found_pos = 0;
	uint32_t cur_pos = 0;

	for (auto run : v_rle_genotypes)
	{
		if (run.first == value)
		{
			if (cur_pos + run.second >= min_pos)
			{
				found_pos = min_pos;
				found = true;
				cur_pos = min_pos;
			}
			else if (cur_pos >= min_pos)
			{
				found_pos = min_pos;
				found = true;
			}
			else
				cur_pos += run.second;
		}
		else
			cur_pos += run.second;

		if (cur_pos >= min_pos && found)
			break;
	}

	return found;
}

// ******************************************************************************
bool CApplication::CompressSample()
{
	CBarrier barrier(3);
	unique_ptr<CSampleFile> sfile(new CSampleFile());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	unique_ptr<CVCF> vfile(new CVCF());
	bool end_of_processing = false;
	vector<run_desc_t> rle_genotypes;

	if (!sfile->OpenForWriting(params.sample_file_name))
	{
		cerr << "Cannot open: " << params.sample_file_name << endl;
		return false;
	}

	if (!vfile->OpenForReading(params.vcf_file_name))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	if (!cfile->OpenForReading(params.db_file_name))
		return false;

	cfile->InitPBWT();
	params.neglect_limit = cfile->GetNeglectLimit();

	uint32_t no_variants = cfile->GetNoVariants();
	uint32_t i_variant = 0;

	string header, v_header;
	vector<string> v_samples, cur_sample;

	cfile->GetHeader(header);
	cfile->GetSamples(v_samples);
	vfile->GetHeader(v_header);
	vfile->GetSamplesList(cur_sample);
	
	if (cur_sample.size() != 1)
	{
		cerr << "File to compress must contain exactly 1 sample\n";
		return false;
	}

	string empty_header;
	sfile->WriteHeaderAndSample(header, params.store_sample_header ? v_header : empty_header, cur_sample.front());

	uint32_t ploidy = cfile->GetPloidy();
	array<uint32_t, 2> sample_pos;
	array<uint32_t, 2> sample_pos_perm;
	uint32_t no_pred_same[2] = { 0, 0 };
	uint32_t no_succ_same[2] = { 0, 0 };

	if (ploidy == 1)
		sample_pos[0] = (uint32_t) v_samples.size();
	else if (ploidy == 2)
	{
		sample_pos[0] = 2u * (uint32_t) v_samples.size();
		sample_pos[1] = sample_pos[0];
	}
	sample_pos_perm = sample_pos;

	// Thread making rev-PBWT and estimating the position of sample to process
	unique_ptr<thread> t_vcf(new thread([&] {
		variant_desc_t desc;
		vector<uint8_t> v_data;

		while (!end_of_processing)
		{
			v_vcf_data_compress.clear();

			for (size_t i = 0; i < no_variants_in_buf && i_variant < no_variants; ++i, ++i_variant)
			{
				cfile->GetVariantGenotypesRaw(rle_genotypes);

				v_data.clear();
				vfile->GetVariant(desc, v_data);

				array<uint8_t, 2> a_sample;

				for (uint32_t j = 0; j < ploidy; ++j)
				{
					run_t runs;
					uint8_t value = (v_data[0] >> (2 * j)) & 0b00000011;
					a_sample[j] = value;

					cfile->EstimateValue(rle_genotypes, sample_pos_perm[j], value, runs, sample_pos_perm[j]);

					v_sample_data_compress.push_back(make_tuple(value, runs, no_pred_same[j], no_succ_same[j]));

					if (runs[0].first == value)
						no_pred_same[j] = min<size_t>(no_pred_same[j]+1, max_tracked_dist);
					else
					{
						uint32_t pos_sample_to_trace;

						if (find_prev_value(rle_genotypes, sample_pos_perm[j], value, pos_sample_to_trace))
						{
							no_pred_same[j] = 1;
							for (uint32_t k = 0; k < d_hist_rle_genotypes.size(); ++k)
							{
								if (cfile->RevertDecode(pos_sample_to_trace, d_hist_rle_genotypes[k], d_hist_tracked_samples[k][j]))
									++no_pred_same[j];
								else
									break;
							}
						}
						else
							no_pred_same[j] = 0;
					}

					if (runs[1].first == value)
						no_succ_same[j] = min<size_t>(no_succ_same[j]+1, max_tracked_dist);
					else
					{
						uint32_t pos_sample_to_trace;

						if (find_next_value(rle_genotypes, sample_pos_perm[j], value, pos_sample_to_trace))
						{
							no_succ_same[j] = 1;
							for (uint32_t k = 0; k < d_hist_rle_genotypes.size(); ++k)
							{
								if (cfile->RevertDecode(pos_sample_to_trace, d_hist_rle_genotypes[k], d_hist_tracked_samples[k][j]))
									++no_succ_same[j];
								else
									break;
							}
						}
						else
							no_succ_same[j] = 0;
					}
				}

				// Docelowo ten bufor cykliczny mo¿na zrobiæ np. na array - powinno byæ szybciej, ale poki co sprawdzam czy sam pomysl jest dobry
				d_hist_rle_genotypes.push_front(rle_genotypes);
				d_hist_tracked_samples.push_front(a_sample);

				if (d_hist_tracked_samples.size() > max_tracked_dist)
				{
					d_hist_rle_genotypes.pop_back();
					d_hist_tracked_samples.pop_back();
				}
			}

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread writing decompressed sample VCF file in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			for (size_t i = 0; i < v_sample_data_io.size(); ++i)
				sfile->Put(get<0>(v_sample_data_io[i]), get<1>(v_sample_data_io[i]), get<2>(v_sample_data_io[i]), get<3>(v_sample_data_io[i]));
			v_sample_data_io.clear();

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		swap(v_sample_data_compress, v_sample_data_io);
		if (v_sample_data_io.empty())
			end_of_processing = true;

		cout << i_variant << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	t_vcf->join();
	t_io->join();

	cout << endl;

	return true;
}

// ******************************************************************************
bool CApplication::DecompressSample()
{
	CBarrier barrier(3);
	unique_ptr<CSampleFile> sfile(new CSampleFile());
	unique_ptr<CCompressedFile> cfile(new CCompressedFile());
	unique_ptr<CVCF> vfile(new CVCF());
	bool end_of_processing = false;
	vector<pair<uint8_t, uint32_t>> rle_genotypes;

	if (!sfile->OpenForReading(params.sample_file_name))
	{
		cerr << "Cannot open: " << params.sample_file_name << endl;
		return false;
	}

	if (!vfile->OpenForWriting(params.vcf_file_name, params.out_type, params.bcf_compression_level))
	{
		cerr << "Cannot open: " << params.vcf_file_name << endl;
		return false;
	}

	if (!cfile->OpenForReading(params.db_file_name))
		return false;

	cfile->InitPBWT();
	params.neglect_limit = cfile->GetNeglectLimit();

	uint32_t no_variants = cfile->GetNoVariants();
	uint32_t i_variant = 0;

	string header;
	vector<string> v_samples;
	string sample_name;
	string v_header;

	cfile->GetHeader(header);
	cfile->GetSamples(v_samples);

	sfile->ReadHeaderAndSample(header, v_header, sample_name);

	vfile->SetHeader(v_header);
	params.id_sample = sample_name;
	vfile->AddSample(params.id_sample);
	vfile->WriteHeader();

	uint32_t ploidy = cfile->GetPloidy();
	array<uint32_t, 2> sample_pos;
	array<uint32_t, 2> sample_pos_perm;
	uint32_t no_pred_same[2] = { 0, 0 };
	uint32_t no_succ_same[2] = { 0, 0 };

	if (ploidy == 1)
		sample_pos[0] = (uint32_t)v_samples.size();
	else if (ploidy == 2)
	{
		sample_pos[0] = 2u * (uint32_t)v_samples.size();
		sample_pos[1] = sample_pos[0];
	}
	sample_pos_perm = sample_pos;

	vfile->SetPloidy(ploidy);

	// Thread making rev-PBWT and estimating position of the sample to process
	unique_ptr<thread> t_vcf(new thread([&] {
		variant_desc_t desc;
		array<uint8_t, 2> a_sample;

		while (!end_of_processing)
		{
			v_sample_d_data_compress.clear();

			for (size_t i = 0; i < no_variants_in_buf && i_variant < no_variants; ++i, ++i_variant)
			{
				cfile->GetVariantGenotypesRawAndDesc(desc, rle_genotypes);

				uint8_t value = ploidy == 2 ? 0b00010000 : 0;		// Data phased

				for (uint32_t j = 0; j < ploidy; ++j)
				{
					run_t runs;
					uint32_t tmp;
					uint8_t v;

					cfile->EstimateValue(rle_genotypes, sample_pos_perm[j], 0, runs, tmp);
					sfile->Get(v, runs, no_pred_same[j], no_succ_same[j]);
					cfile->EstimateValue(rle_genotypes, sample_pos_perm[j], v, runs, sample_pos_perm[j]);

/*					if (runs[0].first == v)
						no_pred_same[j]++;
					else
						no_pred_same[j] = 0;

					if (runs[1].first == v)
						no_succ_same[j]++;
					else
						no_succ_same[j] = 0;*/

					if (runs[0].first == v)
						no_pred_same[j] = min<size_t>(no_pred_same[j] + 1, max_tracked_dist);
					else
					{
						uint32_t pos_sample_to_trace;

						if (find_prev_value(rle_genotypes, sample_pos_perm[j], v, pos_sample_to_trace))
						{
							no_pred_same[j] = 1;
							for (uint32_t k = 0; k < d_hist_rle_genotypes.size(); ++k)
							{
								if (cfile->RevertDecode(pos_sample_to_trace, d_hist_rle_genotypes[k], d_hist_tracked_samples[k][j]))
									++no_pred_same[j];
								else
									break;
							}
						}
						else
							no_pred_same[j] = 0;
					}

					if (runs[1].first == v)
						no_succ_same[j] = min<size_t>(no_succ_same[j] + 1, max_tracked_dist);
					else
					{
						uint32_t pos_sample_to_trace;

						if (find_next_value(rle_genotypes, sample_pos_perm[j], v, pos_sample_to_trace))
						{
							no_succ_same[j] = 1;
							for (uint32_t k = 0; k < d_hist_rle_genotypes.size(); ++k)
							{
								if (cfile->RevertDecode(pos_sample_to_trace, d_hist_rle_genotypes[k], d_hist_tracked_samples[k][j]))
									++no_succ_same[j];
								else
									break;
							}
						}
						else
							no_succ_same[j] = 0;
					}

					value += v << (2 * j);
					a_sample[j] = v;
				}
				v_sample_d_data_compress.push_back(make_pair(desc, value));

				// Docelowo ten bufor cykliczny mo¿na zrobiæ np. na array - powinno byæ szybciej, ale poki co sprawdzam czy sam pomysl jest dobry
				d_hist_rle_genotypes.push_front(rle_genotypes);
				d_hist_tracked_samples.push_front(a_sample);

				if (d_hist_tracked_samples.size() > max_tracked_dist)
				{
					d_hist_rle_genotypes.pop_back();
					d_hist_tracked_samples.pop_back();
				}
			}

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Thread writing compressed sample file file in parts
	unique_ptr<thread> t_io(new thread([&] {
		while (!end_of_processing)
		{
			vector<uint8_t> variants(1);

			for (size_t i = 0; i < v_sample_d_data_io.size(); ++i)
			{
				variants[0] = v_sample_d_data_io[i].second;
				vfile->SetVariant(v_sample_d_data_io[i].first, variants);
			}

			v_sample_d_data_io.clear();

			barrier.count_down_and_wait();
			barrier.count_down_and_wait();
		}
	}));

	// Synchronization
	while (!end_of_processing)
	{
		barrier.count_down_and_wait();
		swap(v_sample_d_data_compress, v_sample_d_data_io);
		if (v_sample_d_data_io.empty())
			end_of_processing = true;

		cout << i_variant << "\r";
		fflush(stdout);
		barrier.count_down_and_wait();
	}

	cout << i_variant << "\r";
	fflush(stdout);

	t_vcf->join();
	t_io->join();

	vfile->Close();
	cout << endl;

	return true;
}

// EOF
