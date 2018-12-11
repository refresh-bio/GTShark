// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include "vcf.h"
#include <iostream>

// ************************************************************************************
CVCF::CVCF()
{
    vcf_file = nullptr;
    vcf_hdr = nullptr;
    rec = nullptr;
    curr_alt_number = 1;
    tmpia = nullptr;
    ploidy = 0; //default
    first_variant = true;
}

// ************************************************************************************
CVCF::~CVCF()
{
}

// ************************************************************************************
bool CVCF::OpenForReading(string & file_name)
{
    if(vcf_file)
        hts_close(vcf_file);
    vcf_file = hts_open(file_name.c_str(), "r"); //  With 'r' opens for reading; any further format mode letters are ignored as the format is detected by checking the first few bytes or BGZF blocks of the file.
    if(!vcf_file)
        return false;
    hts_set_opt(vcf_file, HTS_OPT_CACHE_SIZE, 32000000);
    if(vcf_hdr)
        bcf_hdr_destroy(vcf_hdr);
    vcf_hdr = bcf_hdr_read(vcf_file);
    rec = bcf_init();
    return true;
}

// ************************************************************************************
bool CVCF::OpenForWriting(string & file_name, file_type type, char bcf_compression_level)
{
    if(type == file_type::VCF)
        vcf_file = hts_open(file_name.c_str(), "w");
    else  // file_type::BCF
    {
        char write_mode[5] = "wb";
        write_mode[2] = bcf_compression_level;
        write_mode[3] = '\0';
        vcf_file = hts_open(file_name.c_str(), write_mode);
    }
    if(!vcf_file)
        return false;
    hts_set_opt(vcf_file, HTS_OPT_CACHE_SIZE, 32000000);
    rec = bcf_init();
    return true;
}

// ************************************************************************************
bool CVCF::Close()
{
    if(vcf_file)
    {
        if(hts_close(vcf_file) < 0)
            return false;
        vcf_file = nullptr;
    }
    if(vcf_hdr)
    {
        bcf_hdr_destroy(vcf_hdr);
        vcf_hdr = nullptr;
    }
    if(rec)
    {
        bcf_clear(rec);
        rec = nullptr;
    }
    return true;
    
}

// ************************************************************************************
int CVCF::GetNoSamples()
{
    if(!vcf_file || !vcf_hdr)
        return -1;
    
    return bcf_hdr_nsamples(vcf_hdr);
}

// ************************************************************************************
bool CVCF::GetSamplesList(vector<string> &s_list)
{
    if(!vcf_hdr)
        return false;
    int n = GetNoSamples();
    for (int i = 0; i < n; i++)
        s_list.push_back(vcf_hdr->samples[i]);
    return true;
}

// ************************************************************************************
int CVCF::GetPloidy()
{
    return ploidy;
}

// ************************************************************************************
void CVCF::SetPloidy(int _ploidy)
{
    if(_ploidy != 1 && _ploidy != 2)
    {
        std::cerr << "Unsupported ploidy (" << _ploidy << ")\n";
        exit(1);
    }
    ploidy = _ploidy;
}
/*
bool CVCF::GetMeta(vector<string> &v_meta)
{
    return true;
}

bool CVCF::SetMeta(vector<string> &v_meta)
{
    return true;
}

bool CVCF::GetHeader(vector<string> &v_header)
{
    return true;
}

bool CVCF::SetHeader(vector<string> &v_header)
{
    return true;
}

bool CVCF::Eof()
{
    return true;
}
*/

// ************************************************************************************
bool CVCF::GetHeader(string &v_header)
{
    if(!vcf_hdr)
        return false;
    kstring_t str = {0, 0, nullptr};
    bcf_hdr_format(vcf_hdr, 0, &str);
    char * ptr = strstr(str.s, "#CHROM");
    v_header.assign(str.s, ptr - str.s);
    return true;
}

// ************************************************************************************
bool CVCF::SetHeader(string &v_header)
{
    vcf_hdr = bcf_hdr_init("r");
    string temp = v_header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    bcf_hdr_parse(vcf_hdr, (char *)temp.c_str());
    bcf_hdr_sync(vcf_hdr);
    if(!vcf_hdr)
        return false;
    return true;
}

// ************************************************************************************
bool CVCF::AddSamples(vector<string> &s_list)
{
    if(!vcf_hdr)
        return false;
    for(size_t i = 0; i < s_list.size(); i++)
        bcf_hdr_add_sample(vcf_hdr, s_list[i].c_str());
    
    bcf_hdr_sync(vcf_hdr);
    return true;
}

// ************************************************************************************
bool CVCF::AddSample(string & s_name)
{
    if(!vcf_hdr)
        return false;
    bcf_hdr_add_sample(vcf_hdr, s_name.c_str());
    bcf_hdr_sync(vcf_hdr);
    return true;
}

// ************************************************************************************
bool CVCF::WriteHeader()
{
    if(vcf_hdr && vcf_file)
    {
        if(bcf_hdr_write(vcf_file, vcf_hdr) < 0)
            return false;
        return true;
    }
    return false;
}

// ************************************************************************************
bool CVCF::GetVariant(variant_desc_t &desc, vector<uint8_t> &data)
{
    if(!vcf_file || !vcf_hdr)
        return false;
    if(curr_alt_number == 1) // read new line/record
    {
        bcf_clear(rec);
        if(bcf_read(vcf_file, vcf_hdr, rec) == -1)
            return false;
        
        if(rec->errcode)
        {
            std::cerr << "Error in VCF file\n";
            exit(1);
        }
        bcf_unpack((bcf1_t*)rec, BCF_UN_ALL);
    }
    
    desc.chrom = vcf_hdr->id[BCF_DT_CTG][rec->rid].key; // CHROM
    desc.pos = rec->pos + 1;  // POS
    desc.id = rec->d.id ? rec->d.id : "."; // ID
    
    desc.ref.erase(); //REF
    if (rec->n_allele > 0)
        desc.ref =  rec->d.allele[0];
    else
        desc.ref =  '.';
    
    desc.alt.erase();  // ALT
    if (rec->n_allele > 1) {
        desc.alt += rec->d.allele[curr_alt_number];
    }
    else
        desc.alt = '.';
    if(rec->n_allele > 2) //more than one ALT allele
    {
        desc.alt += ",<M>";
    }
        
    desc.qual.erase(); // QUAL
    if ( bcf_float_is_missing(rec->qual) )
        desc.qual = '.';
    else
    {
        kstring_t s = {0,0,0};
        kputd(rec->qual, &s);
        desc.qual += s.s;
        
    }    
    
    desc.filter.erase();  //FILTER
    if (rec->d.n_flt) {
        for (int i = 0; i < rec->d.n_flt; ++i) {
            if (i) desc.filter += ';';
            desc.filter += vcf_hdr->id[BCF_DT_ID][rec->d.flt[i]].key;
        }
    }
    else
        desc.filter = '.';
    
    desc.info.erase(); // INFO
    if (rec->n_info) {
        int first = 1;
        for (int i = 0; i < rec->n_info; ++i) {
            bcf_info_t *z = &rec->d.info[i];
            if ( !z->vptr )
                continue;
            if ( !first )
                desc.info += ';';
            first = 0;
            if (z->key >= vcf_hdr->n[BCF_DT_ID]) {
                hts_log_error("Invalid BCF, the INFO index is too large");
                // errno = EINVAL;
                return -1;
            }
            desc.info += vcf_hdr->id[BCF_DT_ID][z->key].key;
            if (z->len <= 0) continue;
                desc.info += '=';
            if (z->len == 1)
            {
                switch (z->type)
                {
                    case BCF_BT_INT8:  if ( z->v1.i==bcf_int8_missing ) desc.info += '.'; else desc.info += to_string(z->v1.i); break;
                    case BCF_BT_INT16: if ( z->v1.i==bcf_int16_missing ) desc.info += '.'; else desc.info += to_string(z->v1.i); break;
                    case BCF_BT_INT32: if ( z->v1.i==bcf_int32_missing ) desc.info += '.'; else desc.info += to_string(z->v1.i); break;
                    case BCF_BT_FLOAT: if ( bcf_float_is_missing(z->v1.f) ) desc.info += '.'; else {
                        kstring_t s = {0,0,0};
                        kputd(z->v1.f, &s);
                        desc.info += s.s;}
                        break;
                    case BCF_BT_CHAR:  desc.info += z->v1.i; break;
                    default: hts_log_error("Unexpected type %d", z->type); exit(1); break;
                }
            }
            else
            {
                kstring_t s_info = {0,0, nullptr};
                bcf_fmt_array(&s_info, z->len, z->type, z->vptr);
                desc.info += s_info.s;
            }
        }
        if ( first ) desc.info += '.';
    } else desc.info += '.';

       
    //genotypes
    int *gt_arr = NULL, ngt_arr = 0, ngt;
    ngt = bcf_get_genotypes(vcf_hdr, rec, &gt_arr, &ngt_arr);
    if (ngt <= 0 )
        return false; //genotype not present
    if(first_variant)
    {
        ploidy = ngt/bcf_hdr_nsamples(vcf_hdr);
        if(ploidy != 1 && ploidy != 2)
        {
            std::cerr << "Unsupported ploidy (" << ploidy << ")\n";
            exit(1);
        }
        first_variant = false;
    }
    else
    {
        if((float)ngt/bcf_hdr_nsamples(vcf_hdr) != ploidy)
        {
            std::cerr << "Unsupported ploidy (different in different variants/samples)\n";
            exit(1);
        }
    }
    uint8_t genotype;
    int allele;
    if(ploidy == 2) //ploidy = 2
    {
        for(int i = 0; i < ngt_arr; i+=2)
        {
           
            allele = bcf_gt_allele(gt_arr[i]);
            if(bcf_gt_is_missing(gt_arr[i]))
            {
                genotype = 3;
            }
            else if(allele == 0)
            {
                genotype = 0;
            }
            else if(allele == curr_alt_number)
            {
                genotype = 1;
            }
            else
            {
                genotype = 2;
            }
            
            allele = bcf_gt_allele(gt_arr[i+1]);
            if(bcf_gt_is_missing(gt_arr[i+1]))
            {
                genotype += (uint8_t) (3u << 2);
            }
            else if(allele == 0)
            {
                genotype += (uint8_t) (0u << 2);
            }
            else if(allele == curr_alt_number)
            {
                genotype += (uint8_t) (1u << 2);
            }
            else
            {
                genotype += (uint8_t) (2u << 2);
            }
            
            genotype += (uint8_t) (bcf_gt_is_phased(gt_arr[i+1])  << 4);

            data.push_back(genotype);
        }
        
    }
    else //if (ngt == bcf_hdr_nsamples(vcf_hdr)) //haploid
    {
        for(int i = 0; i < ngt_arr; i++)
        {
            allele = bcf_gt_allele(gt_arr[i]);
            if(bcf_gt_is_missing(gt_arr[i]))
            {
                genotype = 3;
            }
            else if(allele == 0)
            {
                genotype = 0;
            }
            else if(allele == curr_alt_number)
            {
                genotype = 1;
            }
            else
            {
                genotype = 2;
            }
            data.push_back(genotype);
        }
    }

    free(gt_arr);
    
    if(rec->n_allele > 2)
    {
        //if "ALT,<M>", do not add additional line(as VCF was already altered)
        if(rec->n_allele==3 && strcmp(rec->d.allele[2], "<M>") == 0)
            ;
        else
        {
            curr_alt_number++;
            if(curr_alt_number == rec->n_allele)
                curr_alt_number = 1;
        }
    }

    return true;
}

// ************************************************************************************
bool CVCF::SetVariant(variant_desc_t &desc, vector<uint8_t> &data)
{
    bcf_clear(rec);
    
    string record;
   // record = desc.chrom + "\t" + to_string(desc.pos) + "\t" + desc.id + "\t" + desc.ref + "\t" + desc.alt + "\t" + desc.qual +"\t" + desc.filter + "\t" + desc.info;
    record = desc.chrom + "\t0\t" + desc.id + "\t" + desc.ref + "\t" + desc.alt + "\t" + desc.qual +"\t" + desc.filter + "\t" + desc.info;
    kstring_t s;
    s.s = (char*)record.c_str();
    s.m = record.length();
    s.l = 0;
    vcf_parse(&s, vcf_hdr, rec);
    rec->pos = (int32_t) (desc.pos - 1);
  
    // GT
    if(first_variant)
    {
        tmpia = new int[(bcf_hdr_nsamples(vcf_hdr)*sizeof(int))*2];
        first_variant = false;
    }
    if(ploidy == 2) //diploid
    {
        for(size_t i = 0; i < data.size(); i++)
        {
            switch(data[i] & 0x3) //first allele
            {
                case 3:
                    tmpia[i*2] = bcf_gt_missing;
                    break;
                default:
                    tmpia[i*2] = bcf_gt_unphased(data[i] & 0x3);
                    break;
            }
            //phasing info only present in second allele
            if(data[i] & 0x10)
            {
                switch((data[i] & 0xC) >> 2) //second allele
                {
                    case 3:
                       tmpia[i*2+1] = bcf_gt_missing | 1;  //first bit indicates phasing, lack of appropriate function/macro in htslib to do it nicely
                        break;
                    default:
                        tmpia[i*2+1] = bcf_gt_phased((data[i] >> 2) & 0x3);
                        break;
                }
            }
            else
            {
                switch((data[i] & 0xC) >> 2) //second allele
                {
                    case 3:
                        tmpia[i*2+1] = bcf_gt_missing;
                        break;
                    default:
                        tmpia[i*2+1] = bcf_gt_unphased((data[i] >> 2) & 0x3);
                        break;
                }
            }
        }
        bcf_update_genotypes(vcf_hdr, rec, tmpia, bcf_hdr_nsamples(vcf_hdr)*2);
    }
    else  //haploid
    {
        for(size_t i = 0; i < data.size(); i++)
        {
            if(data[i] == 3)
            {
                tmpia[i] = bcf_gt_missing;
            }
            else
            {
                tmpia[i] = bcf_gt_unphased(data[i]);
            }
        }
        bcf_update_genotypes(vcf_hdr, rec, tmpia, bcf_hdr_nsamples(vcf_hdr));
        
    }
   
    bcf_write(vcf_file, vcf_hdr, rec);
    
    return true;

}

// EOF
