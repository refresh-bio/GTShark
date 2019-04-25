# GTShark: Genotype compression in large projects

Requirements
--------------

GTShark requires:

* A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.

Installation
--------------

To download, build and install GTShark use the following commands.
```sh
git clone https://github.com/refresh-bio/GTShark.git
cd GTShark
./install.sh 
make
make install
```
The `install.sh` script downloads and installs the HTSlib library into the `include` and `lib` directories in the `GTShark/htslib` directory. 

By default GTShark is installed in the `bin` directory of the `/usr/local/` directory. A different location prefix can be specified with `prefix` parameter:
```sh
make prefix=/usr/local install
```
---
To uninstall GTShark:
```sh
make uninstall
```
This uninstalls GTShark from the `/usr/local` directory. To uninstall from different location use the `prefix` parameter:
```sh
make prefix=/usr/local uninstall
```
To uninstall the HTSlib library use the provided uninstall script:
```sh
./uninstall.sh 
```
---
To clean the GTShark build use:
```sh
make clean
```
Usage
--------------
* Compress the input VCF/BCF file.
```
Input: <input_vcf> VCF/BCF file. 
Output: <output_db>_gt file with the archive and <output_db>_db file with description of samples and variants.

Usage: gtshark compress-db [options] <input_vcf> <output_db>
Parameters:
  input_vcf - path to input VCF (or VCF.GZ or BCF) file
  output_db - path to output database file
Options:
  -nl <value> - ignore rare variants; value is a limit of alternative alleles (default: 10)
  ```
  
 * Decompress the whole archive.
 ```
Input: <input_db> archive (<input_db>_gt and <input_db>_db). 
Output: <output_vcf> VCF/BCF file.
 
Usage: gtshark decompress-db [options] <input_db> <output_vcf>
Parameters:
  input_db   - path to input database file
  output_vcf - path to output VCF/BCF file
Options:
  -b - output BCF file (VCF file by default)
  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)	
 ```
 
 * Extract a sample from a database (compressed VCF/BCF file).
 ```
Input: <database> archive (<database>_gt and <database>_db) and <sample_id> id of the sample to extract.
Output: <output_sample> VCF/BCF file.

Usage: gtshark extract-sample [options] <database> <sample_id> <output_sample>
Parameters:
  database      - path to database file obtained using `compress-db' command
  sample id     - id of sample to decompress
  output_sample - path to output VCF file containing a single sample
Options:
  -b - output BCF file (VCF file by default)
  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)
 ```
 
* Compress a sample in reference to the existing database (compressed VCF/BCF file).
 ```
Input: <database> archive (<database>_gt and <database>_db) and <input_sample> VCF/BCF file of the sample to compress.
Output: <compressed_sample> sample archive.

Usage: gtshark compress-sample [options] <database> <input_sample> <compressed_sample>
Parameters:
  database          - path to database file obtained using `compress-db' command
  input_sample      - path to input VCF (or VCF.GZ or BCF) file containing a single sample
  compressed_sample - path to output compressed file containing a single sample
Options:
  -sh               - store header of compressed_sample file
 ```


* Decompress a sample in reference to the existing database (compressed VCF/BCF file).
 ```
Input: <database> archive (<database>_gt and <database>_db) and <compressed_sample> sample archive .
Output: <output_sample> VCF/BCF file.

Usage: gtshark decompress-sample [options] <database> <compressed_sample> <output_sample>
Parameters:
  database          - path to database file obtained using `compress-db' command
  compressed_sample - path to compressed file containing a single sample
  output_sample     - path to output VCF file containing a single sample
Options:
  -b - output BCF file (VCF file by default)
  -c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)	
 ```
 
 
Toy example
--------------

There are two example VCF files, `toy_collection.vcf` and `toy_new_sample.vcf`, in the `toy_ex` folder, which can be used to test GTShark. All instructions should be called within `toy_ex` folder.

To compress the example VCF file `toy_collection.vcf` with a collection of samples (6 samples, 11 variant sites including one with two alternating alleles) and store the archive called `toy_archive` in the `toy_ex` folder:
```sh
../gtshark compress-db toy_collection.vcf toy_archive
```
This will create an archive consisting of two files
* `toy_archive_db` - file with description of samples and variants,
* `toy_archive_gt` - main archive with all genotypes compressed.

To decompress the compressed archive `toy_archive` to a VCF file `toy_coll_decomp.vcf`:
```sh
../gtshark  decompress-db toy_archive toy_coll_decomp.vcf
```

To extract a single sample (`s2`) from the compressed archive `toy_archive` to a VCF file `extracted_sample.vcf`:
```sh
../gtshark extract-sample toy_archive s2 extracted_sample.vcf
```

To compress new sample (`toy_new_sample.vcf`) in reference to the existing compressed database/archive `toy_archive` and store it in archive `toy_new_sample_comp`:
```sh
../gtshark compress-sample toy_archive toy_new_sample.vcf toy_new_sample_comp
```

To decompress a sample compressed in reference to database `toy_archive` and stored in archive `toy_new_sample_comp` to a VCF file `toy_new_decomp.vcf`:
```sh
../gtshark decompress-sample toy_archive toy_new_sample_comp toy_new_sample_decomp.vcf
```

For more options see Usage section.

Developers
--------------
The GTShark algorithm was invented by [Sebastian Deorowicz](https://github.com/sebastiandeorowicz) and [Agnieszka Danek](https://github.com/agnieszkadanek).
The implementation is by Sebastian Deorowicz and Agnieszka Danek.

Publication
--------------

<!--- Danek, A., Deorowicz, S., GTC: an attempt to maintenance of huge genome collections compressed, [bioRxiv](http://biorxiv.org/content/early/2017/04/28/131649), 2017; --->

