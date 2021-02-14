# Biotool for annotating variants from a VCF file.

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![codecov](https://codecov.io/gh/ypradat/VariantAnnotator/branch/master/graph/badge.svg?token=H821S3WZHS)](https://codecov.io/gh/ypradat/VariantAnnotator)

The tool is divided in 3 steps
- Manual parsing of the VCF
- Run [vcf2maf](https://github.com/mskcc/vcf2maf) to extract standard information 
- Run [Variant Ensembl's Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) for annotations. 

## 1. What is the tool doing ?

VEP annotates variants with information from multiple external databases and can be configured to answer a lot of
specific needs. For more details, see [VEP's options
page](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html). VEP does not however extract information
like number of reads or somatic status from the VCF file. vcf2maf is supposed to perform these tasks but failed to do on
quite a lot of example VCF and does not always extract all the relevant information. For this reason, manual parsing was
implemented in Python and was tested on TCGA VCF files from the legacy archive portal (see examples) and other VCFs. See
the table for exhaustive details.

### 1.1. Manual parsing

Relies on tags specified by the user to extract relevant info like genotype (GT), somatic status (SS), quality and
filter info (QUAL, INFO) and most importantly reads information (AD, DP, FA, DP4, TAR, TIR). The parser has been tested
on VCF files as produced by
- Mutect v.1 (TCGA GA SNV) and Strelka (TCGA GA Indel)
- sets of callers VarScanSomatic-Strelka-Sniper-Samtools (TCGA HS SNP), GatkSomaticIndel-Pindel-Strelka-VarScanSomatic
  (TCGA HS Indel).
- Mutect v.1.1.7 with no header
- Strelka v.2.9.2

### 1.2. vcf2maf

vcf2maf also runs VEP internally but performs extra work to build some annotation fields including the following
- *Hugo_Symbol* (can be replaced by VEP's *Symbol*)
- *Entrez_Gene_Id*
- *NCBI_Build*
- *Chromosome* (can be extracted from VEP)
- *Start_Position*, *End_Position*
- *Variant_Classification* (rewriting of VEP's *Consequence*)
- *Variant_Type*
- *Reference_Allele*
- *Tumor_Seq_Allele1* (can be discarded)
- *Tumor_Seq_Allele2* (can be discarded)
- *dbSNP_RS*
- *HGVS_p*
- *HGVS_c*
- *HGVS_Short*
- *all_effects*
that are not available from VEP's output.

### 1.3. VEP

Run the VEP annotator on the VCF file from a specific set of options. The options can be changed in the code if required.

## 2. How to run the tool ?

Clone the github repository using

```
git clone --recurse-submodules [path/to/this/repo]
```

or, if you forgot the `--recurse-submodules` option, run `git submodule update --init` after the cloning.

### 2.1 Install VEP

Go to the ensembl-vep/ folder in tools/ and run

```
perl INSTALL.pl
```

The perl script INSTALL.pl may return errors as missing dependencies or other. For instance,the error
`Bio::Root::Version is not installed` may be solved by running `sudo cpanm Bio::Root::Version`. You may have more than
one such library to install. Refer to the github for the details.

The installation from the perl script offers the choice to install cache files (most efficient use of vep) and FASTA
files (to retrieve sequence data for HGVS notations) into `$HOME/.vep`. You may also install plugins for additional
analyses. Download cache files for Homo Sapiens genome 100_GRCh37 (or newer). The total download size is about 12 GB of
data so a stable and fast connection is required here. The uncompressed FASTA file (**DO NOT FORGET** to uncompress this
file or VEP will fail) Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz requires about 3.0 GB of storage.

As detailed in this [page](https://m.ensembl.org/info/docs/tools/vep/script/vep_custom.html), sometimes the Ensembl
latest release does not use the latests release of other databases. If you want the latest release of say, ClinVar, you
should download all relevant files and specify the path to these files along with field names while calling the vep
command. For instance, you can get the latest ClinVar VCF files from their FTP website with

```
# Compressed VCF file
curl -O ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
# Index file
curl -O ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
````

and then add the following to the vep command

```
./vep [...] --custom /path/to/custom_files/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN
```

### 2.2 Install vcf2maf

VEP is required by vcf2maf but you also need the commands from `samtools` and `htslib` available at
[http://www.htslib.org/download/](http://www.htslib.org/download/). Do the following

```
cd samtools-1.x
./configure --prefix=/where/to/install
make
make install
ln -s /where/to/install/bin/samtools /usr/local/bin

cd htslib −1.x
./configure −−prefix=/where/to/install
make
make install

ln −s /where/to/install/bin/htslib /usr/local/bin 
ln −s /where/to/install/bin/tabix /usr/local/bin
ln −s /where/to/install/bin/bgzip /usr/local/bin
```

You may replace `/usr/local/bin` with whatever path where you usually save binaries. Finish the installation of vcf2maf
following the instructions given in the [github](https://github.com/mskcc/vcf2maf).

### 2.3 Example

The main function for annotating a vcf is `run_annotator` in `main` module. Have a look at `run_example_tcga_GA.py` or
`run_example_tcga_HS` to have examples of how to run the tool and at the `run_annotator` documentation for more details
about the options.

## 3. References

McLaren, W., Gil, L., Hunt, S.E. et al. The Ensembl Variant Effect Predictor. Genome Biol 17, 122 (2016).
[https://doi.org/10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4).
