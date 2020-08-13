# Biotool for annotating variants from a VCF file.

The tool is divided in 3 steps
- Manual parsing of the VCF
- Run [vcf2maf](https://github.com/mskcc/vcf2maf) to extract standard information 
- Run [Variant Ensembl's Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) for annotations. 

## 1. What is the tool doing ?

### 1.1. Manual parsing

VEP annotates variants with information from multiple external databases and can be configure for to answer a lot of specific needs. For more details,
see [VEP's options page](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html). VEP does not however extract information like number of reads or somatic status from the VCF file. vcf2maf is supposed to perform these tasks but failed to do on quite a lot of example VCF and does not always extract all the relevant information. For this reason, manual parsing was implemented in Python and was tested on TCGA VCF files from the legacy archive portal (see examples) and other VCFs. See the table for exhaustive details.

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

### 2.1 Install VEP

Clone VEP from the [official github](https://github.com/Ensembl/ensembl-vep). Installation details are provided there. As specified, run

```
git clone https://github.com/Ensembl/ensembl-vep
cd ensembl-vep
perl INSTALL.pl
```

The perl script INSTALL.pl may return errors as missing dependencies or other. For instance,the error `Bio::Root::Version is not installed` may be solved by running `sudo cpanm Bio::Root::Version`. You may have more than one such library to install. Refer to the github for the details.

The installation from the perl script offers the choice to install cache files (most efficient use of vep) and FASTA files (to retrieve sequence data for HGVS notations) into `$HOME/.vep`. You may also install plugins for additional analyses. Download cache files for Homo Sapiens genome 99_GRCh37 (or newer). The total download size is about 12 GB of data so a stable and fast connection is required here. The FASTA file Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz requires 882 MB of storage. *DO NOT FORGET* to unzip this file (or VEP will fail).

### 2.2 Install vcf2maf

VEP is required by vcf2maf may you need also need the set of commands `samtools`, and `htslib` available at [http://www.htslib.org/download/](http://www.htslib.org/download/). Do the following

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

You may replace `/usr/local/bin` with whatever path where you usually save binaries. Finish the installation of vcf2maf following the instructions given in the [github](https://github.com/mskcc/vcf2maf).
