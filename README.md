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
- *dnSNP_RS*
- *HGVS_p*
- *HGVS_c*
- *HGVS_Short*
- *all_effects*

that are not available from VEP's output.

### 1.3. VEP

Run the VEP annotator on the VCF file from a specific set of options. The options can be changed in the code if required.

## 2. How to run the tool ?
