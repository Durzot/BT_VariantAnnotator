# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Example of how to annotate a list of VCF from one project/study.

Example
-----------
python examples/run_example_tcga_HS.py \
    --i_split 1 \
    --n_split 1 \
    --vcf2maf ~/Documents/biotools/informatics/VCF/vcf2maf/vcf2maf.pl \
    --vep_folder ~/Documents/biotools/informatics/VCF/ensembl-vep \
    --vep_data ~/.vep \
    --vep_n_fork 4 \
    --fasta ~/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
"""
import argparse
import os
import pandas as pd
import sys

if "." not in sys.path:
    sys.path.append(".")

from   src.main import run_annotator

#### # SCRIPT PARAMETERS 
#### #####################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--i_split'         , type=int , default=1  , help='the split processed')
parser.add_argument('--n_split'         , type=int , default=1  , help='total number of splits')
parser.add_argument('--vcf2maf'         , type=str , default="" , help='path to the vcf2maf perl script')
parser.add_argument('--vep_folder'      , type=str , default="" , help='path to the folder of the vep command')
parser.add_argument('--vep_data'        , type=str , default="" , help='path to the .vep data folder')
parser.add_argument('--vep_n_fork'      , type=int , default=4  , help='number of forks to be used by VEP')
parser.add_argument('--fasta'           , type=str , default="" , help='path to reference genome FASTA file')
args = parser.parse_args()

print("Parameters", flush=True)
for arg in vars(args):
    print("%s: %s" % (arg, getattr(args, arg)), flush=True)

#### # SCRIPT FUNCTION
#### #####################################################################################################

if __name__ == "__main__":

    vcf_folder = "./examples/data/TCGA_HS/"
    out_folder = "./examples/results/TCGA_HS/"
    vcf_meta_path = os.path.join(vcf_folder, "vcf_meta.txt")

    #### paths to results folders
    dt_folders = {
        'manual_out_folder'  : os.path.join(out_folder, "tmp/out_manual"),
        'vcf2maf_tmp_folder' : os.path.join(out_folder, "tmp/tmp_vcf2maf"),
        'vcf2maf_out_folder' : os.path.join(out_folder, "tmp/out_vcf2maf"),
        'vep_out_folder'     : os.path.join(out_folder, "tmp/out_vep"),
        'maf_folder'         : os.path.join(out_folder, "maf"),
    }

    #### # 1. LOAD
    #### # ##################################################################################################

    for k, v in dt_folders.items():
        if "folder" in k:
            os.makedirs(v, exist_ok=True)

    #### load meta data
    df_meta = pd.read_csv(
        filepath_or_buffer = vcf_meta_path,
        sep                = "\t"
    )

    vcf_files = [x for x in os.listdir(vcf_folder) if x.endswith(".vcf")]

    #### # 2. SPLIT
    #### # ##################################################################################################

    count_one_split = len(vcf_files)//args.n_split

    if args.i_split == args.n_split:
        vcf_files  = vcf_files[(args.i_split-1)*count_one_split:]
    else:
        vcf_files  = vcf_files[(args.i_split-1)*count_one_split:args.i_split*count_one_split]

    count = 0
    count_total = len(vcf_files)

    #### # 3. PROCESS
    #### # ##################################################################################################

    #### loop over the list
    for vcf_file in vcf_files:
        count += 1
        print("="*80, flush=True)
        print("vcf %d/%d" % (count, count_total), flush=True)
        print("processing %s\n" % vcf_file, flush=True)

        #### get vcf identifiers
        mask_vcf_file  = df_meta["file_name_HiSeq"] == vcf_file
        index_vcf_file = mask_vcf_file[mask_vcf_file].index[0]

        dt_identifiers = {
            "Tumor_Sample"                : df_meta.loc[index_vcf_file, "tumor_sample"],
            "Tumor_Sample_Barcode"        : df_meta.loc[index_vcf_file, "tumor_sample_barcode"],
            "Matched_Norm_Sample_Barcode" : df_meta.loc[index_vcf_file, "normal_sample_barcode"],
            "Tumor_Sample_Site"           : df_meta.loc[index_vcf_file, "tumor_sample_barcode"].split("-")[3][:2],
        }

        #### get parameter values
        col_normal = dt_identifiers["Matched_Norm_Sample_Barcode"]
        col_tumor  = dt_identifiers["Tumor_Sample_Barcode"]
        normal_id  = dt_identifiers["Matched_Norm_Sample_Barcode"]
        tumor_id   = dt_identifiers["Tumor_Sample_Barcode"]

        vcf_type = df_meta.loc[index_vcf_file, "vcf_type"]

        if vcf_type == "indel":
            infos_n_reads = ["AD", "DP4", "DP", "TAR", "TIR"]
        else:
            infos_n_reads = ["AD", "DP4", "DP", "FA"]

        infos_other   = ["SS", "GT"]

        run_annotator(
            vcf_folder        = vcf_folder,
            vcf_file          = vcf_file,
            col_normal        = col_normal,
            col_tumor         = col_tumor,
            normal_id         = normal_id,
            tumor_id          = tumor_id,
            infos_n_reads     = infos_n_reads,
            infos_other       = infos_other,
            vcf2maf           = args.vcf2maf,
            vep_folder        = args.vep_folder,
            vep_data          = args.vep_data,
            # vep_custom      = "~/.vep/custom/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN",
            vep_n_fork        = args.vep_n_fork,
            vep_overwrite     = True,
            vcf2maf_overwrite = True,
            fasta             = args.fasta,
            dt_folders        = dt_folders,
            dt_identifiers    = dt_identifiers
        )
