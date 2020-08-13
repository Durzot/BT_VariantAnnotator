# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Main functions for running each step and the assembling step.
"""
import os
import numpy  as np
import pandas as pd

from   .manual  import run_manual_annotator
from   .vcf2maf import run_vcf2maf_annotator
from   .vep     import run_vep_annotator

DataFrame = pd.core.frame.DataFrame

#### # FUNCTION FOR ONE VCF
#### # #######################################################################################################

def run_annotator_one(vcf_folder: str, vcf_file: str, col_normal: str, col_tumor: str, tumor_id: str, normal_id: str,
infos_n_reads: list, infos_other: list, vcf2maf: str, vep_folder: str, vep_data: str, fasta: str, dt_folders: dict):
    """
    Run the manual, vcf2maf and vep annotations on one VCF file and assemble.

    Parameters
    --------
    vcf_file: str
        name of the vcf file
    vcf_folder: str
        path to the folder where the vcf is
    col_normal: str
        name of the column in the vcf for the normal sample
    col_tumor: str
        name of the column in the vcf for the tumor sample
    infos_n_reads: list
        list of sigles that contain read info
    infos_other: list
        list of sigles that need extraction
    vcf2maf: str
        path to the vcf2maf perl script
    vep_folder: str
        path to the folder where the vep command is
    vep_data: str
        path to the .vep data where the reference genome is located
    fasta: str
        relative path to fasta file from vep_folder
    vcf_folder: str
        path to the folder where the vcf files are
    dt_folders: dict
        dict with the following keys:
            * manual_out_folder
            * vcf2maf_tmp_folder
            * vcf2maf_out_folder
            * vep_out_folder
            * maf_folder
    """
    vcf_path = os.path.join(vcf_folder, vcf_file)

    out_file         = vcf_file.replace(".vcf", ".txt")
    manual_out_path  = os.path.join(dt_folders["manual_out_folder"], out_file)
    vcf2maf_out_path = os.path.join(dt_folders["vcf2maf_out_folder"], out_file)
    vep_out_path     = os.path.join(dt_folders["vep_out_folder"], out_file)

    #### # 1. RUN EACH ANNOTATOR
    #### # ###################################################################################################

    run_manual_annotator(
        vcf_path      = vcf_path,
        out_path      = manual_out_path,
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )

    run_vcf2maf_annotator(
        vcf2maf    = vcf2maf,
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = vcf_path,
        out_path   = vcf2maf_out_path,
        tmp_folder = dt_folders["vcf2maf_tmp_folder"],
        tumor_id   = tumor_id,
        normal_id  = normal_id,
        fasta      = fasta
    )

    run_vep_annotator(
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = vcf_path,
        out_path   = vep_out_path,
        fasta      = fasta
    )

    #### # 2. ASSEMBLE ANNOTATIONS
    #### ######################################################################################################

    ddf_vcf = {}

    #### vep manual
    ddf_vcf["manual"] = pd.read_csv(
        filepath_or_buffer = manual_out_path,
        sep                = "\t"
    )

    #### vep alone
    skipsymbol = "##"
    with open(vep_out_path) as file:
        skiprows = sum(line.startswith(skipsymbol) for line in file.readlines())

    df_alone = pd.read_csv(
        filepath_or_buffer = vep_out_path,
        sep                = "\t",
        skiprows           = skiprows
    )

    df_alone_extra = df_alone.Extra.apply(lambda x: dict(item.split("=") for item in x.split(";")))
    df_alone_extra = df_alone_extra.apply(pd.Series)
    df_alone = pd.concat((df_alone, df_alone_extra), axis=1)

    ddf_vcf["alone"] = df_alone

    #### maf vcf2maf
    skipsymbol = "#"
    with open(vcf2maf_out_path) as file:
        skiprows = sum(line.startswith(skipsymbol) for line in file.readlines())

    ddf_vcf["vcf2maf"] = pd.read_csv(
        filepath_or_buffer = vcf2maf_out_path,
        sep                = "\t",
        skiprows           = skiprows
    )

    #### checked that all scripts ran ok the vcf file
    n_manual  = ddf_vcf["manual"].shape[0]
    n_vcf2maf = ddf_vcf["vcf2maf"].shape[0]
    n_alone   = ddf_vcf["alone"].shape[0]

    if n_manual != n_vcf2maf:
        print("SHAPE error: manual %d vs vcf2maf %d " % (n_manual, n_vcf2maf), flush=True)

    if n_manual != n_alone:
        print("SHAPE error: manual %d vs alone %d " % (n_manual, n_alone), flush=True)

    #### choose column per source
    choice_columns = {
        "Hugo_Symbol"                 : "vcf2maf",
        "Entrez_Gene_Id"              : "vcf2maf",
        "NCBI_Build"                  : "vcf2maf",
        "Chromosome"                  : "vcf2maf",
        "Start_Position"              : "vcf2maf",
        "End_Position"                : "vcf2maf",
        "Variant_Quality"             : "manual",
        "Filter"                      : "manual",
        "Variant_Classification"      : "vcf2maf",
        "Variant_Type"                : "vcf2maf",
        "Reference_Allele"            : "vcf2maf",
        "Tumor_Seq_Allele1"           : "vcf2maf",
        "Tumor_Seq_Allele2"           : "vcf2maf",
        "dbSNP_RS"                    : "vcf2maf",
        "Tumor_Sample"                : "manual",
        "Tumor_Sample_Barcode"        : "manual",
        "Tumor_Sample_Site"           : "manual",
        "Matched_Norm_Sample_Barcode" : "manual",
        "Match_Norm_Seq_Allele1"      : "vcf2maf",
        "Match_Norm_Seq_Allele2"      : "vcf2maf",
        "Tumor_Sample_UUID"           : "manual",
        "Matched_Norm_Sample_UUID"    : "manual",
        "HGVSc"                       : "vcf2maf",
        "HGVSp"                       : "vcf2maf",
        "HGVSp_Short"                 : "vcf2maf",
        "all_effects"                 : "vcf2maf",
        "Location"                    : "alone",
        "Gene"                        : "alone",
        "Feature"                     : "alone",
        "Feature_type"                : "alone",
        "cDNA_position"               : "alone",
        "CDS_position"                : "alone",
        "Protein_position"            : "alone",
        "Amino_acids"                 : "alone",
        "Codons"                      : "alone",
        "Existing_variation"          : "alone",
        "Consequence"                 : "alone",
        "IMPACT"                      : "alone",
        "SIFT"                        : "alone",
        "PolyPhen"                    : "alone",
        "CLIN_SIG"                    : "alone",
        "STRAND"                      : "alone",
        "SYMBOL_SOURCE"               : "alone",
        "HGNC_ID"                     : "alone",
        "BIOTYPE"                     : "alone",
        "CCDS"                        : "alone",
        "ENSP"                        : "alone",
        "SWISSPROT"                   : "alone",
        "TREMBL"                      : "alone",
        "UNIPARC"                     : "alone",
        "EXON"                        : "alone",
        "INTRON"                      : "alone",
        "AF"                          : "alone",  #### AF from 1000 Genomes Phase 3
        "AA_AF"                       : "alone",  #### NHLBI-ESP populations
        "EA_AF"                       : "alone",  #### NHLBI-ESP populations
        "gnomAD_AF"                   : "alone",  #### gnomAD exome populations
        "MAX_AF"                      : "alone",  #### highest AF from any of 1000 genomes, ESP or gnomAD
        "MAX_AF_POPS"                 : "alone",  #### highest AF from any of 1000 genomes, ESP or gnomAD
    }

    maf_columns = []
    cols_reads = [x for x in ddf_vcf["manual"].columns if x.startswith(("n_", "t_"))]

    #### add column from dict
    for column,source in choice_columns.items():
        if column not in ddf_vcf[source].columns:
            print("WARNING: %s is not in %s" % (column, source), flush=True)
        else:
            maf_columns.append(ddf_vcf[source][column])

    #### add read columns from manual extraction of reads
    for column in cols_reads:
        if column not in ddf_vcf["manual"].columns:
            print("WARNING: %s is not in %s" % (column, source), flush=True)
        else:
            maf_columns.append(ddf_vcf["manual"][column])

    #### # 3. SAVE
    #### ######################################################################################################

    maf_out_path = os.path.join(dt_folders["maf_folder"], out_file)
    df_maf = pd.concat(maf_columns, axis=1)
    df_maf.to_csv(
        path_or_buf = maf_out_path,
        sep         = "\t",
        header      = True,
        index       = False
    )
    print("maf file saved at %s" % maf_out_path, flush=True)


def run_annotator(i_split: int, n_split: int, vcf2maf: str, vep_folder: str, vep_data: str, fasta: str, vcf_folder: str,
out_folder: str, vcf_list_path: str=None, vcf_meta_path: str=None):
    """
    Run the manual, vcf2maf and vep annotations and assemble.


    Parameters
    --------
    i_split: int
        the split processed.
    n_split: int
        the number of splits for processing a set of vcf files
    vcf2maf: str
        path to the vcf2maf perl script
    vep_folder: str
        path to the folder where the vep command is
    vep_data: str
        path to the .vep data where the reference genome is located
    fasta: str
        relative path to fasta file from vep_folder
    vcf_folder: str
        path to the folder where the vcf files are
    out_folder: str
        path to the folder where subfolders with output results will be saved
    vcf_list_path: str, optional
        path to the file containing the list of vcf files to be processed. If not specified, all vcf files found in the
        vcf_folder are processed.
    vcf_meta_path: str, optional
        path to the file that contain the names of the tumor and normal columns for each vcf file. The file must be a
        .txt file containing 3 columns: "vcf_name", "normal_column", "tumor_column". If not specified,
        all vcf files must have columns with the names NORMAL and PRIMARY.
    """

    #### paths to results folders
    dt_folders = {
        'manual_out_folder'  : os.path.join(out_folder, "tmp/manual/out"),
        'vcf2maf_tmp_folder' : os.path.join(out_folder, "tmp/vcf2maf/tmp"),
        'vcf2maf_out_folder' : os.path.join(out_folder, "tmp/vcf2maf/out"),
        'vep_out_folder'     : os.path.join(out_folder, "tmp/vep/out"),
        'maf_folder'         : os.path.join(out_folder, "tmp/maf"),
    }

    #### make folders if they do not exist already
    for k, v in dt_folders.items():
        os.makedirs(v, exist_ok=True)

    #### load meta data
    df_meta = pd.read_csv(
        filepath_or_buffer = dt_paths["vcf_meta"],
        sep                = "\t"
    )

    if os.path.exists(dt_paths["vcf_list"]):
        with open(dt_paths["vcf_list"]) as file:
            vcf_files = file.read().splitlines()
    else:
        #### what if compressed ?
        vcf_files = [x for x in os.listdir(vcf_folder) if x.endswith(".vcf")]

    #### get list of vcfs for the split
    count_one_split = len(vcf_files)//args.n_split

    if args.i_split == args.n_split:
        vcf_files  = vcf_files[(args.i_split-1)*count_one_split:]
    else:
        vcf_files  = vcf_files[(args.i_split-1)*count_one_split:args.i_split*count_one_split]

    count = 0
    count_total = len(vcf_files)

    #### loop over the list
    for vcf_file in vcf_files:

        col_normal    = "NORMAL"
        col_tumor     = "PRIMARY"
        normal_id     = "TCGA-A1-A0SD-10A-01D-A110-09"
        tumor_id      = "TCGA-A1-A0SD-01A-11D-A10Y-09"
        infos_n_reads = ["AD", "DP4", "DP", "TAR", "TIR"]
        infos_other   = ["SS", "GT"]

        run_annotator_one(
            vcf_folder    = vcf_folder,
            vcf_file      = vcf_file,
            col_normal    = col_normal,
            col_tumor     = col_tumor,
            normal_id     = normal_id,
            tumor_id      = tumor_id,
            infos_n_reads = infos_n_reads,
            infos_other   = infos_other,
            vcf2maf       = vcf2maf,
            vep_folder    = vep_folder,
            vep_data      = vep_data,
            fasta         = fasta,
            dt_folders    = dt_folders
        )
