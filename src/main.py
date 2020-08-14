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

def run_annotator(vcf_folder: str, vcf_file: str, col_normal: str, col_tumor: str, tumor_id: str, normal_id: str,
infos_n_reads: list, infos_other: list, vcf2maf: str, vep_folder: str, vep_data: str, fasta: str, dt_folders: dict,
dt_identifiers: dict=None):
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
    dt_identifiers: dict, optional
        dict with key, value pairs that will be added as single-value columns in the maf file
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

    if dt_identifiers is not None:
        for k,v in dt_identifiers.items():
            ddf_vcf["manual"].insert(0, k, v)

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

    #### add identifier columns if any
    if dt_identifiers is not None:
        for column in dt_identifiers.keys():
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
