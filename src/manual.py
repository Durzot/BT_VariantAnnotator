# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Functions for manually parsing VCF files.
"""

import os
import numpy  as np
import pandas as pd
import re

from .util import load_vcf

DataFrame = pd.core.frame.DataFrame

#### # SUBFUNCTIONS
#### #####################################################################################################

def get_matches(regex: str, string: str) -> list:
    """
    Get the  matches from applying a regex on a string.
    """
    matches = re.finditer(regex, string, re.MULTILINE)
    lt_matches = []

    for match in matches:
        lt_matches += [match.group()]

    return lt_matches

def get_format_info(df_vcf, col_caller, info="FT") -> pd.Series:
    mask_caller = ~(df_vcf[col_caller].str.startswith("."))
    col_position = "position_%s" % (info)
    df_vcf.loc[:, col_position] = np.nan

    def _get_pos(x):
        if info not in x:
            return np.nan
        else:
            return x.split(":").index(info)

    def _get_val(x):
        if np.isnan(x[1]):
            return np.nan
        else:
            return x[0].split(":")[int(x[1])]

    df_vcf.loc[mask_caller, col_position] = df_vcf.loc[mask_caller, "FORMAT"].apply(_get_pos)
    s_filter = df_vcf.loc[mask_caller, [col_caller, col_position]].apply(_get_val, axis=1)
    del df_vcf[col_position]

    return s_filter

def add_chr_prefix(filepath_orig, filepath_dest):
    with open(filepath_orig, 'r') as vcf_file_orig, open(filepath_dest, 'w') as vcf_file_dest:
        vcf_lines = vcf_file_orig.readlines()

        for vcf_line in vcf_lines:
            if vcf_line.startswith("#"):
                vcf_file_dest.write(vcf_line)
            elif any(vcf_line.startswith(str(x)) for x in list(range(23)) + ["X", "Y"]):
                vcf_line = "chr" + vcf_line
                vcf_file_dest.write(vcf_line)
            else:
                vcf_file_dest.write(vcf_line)

def remove_chr_prefix(filepath_orig, filepath_dest):
    with open(filepath_orig, 'r') as vcf_file_orig, open(filepath_dest, 'w') as vcf_file_dest:
        vcf_lines = vcf_file_orig.readlines()

        for vcf_line in vcf_lines:
            if vcf_line.startswith("#"):
                vcf_file_dest.write(vcf_line)
            elif vcf_line.startswith("chr"):
                vcf_file_dest.write(vcf_line[3:])
            else:
                vcf_file_dest.write(vcf_line)



def process_info(df_vcf: DataFrame) -> DataFrame:
    """
    Processes the INFO field with contains the key=value fields like
        DB: dbSNP membership
        DP: read depth
        Gene: Hugo gene symbol
        MQ0: Number of mapping quality zero reads
        SS: Somatic status
        VC: variant classification
        VT: variant type
        TID: ensembl transcript identifier
        VLSC: final somatic score between 0 and 255
    """
    def keyval_to_dict(x):
        if x == ".":
            return np.nan
        else:
            return dict(item.split("=") if "=" in item else [item, item] for item in x.split(";"))

    df_vcf_info = df_vcf["INFO"].apply(keyval_to_dict).apply(pd.Series)
    df_vcf_info = df_vcf_info.loc[:, df_vcf_info.isnull().mean(axis=0) < 1]

    return df_vcf_info

def extract_n_reads(df_vcf_reads: DataFrame, infos_n_reads: list) -> DataFrame:
    """
    TODO: clarify the rules/meaning of DP, AD, DP4 and other.
    """
    df_n_reads = pd.DataFrame(index = df_vcf_reads.index)

    if set(["DP", "AD", "TIR", "TAR", "DP4"]).issubset(set(infos_n_reads)):
        #### depth
        df_n_reads.loc[:, "depth"] = df_vcf_reads["DP"].astype(float)

        #### ref_count, alt_count in this order: TIR/TAR, DP4
        if df_vcf_reads["TAR"].isnull().mean() < 1:
            df_tar_counts = df_vcf_reads["TAR"].str.split(",").apply(pd.Series)
            df_n_reads.loc[:,"ref_count"] = df_tar_counts[0]
        else:
            df_n_reads.loc[:,"ref_count"] = np.nan


        if df_vcf_reads["TIR"].isnull().mean() < 1:
            df_tir_counts = df_vcf_reads["TIR"].str.split(",").apply(pd.Series)
            df_n_reads.loc[:,"alt_count"] = df_tir_counts[0]
        else:
            df_n_reads.loc[:,"alt_count"] = np.nan

        if df_vcf_reads["DP4"].isnull().mean() < 1:
            df_dp4_counts = df_vcf_reads["DP4"].str.split(",").apply(pd.Series)

            mask_null = df_n_reads["ref_count"].isnull()
            mask = (mask_null) & (df_dp4_counts[[0,1]].isnull().sum(axis=1) == 0)
            df_n_reads.loc[mask,"ref_count"] = df_dp4_counts.loc[mask,[0,1]].astype(float).agg(sum, 1)

            mask_null = df_n_reads["alt_count"].isnull()
            mask = (mask_null) & (df_dp4_counts[[2,3]].isnull().sum(axis=1) == 0)
            df_n_reads.loc[mask,"alt_count"] = df_dp4_counts.loc[mask,[2,3]].astype(float).agg(sum, 1)

        #### set alt count to 0 where GT is 0/0
        #### and fill ref count by depth where applicable
        mask_gt = df_vcf_reads.GT == "0/0"
        mask_null = df_n_reads["alt_count"].isnull()
        mask = mask_gt & mask_null
        df_n_reads.loc[mask,"alt_count"] = 0

        mask_null = df_n_reads["ref_count"].isnull()
        mask = mask_gt & mask_null
        df_n_reads.loc[mask,"ref_count"] = df_n_reads.loc[mask,"depth"]

    elif set(["DP", "AD", "DP4"]).issubset(set(infos_n_reads)):
        #### depth
        df_n_reads.loc[:, "depth"] = df_vcf_reads["DP"].astype(float)

        #### ref_count, alt_count in this order: AD, DP4
        if df_vcf_reads["AD"].isnull().mean() < 1:
            df_ad_counts = df_vcf_reads["AD"].str.split(",").apply(pd.Series)
            df_n_reads.loc[:,"ref_count"] = df_ad_counts[0]
            df_n_reads.loc[:,"alt_count"] = df_ad_counts[1]
        else:
            df_n_reads.loc[:,"ref_count"] = np.nan
            df_n_reads.loc[:,"alt_count"] = np.nan

        if df_vcf_reads["DP4"].isnull().mean() < 1:
            df_dp4_counts = df_vcf_reads["DP4"].str.split(",").apply(pd.Series)

            mask_null = df_n_reads["ref_count"].isnull()
            mask = (mask_null) & (df_dp4_counts[[0,1]].isnull().sum(axis=1) == 0)
            df_n_reads.loc[mask,"ref_count"] = df_dp4_counts.loc[mask,[0,1]].astype(float).agg(sum, 1)

            mask_null = df_n_reads["alt_count"].isnull()
            mask = (mask_null) & (df_dp4_counts[[2,3]].isnull().sum(axis=1) == 0)
            df_n_reads.loc[mask,"alt_count"] = df_dp4_counts.loc[mask,[2,3]].astype(float).agg(sum, 1)

        #### set alt count to 0 where GT is 0/0
        #### and fill ref count by depth where applicable
        mask_gt = df_vcf_reads.GT == "0/0"
        mask_null = df_n_reads["alt_count"].isnull()
        mask = mask_gt & mask_null
        df_n_reads.loc[mask,"alt_count"] = 0

        mask_null = df_n_reads["ref_count"].isnull()
        mask = mask_gt & mask_null
        df_n_reads.loc[mask,"ref_count"] = df_n_reads.loc[mask,"depth"]

    elif set(["DP", "TAR", "TIR"]).issubset(set(infos_n_reads)):
        #### extracts depth=DP, ref_count=TAR[0], alt_count=TIR[1]
        df_n_reads.loc[:,"depth"]     = df_vcf_reads["DP"]

        if df_vcf_reads["TAR"].isnull().mean() < 1:
            df_tar_counts = df_vcf_reads["TAR"].str.split(",").apply(pd.Series)
            df_n_reads.loc[:,"ref_count"] = df_tar_counts[0]
        else:
            df_n_reads.loc[:,"ref_count"] = np.nan

        if df_vcf_reads["TIR"].isnull().mean() < 1:
            df_tir_counts = df_vcf_reads["TIR"].str.split(",").apply(pd.Series)
            df_n_reads.loc[:,"alt_count"] = df_tir_counts[0]
        else:
            df_n_reads.loc[:,"alt_count"] = np.nan

    elif set(["RC", "AC"]).issubset(set(infos_n_reads)):
        #### extracts depth=RC+AC, ref_count=RC, alt_count=AC

        df_n_reads.loc[:,"depth"] = df_vcf_reads[["RC", "AC"]].astype(float).agg(sum, axis=1)
        df_n_reads.loc[:,"ref_count"] = df_vcf_reads["RC"].astype(float)
        df_n_reads.loc[:,"alt_count"] = df_vcf_reads["AC"].astype(float)

    df_n_reads = df_n_reads.astype(float)
    return df_n_reads

def process_reads(df_vcf: DataFrame, col_normal: str, col_tumor: str, infos_n_reads: list, infos_other: list) -> DataFrame:
    """
    Extract read counts and somatic status in normal and tumor samples from the following info fields:
        SS
        AD or DP4
        DP
        GT
        RC
        AC
        RS
    """
    sources = [col_normal, col_tumor]
    prefixs = ["n", "t"]
    infos   = infos_n_reads + infos_other

    ddf_vcf_reads = {}
    for source in sources:
        ddf_vcf_reads[source] = pd.DataFrame(index = df_vcf.index)

        for info in infos:
            s_info = get_format_info(df_vcf, source, info)
            ddf_vcf_reads[source].insert(0, info, get_format_info(df_vcf, source, info))
            ddf_vcf_reads[source] = ddf_vcf_reads[source].replace([".", ".,.", ".,.,.", ".,.,.,."], np.nan)

    for source, prefix in zip(sources, prefixs):
        df_vcf_reads = ddf_vcf_reads[source]

        #### extract read counts
        df_n_reads = extract_n_reads(
            df_vcf_reads  = df_vcf_reads,
            infos_n_reads = infos_n_reads
        )

        df_vcf_reads = pd.concat((df_vcf_reads, df_n_reads), axis=1)
        df_vcf_reads = df_vcf_reads.rename(
            columns = {x: "%s_%s" % (prefix, x) for x in df_vcf_reads.columns}
        )
        ddf_vcf_reads[source] = df_vcf_reads

    df_reads = pd.concat(ddf_vcf_reads, axis=1)
    df_reads.columns = df_reads.columns.droplevel()

    return df_reads

def process_assemble(df_vcf: DataFrame, df_vcf_info: DataFrame, df_vcf_reads: DataFrame) -> DataFrame:
    df_vcf_post = pd.concat((df_vcf, df_vcf_info, df_vcf_reads), axis=1)

    rename_columns = {
        "Gene"   : "Hugo_Symbol",
        "#CHROM" : "Chromosome",
        "POS"    : "Position",
        "ID"     : "dbSNP_RS",
        "REF"    : "Reference_Allele",
        "REF"    : "Tumor_Seq_Allele1",
        "ALT"    : "Tumor_Seq_Allele2",
        "QUAL"   : "Variant_Quality",
        "FILTER" : "Filter",
        "VC"     : "Variant_Classification",
        "VT"     : "Variant_Type",
        "TID"    : "Transcript_ID",
    }

    df_vcf_post = df_vcf_post.rename(columns=rename_columns)

    count_columns = [x for x in df_vcf_post.columns if x.startswith("n_") or x.startswith("t_")]
    other_columns = [x for x in df_vcf_post.columns if x in rename_columns.values()]

    df_vcf_post = df_vcf_post[other_columns + count_columns]
    df_vcf_post = df_vcf_post.replace(to_replace=".", value=np.nan)

    return df_vcf_post

#### # MAIN FUNCTION
#### #####################################################################################################

def run_manual_annotator(vcf_path: str, out_path:str,  col_normal: str, col_tumor: str, infos_n_reads: list, infos_other: list):
    """
    Manually parse VCF file and save at the path specified.

    Paramters
    ---------
    vcf_path: str
        path to the vcf file
    out_path: str
        path where output should be saved
    col_normal: str
        name of the column in the vcf for the normal sample
    col_tumor: str
        name of the column in the vcf for the tumor sample
    infos_n_reads: list
        list of sigles that contain read info
    infos_other: list
        list of sigles that need extraction
    """

    #### read the filtered vcf
    df_vcf = load_vcf(
        filepath = vcf_path,
    )

    #### process variant info
    df_vcf_info = process_info(df_vcf)

    #### process normal and primary/metastatic reads details
    df_vcf_reads = process_reads(
        df_vcf        = df_vcf,
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )
    cols_reads = df_vcf_reads.columns.tolist()

    #### assemble the dataframes
    df_vcf_post = process_assemble(df_vcf, df_vcf_info, df_vcf_reads)

    #### SAVE
    df_vcf_post.to_csv(
        path_or_buf = out_path,
        sep         = "\t",
        header      = True,
        index       = False
    )
    print("manually extracted vep saved at %s\n" % out_path, flush=True)
