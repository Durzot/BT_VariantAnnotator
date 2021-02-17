# -*- coding: utf-8 -*-
"""
@created: Aug 13 2020
@modified: Oct 30 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Main functions for running each step and the assembling step.
"""
import os
import numpy  as np
import pandas as pd
from   typing import Union

from   ._manual  import run_manual_annotator
from   ._vcf2maf import run_vcf2maf_annotator
from   ._vep     import run_vep_annotator

DataFrame = pd.core.frame.DataFrame

#### # FUNCTION FOR ONE VCF
#### # #######################################################################################################

from dataclasses import dataclass, field

@dataclass
class VepConfig:
    """
    Config for running VEP inside VCF2MAF and separately (custom options, optional).

    Parameters
    --------
    data: str
        path to the .vep data where the reference genome is located. Default: $HOME/.vep
    fasta: str
        relative path to fasta file from folder. Default
        "$HOME/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
    n_fork: int, optional.
        number of forks to be used when running VEP. Use at least 2.
    custom_run: bool, optional
        set to True to run VEP separately from vcf2maf.
    custom_opt: str or list, optional.
        additional options to add to the vep cmd. For instance
        '--custom ~/.vep/custom/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN'
    custom_overwrite: bool, optional.
        set to True to overwrite any existing previous custom run of VEP.
    """
    data: str=os.path.expanduser("~/.vep")
    fasta: str=os.path.expanduser("~/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa")
    n_fork: int=4
    custom_run: bool=False
    custom_opt: Union[str, list]=None
    custom_overwrite: bool=False

@dataclass
class Vcf2mafConfig:
    """
    Run vcf2maf. For VEP-related options, see VepConfig class.

    Parameters
    --------
    run: bool, optional
        set to False to not use vcf2maf.
    overwrite: bool, optional.
        set to True to overwrite any existing previous run of vcf2maf.
    """
    run: bool=True
    overwrite: bool=False

def run_annotator(vcf_folder: str, vcf_file: str, col_normal: str, col_tumor: str, infos_n_reads: list,
                  infos_other: list, dt_folders: dict, vcf2maf_config: Vcf2mafConfig, vep_config: VepConfig,
                  dt_identifiers: dict=None,  tumor_id: str=None, normal_id: str=None) -> None:
    """
    Run the manual, vcf2maf and/or vep annotations on one VCF file and assemble.

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
    dt_folders: dict
        dict with the following keys:
            * manual_out_folder
            * vcf2maf_tmp_folder
            * vcf2maf_out_folder
            * vep_out_folder
            * maf_folder
    vcf2maf_config: object
        See Vcf2mafConfig class.
    vep_config: object
        See VepConfig class.
    dt_identifiers: dict, optional
        dict with key, value pairs that will be added as single-value columns in the maf file
    tumor_id: str, optional.
        Tumor identifier if different from col_tumor (e.g if col_tumor = "PRIMARY").
    normal_id: str, optional.
        Normal identifier if different from col_normal (e.g if col_normal = "NORMAL").
    """
    vcf_path = os.path.join(vcf_folder, vcf_file)
    out_file = vcf_file.replace(".vcf", ".txt")

    if tumor_id is None:
        tumor_id = col_tumor

    if normal_id is None:
        normal_id = col_normal

    #### # 1. RUN EACH ANNOTATOR
    #### # ###################################################################################################

    manual_out_path  = os.path.join(dt_folders["manual_out_folder"], out_file)
    run_manual_annotator(
        vcf_path      = vcf_path,
        out_path      = manual_out_path,
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )

    if vcf2maf_config.run:
        vcf2maf_out_path = os.path.join(dt_folders["vcf2maf_out_folder"], out_file)
        run_vcf2maf_annotator(
            vep_data     = vep_config.data,
            vep_n_fork   = vep_config.n_fork,
            vcf_path     = vcf_path,
            out_path     = vcf2maf_out_path,
            tmp_folder   = dt_folders["vcf2maf_tmp_folder"],
            tumor_id     = tumor_id,
            normal_id    = normal_id,
            fasta        = vep_config.fasta,
            overwrite    = vcf2maf_config.overwrite
        )

    if vep_config.custom_run:
        vep_out_path     = os.path.join(dt_folders["vep_out_folder"], out_file)
        run_vep_annotator(
            vep_data   = vep_config.data,
            vep_n_fork = vep_config.n_fork,
            vcf_path   = vcf_path,
            out_path   = vep_out_path,
            fasta      = vep_config.fasta,
            vep_custom = vep_config.custom_opt,
            overwrite  = vep_config.custom_overwrite,
        )

    #### # 2. ASSEMBLE ANNOTATIONS
    #### ######################################################################################################

    ddf_maf = {}

    #### vep manual
    ddf_maf["manual"] = pd.read_csv(
        filepath_or_buffer = manual_out_path,
        sep                = "\t"
    )
    if dt_identifiers is not None:
        for k,v in dt_identifiers.items():
            ddf_maf["manual"].insert(0, k, v)

    #### maf vcf2maf
    if vcf2maf_config.run:
        skipsymbol = "#"
        with open(vcf2maf_out_path) as file:
            skiprows = sum(line.startswith(skipsymbol) for line in file.readlines())

        ddf_maf["vcf2maf"] = pd.read_csv(
            filepath_or_buffer = vcf2maf_out_path,
            sep                = "\t",
            skiprows           = skiprows
        )

    #### vep custom
    if vep_config.custom_run:
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

        ddf_maf["alone"] = df_alone

    #### check that all scripts ran ok the vcf file
    #### by check that all maf files have the same number of lines
    list_n_rows = [df_maf.shape[0] for df_maf in ddf_maf.values()]
    if len(set(list_n_rows)) > 1:
        print("SHAPE error: not all MAF have the same number of rows")

    #### choose column per source
    #### take all columns from VCF2MAF except for some that are from the manual extraction
    #### for instance: reads number are sometimes not extracted by VCF2MAF 
    maf_columns = []

    cols_manual = [x for x in ddf_maf["manual"].columns if x.startswith(("n_", "t_"))] + [
        "Variant_Quality",
        "Filter_VCF",
        "Tumor_Sample_UUID",
        "Matched_Norm_Sample_UUID",
    ]

    for column in cols_manual:
        if column not in ddf_maf["manual"].columns:
            #### UUID is for instance only available for TCGA projects
            print("WARNING: %s is not in manual" % column, flush=True)
        else:
            maf_columns.append(ddf_maf["manual"][column])

    if vcf2maf_config.run and vep_config.custom_run:
        #### add flags column source

        #### vcf2maf
        ddf_maf["vcf2maf"].columns = ["%s_VCF2MAF" % x for x in ddf_maf["vcf2maf"].columns]
        for column in ddf_maf["vcf2maf"].columns:
            if column in dt_identifiers.keys() or column in cols_manual:
                #### prioritize identifiers from input
                pass
            else:
                maf_columns.append(ddf_maf["vcf2maf"][column])

        #### vep
        ddf_maf["alone"].columns = ["%s_VEP" % x for x in ddf_maf["alone"].columns]
        for column in ddf_maf["alone"].columns:
            maf_columns.append(ddf_maf["alone"][column])

    elif vcf2maf_config.run:
        #### vcf2maf
        for column in ddf_maf["vcf2maf"].columns:
            if column in dt_identifiers.keys() or column in cols_manual:
                #### prioritize identifiers from input
                pass
            else:
                maf_columns.append(ddf_maf["vcf2maf"][column])

    elif vep_config.custom_run:
        #### vep
        for column in ddf_maf["alone"].columns:
            maf_columns.append(ddf_maf["alone"][column])

    #### add identifier columns if any
    if dt_identifiers is not None:
        for column in dt_identifiers.keys():
            if column not in ddf_maf["manual"].columns:
                print("WARNING: %s is not in manual" % column, flush=True)
            else:
                maf_columns.append(ddf_maf["manual"][column])

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
