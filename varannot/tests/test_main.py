# -*- coding: utf-8 -*-
"""
@created: Aug 13 2020
@modified: Oct 30 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Test functions from vep module.
"""

import os
from .._util import set_wd_to_repo
from .._main import run_annotator
from .._main import VepConfig
from .._main import Vcf2mafConfig

def test_main():
    vep_config = VepConfig(
        data             = "~/.vep",
        n_fork           = 4,
        fasta            = "~/.vep/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa",
        custom_run       = True,
        custom_overwrite = True
    )
    vcf2maf_config = Vcf2mafConfig(
        run       = True,
        overwrite = True
    )

    current_wd = set_wd_to_repo()

    #### # 1. TCGA GA
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_GA/"
    out_folder = "./examples/results/TCGA_GA/"
    os.makedirs(out_folder, exist_ok=True)

    #### paths to results folders
    dt_folders = {
        'manual_out_folder'  : os.path.join(out_folder, "tmp/out_manual"),
        'vcf2maf_tmp_folder' : os.path.join(out_folder, "tmp/tmp_vcf2maf"),
        'vcf2maf_out_folder' : os.path.join(out_folder, "tmp/out_vcf2maf"),
        'vep_out_folder'     : os.path.join(out_folder, "tmp/out_vep"),
        'maf_folder'         : os.path.join(out_folder, "maf"),
    }

    #### make folders if they do not exist already
    for k, v in dt_folders.items():
        os.makedirs(v, exist_ok=True)

    #### Indel TCGA_GA
    vcf_file      = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.indel.capture.tcga.vcf"
    col_normal    = "NORMAL"
    col_tumor     = "PRIMARY"
    normal_id     = "TCGA-A1-A0SD-10A-01D-A110-09"
    tumor_id      = "TCGA-A1-A0SD-01A-11D-A10Y-09"
    infos_n_reads = ["AD", "DP4", "DP", "TAR", "TIR"]
    infos_other   = ["SS", "GT"]

    dt_identifiers = {
        "Tumor_Sample"                : "TCGA-A1-A0SD",
        "Tumor_Sample_Barcode"        : "TCGA-A1-A0SD-01A-11D-A10Y-09",
        "Matched_Norm_Sample_Barcode" : "TCGA-A1-A0SD-10A-01D-A110-09",
        "Tumor_Sample_Site"           : "01",
    }

    run_annotator(
        vcf_folder     = vcf_folder,
        vcf_file       = vcf_file,
        col_normal     = col_normal,
        col_tumor      = col_tumor,
        normal_id      = normal_id,
        tumor_id       = tumor_id,
        infos_n_reads  = infos_n_reads,
        infos_other    = infos_other,
        dt_folders     = dt_folders,
        dt_identifiers = dt_identifiers,
        vcf2maf_config = vcf2maf_config,
        vep_config     = vep_config
    )

    #### SNP TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.oxoG.snp.capture.tcga.vcf"

    #### # 2. TCGA HS
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_HS/"

    #### Indel TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.indel.0e81f9c986154ce89e59240c3f09534f.vcf"

    #### SNP TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.snv.0e81f9c986154ce89e59240c3f09534f.vcf"
