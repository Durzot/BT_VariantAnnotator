# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Test functions from manual module.
"""

import os
from ..manual import run_manual_annotator

def test_manual():

    #### # 1. TCGA GA
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_GA/"
    out_folder = "./examples/results/TCGA_GA/out_manual"
    os.makedirs(out_folder, exist_ok=True)

    #### Indel TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.indel.capture.tcga.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    col_normal = "NORMAL"
    col_tumor  = "PRIMARY"
    infos_n_reads = ["AD", "DP", "FA"]
    infos_other = ["SS", "GT"]

    run_manual_annotator(
        vcf_path      = os.path.join(vcf_folder, vcf_file),
        out_path      = os.path.join(out_folder, out_file),
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )

    #### SNP TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.oxoG.snp.capture.tcga.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    col_normal = "NORMAL"
    col_tumor  = "PRIMARY"
    infos_n_reads = ["AD", "DP", "FA"]
    infos_other = ["SS", "GT"]

    run_manual_annotator(
        vcf_path      = os.path.join(vcf_folder, vcf_file),
        out_path      = os.path.join(out_folder, out_file),
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )

    #### # 2. TCGA HS
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_HS/"
    out_folder = "./examples/results/TCGA_HS/out_manual"
    os.makedirs(out_folder, exist_ok=True)

    #### Indel TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.indel.0e81f9c986154ce89e59240c3f09534f.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    col_normal = "TCGA-A1-A0SD-10A-01D-A110-09"
    col_tumor  = "TCGA-A1-A0SD-01A-11D-A10Y-09"
    infos_n_reads = ["AD", "DP4", "DP", "TAR", "TIR"]
    infos_other = ["SS", "GT"]

    run_manual_annotator(
        vcf_path      = os.path.join(vcf_folder, vcf_file),
        out_path      = os.path.join(out_folder, out_file),
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )

    #### SNP TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.snv.0e81f9c986154ce89e59240c3f09534f.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    col_normal = "TCGA-A1-A0SD-10A-01D-A110-09"
    col_tumor  = "TCGA-A1-A0SD-01A-11D-A10Y-09"
    infos_n_reads = ["AD", "DP4", "DP", "FA"]
    infos_other = ["SS", "GT"]

    run_manual_annotator(
        vcf_path      = os.path.join(vcf_folder, vcf_file),
        out_path      = os.path.join(out_folder, out_file),
        col_normal    = col_normal,
        col_tumor     = col_tumor,
        infos_n_reads = infos_n_reads,
        infos_other   = infos_other
    )
