# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Test functions from vcf2maf module.
"""

import os
from ..vcf2maf import run_vcf2maf_annotator

def test_vcf2maf():
    vcf2maf    = "~/Documents/biotools/informatics/VCF/mskcc-vcf2maf-5453f80/vcf2maf.pl"
    vep_folder = "~/Documents/biotools/informatics/VCF/ensembl-vep"
    vep_data   = "~/.vep"
    fasta      = "~/.vep/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

    #### # 1. TCGA GA
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_GA/"
    tmp_folder = "./examples/results/TCGA_GA/tmp_vcf2maf"
    out_folder = "./examples/results/TCGA_GA/out_vcf2maf"
    os.makedirs(out_folder, exist_ok=True)

    #### Indel TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.indel.capture.tcga.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    normal_id = "TCGA-A1-A0SD-10A-01D-A110-09"
    tumor_id  = "TCGA-A1-A0SD-01A-11D-A10Y-09"

    run_vcf2maf_annotator(
        vcf2maf    = vcf2maf,
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        tmp_folder = tmp_folder,
        tumor_id   = tumor_id,
        normal_id  = normal_id,
        fasta      = fasta
    )

    #### SNP TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.oxoG.snp.capture.tcga.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    normal_id = "TCGA-A1-A0SD-10A-01D-A110-09"
    tumor_id  = "TCGA-A1-A0SD-01A-11D-A10Y-09"

    run_vcf2maf_annotator(
        vcf2maf    = vcf2maf,
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        tmp_folder = tmp_folder,
        tumor_id   = tumor_id,
        normal_id  = normal_id,
        fasta      = fasta
    )

    #### # 2. TCGA HS
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_HS/"
    tmp_folder = "./examples/results/TCGA_HS/tmp_vcf2maf"
    out_folder = "./examples/results/TCGA_HS/out_vcf2maf"
    os.makedirs(out_folder, exist_ok=True)

    #### Indel TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.indel.0e81f9c986154ce89e59240c3f09534f.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    normal_id = "TCGA-A1-A0SD-10A-01D-A110-09"
    tumor_id  = "TCGA-A1-A0SD-01A-11D-A10Y-09"

    run_vcf2maf_annotator(
        vcf2maf    = vcf2maf,
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        tmp_folder = tmp_folder,
        tumor_id   = tumor_id,
        normal_id  = normal_id,
        fasta      = fasta
    )

    #### SNP TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.snv.0e81f9c986154ce89e59240c3f09534f.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")
    normal_id = "TCGA-A1-A0SD-10A-01D-A110-09"
    tumor_id  = "TCGA-A1-A0SD-01A-11D-A10Y-09"

    run_vcf2maf_annotator(
        vcf2maf    = vcf2maf,
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        tmp_folder = tmp_folder,
        tumor_id   = tumor_id,
        normal_id  = normal_id,
        fasta      = fasta
    )