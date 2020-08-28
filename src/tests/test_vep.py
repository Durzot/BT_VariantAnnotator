# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Test functions from vep module.
"""

import os
from ..vep import run_vep_annotator

def test_vep():
    vep_folder = "~/Documents/biotools/informatics/VCF/ensembl-vep"
    vep_data   = "~/.vep"
    vep_n_fork = 4
    fasta      = "~/.vep/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

    #### # 1. TCGA GA
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_GA/"
    out_folder = "./examples/results/TCGA_GA/tmp/out_vep"
    os.makedirs(out_folder, exist_ok=True)

    #### Indel TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.indel.capture.tcga.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")

    run_vep_annotator(
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        fasta      = fasta,
        vep_n_fork = vep_n_fork
    )

    #### SNP TCGA_GA
    vcf_file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.oxoG.snp.capture.tcga.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")

    run_vep_annotator(
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        fasta      = fasta,
        vep_n_fork = vep_n_fork
    )

    #### # 2. TCGA HS
    #### # ########################################################################################################

    vcf_folder = "./examples/data/TCGA_HS/"
    out_folder = "./examples/results/TCGA_HS/tmp/out_vep"
    os.makedirs(out_folder, exist_ok=True)

    #### Indel TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.indel.0e81f9c986154ce89e59240c3f09534f.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")

    run_vep_annotator(
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        fasta      = fasta,
        vep_n_fork = vep_n_fork
    )

    #### SNP TCGA_HS
    vcf_file =  "genome.wustl.edu.TCGA-A1-A0SD.snv.0e81f9c986154ce89e59240c3f09534f.vcf"
    out_file = vcf_file.replace(".vcf", ".txt")

    run_vep_annotator(
        vep_folder = vep_folder,
        vep_data   = vep_data,
        vcf_path   = os.path.join(vcf_folder, vcf_file),
        out_path   = os.path.join(out_folder, out_file),
        fasta      = fasta,
        vep_n_fork = vep_n_fork
    )
