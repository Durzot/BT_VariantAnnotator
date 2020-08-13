# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Python wrapper around vcf2maf perl script.
"""

import os

def run_vcf2maf_annotator(vcf2maf: str, vep_folder: str, vep_data: str, vcf_path: str, out_path: str, tmp_folder: str,
tumor_id: str, normal_id: str, fasta: str, overwrite: bool=False):
    """
    Run vcf2maf reannotator. Details may found at https://github.com/mskcc/vcf2maf.

    Parameters
    ----------
    vcf2maf: str
        path to the vcf2maf perl script
    vep_folder: str
        path to the folder where the vep command is
    vep_data: str
        path to the .vep data where the reference genome is located
    vcf_path: str
        path to the vcf file
    out_path: str
        path where output should be saved
    tmp_folder: str
        dir where intermediary vep files are saved
    tumor_id: str
        identifier of tumor sample
    tumor_id: str
        identifier of normal sample
    fasta: str
        relative path to fasta file from vep_folder
    overwrite: bool
        if the output file already exists (from previous run), should it be overwritten?
    """
    need_run = True

    if os.path.exists(out_path) and not overwrite:
        need_run = False
    elif os.path.exists(out_path):
        os.remove(out_path)

    if need_run:
        os.system('perl %s \
            --input-vcf %s \
            --output-maf %s \
            --tmp-dir %s \
            --tumor-id %s \
            --normal-id %s \
            --vep-path %s \
            --vep-data %s \
            --buffer-size 5000 \
            --vep-forks 4 \
            --ncbi-build GRCh37 \
            --ref-fasta %s \
            --filter-vcf 0' % (vcf2maf, vcf_path, out_path, tmp_folder, tumor_id, normal_id, vep_folder, vep_data, fasta)
        )
    else:
        print("output file %s already exists and overwrite is set to False" % out_path)
