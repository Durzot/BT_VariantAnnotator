# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Python wrapper around VEP command.
"""

import os

def run_vep_annotator(vep_folder: str, vep_data: str, vcf_path: str, out_path: str, fasta: str, overwrite: bool=False):
    """
    Run variant ensembl predictor alone with custom options. See options details at
    https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_af

    Parameters
    ---------
    vep_folder: str
        path to the folder where the vep command is
    vep_data: str
        path to the .vep data where the reference genome is located
    vcf_path: str
        path to the vcf file
    out_path: str
        path where output should be saved
    fasta: str
        relative path to fasta file from vep_folder
    overwrite: bool
        if the output file already exists (from previous run), should it be overwritten?
    """
    vep = os.path.join(vep_folder, "vep")
    need_run = True

    if os.path.exists(out_path) and not overwrite:
        need_run = False
    elif os.path.exists(out_path):
        os.remove(out_path)

    if need_run:
        os.system('%s \
            --dir %s \
            --af \
            --af_gnomad \
            --af_esp \
            --clin_sig_allele 0 \
            --max_af \
            --af_1k \
            --no_progress \
            --no_stats \
            --appris \
            --biotype \
            --buffer_size 500 \
            --canonical \
            --ccds \
            --check_existing \
            --distance 5000 \
            --hgvs \
            --fork 4 \
            --numbers \
            --mane \
            --pick \
            --polyphen b \
            --protein \
            --pubmed \
            --regulatory \
            --sift b \
            --species homo_sapiens \
            --symbol \
            --transcript_version \
            --tsl \
            --uniprot \
            --input_file %s \
            --output_file %s \
            --fasta %s \
            --offline ' % (vep, vep_data, vcf_path, out_path, fasta)
        )
    else:
        print("output file %s already exists and overwrite is set to False" % out_path)
