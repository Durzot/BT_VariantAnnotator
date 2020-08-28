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
from   typing import Union

def run_vep_annotator(vep_folder: str, vep_data: str, vcf_path: str, out_path: str, fasta: str, vep_custom: Union[str,list]=None, overwrite: bool=False, vep_n_fork: int=4):
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
    vep_custom: str or list
        additional options to add to the vep cmd. For instance
        '~/.vep/custom/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN'
    overwrite: bool, optional.
        if the output file already exists (from previous run), should it be overwritten?
    vep_n_fork: int, optional.
        number of forks to be used when running VEP.
    """
    vep = os.path.join(vep_folder, "vep")
    need_run = True

    if os.path.exists(out_path) and not overwrite:
        need_run = False

    if need_run:
        print("STATUS: RUNNING VEP")

        if os.path.exists(out_path):
            os.remove(out_path)
            print("removed existing file: %s" % out_path)

        cmd = """%s \
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
            --fork %s \
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
            --cache \
            --offline """ % (vep, vep_data, vep_n_fork, vcf_path, out_path, fasta)

        if vep_custom is not None:
            if type(vep_custom) == list:
                for v_custom in vep_custom:
                    cmd += "--custom %s " % v_custom
            elif type(vep_custom) == str:
                cmd += "--custom %s " % vep_custom
            else:
                raise ValueError("vep_custom should be of type list or str")

        os.system(cmd)
    else:
        print("output file %s already exists and overwrite is set to False" % out_path)
