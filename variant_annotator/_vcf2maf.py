# -*- coding: utf-8 -*-
"""
@created: Aug 13 2020
@modified: Oct 30 2020
@author: Yoann Pradat
    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Python wrapper around vcf2maf perl script.
"""

import os

def run_vcf2maf_annotator(vep_data: str, vep_n_fork: int, vcf_path: str, out_path: str, tmp_folder: str, tumor_id: str, normal_id: str, fasta: str, overwrite: bool=False):
    """
    Run vcf2maf reannotator. Details may found at https://github.com/mskcc/vcf2maf.
    Parameters
    ----------
    vep_data: str
        path to the .vep data where the reference genome is located
    vep_n_fork: int
        number of forks to be used by VEP.
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
    vcf2maf_path = os.path.normpath(os.path.join(__file__, "../../tools/vcf2maf/vcf2maf.pl"))
    vep_path     = os.path.normpath(os.path.join(__file__, "../../tools/ensembl-vep"))

    need_run = True
    vcf_file = out_path.split("/")[-1]
    tmp_file = vcf_file.replace(".txt", ".vep.vcf")
    tmp_path = os.path.join(tmp_folder, tmp_file)

    if os.path.exists(out_path) and not overwrite:
        need_run = False

    if need_run:
        print("STATUS: RUNNING VCF2MAF")

        if os.path.exists(tmp_path):
            os.remove(tmp_path)
            print("removed existing file: %s" % tmp_path)

        if os.path.exists(out_path):
            os.remove(out_path)
            print("removed existing file: %s" % out_path)

        os.system('perl %s \
            --input-vcf %s \
            --output-maf %s \
            --tmp-dir %s \
            --tumor-id %s \
            --normal-id %s \
            --vep-path %s \
            --vep-data %s \
            --buffer-size 5000 \
            --vep-forks %d \
            --ncbi-build GRCh37 \
            --ref-fasta %s \
            --filter-vcf 0' % \
            (vcf2maf_path, vcf_path, out_path, tmp_folder, tumor_id, normal_id, vep_path, vep_data, vep_n_fork, fasta)
        )
    else:
        print("output file %s already exists and overwrite is set to False" % out_path)
