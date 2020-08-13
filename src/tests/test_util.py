# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Test functions in util.py module.
"""

import os
from ..util import load_vcf

def test_load_vcf():
    folder = "./examples/data/TCGA_GA/"

    #### Indel with header
    file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.indel.capture.tcga.vcf"
    df_vcf = load_vcf(
        filepath = os.path.join(folder, file)
    )

    #### SNV with header
    file = "TCGA-A1-A0SB_db9d40fb-bfce-4c3b-a6c2-41c5c88982f1_a3254f8e-3bbd-42fc-abea-a5f25b7648b3.oxoG.snp.capture.tcga.vcf"
    df_vcf = load_vcf(
        filepath = os.path.join(folder, file)
    )
