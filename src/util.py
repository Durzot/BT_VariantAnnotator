# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Useful functions
"""


import numpy  as np
import os
import pandas as pd
import re

DataFrame = pd.core.frame.DataFrame

def load_vcf(filepath: str, no_header: bool=False) -> DataFrame:
    """
    Load VCF file from the specified filepath into a pandas DataFrame.

    Parameters
    ----------
    filepath:  str
        Path to the file.
    no_header:  bool
        If True, set column names to default names.

    Returns
    -------
    df: DataFrame
        The data loaded in a DataFrame
    """

    if not os.path.exists(filepath):
        raise ValueError("The file %s does not exist.")
    else:
        if no_header:
            df_vcf = pd.read_csv(
                filepath_or_buffer = filepath,
                sep                = "\t",
                skiprows           = 0,
                low_memory         = False,
            )

        else:
            skipsymbol = "##"
            with open(filepath, "r") as file:
                skiprows = sum(line.startswith(skipsymbol) for line in file.readlines())

            df_vcf = pd.read_csv(
                filepath_or_buffer = filepath,
                sep                = "\t",
                skiprows           = skiprows,
                low_memory         = False,
            )

    return df_vcf

def write_vcf(filepath_orig: str, filepath_dest: str, df_vcf: DataFrame) -> None:
    headersymbol = "##"
    headerrows = []

    with open(filepath_orig, "r") as file:
        while True:
            line = file.readline()
            if line.startswith(headersymbol):
                headerrows.append(line)
            else:
                break

    with open(filepath_dest, "w") as file:
        for line in headerrows:
            file.write(line)
        file.write(df_vcf.to_csv(sep="\t", index=False))
