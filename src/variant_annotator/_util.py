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

#### modify if the repository was cloned under a different name
REPO_FOLDER = "VariantAnnotator"

def set_wd_to_repo():
    current_wd = os.getcwd()
    if REPO_FOLDER not in os.getcwd():
        raise ValueError("Please set the working directory to a location in the repository %s" % REPO_FOLDER)
    else:
        while not os.getcwd().endswith(REPO_FOLDER):
            os.chdir("..")
    return current_wd

def get_path_to_repo() -> str:
    current_wd = os.getcwd()
    if REPO_FOLDER not in os.getcwd():
        raise ValueError("Please set the working directory to a location in the repository %s" % REPO_FOLDER)
    else:
        while not os.getcwd().endswith(REPO_FOLDER):
            os.chdir("..")
    repo_path = os.getcwd()
    os.chdir(current_wd)
    return repo_path


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
        raise ValueError("The file %s does not exist." % filepath)
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
