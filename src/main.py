# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Main functions for running each step and the assembling step.
"""

def run_annotator(i_split: int, n_split: int, vcf2maf: str, vep_folder: str, vep_data: str, fasta: str):
    """
    Run the manual, vcf2maf and vep annotations and assemble.


    Paramters
    --------
    i_split: int
        the split processed.
    n_split: int
        the number of splits for processing a set of vcf files
    vcf2maf: str
        path to the vcf2maf perl script
    vep_folder: str
        path to the folder where the vep command is
    vep_data: str
        path to the .vep data where the reference genome is located
    fasta: str
        relative path to fasta file from vep_folder
    """


#### # SCRIPT PARAMETERS 
#### #####################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--i_split'    , type=int , default=1  , help='the split processed')
parser.add_argument('--n_split'    , type=int , default=1  , help='total number of splits')
parser.add_argument('--vcf2maf'    , type=str , default="" , help='path to the vcf2maf perl script')
parser.add_argument('--vep_folder' , type=str , default="" , help='path to the folder of the vep command')
parser.add_argument('--vep_data'   , type=str , default="" , help='path to the .vep data folder')
parser.add_argument('--fasta'      , type=str , default="" , help='path to reference genome FASTA file')
args = parser.parse_args()

print("Parameters", flush=True)
for arg in vars(args):
    print("%s: %s" % (arg, getattr(args, arg)), flush=True)


if __name__ == "__main__":

    #### change working dir
    current_wd = setwd_to_data()

    #### Parameters
    dt_paths = {
        'vcf_meta'        : './raw/vcf/mbc/info/mbc_vcf_identifiers.txt',
        'vcf_dir'         : './raw/vcf/mbc/files_filtered/',
        'vcf_list'        : './raw/vcf/mbc/info/list_vcf_mbc.txt',
        'vep_manual_dir'  : './raw/vep/vep_manual/mbc',
        'vep_alone_dir'   : './raw/vep/vep_alone/mbc',
        'vep_vcf2maf_dir' : './raw/vep/vep_vcf2maf/mbc/vep',
        'maf_vcf2maf_dir' : './raw/vep/vep_vcf2maf/mbc/maf',
        'maf_one_dir'     : './raw/maf/mbc',
    }

    #### make dir if does not exist already
    for k, v in dt_paths.items():
        if "dir" in k:
            os.makedirs(v, exist_ok=True)

    #### load meta data
    df_meta = pd.read_csv(
        filepath_or_buffer = dt_paths["vcf_meta"],
        sep                = "\t"
    )

    if os.path.exists(dt_paths["vcf_list"]):
        with open(dt_paths["vcf_list"]) as file:
            vcf_names_orig = file.read().splitlines()
    else:
        vcf_names_orig = df_meta.loc[df_meta.keep == True, "vcf"].tolist()

    #### get list of vcfs for the split
    count_one_split = len(vcf_names_orig)//args.n_split

    if args.i_split == args.n_split:
        vcf_names_orig  = vcf_names_orig[(args.i_split-1)*count_one_split:]
    else:
        vcf_names_orig  = vcf_names_orig[(args.i_split-1)*count_one_split:args.i_split*count_one_split]

    count = 0
    count_total = len(vcf_names_orig)

    #### loop over the list
    for vcf_name_orig in vcf_names_orig:
        vcf_name = vcf_name_orig
        vcf_path = os.path.join(dt_paths['vcf_dir'], vcf_name)

        count += 1

        print("="*80, flush=True)
        print("vcf %d/%d" % (count, count_total), flush=True)
        print("processing %s\n" % vcf_name, flush=True)

        #### get vcf type
        if "mutect" in vcf_name_orig:
            vcf_type = "snv"
        elif "indel" in vcf_name_orig:
            vcf_type = "indel"

        #### get vcf identifiers
        mask_vcf_name = df_meta.vcf == vcf_name_orig
        index_vcf_name = mask_vcf_name[mask_vcf_name].index[0]

        vcf_identifiers = {
            "Tumor_Sample"                : df_meta.loc[index_vcf_name, "Subject.ID"],
            "Tumor_Sample_Barcode"        : df_meta.loc[index_vcf_name, "tumor_sample"],
            "Matched_Norm_Sample_Barcode" : df_meta.loc[index_vcf_name, "normal_sample"],
        }

        if vcf_type == "snv":
            infos_n_reads = ["RC", "AC"]
            infos_other = ["RS", "AS", "GT"]

            col_tumor  = vcf_identifiers["Tumor_Sample_Barcode"]
            col_normal = vcf_identifiers["Matched_Norm_Sample_Barcode"]

        elif vcf_type == "indel":
            infos_n_reads = ["DP", "TAR", "TIR"]
            infos_other = ["SGT"]

            col_tumor  = "TUMOR"
            col_normal = "NORMAL"

        process_vcf(
            vcf_identifiers = vcf_identifiers,
            vcf_path        = vcf_path,
            vcf_name        = vcf_name,
            vcf2maf         = args.vcf2maf,
            vep_folder      = args.vep_folder,
            vep_data        = args.vep_data,
            fasta           = args.fasta,
            dt_paths        = dt_paths,
            col_tumor       = col_tumor,
            col_normal      = col_normal,
            infos_n_reads   = infos_n_reads,
            infos_other     = infos_other
        )

    #### revert to previous work dir
    os.chdir(current_wd)
