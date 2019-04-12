#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import numpy as np
import pysam
from pysam import VariantFile

def get_args_snps():
    """
    Get commandline arguments
    Args: None
    Returns:
    args (Namespace): the parsed arguments.
    """
    parser = argparse.ArgumentParser(description='Command Line arguments for snp detection in primers.')
    parser.add_argument("-vcf", "--variant-call-file", dest="vcf",
                        help="Tabix indexed VCF.")
    parser.add_argument("-in", "--primer-input-file", dest="p_info",
                        help="The output of the primer pipeline.")
    parser.add_argument("-o", "--output", dest="output",
                        help="The name of the output file")
    args = parser.parse_args()
    return args

def create_df(primer_pipeline_output):
    """
    Takes output of primer pipeline and generates dataframe.

    Args:
    total_primers (file): the output of the primer pipeline
    Returns:
    dataframe (pd.DataFrame): a pandas dataframe
    """
    primer_df = pd.read_csv(primer_pipeline_output, header=0)
    return primer_df

def primer_range_left(seqid, rank, chrm, p_left, position1):
    """
    Takes the chromosome, primer sequence, and position and creates a
    query range for tabix to search by.

    Args:
        seqid (pd.Series): pandas column containing seqids
        rank (pd.Series): pandas column containing primer rank
        chrm (pd.Series): pandas column of chrm info
        p_left (pd.Series): pandas column of the left primer seq
        position1 (pd.Series): pandas column of the left primer position.
    Returns:
        p_left_info (pd.DataFrame): the positional info in a frame format.
    """
    p_left_info = pd.DataFrame()
    p_left_info['Sequence ID'] = seqid
    p_left_info['Primer Rank'] = rank
    p_left_info['Chromosome'] = chrm.apply(str)
    p_left_info['P_Len'] = p_left.apply(lambda x: len(x))
    p_left_info['Position1'] = position1
    p_left_info['Position2'] = p_left_info['Position1']\
        + p_left_info['P_Len']
    return p_left_info

def primer_range_right(seqid, rank, chrm, p_right, position2):
    """
    Takes the chromosome, primer sequence, and position and creates a
    query range for tabix to search by.

    Args:
        seqid (pd.Series): pandas column containing seqids
        rank (pd.Series): pandas column containing primer rank
        chrm (pd.Series): pandas column of chrm info
        p_right (pd.Series): pandas column of the right primer seq
        position1 (pd.Series): pandas column of the right primer position.
    Returns:
        p_right_info (pd.DataFrame): the positional info in a frame format.
    """
    p_right_info = pd.DataFrame()
    p_right_info['Sequence ID'] = seqid
    p_right_info['Primer Rank'] = rank
    p_right_info['Chromosome'] = chrm.apply(str)
    p_right_info['P_Len'] = p_right.apply(lambda x: len(x))
    p_right_info['Position2'] = position2
    p_right_info['Position1'] = p_right_info['Position2']\
        - p_right_info['P_Len']
    return p_right_info


def match_pinfo_to_vcf(p_info, vcf):
    """
    Normalizes chromosome column to reference VCF info.
    Args:
        p_info (pd.DataFrame): the primer information df
        vcf (file): the tabix indexed vcf input
    Returns:
        p_info (pd.DataFrame): the primer information df normalized to genome
    """
    switch = 0
    try:
        for rec in vcf.fetch('chr1', 1000000, 1000050):
            print("Updating switch: 1")
        switch = 1
    except:
        for rec in vcf.fetch('1', 1000000, 1000050):
            print("Updating switch: 2")
        switch = 2

    if switch == 1 and not p_info['Chromosome']\
            .str.contains("chr").any():
        p_info['Chromosome'] = 'chr' + p_info['Chromosome']
    elif switch == 2 and p_info['Chromosome']\
            .str.contains("chr").any():
        p_info['Chromosome'] = p_info['Chromosome']\
            .str.replace('chr', '')
    else:
        pass
    return p_info


def tabix_fetch(seqids, ranks, chrom, position1, position2, vcf_in):
    """
    Takes p_info positions and fetches SNPs from tabix indexed VCF.
    Args:
        seqids (pd.Series): pandas column of seqids
        ranks (pd.Series): pandas column of primer ranks
        chrom (pd.Series): pandas column of chrm info
        position1 (pd.Series): pandas column of pos1 info
        position2 (pd.Series): pandas column of pos2 info
    Returns:
        snp_list (list): list containing vcf info for primer positions.
    """
    snp_list = []
    for seqid, rank, chrm, pos1, pos2 in zip(seqids, ranks, chrom, position1, position2):
        for row in vcf_in.fetch(chrm, pos1, pos2):
            snp_list.append((seqid, rank, str(row).strip('\n').split('\t')))
    select_list = []
    for item in snp_list:
        select_list.append((item[0], item[1], item[2][0], item[2][1], item[2][2], item[2][7].split(';')))
    final_list = []
    for item in select_list:
        geneinfo = [j for j in item[5] if "GENEINFO=" in j]
        caf = [j for j in item[5] if "CAF=" in j]
        topmed = [j for j in item[5] if "TOPMED=" in j]
        final_list.append((item[0], item[1], item[2], item[3], item[4], geneinfo, caf, topmed))
    return final_list

def tabix_results_to_df(tabix_list, which_primer, column_name):
    """
    Takes the results from the tabix search and creates a dataframe.
    Args:
        tabix_list (list): tabix info from primer_left
        which_primer (str): should be L or R to denote left or right in names
        column_name (str): column name of snp count (should specify left or right)
    Returns:
        tabix_frame (pd.DataFrame): list organized into dataframe
    """
    tabix_frame = pd.DataFrame(tabix_list)
    tabix_frame.columns = ["Sequence ID", "Primer Rank", "Chromosome", "SNPPosition", "rs_id",
                           "GeneInfo", "CommonAlleleFreq", "TopMedFreq"]

    tabix_frame["GeneInfo"] = tabix_frame["GeneInfo"].apply(lambda x: "NA" if len(x)==0 else x[0].split("=")[1])
    tabix_frame["CommonAlleleFreq"] = tabix_frame["CommonAlleleFreq"].apply(lambda x: "NA" if len(x)==0 else x[0].split('=')[1].split(',')[1])
    tabix_frame["TopMedFreq"] = tabix_frame["TopMedFreq"].apply(lambda x: "NA" if len(x)==0 else x[0].split('=')[1].split(',')[1])
    tabix_frame = tabix_frame.groupby(["Sequence ID", "Primer Rank", "Chromosome"]).agg({'SNPPosition': ';'.join,
                                                                                         'rs_id': ';'.join,
                                                                                         'GeneInfo': 'first',
                                                                                         'CommonAlleleFreq': ';'.join,
                                                                                         'TopMedFreq': ';'.join}).reset_index()
    tabix_frame[column_name] = tabix_frame["rs_id"].apply(lambda x: len(x.split(';')))
    tabix_frame.columns = ["Sequence ID", "Primer Rank", "Chromosome", "%s_SNPPosition" %which_primer,
                           "%s_rs_id" %which_primer, "%s_GeneInfo" %which_primer,
                           "%s_CommonAlleleFreq" %which_primer, "%s_TopMedFreq" %which_primer, column_name] 
    return tabix_frame

def merge_left_right(left_df, right_df, total):
    """
    Merge left and right primer tabix info dataframes.
    Args:
        left_df (pd.DataFrame): primer left tabix dataframe
        right_df (pd.DataFrame): primer right tabix dataframe
    returns:
        merged_tabix_df (pd.DataFrame): left and right merged df
    """
    merged_tabix_df = pd.merge(total, left_df, on=['Sequence ID', 'Primer Rank'], how='left')\
        .merge(right_df, on=['Sequence ID', 'Primer Rank'], how='left')
    return merged_tabix_df

def main():
    """
    One func to rule them all!
    """
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
    args = get_args_snps()
    vcf_in = VariantFile(args.vcf)
    p_info = create_df(args.p_info)
    p_left = primer_range_left(p_info["Sequence ID"],
                               p_info["Primer Rank"],
                               p_info["Chromosome"],
                               p_info["Primer Left Seq"],
                               p_info["Position1"])
    p_right = primer_range_right(p_info["Sequence ID"],
                                 p_info["Primer Rank"],
                                 p_info["Chromosome"],
                                 p_info["Primer Right Seq"],
                                 p_info["Position2"])
    # for testing
    p_left.to_csv('p_left.csv', index=False)
    pn_left = match_pinfo_to_vcf(p_left, vcf_in)
    pn_right = match_pinfo_to_vcf(p_right, vcf_in)
    left_snps = tabix_fetch(pn_left["Sequence ID"],
                            pn_left["Primer Rank"],
                            pn_left["Chromosome"],
                            pn_left["Position1"],
                            pn_left["Position2"],
                            vcf_in)
    right_snps = tabix_fetch(pn_right["Sequence ID"],
                             pn_right["Primer Rank"],
                             pn_right["Chromosome"],
                             pn_right["Position1"],
                             pn_right["Position2"],
                             vcf_in)
    left_df = tabix_results_to_df(left_snps, "L", "Left SNP Count")
    right_df = tabix_results_to_df(right_snps, "R", "Right SNP Count")
    merged_df = merge_left_right(left_df, right_df, p_info)
    merged_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
