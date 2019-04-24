#!/usr/bin/env python3

"""
Dependencies required to run program:
 - python3.6+
    - pandas>=0.22.0
    - argparse
    - sequence_info
A subset of functions used after PCR.
"""

import sys
import os
import argparse
import pandas as pd
from primer_tk import sequence_info as seqinf

def add_post_subparser(subparser):
    """ Add subparser for postprocessing step.

    Args:
        subparser (?): Subparser object.

    Returns: None
    """
    parser = subparser.add_parser("post", help='Parses info from pcr_output')
    parser.add_argument("-i", "--pcr_output",
                        dest="pcrfile", default="pcr_output.txt",
                        help="use output of isPCR")

    parser.add_argument("-tp", "--total_primers",
                        help="the pre-PCR master primer file that\
                              contains all sample + primer info.")


def fasta_parser(pcrfile):
    """
    Iterates through sudo-fasta output from in silico pcr
    and extracts all useful info related to samples.
    Args:
        pcrfile (file): file output from in silico pcr of primers.
    Returns:
        seqs (list): list of seqs from pcr output
        headers (list): list of headers from pcr output
    """
    seqs = []
    headers = []
    if not os.path.exists(pcrfile):
        sys.exit("ERROR: pcr_file '%s' DNE"%(pcrfile))
    with open(pcrfile) as pcr_file:
        sequence = ""
        header = None
        for line in pcr_file:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append([sequence])
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append([sequence])
    return seqs, headers

def gc_percent_seqs(seqs_list):
    """
    Use seqinf to calculate GC content for output from isPCR.
    Args:
        seqs_list (list): list of seqs output from isPCR
    Returns:
        gc_list (list): list of gc content of seqs.
    """
    upper_seqs = [[nuc.upper() for nuc in seq] for seq in seqs_list]
    gc_list = []
    for seqs in upper_seqs:
        for seq in seqs:
            gc_list.append(seqinf.Sequence(seq).gc_percent())
    return gc_list

# Below is a plethora of list modification functions. perhaps this can be improved.
def split_headers_list(headers):
    """
    Split the headers list into unique elements instead of 1 string.
    Args:
        headers (list): headers output from fasta_parser
    Returns:
        split_headers (list): headers list split into components.
        no_chrom (list): all things that are not the chromosome.
    """
    split_headers = [item.split(" ") for item in headers]
    no_chrom = [item[1:] for item in split_headers]
    return split_headers, no_chrom

def chr_split_list(split_headers):
    """
    Gets chromosome info from split_headers.
    Args:
        split_headers (list): header list split into parts
    Returns:
        chr_list (list): list of chromosome values
    """
    chr_list = []
    for item in split_headers:
        if item[0]:
            new_item = item[0].split(":")
            chr_list.append(new_item)
    return chr_list

def pos_split_list(chr_list):
    """
    Gets position info from chr_list.
    Args:
        chr_list (list): list of chr info with position
    Returns:
        pos_list (list): list with position info
    """
    pos_list = []
    for pos in chr_list:
        if "+" in pos[1]:
            for_pos = pos[1].split("+")
            pos_list.append(for_pos)
        if "-" in pos[1]:
            rev_pos = pos[1].split("-")
            pos_list.append(rev_pos)
    return pos_list


def split_name_pos(no_chrom_list):
    """
    Splits the name and position in no_chrom_list.
    Args:
        no_chrom_list (list): list with header info not containing chromosome
    Returns:
        name_pos_split_list: split name and position list
    """
    name_pos_split_list = []
    for item in no_chrom_list:
        new_split = item[0].replace("_", "").split(":")
        name_pos_split_list.append([new_split[1]])
    return name_pos_split_list

def merge_info(chr_list, pos_split, name_pos_split_list, no_chrom_list):
    """
    Merges all lists created above to be used for DataFrame.
    Args:
        chr_list (list): the list with the chromosome info of product
        pos_split_list (list): the list with the split positions
        name_pos_split_list (list): the list with the names and positions split
        no_chrom_list (list): the list without chrom info
    Returns:
        merged_list (list): list with all info to be taken into dataframe
    """
    merged_list = []
    for chrom, pos, name, prims in zip(chr_list, pos_split, name_pos_split_list,
                                       no_chrom_list):
        merged_items = chrom + pos + name + prims
        merged_list.append(merged_items)
    return merged_list

def generate_pcr_df(merged_list, gc_list):
    """
    Generate a DataFrame from all the relevant merged_list information.
    Args:
        merged_list (list): list of all information about pcr products
    Returns:
        pcr_df (DataFrame): dataframe with all PCR product information
        good_primers_df (DataFrame): dataframe with primers that are on target
        bad_primers_df (DataFrame): dataframe with primers that failed
    """
    pcr_df = pd.DataFrame.from_records(merged_list)
    pcr_df.columns = ['Chromosome', 'Position1', 'Position2', 'Sample_Pos',
                      'Sequence ID', 'PCR_Prod_Len', 'Forward', 'Reverse']
    pcr_df['Product GC%'] = gc_list
    pcr_df['ForwardPrimLen'] = [len(primer) for primer in pcr_df['Forward']]
    pcr_df['ReversePrimLen'] = [len(primer) for primer in pcr_df['Reverse']]
    pcr_df['Position1'] = pcr_df['Position1'].astype(int)
    pcr_df['Position2'] = pcr_df['Position2'].astype(int)
    pcr_df['Sample_Pos'] = pcr_df['Sample_Pos'].astype(int)
    # Filter primers based on if they amplified in correct position
    good_primers_df = pcr_df.loc[(pcr_df['Position1'] <= pcr_df['Sample_Pos']) &
                                 (pcr_df['Sample_Pos'] <= pcr_df['Position2'])]
    # These primers failed filtering
    bad_primers_df = pcr_df.loc[~((pcr_df['Position1'] <= pcr_df['Sample_Pos']) &
                                  (pcr_df['Sample_Pos'] <= pcr_df['Position2']))]

    return pcr_df, good_primers_df, bad_primers_df

def merge_good_total(good_primers, total_primers):
    """
    Merges good primers and total primers on Sequence ID.
    Picks highest ranking primers.
    Args:
        good_primers (file): dataframe generated by pcr analysis
        total_primers (file): dataframe generated by Primer3
    Returns:
        merged_df (DataFrame): good_total merged dataframe on Sequence ID
    """
    total_primers_df = pd.read_csv(total_primers)
    merged_df = total_primers_df.merge(good_primers, on='Sequence ID', how='left')
    merged_df['FwdPrimerPos'] = merged_df['Position1'].apply(str) + '-' +\
        (merged_df['Position1'] + merged_df['Primer Left Len']).apply(str)
    merged_df['RvsPrimerPos'] = (merged_df['Position2'] - merged_df['Primer Right Len']).apply(str)\
        + '-' + merged_df['Position2'].apply(str)
    return merged_df

def filter_merged(merged_df):
    """
    Filters the merged df to find only primers which align with initial primer pairs.
    Args:
        merged_df (DataFrame): output df from merge_good_total()
    Returns:
        useful_filtered (DataFrame): the filtered dataframe
    """
    filtered_df = merged_df.loc[(merged_df['Primer Left Seq'] == merged_df['Forward']) &
                                (merged_df['Primer Right Seq'] == merged_df['Reverse'])]
    useful_filtered = filtered_df[['Sequence ID', 'Primer Rank', 'Primer Left Seq',
                                   'Primer Right Seq', 'Primer Left Len', 'Primer Right Len',
                                   'Primer Left TM', 'Primer Right TM', 'Primer Left GC %',
                                   'Primer Right GC %', 'Chromosome', 'Position1',
                                   'Position2', 'PCR_Prod_Len',
                                   'Product GC%', 'FwdPrimerPos', 'RvsPrimerPos']]
    return useful_filtered

def top_ranked_final_primers(filter_merged_df):
    """
    Drops duplicate sequence ids and keeps first (which also corresponds)
    to the highest ranking primer pair for each sample.
    Args:
        filter_merged_df (DataFrame): input from filter_merged, where primers are only equal
                                      to on target primers from initial primer generation.
    Returns:
        top_ranked_df (DataFrame): outputs only the highest scoring primer pair
                                   at each position
    """
    top_ranked_df = filter_merged_df.drop_duplicates('Sequence ID', keep='first')
    return top_ranked_df

def to_order_plate(top_ranking_final_primers):
    """
    Takes output from top_ranked_final_primers and organizes to easy order IDT plates
    (Forward and Reverse).
    Args:
        top_ranked_final_primers (DataFrame): filtered, top ranked dataframe primers
    Returns:
        idt_order_sheet_plate_f (DataFrame): 96 well plate order format forward primers
        idt_order_sheet_plate_r (DataFrame): 96 well plate order format reverse primers
    """
    # Generate well numbers
    well_and_nums = []
    wells = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    num_wells = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    for well in wells:
        for num in num_wells:
            well_and_nums.append(well+num)

    filt_ranked_df = top_ranking_final_primers
    plate_order_sheet_f = pd.DataFrame(columns=['Sequence Name', 'Sequence'])
    plate_order_sheet_r = pd.DataFrame(columns=['Sequence Name', 'Sequence'])
    plate_order_sheet_f['Sequence Name'] = filt_ranked_df['Sequence ID']+'_F'
    plate_order_sheet_f['Sequence'] = filt_ranked_df['Primer Left Seq']
    plate_order_sheet_f.insert(0, 'Well Position', well_and_nums[:len(plate_order_sheet_f['Sequence'])])
    plate_order_sheet_r['Sequence Name'] = filt_ranked_df['Sequence ID']+'_R'
    plate_order_sheet_r['Sequence'] = filt_ranked_df['Primer Right Seq']
    plate_order_sheet_r.insert(0, 'Well Position', well_and_nums[:len(plate_order_sheet_r['Sequence'])])
    return plate_order_sheet_f, plate_order_sheet_r
