#!/usr/bin/env python

"""
Dependencies required to run program:
 - python3.6+
    - pandas>=0.22.0
    - argparse
    - sequence_info
"""

import sys
import argparse
import pandas as pd
from primer_tk import sequence_info as seqinf

def add_post_subparser(subparser):
    """ Add subparser for postprocessing sv step.

    Args:
        subparser (?): Subparser object.

    Returns: None
    """
    parser = subparser.add_parser("post_sv", help='Parses info from pcr_output')
    parser.add_argument("-f", "--flank_file", dest="flank_file",
                        help="use flanking_regions file from output of genome_iterator_sv.py")
    parser.add_argument("-tp", "--total_primers", dest="total_primers",
                        help="the pre-PCR master primer file that\
                              contains all sample + primer info")
    parser.add_argument("-all", "--all_final_primers", default="all_final_primers_sv.csv",
                        help="all primers generated for targets")
    parser.add_argument("-top", "--top_final_primers", default="top_final_primers_sv.csv",
                        help="top primers generated for targets")
    parser.add_argument("-plate", "--plate_basename", default="plated_primers",
                        help="the basename of the primers in plate format ready to order.")

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
    with open(pcrfile) as pcr_file:
        sequence = ""
        header = None
        for line in pcr_file:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append(sequence)
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append(sequence)
    return seqs, headers

def amp_header_region(infile):
    """
    Extract relevant info from total_primer_list.csv to use for GC.
    Args:
        infile (file): total_list file to parse primer info
    Returns:
        header_start_stop (list): list of information to use when parsing fasta
    """
    primer_df = pd.read_csv(infile)
    header_start_stop = [[header, start, stop] for header, start, stop\
                             in zip(primer_df['Sequence ID'],
                                    primer_df['Primer Left Pos'],
                                    primer_df['Primer Right Pos'])]
    return header_start_stop

def get_gc_region(seqs, headers, positions):
    """
    Slice out the exact region of the fasta the primer amplified.
    Args:
        seqs (list): the list of sequences from flank_file
        headers (list): the list of headers from flank_file
        positions (list): list of tuples with header, p_len, p_seq
    Returns:
        sliced_seq (list): the exact amplified sequence.
    """
    sliced_seq = []
    for header, seq in zip(headers, seqs):
        for pos in positions:
            if pos[0] == header:
                sliced_seq.append((header, len(seq[pos[1]:pos[2]+1]), seq[pos[1]:pos[2]+1]))
    return sliced_seq

def calc_gc(sliced_seq):
    """
    Takes sliced seq list with header, product len, and product seq and returns
    the header, product len, and product GC%. A bit redundant but I might find
    it useful to have this exact sequence in the future.
    Args:
        sliced_seq (list): contains (header, product_len, product_seq)
    Returns:
        sample_gc (list): contains (header, product_len, product GC%)
    """
    sample_gc = []
    for header, p_len, p_seq in sliced_seq:
        sample_gc.append((header, p_len, seqinf.Sequence(p_seq).gc_percent()))
    return sample_gc

def merge_dfs(gc_calc, total_df, seqs):
    """
    Merges the product information onto the total primers sheet.
    Args:
        gc_calc (list): list containing header, product len, gc info
        total_df (pandas object): dataframe containing primer information
    Returns:
        merged_df (pandas object): dataframe containing updated primer information
    """

    gc = pd.DataFrame(gc_calc, columns=['Sequence Name', 'Product Len', 'Product GC%'])
    total = pd.read_csv(total_df)
    merged = pd.concat([total, gc], axis=1) # side-by-side merge
    merged['Chromosome'] = [seqid.split(':')[0].split('_')[2] for seqid in merged['Sequence ID']]
    merged['Position1'] = [seqid.split(':')[1].strip('_').split('-')[0] for seqid in merged['Sequence ID']]
    merged['Position2'] = [seqid.split(':')[1].strip('_').strip('__BP1').strip('__BP2').\
                           split('-')[1] for seqid in merged['Sequence ID']]
    seqlen = []
    for seq in seqs:
        seqlen.append(len(seq[:-1]))
    flanksize = (seqlen[0]/2)
    merged['FlankSize'] = int(flanksize)
    merged['NucsBeforeBP'] = (merged['Position1']).apply(int)\
        - ((merged['Position1'].apply(int) - merged['FlankSize'])\
        + (merged['Primer Left Pos']).apply(int))
    merged['NucsAfterBP'] = merged['Product Len'].apply(int) -\
        (merged['FlankSize'] - merged['Primer Left Pos'].apply(int))
    merged['Position1'] = ((merged['Position1'].apply(int) - merged['FlankSize'])\
                               + merged['Primer Left Pos'].apply(int))
    merged['Position2'] = (merged['Position2'].apply(int) + merged['NucsAfterBP'].apply(int))
    merged['FwdPrimPos'] = (merged['Position1']).apply(str) + '-'\
        + (merged['Position1'] + merged['Primer Left Len']).apply(str)
    merged['RvsPrimPos'] = (merged['Position2'] - merged['Primer Right Len']).apply(str)\
        + '-' + (merged['Position2']).apply(str) 
    merged_df = merged[['Sequence ID', 'Primer Rank', 'Primer Left Seq', 'Primer Right Seq',
                        'Primer Left Len', 'Primer Right Len', 'Primer Left TM',
                        'Primer Right TM', 'Primer Left GC %', 'Primer Right GC %',
                        'Product Len', 'Product GC%', 'Chromosome', 'Position1', 'Position2',
                        'NucsBeforeBP', 'NucsAfterBP', 'FwdPrimPos', 'RvsPrimPos']]
    return merged_df

def to_order_plate(top_ranking_final_primers):
    """
    Takes output from top_ranked_final_primers and organizes to easy order plates
    (Forward and Reverse).
    Args:
        top_ranked_final_primers (DataFrame): filtered, top ranked dataframe primers
    Returns:
        idt_order_sheet_plate_f (DataFrame): 96 well plate order format forward primers
        idt_order_sheet_plate_r (DataFrame): 96 well plate order format reverse primers
    """
    # Generate well numbers
    if len(top_ranking_final_primers) == 0:
        sys.exit("No primers were kept for any target. Try starting over with relaxed parameters.")
    well_and_nums = []
    wells = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    num_wells = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    for well in wells:
        for num in num_wells:
            well_and_nums.append(well+num)

    # this is required because sometimes more than 96 primers are generated and the list is consumed.
    big_well_nums = well_and_nums * 50
    top_ranking_final_primers.set_index('Sequence ID', inplace=True)
    plate_order_sheet_f = pd.DataFrame(columns=['Sequence Name', 'Sequence'])
    plate_order_sheet_r = pd.DataFrame(columns=['Sequence Name', 'Sequence'])
    names = top_ranking_final_primers.index.tolist()
    f_seq = top_ranking_final_primers['Primer Left Seq'].tolist()
    r_seq = top_ranking_final_primers['Primer Right Seq'].tolist()
    plate_order_sheet_f['Sequence Name'] = names
    plate_order_sheet_f['Sequence Name'] = plate_order_sheet_f['Sequence Name'] + '_F'
    plate_order_sheet_f['Sequence'] = f_seq
    plate_order_sheet_f.insert(0, 'Well Position', big_well_nums[:len(plate_order_sheet_f['Sequence'])])
    plate_order_sheet_r['Sequence Name'] = names
    plate_order_sheet_r['Sequence Name'] = plate_order_sheet_r['Sequence Name'] + '_R'
    plate_order_sheet_r['Sequence'] = r_seq
    plate_order_sheet_r.insert(0, 'Well Position', big_well_nums[:len(plate_order_sheet_r['Sequence'])])
    return plate_order_sheet_f, plate_order_sheet_r

