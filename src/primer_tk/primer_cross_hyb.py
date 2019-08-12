#!/usr/bin/env python3
"""
Modules required for program.
 - python3.6+
    - pandas>=0.22.0
    - argparse
    - sequence_info
"""

import argparse
from primer_tk import sequence_info as seqinf

def add_pre_subparser(subparser):
    """ Add subparser for 'pre' step.

    Args:
        subparser (?): Subparser object.

    Returns: None
    """
    parser = subparser.add_parser("pre", help="Preprocessing for snv/indel",
                                  description="Command Line argument for total primer"
                                              "input file to check if primers have a degree"
                                              "of complementarity with each other as defined"
                                              "by the user. Default is 60% (fairly strict).")

    parser.add_argument("-d", "--primer3_dump", dest="dump", required=True,
                        help="Primer3 stdout passed into a 'dump' file to be used as input")

    parser.add_argument("-o", "--outfile_name", dest="outfile", required=True,
                        help="The output filename for all primer information.")
    parser.add_argument("-nd", "--no_dimer", default="no_dimer_df.csv",
                        help="The primers left after dimers removed.")
    parser.add_argument("-spcr", "--standard_pcr_file", default="standard_pcr.txt",
                        help="The file to be used for standard pcr input")
    parser.add_argument("-mpcr", "--multiplex_pcr_file", default="multiplex_pcr.txt",
                        help="The file to be used for multiplex pcr input")
    parser.add_argument("-pa", "--percent_alignment", dest="percent_alignment",
                        default="60", help="Percent match between 2 primers for pair to be\
                                            discarded. EX: primer_len = 22, percent_aln = 60\
                                            dimer_len = (60/100) * 22 = 13.2 -> 13.")
    parser.add_argument("-pcr", "--pcr_type", dest="pcr", required=True,
                        choices=['standard', 'multiplex'],
                        help="perform standard or multiplex pcr on given inputs.")

def add_pre_sv_subparser(subparser):
    """ Add subparser for 'pre' step.

    Args:
        subparser (?): Subparser object.

    Returns: None
    """
    parser = subparser.add_parser("pre_sv", help="Preprocessing for SV's")

    parser.add_argument("-d", "--primer3_dump", dest="dump", required=True,
                        help="Primer3 stdout passed into a 'dump' file to be used as input")

    parser.add_argument("-o", "--outfile_name", dest="outfile", required=True,
                        help="The output filename for all primer information.")
    parser.add_argument("-pcr", "--pcrfile", default="standard_pcr.txt",
                        help="The pseudopcr file")

def get_fprimer_percent_aln(fprimer, percent_alignment):
    """
    Gets the len of the fprimer and calculates minimum percent alignment
    based on user input of maximum alignment between 2 primers allowed.
    See percent alignment help for info.
    Args:
        fprimer (str): the forward primer
    Returns:
        fp_len (int): user defined total number of bases to allowed to match.
    """
    fp_len = []
    for fseq in fprimer:
        f_len = len(fseq)
        min_dimer_alignment = int(f_len*(percent_alignment/100))
        fp_len.append(min_dimer_alignment)
    return fp_len

def primer_dimer_local(aln_len_list, names, seq1, seq2):
    """
    Looks for complementarity between all primers in list,
    and returns alignment + score. No return, prints visualization
    between aligned primers (forming dimer) if they meet the min score set in the func.
    Args:
        min_dimer_alignment (int): the minimum alignment score between 2 primers before
            it is thrown out
        names (string): sequence ID
        seq1 (string): primer1 sequence
        seq2 (string): primer2 sequence
    Returns:
        name, seq1, seq2, score (generator): yields any primer pair that has score
        greater than or equal to min_dimer_alignment.
    """
    for name, fseq, aln_len in zip(names, seq1, aln_len_list):
        for rseq in seq2:
            pd_test = seqinf.PrimerDimer(fseq, seqinf.Sequence(rseq).complement(), aln_len)
            pd_output = pd_test.pd_local()
            for item in pd_output:
                yield(name+"\n", seqinf.PrimerDimer.format_alignment_compl(item[0],
                                                                           item[1],
                                                                           item[2],
                                                                           item[3],
                                                                           item[4]))

def list_from_gen(gen_expression):
    """
    Creates a list from the generator function for dimer alignment.
    Items in each_item are name, seq1, seq2, score.
    Args:
        gen_expression (generator object): the output from primer_dimer_local
    Returns:
        primer_compare_list (list): generaor converted to list
    """
    primer_compare_list = []
    for each_item in gen_expression:
        primer_compare_list.append(str(each_item))
    return primer_compare_list

def p_list_formatter(primer_list):
    """
    Reformat the primer list (remove unnecessary characters from biopython2 output).
    Args:
        primer_list (list): list from list_from_gen output
    Returns:
        primer_dimers (list): list with unnecessary chars removed.
    """
    reformat_p_list = []
    primer_dimers = []
    reformat_p_list = [each_item.replace('\\n', ' ').split() for each_item in primer_list]
    for each_item in reformat_p_list:
        primer_dimers.append((each_item[0].replace('(', '').replace('\'', ''),
                              each_item[2].replace('\'', ''), each_item[4],
                              each_item[5].replace('\')', '')))
    return primer_dimers

def dimer_true(dataframe, col_num, dimer_list):
    """
    Boolean masks to let us know which primers from original df
    form dimers. If so they are dropped.
    Args:
        dataframe (pd.DataFrame): the primer dataframe
        col_num (int): the column number to check for primer match
        dimer_list (list): the list containing primer dimer info
    Returns:
        out_series (pd.Series): boolean masked series, True if primer is dimer, else False
    """
    out_series = dataframe.iloc[:, col_num].isin([seq[1] or seq[2] for seq in dimer_list])
    return out_series

def all_vs_all_pcr(df_boolean, outfile):
    """
    Creates all vs all pcr input to check for off target PCR amplification.
    This function assumes you used the output created from primer_cross_hyb.py
    Args:
        df_boolean (pd.DataFrame): dataframe with removed cross hybridizing primers
    Returns:
        nothing: writes a file output
    """
    no_dimer_df = df_boolean
    all_vs_all = open(outfile, 'w')
    for seqid, primer_left in zip(no_dimer_df['Sequence ID'], no_dimer_df['Primer Left Seq']):
        for primer_right in no_dimer_df['Primer Right Seq']:
            all_vs_all.write(str(seqid) + '\t' + str(primer_left) + '\t' + str(primer_right) + '\n')

    for seqid, primer_right in zip(no_dimer_df['Sequence ID'], no_dimer_df['Primer Right Seq']):
        for primer_left in no_dimer_df['Primer Left Seq']:
            all_vs_all.write(str(seqid) + '\t' + str(primer_right) + '\t' + str(primer_left) + '\n')

    for seqid, primer_left in zip(no_dimer_df['Sequence ID'], no_dimer_df['Primer Left Seq']):
        for primer_left2 in no_dimer_df['Primer Left Seq']:
            all_vs_all.write(str(seqid) + '\t' + str(primer_left) + '\t' + str(primer_left2) + '\n')

    for seqid, primer_right in zip(no_dimer_df['Sequence ID'], no_dimer_df['Primer Right Seq']):
        for primer_right2 in no_dimer_df['Primer Right Seq']:
            all_vs_all.write(str(seqid) + '\t' + str(primer_right) +\
                                 '\t' + str(primer_right2) + '\n')

    all_vs_all.close()

def standard_pcr(primer_df, outfile):
    """
    Creates standard PCR input file (one primer with one other) to run in silico PCR.
    Args:
        primer_df (pd.DataFrame): dataframe with primer info
    Returns:
        nothing: writes a file output
    """
    standard = open(outfile, 'w')
    for seqid, primer_left, primer_right in zip(primer_df['Sequence ID'],
                                                primer_df['Primer Left Seq'],
                                                primer_df['Primer Right Seq']):
        standard.write(str(seqid) + '\t' + str(primer_left) + '\t' + str(primer_right) + '\n')
    standard.close()
