#!/usr/bin/env python3

"""
Modules required for program.
 - python3.6+
    - pandas>=0.22.0
    - argparse
"""
import argparse

def get_args_hyb():
    """
    get commandline arguments.
    Args: None
    Returns:
    args (NameSpace): the commandline arguments
    """

    parser = argparse.ArgumentParser(description='Command Line argument for total primers\
                                                  generation.')
    parser.add_argument("-d", "--primer3_dump",
                        dest="dump",
                        help="Primer3 stdout passed into a 'dump' file to be used as input")

    parser.add_argument("-o", "--outfile_name",
                        dest="outfile",
                        help="The output filename for all primer information.")
    args = parser.parse_args()
    return args


def standard_pcr(primer_df):
    """
    Creates standard PCR input file (one primer with one other) to run in silico PCR.
    Args:
        input_file (file): no_dimer_df.csv generated from previous step in program
    Returns:
        nothing: writes a file output
    """
    standard = open('standard_pcr.txt', 'w')
    for seqid, primer_left, primer_right in zip(primer_df['Sequence ID'],
                                                primer_df['Primer Left Seq'],
                                                primer_df['Primer Right Seq']):
        standard.write(str(seqid) + '\t' + str(primer_left) + '\t' + str(primer_right) + '\n')
    standard.close()
