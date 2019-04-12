#!/usr/bin/env python3

""" Modules required for program.
 - python3.6+
    - pandas>=0.22.0
    - argparse
    - mp_class
    - sequence_info
"""
import sys
from mp_class import MissingPrimers
from mp_class import create_df
import standard_pcr_setup as sps

def main():
    """
    Function for all steps leading up to PCR.
    """
    args = sps.get_args_hyb()
    # 1) Initialize primer lists by rank for each sample
    prim_list_0 = MissingPrimers(args.dump, 0).samp_primer_info
    prim_list_1 = MissingPrimers(args.dump, 1).samp_primer_info
    prim_list_2 = MissingPrimers(args.dump, 2).samp_primer_info
    prim_list_3 = MissingPrimers(args.dump, 3).samp_primer_info
    prim_list_4 = MissingPrimers(args.dump, 4).samp_primer_info
    # 2) Generate the output df
    primer_df = create_df([prim_list_0, prim_list_1, prim_list_2,
                           prim_list_3, prim_list_4])
    # 3) Generate csv output
    primer_df = primer_df.loc[~(primer_df['Primer Left Seq'] == 'NA')]
    primer_df.to_csv(args.outfile, index=False)
    # 4) create standard pcr input
    sps.standard_pcr(primer_df)

if __name__ == "__main__":
    main()
