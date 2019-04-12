#!/usr/bin/env python3

"""
Dependencies required to run program:
 - python3.6+
    - pandas>=0.22.0
    - argparse
    - sequence_info
"""

import analyze_pcr_output as ap

def main():
    """
    Generates pcr analysis dataframes and applies primer filtering based on
    off-target amplification. Then compares good primers to initial primer
    list to find which primer pair was generated and top ranking.
    Finally, produces easy to use IDT order sheet in plate format (standard PCR only).
    """
    # 1) Get arguments
    args = ap.get_args_pcr_analysis()
    # 2) Generate seqs and headers lists
    seqs, headers = ap.fasta_parser(args.flank_file)
    # 3) Calculate GC of each PCR product and store in list
    positions_to_compare = ap.amp_header_region(args.total_primers)
    sliced_seqs = ap.get_gc_region(seqs, headers, positions_to_compare)
    gc_calc = ap.calc_gc(sliced_seqs)
    merged_df = ap.merge_dfs(gc_calc, args.total_primers, seqs)
    merged_df.to_csv('total_list_gc.csv', index=False)
    merged_df.drop_duplicates('Sequence ID', keep='first', inplace=True)
    merged_df.to_csv('top_ranked_final_primers.csv', index=False)

if __name__ == "__main__":
    main()
