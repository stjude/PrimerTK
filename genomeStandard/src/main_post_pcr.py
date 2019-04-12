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
    seqs, headers = ap.fasta_parser(args.pcrfile)
    # 3) Calculate GC of each PCR product and store in list
    gc_list = ap.gc_percent_seqs(seqs)
    # 4) Split up the header line and get no_chrom_list
    split_headers, no_chrom = ap.split_headers_list(headers)
    # 5) Generate chrom list
    chrom_list = ap.chr_split_list(split_headers)
    # 6) Get split positions list
    pos_split = ap.pos_split_list(chrom_list)
    # 7) Need to delete positions after extracting positions (ugly)
    for item in chrom_list:
        del item[-1]
    # 7) Get name and pos split list
    name_pos_list = ap.split_name_pos(no_chrom)
    # 8) Merge all these lists for dataframe
    merged_list = ap.merge_info(chrom_list, pos_split, name_pos_list, no_chrom)
    all_pcr_df, good_primers_df, bad_primers_df = ap.generate_pcr_df(merged_list, gc_list)
    # 9) Output file generation
    all_pcr_df.to_csv('pcr_product_info.csv', index=False)
    # 10) Merge good primers df with toal primers df
    merged_df = ap.merge_good_total(good_primers_df, args.total_primers)
    # 11) Keep only primers which match bw good and total primers
    filtered_df = ap.filter_merged(merged_df)
    filtered_df.to_csv('all_final_primers.csv', index=False)
    # 12) Output only top ranked final primers after filter
    top_ranked_df = ap.top_ranked_final_primers(filtered_df)
    top_ranked_df.to_csv('top_final_primers.csv', index=False)
    # 13) generate easy order plate (only for standard PCR atm)
    plate_forward_primers, plate_reverse_primers = ap.to_order_plate(top_ranked_df)
    plate_forward_primers.to_csv('plate_forward_primers.csv', index=False)
    plate_reverse_primers.to_csv('plate_reverse_primers.csv', index=False)

if __name__ == "__main__":
    main()
