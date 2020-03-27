#!/usr/bin/env python3
""" Core modules """

import os
from pysam import VariantFile
from primer_tk import genome_iterator as gi
from primer_tk import analyze_pcr_output as ap
from primer_tk import analyze_pcr_output_sv as apv
from primer_tk.mp_class import MissingPrimers as mp
from primer_tk.mp_class import create_df as cd
from primer_tk.mp_class_sv import MissingPrimers as mpv
from primer_tk.mp_class_sv import create_df as cdv
from primer_tk import utils
import primer_tk.primer_cross_hyb as pch
from primer_tk import primer_tabix as pt
from primer_tk import genome_iterator_sv as giv


def iterator_sv(args):
    """ Use an input regions file with SV positions to pull
        down flanking sequence on both sides of SV to generate
        primers upstream and downstream of the SV.

    Args:
        args (Namespace): Argparse results.

    Returns: None
    """

    dataset_name = os.path.splitext(str(args.regions_file))[0]
    genome = giv.genome_iterator(args.ref_genome)
    # 2) create dataframe from input regions file
    if args.sv in ('deletion', 'inversion'):
        small_regions = giv.file_extension(args.regions_file, args.sv)
    # 3) ensure proper proper number of columns in dataframe
        assert len(list(small_regions)) == 5, "DataFrame contains more/less than 5 columns...\
                                           Improper format."
    # 4) format dataframe "chr" column to match reference genome
        small_regions = giv.match_chr_to_genome(small_regions, genome, args.sv)
    # 5) generate flanking regions fasta based on position in input file
        flanking = open("flanking_regions.%s.fasta" % dataset_name, 'w')
        if args.sv == 'deletion':
            flank_data = giv.flanking_regions_fasta_deletion(genome, small_regions, args.flanking_region_size)
            primer3_in = open("primer3_input.%s.txt" % dataset_name, 'w')
            for head, seq in flank_data:
                flanking.write(">"+head+'\n'+seq+'\n')
            # 6) generate primer3 input file
                primer3_in.write(utils.primer3_input(head, seq, args))
        elif args.sv == 'inversion':
            flank_data = giv.flanking_regions_fasta_inversion(genome, small_regions, args.flanking_region_size)
            primer3_in = open("primer3_input.%s.txt" % dataset_name, 'w')
            for head, seq in flank_data:
                flanking.write(">"+head+'\n'+seq+'\n')
            # 6) generate primer3 input file
                primer3_in.write(utils.primer3_input(head, seq, args))
        flanking.close()
        primer3_in.close()

    elif args.sv == 'insertion':
        small_regions = giv.file_extension(args.regions_file, args.sv)
        assert len(list(small_regions)) == 10, "DataFrame contains more/less than 10 columns... Exiting."
        small_regions = giv.match_chr_to_genome(small_regions, genome, args.sv)
        flanking = open("flanking_regions.%s.fasta" %dataset_name, 'w')
        flank_data = giv.flanking_region_fasta_insertion(genome, small_regions, args.flanking_region_size)
        primer3_in = open("primer3_input.%s.txt" % dataset_name, 'w')
        for head, seq in flank_data.items():
            flanking.write(">"+head+'\n'+seq+'\n')
            primer3_in.write(utils.primer3_input(head, seq, args))
    elif args.sv == 'translocation':
        small_regions = giv.file_extension(args.regions_file, args.sv)
        assert len(list(small_regions)) == 8, "DataFrame contains more/less than 8 columns... Exiting."
        small_regions = giv.match_chr_to_genome(small_regions, genome, args.sv)
        flanking = open("flanking_regions.%s.fasta" %dataset_name, 'w')
        flank_data = giv.flanking_region_fasta_translocation(genome, small_regions, args.flanking_region_size)
        primer3_in = open("primer3_input.%s.txt" %dataset_name, 'w')
        for head, seq in flank_data.items():
            flanking.write(">"+head+'\n'+seq+'\n')
            primer3_in.write(utils.primer3_input(head, seq, args))
    flanking.close()
    primer3_in.close()

def iterator(args):
    """
    Use an input regions file with specific region of interest\
    to design primers around, then run primer3.

    Args:
        args (Namespace): Argparse object or None.

    Returns: None
    """

    dataset_name = os.path.splitext(str(args.regions_file))[0]
    # 1) create genome tuple from provided reference
    genome = gi.genome_iterator(args.ref_genome)
    # 2) create dataframe from input regions file
    small_regions = gi.file_extension(args.regions_file)
    # 3) ensure proper proper number of columns in dataframe
    assert len(list(small_regions)) == 4, "DataFrame contains more than 4 columns...\
                                           Improper format."
    # 4) format dataframe "chr" column to match reference genome
    small_regions = gi.match_chr_to_genome(small_regions, genome)
    # 5) generate flanking regions fasta based on position in input file
    flanking = open("flanking_regions.%s.fasta" % dataset_name, 'w')
    flank_data = gi.create_flanking_regions_fasta(genome, small_regions, args.flanking_region_size)
    primer3_in = open("primer3_input.%s.txt" % dataset_name, 'w')
    for head, seq in flank_data:
        flanking.write(">"+head+'\n'+seq+'\n')
    # 6) generate primer3 input file
        primer3_in.write(utils.primer3_input(head, seq, args))
    flanking.close()
    primer3_in.close()

def pre(args):
    """ Function for all steps leading up to PCR. """
    # 1) Initialize primer lists by rank for each sample
    prim_list_0 = mp(args.dump, 0).samp_primer_info
    prim_list_1 = mp(args.dump, 1).samp_primer_info
    prim_list_2 = mp(args.dump, 2).samp_primer_info
    prim_list_3 = mp(args.dump, 3).samp_primer_info
    prim_list_4 = mp(args.dump, 4).samp_primer_info
    # 2) Generate the output df
    primer_df = cd([prim_list_0, prim_list_1, prim_list_2,
                    prim_list_3, prim_list_4])
    # 3) Generate csv output
    primer_df = primer_df.loc[~(primer_df['Primer Left Seq'] == 'NA')]
    primer_df.to_csv(args.outfile, index=False)
    primer_df_standard = primer_df.copy()
    # 5) Get length of forward primers for percent alignment check
    if args.pcr == 'multiplex':
        fp_len = pch.get_fprimer_percent_aln(primer_df['Primer Left Seq'], int(args.percent_alignment))
    # 6) Generate primer dimer pairs for all vs all input
        primer1_2_compare = pch.list_from_gen(pch.primer_dimer_local(fp_len, primer_df['Sequence ID'],
                                                                     primer_df['Primer Left Seq'],
                                                                     primer_df['Primer Right Seq']))
        primer2_1_compare = pch.list_from_gen(pch.primer_dimer_local(fp_len, primer_df['Sequence ID'],
                                                                     primer_df['Primer Right Seq'],
                                                                     primer_df['Primer Left Seq']))
        primer1_1_compare = pch.list_from_gen(pch.primer_dimer_local(fp_len, primer_df['Sequence ID'],
                                                                     primer_df['Primer Left Seq'],
                                                                     primer_df['Primer Left Seq']))
        primer2_2_compare = pch.list_from_gen(pch.primer_dimer_local(fp_len, primer_df['Sequence ID'],
                                                                     primer_df['Primer Right Seq'],
                                                                     primer_df['Primer Right Seq']))
    # 7) Reformat output, Biopython + Class writes ugly output.
        primer1_2_dimers = pch.p_list_formatter(primer1_2_compare)
        primer2_1_dimers = pch.p_list_formatter(primer2_1_compare)
        primer1_1_dimers = pch.p_list_formatter(primer1_1_compare)
        primer2_2_dimers = pch.p_list_formatter(primer2_2_compare)
    # 8) Write output of all primer pairs that form dimers.
        pd_file = open('Primer_Dimers.txt', 'w')
        pd_file.write('#This is a list of possible primer dimers based on the input complementarity.\n')
        pd_file.write('#Comparison of forward primers with all reverse primers...\
                       \n#Sequence ID \t \t Complementarity Score\n')
        for seq in primer1_2_dimers:
            pd_file.write(str(seq) + '\n')
            pd_file.write('#Comparison of reverse primers with all forward primers...\n')
        for seq in primer2_1_dimers:
            pd_file.write(str(seq) + '\n')
            pd_file.write('#Comparison of forward primers with all other forward primers...\n')
        for seq in primer1_1_dimers:
            pd_file.write(str(seq) + '\n')
            pd_file.write('#Comparison of revers primers with all other reverse primers...\n')
        for seq in primer2_2_dimers:
            pd_file.write(str(seq) + '\n')
        pd_file.close()
    # 9) Searche dataframe for dimers in list. If present, marked with True boolean in new column.
    # Else, False. Used for filtering in next step.
        primer_df['Dimers1_2F'] = pch.dimer_true(primer_df, 2, primer1_2_dimers)
        primer_df['Dimers2_1F'] = pch.dimer_true(primer_df, 2, primer2_1_dimers)
        primer_df['Dimers1_1F'] = pch.dimer_true(primer_df, 2, primer1_1_dimers)
        primer_df['Dimers2_2F'] = pch.dimer_true(primer_df, 2, primer2_2_dimers)
        primer_df['Dimers1_2R'] = pch.dimer_true(primer_df, 3, primer1_2_dimers)
        primer_df['Dimers2_1R'] = pch.dimer_true(primer_df, 3, primer2_1_dimers)
        primer_df['Dimers1_1R'] = pch.dimer_true(primer_df, 3, primer1_1_dimers)
        primer_df['Dimers2_2R'] = pch.dimer_true(primer_df, 3, primer2_2_dimers)
    # 10) Search for true statements, True statements indicate dimer duo. Only keep non-dimers
        df_bool = (primer_df.loc[~(primer_df['Dimers1_2F'] | primer_df['Dimers2_1F']\
                                   | primer_df['Dimers1_1F'] | primer_df['Dimers2_2F']\
                                   | primer_df['Dimers1_2R'] | primer_df['Dimers2_1R']\
                                   | primer_df['Dimers1_1R'] | primer_df['Dimers2_2R'] == True)])
    # 11) Write output to csv
        df_bool.to_csv(args.no_dimer, index=False)
        pch.all_vs_all_pcr(df_bool, args.multiplex_pcr_file)
    elif args.pcr == 'standard':
        pch.standard_pcr(primer_df_standard, args.standard_pcr_file)
    else:
        print("Please select pcr setup")

def pre_sv(args):
    """
    Function for all steps leading up to PCR.
    """
    # 1) Initialize primer lists by rank for each sample
    prim_list_0 = mpv(args.dump, 0).samp_primer_info
    prim_list_1 = mpv(args.dump, 1).samp_primer_info
    prim_list_2 = mpv(args.dump, 2).samp_primer_info
    prim_list_3 = mpv(args.dump, 3).samp_primer_info
    prim_list_4 = mpv(args.dump, 4).samp_primer_info
    # 2) Generate the output df
    primer_df = cdv([prim_list_0, prim_list_1, prim_list_2,
                     prim_list_3, prim_list_4])
    # 3) Generate csv output
    primer_df = primer_df.loc[~(primer_df['Primer Left Seq'] == 'NA')]
    primer_df.to_csv(args.outfile, index=False)
    # 4) create standard pcr input
    pch.standard_pcr(primer_df, args.pcrfile)

def post(args):
    """
    Generates pcr analysis dataframes and applies primer filtering based on
    off-target amplification. Then compares good primers to initial primer
    list to find which primer pair was generated and top ranking.
    Finally, produces easy to use IDT order sheet in plate format (standard PCR only).
    """
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
    all_pcr_df.to_csv(args.pcr_product_info, index=False)
    # 10) Merge good primers df with toal primers df
    merged_df = ap.merge_good_total(good_primers_df, args.total_primers)
    # 11) Keep only primers which match bw good and total primers
    filtered_df = ap.filter_merged(merged_df, args.off_target_max)
    filtered_df.to_csv(args.all_primer_info, index=False)
    # 12) Output only top ranked final primers after filter
    top_ranked_df = ap.top_ranked_final_primers(filtered_df)
    top_ranked_df.to_csv(args.top_final_primers, index=False)
    plate_order_f, plate_order_r = ap.to_order_plate(top_ranked_df)
    plate_order_f.to_csv(args.plate_basename+'_F.csv', index=False)
    plate_order_r.to_csv(args.plate_basename+'_R.csv', index=False)

def post_sv(args):
    """
    Generates pcr analysis dataframes and applies primer filtering based on
    off-target amplification. Then compares good primers to initial primer
    list to find which primer pair was generated and top ranking.
    Finally, produces easy to use IDT order sheet in plate format (standard PCR only).
    """
    # 2) Generate seqs and headers lists
    seqs, headers = apv.fasta_parser(args.flank_file)
    # 3) Calculate GC of each PCR product and store in list
    positions_to_compare = apv.amp_header_region(args.total_primers)
    sliced_seqs = apv.get_gc_region(seqs, headers, positions_to_compare)
    gc_calc = apv.calc_gc(sliced_seqs)
    merged_df = apv.merge_dfs(gc_calc, args.total_primers, seqs)
    merged_df.to_csv(args.all_final_primers, index=False)
    merged_df.drop_duplicates('Sequence ID', keep='first', inplace=True)
    merged_df.to_csv(args.top_final_primers, index=False)
    plate_order_f, plate_order_r = apv.to_order_plate(merged_df)
    plate_order_f.to_csv(args.plate_basename+'_F.csv', index=False)
    plate_order_r.to_csv(args.plate_basename+'_R.csv', index=False)

def tabix(args):
    """
    Annotates primers with SNP information.
    """
    vcf_in = VariantFile(args.vcf)
    p_info = pt.create_tabix_df(args.p_info)
    p_left = pt.primer_range_left(p_info["Sequence ID"],
                                  p_info["Primer Rank"],
                                  p_info["Chromosome"],
                                  p_info["Primer Left Seq"],
                                  p_info["Position1"])
    p_right = pt.primer_range_right(p_info["Sequence ID"],
                                    p_info["Primer Rank"],
                                    p_info["Chromosome"],
                                    p_info["Primer Right Seq"],
                                    p_info["Position2"])
    pn_left = pt.match_pinfo_to_vcf(p_left, vcf_in)
    pn_right = pt.match_pinfo_to_vcf(p_right, vcf_in)
    left_snps = pt.tabix_fetch(pn_left["Sequence ID"],
                               pn_left["Primer Rank"],
                               pn_left["Chromosome"],
                               pn_left["Position1"],
                               pn_left["Position2"],
                               vcf_in)
    right_snps = pt.tabix_fetch(pn_right["Sequence ID"],
                                pn_right["Primer Rank"],
                                pn_right["Chromosome"],
                                pn_right["Position1"],
                                pn_right["Position2"],
                                vcf_in)
    left_df = pt.tabix_results_to_df(left_snps, "L", "Left SNP Count")
    right_df = pt.tabix_results_to_df(right_snps, "R", "Right SNP Count")
    merged_df = pt.merge_left_right(left_df, right_df, p_info)
    merged_df.to_csv(args.output, index=False)
