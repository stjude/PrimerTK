#!/usr/bin/env python3
""" Core modules """

from primer_tk import genome_iterator
from primer_tk import analyze_pcr_output as ap
from primer_tk.mp_class import MissingPrimers

def iterator_sv(args):
    """ Use an input regions file with SV positions to pull
        down flanking sequence on both sides of SV to generate
        primers upstream and downstream of the SV.

    Args:
        args (Namespace): Argparse results.

    Returns: None
    """

    timestr = time.strftime("%Y%m%d-%H%M%S")

    args = get_args_sv()
    # 1) create genome tuple from provided reference
    genome = genome_iterator(args.ref_genome)
    # 2) create dataframe from input regions file
    small_regions = file_extension(args.regions_file)
    # 3) ensure proper proper number of columns in dataframe
    assert len(list(small_regions)) == 5, "DataFrame contains more/less than 5 columns...\
                                           Improper format."
    # 4) format dataframe "chr" column to match reference genome
    small_regions = match_chr_to_genome(small_regions, genome)
    # 5) generate flanking regions fasta based on position in input file
    flanking = open("flanking_regions.%s.fasta" % timestr, 'w')
    if args.sv == 'deletion':
        flank_data = flanking_regions_fasta_deletion(genome, small_regions, args.flanking_region_size)
        primer3_in = open("primer3_input.%s.txt" % timestr, 'w')
        for head, seq in flank_data:
            flanking.write(">"+head+'\n'+seq+'\n')
            # 6) generate primer3 input file
            primer3_in.write("SEQUENCE_ID="+head+'\n'
                             +"SEQUENCE_TEMPLATE="+seq+'\n'
                             +"SEQUENCE_TARGET=%s" %args.sequence_target+'\n'
                             +"PRIMER_FIRST_BASE_INDEX=1"+'\n'
                             +"PRIMER_TASK=pick_detection_primers"+'\n'
                             +"PRIMER_MIN_THREE_PRIME_DISTANCE=3"+'\n'
                             +"PRIMER_MAX_LIBRARY_MISPRIMING=12.00"+'\n'
                             +"PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00"+'\n'
                             +"PRIMER_PRODUCT_SIZE_RANGE=%s" %args.product_size_range+'\n'
                             +"PRIMER_MAX_END_STABILITY=9.0"+'\n'
                             +"PRIMER_MAX_SELF_ANY_TH=45.00"+'\n'
                             +"PRIMER_MAX_SELF_END_TH=35.00"+'\n'
                             +"PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00"+'\n'
                             +"PRIMER_PAIR_MAX_COMPL_END_TH=35.00"+'\n'
                             +"PRIMER_MAX_HAIRPIN_TH=24.00"+'\n'
                             +"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00"+'\n'
                             +"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00"+'\n'
                             +"PRIMER_TM_FORMULA=1"+'\n' # use SantaLucia parameters
                             +"PRIMER_SALT_CORRECTIONS=1"+'\n' # SantaLucia 1998 paper
                             +"PRIMER_SALT_MONOVALENT=50.0"+'\n' # mM conc of monovalent salt cations
                             +"PRIMER_INTERNAL_SALT_MONOVALENT=50.0"+'\n' # same as above
                             +"PRIMER_SALT_DIVALENT=1.5"+'\n'
                             +"PRIMER_INTERNAL_SALT_DIVALENT=1.5"+'\n'
                             +"PRIMER_DNTP_CONC=0.6"+'\n' # assume no dntps are present when hybridizing
                             +"PRIMER_INTERNAL_DNTP_CONC=0.6"+'\n'
                             +"PRIMER_DNA_CONC=50.0"+'\n'
                             +"PRIMER_INTERNAL_DNA_CONC=50.0"+'\n'
                             +"PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1"+'\n'
                             +"PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1"+'\n'
                             +"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s" %args.thermopath+'\n'
                             +"PRIMER_PICK_LEFT_PRIMER=1"+'\n'
                             +"PRIMER_PICK_RIGHT_PRIMER=1"+'\n'
                             +"PRIMER_PICK_INTERNAL_OLIGO=1"+'\n'
                             +"PRIMER_MAX_POLY_X=3"+'\n'
                             +"PRIMER_LEFT_NUM_RETURNED=5"+'\n'
                             +"PRIMER_RIGHT_NUM_RETURNED=5"+'\n'
                             +"PRIMER_OPT_SIZE=%s" %args.primer_opt_size+'\n'
                             +"PRIMER_MIN_SIZE=%s" %args.primer_min_size+'\n'
                             +"PRIMER_MAX_SIZE=%s" %args.primer_max_size+'\n'
                             +"PRIMER_MIN_TM=%s" %args.primer_min_tm+'\n'
                             +"PRIMER_OPT_TM=%s" %args.primer_opt_tm+'\n'
                             +"PRIMER_MAX_TM=%s" %args.primer_max_tm+'\n'
                             +"PRIMER_MAX_NS_ACCEPTED=1"+'\n'
                             +"PRIMER_NUM_RETURN=5"+'\n'
                             +"P3_FILE_FLAG=1"+'\n'
                             +"PRIMER_EXPLAIN_FLAG=1"+'\n'
                             +"PRIMER_MISPRIMING_LIBRARY=%s" %args.mispriming+'\n'
                             +"PRIMER_MIN_GC=%s" %args.primer_min_gc+'\n'
                             +"PRIMER_OPT_GC_PERCENT=%s" %args.primer_opt_gc+'\n'
                             +"PRIMER_MAX_GC=%s" %args.primer_max_gc+'\n'
                             +"PRIMER_PAIR_MAX_DIFF_TM=3"+'\n'
                             +"="+'\n')
    elif args.sv == 'inversion':
        flank_data = flanking_regions_fasta_inversion(genome, small_regions, args.flanking_region_size)
        primer3_in = open("primer3_input.%s.txt" % timestr, 'w')
        for head, seq in flank_data:
            flanking.write(">"+head+'\n'+seq+'\n')
            # 6) generate primer3 input file
            primer3_in.write("SEQUENCE_ID="+head+'\n'
                             +"SEQUENCE_TEMPLATE="+seq+'\n'
                             +"SEQUENCE_TARGET=%s" %args.sequence_target+'\n'
                             +"PRIMER_FIRST_BASE_INDEX=1"+'\n'
                             +"PRIMER_TASK=pick_detection_primers"+'\n'
                             +"PRIMER_MIN_THREE_PRIME_DISTANCE=3"+'\n'
                             +"PRIMER_MAX_LIBRARY_MISPRIMING=12.00"+'\n'
                             +"PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00"+'\n'
                             +"PRIMER_PRODUCT_SIZE_RANGE=%s" %args.product_size_range+'\n'
                             +"PRIMER_MAX_END_STABILITY=9.0"+'\n'
                             +"PRIMER_MAX_SELF_ANY_TH=45.00"+'\n'
                             +"PRIMER_MAX_SELF_END_TH=35.00"+'\n'
                             +"PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00"+'\n'
                             +"PRIMER_PAIR_MAX_COMPL_END_TH=35.00"+'\n'
                             +"PRIMER_MAX_HAIRPIN_TH=24.00"+'\n'
                             +"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00"+'\n'
                             +"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00"+'\n'
                             +"PRIMER_TM_FORMULA=1"+'\n' # use SantaLucia parameters
                             +"PRIMER_SALT_CORRECTIONS=1"+'\n' # SantaLucia 1998 paper
                             +"PRIMER_SALT_MONOVALENT=50.0"+'\n' # mM conc of monovalent salt cations
                             +"PRIMER_INTERNAL_SALT_MONOVALENT=50.0"+'\n' # same as above
                             +"PRIMER_SALT_DIVALENT=1.5"+'\n'
                             +"PRIMER_INTERNAL_SALT_DIVALENT=1.5"+'\n'
                             +"PRIMER_DNTP_CONC=0.6"+'\n' # assume no dntps are present when hybridizing
                             +"PRIMER_INTERNAL_DNTP_CONC=0.6"+'\n'
                             +"PRIMER_DNA_CONC=50.0"+'\n'
                             +"PRIMER_INTERNAL_DNA_CONC=50.0"+'\n'
                             +"PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1"+'\n'
                             +"PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1"+'\n'
                             +"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s" %args.thermopath+'\n'
                             +"PRIMER_PICK_LEFT_PRIMER=1"+'\n'
                             +"PRIMER_PICK_RIGHT_PRIMER=1"+'\n'
                             +"PRIMER_PICK_INTERNAL_OLIGO=1"+'\n'
                             +"PRIMER_MAX_POLY_X=3"+'\n'
                             +"PRIMER_LEFT_NUM_RETURNED=5"+'\n'
                             +"PRIMER_RIGHT_NUM_RETURNED=5"+'\n'
                             +"PRIMER_OPT_SIZE=%s" %args.primer_opt_size+'\n'
                             +"PRIMER_MIN_SIZE=%s" %args.primer_min_size+'\n'
                             +"PRIMER_MAX_SIZE=%s" %args.primer_max_size+'\n'
                             +"PRIMER_MIN_TM=%s" %args.primer_min_tm+'\n'
                             +"PRIMER_OPT_TM=%s" %args.primer_opt_tm+'\n'
                             +"PRIMER_MAX_TM=%s" %args.primer_max_tm+'\n'
                             +"PRIMER_MAX_NS_ACCEPTED=1"+'\n'
                             +"PRIMER_NUM_RETURN=5"+'\n'
                             +"P3_FILE_FLAG=1"+'\n'
                             +"PRIMER_EXPLAIN_FLAG=1"+'\n'
                             +"PRIMER_MISPRIMING_LIBRARY=%s" %args.mispriming+'\n'
                             +"PRIMER_MIN_GC=%s" %args.primer_min_gc+'\n'
                             +"PRIMER_OPT_GC_PERCENT=%s" %args.primer_opt_gc+'\n'
                             +"PRIMER_MAX_GC=%s" %args.primer_max_gc+'\n'
                             +"PRIMER_PAIR_MAX_DIFF_TM=3"+'\n'
                             +"="+'\n')
    flanking.close()
    primer3_in.close()


def iterator(args):
    """ Use an input regions file with specific region of interest\
        to design primers around, then run primer3.

    Args:
        args (Namespace): Argparse object or None.

    Returns: None
    """

    # dedicated string of time for filename output.
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # 1) create genome tuple from provided reference
    genome = genome_iterator(args.ref_genome)
    # 2) create dataframe from input regions file
    small_regions = file_extension(args.regions_file)
    # 3) ensure proper proper number of columns in dataframe
    assert len(list(small_regions)) == 4, "DataFrame contains more than 4 columns...\
                                           Improper format."
    # 4) format dataframe "chr" column to match reference genome
    small_regions = match_chr_to_genome(small_regions, genome)
    # 5) generate flanking regions fasta based on position in input file
    flanking = open("flanking_regions.%s.fasta" % timestr, 'w')
    flank_data = create_flanking_regions_fasta(genome, small_regions, args.flanking_region_size)
    primer3_in = open("primer3_input.%s.txt" % timestr, 'w')
    for head, seq in flank_data:
        flanking.write(">"+head+'\n'+seq+'\n')
    # 6) generate primer3 input file
        primer3_in.write("SEQUENCE_ID="+head+'\n'
                         +"SEQUENCE_TEMPLATE="+seq+'\n'
                         +"SEQUENCE_TARGET=%s" %args.sequence_target+'\n'
                         +"PRIMER_FIRST_BASE_INDEX=1"+'\n'
                         +"PRIMER_TASK=pick_detection_primers"+'\n'
                         +"PRIMER_MIN_THREE_PRIME_DISTANCE=3"+'\n'
                         +"PRIMER_MAX_LIBRARY_MISPRIMING=12.00"+'\n'
                         +"PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00"+'\n'
                         +"PRIMER_PRODUCT_SIZE_RANGE=%s" %args.product_size_range+'\n'
                         +"PRIMER_MAX_END_STABILITY=9.0"+'\n'
                         +"PRIMER_MAX_SELF_ANY_TH=45.00"+'\n'
                         +"PRIMER_MAX_SELF_END_TH=35.00"+'\n'
                         +"PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00"+'\n'
                         +"PRIMER_PAIR_MAX_COMPL_END_TH=35.00"+'\n'
                         +"PRIMER_MAX_HAIRPIN_TH=24.00"+'\n'
                         +"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00"+'\n'
                         +"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00"+'\n'
                         +"PRIMER_TM_FORMULA=1"+'\n' # use SantaLucia parameters
                         +"PRIMER_SALT_CORRECTIONS=1"+'\n' # SantaLucia 1998 paper
                         +"PRIMER_SALT_MONOVALENT=50.0"+'\n' # mM conc of monovalent salt cations
                         +"PRIMER_INTERNAL_SALT_MONOVALENT=50.0"+'\n' # same as above
                         +"PRIMER_SALT_DIVALENT=1.5"+'\n'
                         +"PRIMER_INTERNAL_SALT_DIVALENT=1.5"+'\n'
                         +"PRIMER_DNTP_CONC=0.6"+'\n'
                         +"PRIMER_INTERNAL_DNTP_CONC=0.6"+'\n'
                         +"PRIMER_DNA_CONC=50.0"+'\n'
                         +"PRIMER_INTERNAL_DNA_CONC=50.0"+'\n'
                         +"PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1"+'\n'
                         +"PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1"+'\n'
                         +"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s" %args.thermopath+'\n'
                         +"PRIMER_PICK_LEFT_PRIMER=1"+'\n'
                         +"PRIMER_PICK_RIGHT_PRIMER=1"+'\n'
                         +"PRIMER_PICK_INTERNAL_OLIGO=1"+'\n'
                         +"PRIMER_MAX_POLY_X=3"+'\n'
                         +"PRIMER_LEFT_NUM_RETURNED=5"+'\n'
                         +"PRIMER_RIGHT_NUM_RETURNED=5"+'\n'
                         +"PRIMER_OPT_SIZE=%s" %args.primer_opt_size+'\n'
                         +"PRIMER_MIN_SIZE=%s" %args.primer_min_size+'\n'
                         +"PRIMER_MAX_SIZE=%s" %args.primer_max_size+'\n'
                         +"PRIMER_MIN_TM=%s" %args.primer_min_tm+'\n'
                         +"PRIMER_OPT_TM=%s" %args.primer_opt_tm+'\n'
                         +"PRIMER_MAX_TM=%s" %args.primer_max_tm+'\n'
                         +"PRIMER_MAX_NS_ACCEPTED=1"+'\n'
                         +"PRIMER_NUM_RETURN=5"+'\n'
                         +"P3_FILE_FLAG=1"+'\n'
                         +"PRIMER_EXPLAIN_FLAG=1"+'\n'
                         +"PRIMER_MISPRIMING_LIBRARY=%s" %args.mispriming+'\n'
                         +"PRIMER_MIN_GC=%s" %args.primer_min_gc+'\n'
                         +"PRIMER_OPT_GC_PERCENT=%s" %args.primer_opt_gc+'\n'
                         +"PRIMER_MAX_GC=%s" %args.primer_max_gc+'\n'
                         +"PRIMER_PAIR_MAX_DIFF_TM=3"+'\n'
                         +"="+'\n')
    flanking.close()
    primer3_in.close()

def pre(args):
    """ Function for all steps leading up to PCR. """
    # 1) Initialize primer lists by rank for each sample
    print(args.dump)
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
    primer_df_standard = primer_df.copy()
    # 5) Get length of forward primers for percent alignment check
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
    df_bool.to_csv('no_dimer_df.csv', index=False)
    # 12) create allvsall pcr, standard pcr, or both
    if args.pcr == 'multiplex':
        pch.all_vs_all_pcr(df_bool)
    elif args.pcr == 'standard':
        pch.standard_pcr(primer_df_standard)
    else:
        print("Please select pcr setup")

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



