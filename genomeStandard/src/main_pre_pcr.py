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
import primer_cross_hyb as pch

def main():
    """ Function for all steps leading up to PCR. """
    args = pch.get_args_hyb()
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

if __name__ == "__main__":
    main()
