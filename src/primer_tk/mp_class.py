#!/usr/bin/env python3

"""
Dependencies required to run program:
 - python3.6+
    - pandas>=0.22.0
    - numpy>=1.16.0
"""

import sys
import pandas as pd
import numpy as np

class MissingPrimers:
    """
    Groups Primer3 output by sample id into a list of lists,
    adds NA to all groups that do not have 5 primers generated,
    and then outputs a list corresponding to each primer rank.
    Args:
        file_name (string): name of file to be passed to program
        value (int): rank of primer pair to store as list [0-4] only
    Returns:
        samp_primer_info (list): list of lists with all values required for
                                 downstream analysis (id, lseq, rseq, l_tm, r_tm,
                                 l_gc, r_gc, product_size)
    """
    def __init__(self, file_name, value):
        """
        Initialize Values
        """
        self.file_name = file_name
        if value >= 0 and value <=4:
            self.value = str(value)
        else:
            sys.exit("Improper value selected! Acceptable values: 0,1,2,3,4.")
        self.dump_list = self.__group_seqids()
        self.filled_primers = self.__fill_empty_values()
        self.samp_primer_info = self.__gather_primer_info()

    def __group_seqids(self):
        """
        Group by sample_id
        """
        sequence = ""
        with open(self.file_name) as dump:
            for line in dump:
                if not line.startswith('='):
                    sequence += line
                if line.startswith('='):
                    fixed_line = line.replace('=', '@@@')
                    sequence += fixed_line
                full_string = ''.join([line for line in sequence])
                primer_split = [[string] for string in full_string.\
                                    split('@@@') if string is not ' ']
        primers_info = primer_split[:-1]
        return primers_info

    def __fill_empty_values(self):
        """
        Fill in missing required values with NA
        """
        sequence = []
        for sample in self.dump_list:
            for string in sample:
                if "PRIMER_LEFT_0_SEQUENCE=" not in string:
                    string = string + "PRIMER_LEFT_0_SEQUENCE=NA\n"\
                        + "PRIMER_RIGHT_0_SEQUENCE=NA\n"\
                        + "PRIMER_LEFT_0_TM=NA\n" + "PRIMER_RIGHT_0_TM=NA\n"\
                        + "PRIMER_LEFT_0_GC_PERCENT=NA\n"\
                        + "PRIMER_RIGHT_0_GC_PERCENT=NA\n"\
                        + "PRIMER_PAIR_0_PRODUCT_SIZE=NA\n"
                if "PRIMER_LEFT_1_SEQUENCE=" not in string:
                    string = string + "PRIMER_LEFT_1_SEQUENCE=NA\n"\
                        + "PRIMER_RIGHT_1_SEQUENCE=NA\n"\
                        + "PRIMER_LEFT_1_TM=NA\n" + "PRIMER_RIGHT_1_TM=NA\n"\
                        + "PRIMER_LEFT_1_GC_PERCENT=NA\n"\
                        + "PRIMER_RIGHT_1_GC_PERCENT=NA\n"\
                        + "PRIMER_PAIR_1_PRODUCT_SIZE=NA\n"
                if "PRIMER_LEFT_2_SEQUENCE=" not in string:
                    string = string + "PRIMER_LEFT_2_SEQUENCE=NA\n"\
                        + "PRIMER_RIGHT_2_SEQUENCE=NA\n"\
                        + "PRIMER_LEFT_2_TM=NA\n" + "PRIMER_RIGHT_2_TM=NA\n"\
                        + "PRIMER_LEFT_2_GC_PERCENT=NA\n"\
                        + "PRIMER_RIGHT_2_GC_PERCENT=NA\n"\
                        + "PRIMER_PAIR_2_PRODUCT_SIZE=NA\n"
                if "PRIMER_LEFT_3_SEQUENCE=" not in string:
                    string = string + "PRIMER_LEFT_3_SEQUENCE=NA\n"\
                        + "PRIMER_RIGHT_3_SEQUENCE=NA\n"\
                        + "PRIMER_LEFT_3_TM=NA\n" + "PRIMER_RIGHT_3_TM=NA\n"\
                        + "PRIMER_LEFT_3_GC_PERCENT=NA\n"\
                        + "PRIMER_RIGHT_3_GC_PERCENT=NA\n"\
                        + "PRIMER_PAIR_3_PRODUCT_SIZE=NA\n"
                if "PRIMER_LEFT_4_SEQUENCE=" not in string:
                    string = string + "PRIMER_LEFT_4_SEQUENCE=NA\n"\
                        + "PRIMER_RIGHT_4_SEQUENCE=NA\n"\
                        + "PRIMER_LEFT_4_TM=NA\n" + "PRIMER_RIGHT_4_TM=NA\n"\
                        + "PRIMER_LEFT_4_GC_PERCENT=NA\n"\
                        + "PRIMER_RIGHT_4_GC_PERCENT=NA\n"\
                        + "PRIMER_PAIR_4_PRODUCT_SIZE=NA\n"
                sequence.append(string.lstrip('\n').split('\n'))

        return sequence

    def __gather_primer_info(self):
        """
        Return the final list of lists with all filled values by rank
        """
        sample_info = []
        for item in self.filled_primers:
            for p_info in item:
                if "SEQUENCE_ID=" in p_info:
                    seq_id = p_info[12:]+self.value
                if "PRIMER_LEFT_%s_SEQUENCE=" %self.value in p_info:
                    p_left = p_info[23:]
                if "PRIMER_RIGHT_%s_SEQUENCE=" %self.value in p_info:
                    p_right = p_info[24:]
                if "PRIMER_LEFT_%s_TM=" %self.value in p_info:
                    left_tm = p_info[17:]
                if "PRIMER_RIGHT_%s_TM=" %self.value in p_info:
                    right_tm = p_info[18:]
                if "PRIMER_LEFT_%s_GC_PERCENT=" %self.value in p_info:
                    left_gc = p_info[25:]
                if "PRIMER_RIGHT_%s_GC_PERCENT=" %self.value in p_info:
                    right_gc = p_info[26:]
                if "PRIMER_PAIR_%s_PRODUCT_SIZE=" %self.value in p_info:
                    product_size = p_info[27:]
                    sample_info.append([seq_id, p_left, p_right, left_tm,
                                        right_tm, left_gc, right_gc,
                                        product_size])

        return sample_info

def create_df(primer_lists):
    """
    Creates a dataframe from gather_primer_info
    Args:
        primer_lists (list): primer information lists to be passed for dataframe.
    Returns:
        primer_df (pd.DataFrame): Dataframe with all primer information.
    """
    primer_df = pd.DataFrame(np.ma.row_stack(primer_lists),
                             columns=['Sequence ID', 'Primer Left Seq', 'Primer Right Seq',
                                      'Primer Left TM', 'Primer Right TM',
                                      'Primer Left GC %', 'Primer Right GC %',
                                      'Primer3 Predicted Product'])
    primer_df = primer_df.sort_values('Sequence ID').reset_index().drop(labels='index', axis=1)
    rank = [num[-1] for num in primer_df['Sequence ID']]
    primer_df.insert(1, 'Primer Rank', rank)
    left_len = [len(length) for length in primer_df['Primer Left Seq']]
    right_len = [len(length) for length in primer_df['Primer Right Seq']]
    primer_df.insert(4, 'Primer Left Len', left_len)
    primer_df.insert(5, 'Primer Right Len', right_len)
    primer_df['Sequence ID'] = [seqid[:-1] for seqid in primer_df['Sequence ID']]
    return primer_df
