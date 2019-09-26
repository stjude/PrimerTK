#!/usr/bin/env python3

"""
Modules required for program.
 - python3.6+
   - pandas >= 0.22.0
   - sequence_info
"""

import unittest
import pandas as pd
from primer_tk import analyze_pcr_output as ap

class TestAnalyzePcrOutput(unittest.TestCase):
    """
    Subclass of unittest to test analyze_pcr_output.py
    """
    def setUp(self):
        self.pcr_file = 'test/data/pcr_standard.fa'
        self.primer_file = 'test/data/total_primers_standard.csv'
        self.seqs, self.headers = ap.fasta_parser(self.pcr_file)
        self.gc_list = ap.gc_percent_seqs(self.seqs)
        self.split_headers, self.no_chrom = ap.split_headers_list(self.headers)
        self.chr_list = ap.chr_split_list(self.split_headers)
        self.pos_list = ap.pos_split_list(self.chr_list)
        self.name_pos_split_list = ap.split_name_pos(self.no_chrom)
        for item in self.chr_list:
            del item[-1]
        self.merged = ap.merge_info(self.chr_list, self.pos_list,
                                    self.name_pos_split_list,
                                    self.no_chrom)
        self.pcr_df, self.good_primers, self.bad_primers = ap.generate_pcr_df(self.merged,
                                                                              self.gc_list)
        self.merged_df = ap.merge_good_total(self.good_primers, self.primer_file)
        self.filtered = ap.filter_merged(self.merged_df, 1)

    def test_fasta_parser(self):
        """
        Ensure that seq and header are proper length.
        """
        seqs, headers = ap.fasta_parser(self.pcr_file)
        self.assertTrue(len(seqs[0][0]) == 306)
        self.assertTrue(len(headers[0]) == 98)

    def test_gc_percent_seqs(self):
        """
        Calcs gc percent of list of seqs
        """
        self.assertEqual(int(ap.gc_percent_seqs(self.seqs)[0]), 49)

    def test_split_headers_list(self):
        """
        Split header information and return tuple of parsed header info.
        """
        self.assertTrue(len(ap.split_headers_list(self.headers)[0][0]) >
                        len(ap.split_headers_list(self.headers)[1][0]))

    def test_chr_split_list(self):
        """
        Parses chrm and position out of headers list.
        """
        self.assertEqual(ap.chr_split_list(self.split_headers)[0][0], '1')

    def test_split_name_pos(self):
        """
        Gets the expected position of the snp.
        """
        self.assertEqual(ap.split_name_pos(self.no_chrom)[0][0], '15000000')

    def test_merge_info(self):
        """
        Merges all the info extracted from the fasta files.
        """
        self.assertEqual(len(ap.merge_info(self.chr_list, self.pos_list,
                                       self.name_pos_split_list,
                                           self.no_chrom)[0]), 8)

    def test_generate_pcr_df_off_target(self):
        """
        Creates a df from all info.
        """
        pcr_df, good_primers, bad_primers = ap.generate_pcr_df(self.merged, self.gc_list)
        self.assertEqual(len(pcr_df), 2)
        self.assertEqual(len(good_primers), 1)
        self.assertEqual(len(bad_primers), 1)

    def test_merge_good_total(self):
        """
        Merges pcr df with total primers df.
        """
        self.assertEqual(len(ap.merge_good_total(self.good_primers, self.primer_file)), 12)

    def test_filter_merged(self):
        """
        Drops off target products.
        """
        self.assertEqual(len(ap.filter_merged(self.merged_df, 0)), 0)

    def top_ranked_final_primers(self):
        """
        Selects the top ranking primer pair from each target.
        """
        self.assertEqual(len(ap.top_ranked_final_primers(self.filtered)), 1)

    def test_to_order_plate(self):
        """
        Generates a ready made primer ordering sheet.
        """
        forward, reverse = ap.to_order_plate(self.filtered)
        self.assertEqual(len(forward.columns), 3)
        self.assertEqual(len(forward), 1)

    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
