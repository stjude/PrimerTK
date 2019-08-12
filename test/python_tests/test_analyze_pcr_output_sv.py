#!/usr/bin/env python3

"""
Module required for program.
 - python3.6+
    - pandas >= 0.22.0
"""

import unittest
from primer_tk import analyze_pcr_output_sv as ap

class TestAnalyzePcrOutputSV(unittest.TestCase):
    """
    Subclass of unittest to test analyze_pcr_output_sv.py
    """
    def setUp(self):
        self.flanking = 'test/data/flanking_regions_sv.fasta'
        self.primer_file = 'test/data/total_list_sv.csv'
        self.seqs, self.headers = ap.fasta_parser(self.flanking)
        self.positions = ap.amp_header_region(self.primer_file)

    def test_fasta_parser(self):
        """
        Ensure that seq and header are appropriately parsed.
        """
        self.assertEqual(len(self.seqs), 3)
        self.assertEqual(len(self.headers), 3)

    def test_amp_header_region(self):
        """
        Extract header region from primer info file.
        """
        self.assertTrue(all(len(item) == 3 for item in self.positions))

    def test_get_gc_region(self):
        """
        Extract gc actual sequence region between primers based on positions.
        """
        sliced_seqs = ap.get_gc_region(self.seqs, self.headers, self.positions)
        self.assertTrue(all(len(item[2]) < 401 for item in sliced_seqs))
    def test_calc_gc(self):
        """
        Calculate the sample GC content for each sliced sequence.
        """
        sliced_seqs = ap.get_gc_region(self.seqs, self.headers, self.positions)
        gc_calc = ap.calc_gc(sliced_seqs)
        self.assertTrue(all(isinstance(item[2], float) for item in gc_calc))

    def test_merge_dfs(self):
        """
        Merges calc_gc into the original dataframe for more
        information about product length and gc content.
        """
        sliced_seqs = ap.get_gc_region(self.seqs, self.headers, self.positions)
        gc_calc = ap.calc_gc(sliced_seqs)
        self.assertEqual(len(ap.merge_dfs(gc_calc, self.primer_file, self.seqs).columns), 19)
        self.assertEqual(len(ap.merge_dfs(gc_calc, self.primer_file, self.seqs)['Sequence ID']), 13)

    def test_to_order_plate(self):
        """
        Generates a plate order format for generated primers.
        """
        sliced_seqs = ap.get_gc_region(self.seqs, self.headers, self.positions)
        gc_calc = ap.calc_gc(sliced_seqs)
        merged = ap.merge_dfs(gc_calc, self.primer_file, self.seqs)
        plate_order_sheet_f, plate_order_sheet_r = ap.to_order_plate(merged)
        self.assertEqual(len(plate_order_sheet_f), 13)
        self.assertEqual(len(plate_order_sheet_r), 13)

    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
