#!/usr/bin/env python3
"""
Unit tester for mp_class_sv.py

Python version: Python 3.6.8 :: Anaconda, Inc.
Dependencies:
 - unittest
 - pandas >= 0.22.0
"""

import unittest
import pandas as pd
from primer_tk.mp_class_sv import MissingPrimers
from primer_tk.mp_class_sv import create_df

class TestMpClassSV(unittest.TestCase):
    """
    Subclass of unittest to test mp_class_sv.py
    """
    def setUp(self):
        self.mp = MissingPrimers("test/data/primer_dump_sv.txt", 1)

    def test_group_seqids(self):
        """
        Asserts proper list is  generated for given number of sample values.
        """
        self.assertEqual(len(self.mp._MissingPrimers__group_seqids()), 3)

    def test_group_seqids_unfilled(self):
        """
        Tests to ensure NA values have not been inserted into primer3 string.
        """
        no_missing_vals = self.mp._MissingPrimers__group_seqids()[1]
        self.assertFalse("NA" in no_missing_vals)

    def test_fill_empty_values(self):
        """
        Tests to ensure that NA vals have been inserted into samples with missing primers.
        """
        all_params = self.mp._MissingPrimers__fill_empty_values()
        missing_params = all_params[1]
        self.assertTrue("=NA" in missing_params[-2])

    def test_gather_primer_info(self):
        """
        Test to ensure that if 1 value is passed, all sample names
        end with 1. Ex: SampleABC1
        """
        primer_info = self.mp._MissingPrimers__gather_primer_info()
        self.assertTrue(all(item[0].endswith('1') for item in primer_info))

    def test_len_list_output(self):
        """
        Tests to ensure the proper number of params have been kept from the p3 output.
        Slightly longer than normal output because coordinate positions are kept for downstream.
        """
        prim_list_1 = self.mp.samp_primer_info
        for sample_data in prim_list_1:
            self.assertEqual(len(sample_data), 10)

    def test_create_df(self):
        """
        Tests to ensure a df is created from primer list.
        """
        prim_list_1 = self.mp.samp_primer_info
        primer_df = create_df([prim_list_1])
        self.assertTrue(isinstance(primer_df, pd.DataFrame))

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
