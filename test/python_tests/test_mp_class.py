#!/usr/bin/env python3
"""
Unit tester for mp_class.py

Python version: Python 3.6.8 :: Anaconda, Inc.
Date: 03/20/2019
Dependencies:
 - unittest
 - pandas>=0.22.0
 - numpy>=1.16.0
"""

import unittest
import pandas as pd
import numpy as np
from primer_tk.mp_class import MissingPrimers
from primer_tk.mp_class import create_df

class TestMissingPrimers(unittest.TestCase):
    """
    Subclass of unittest to test mp_class.py
    """
    def setUp(self):
        self.mp = MissingPrimers("test/data/primer_dump_standard.txt", 1)

    def test_group_seqids(self):
        """
        Creates a list of lists, each sample info is stored in its own list
        within a list. Since there were 3 samples, the length of the list
        should be 3.
        """
        self.assertEqual(len(self.mp._MissingPrimers__group_seqids()), 3)

    def test_group_seqids_unfilled(self):
        """
        Tests to ensure that NA values have not been inserted into
        primer3 string.
        """
        no_missing_vals = self.mp._MissingPrimers__group_seqids()[1]
        self.assertFalse("NA" in no_missing_vals)

    def test_fill_empty_values(self):
        """
        Tests to ensure that NA values have been inserted into the 
        sample with missing primers.
        """
        all_params = self.mp._MissingPrimers__fill_empty_values()
        missing_params = all_params[1]
        self.assertTrue("=NA" in missing_params[-2])#-1 is empty space
                
    def test_gather_primer_info(self):
        """
        Test to ensure that if value 1 is passed, all sample names
        end with 1.
        """
        primer_info = self.mp._MissingPrimers__gather_primer_info()
        self.assertTrue(all(item[0].endswith('1') for item in primer_info))

    def test_invalid_value(self):
        """
        Tests to ensure class does not succeed if invalid value is passed.
        Value options are currently [0-4].
        """
        with self.assertRaises(SystemExit) as code_message:
            MissingPrimers("primer_dump_standard.txt", 6)
        self.assertEqual(code_message.exception.code,
                         "Improper value selected! Acceptable values: 0,1,2,3,4.")

    
    def test_len_list_output(self):
        """
        Tests to ensure that the list created for each sample contains 8 values.
        This corresponds to the number of df columns.
        """
        prim_list_1 = self.mp.samp_primer_info
        for sample_data in prim_list_1:
            self.assertEqual(len(sample_data), 8)

    def test_create_df(self):
        """
        Tests to ensure that a df is created from primer list
        and contains the proper number of rows.
        """
        prim_list_1 = self.mp.samp_primer_info
        primer_df = create_df([prim_list_1])
        self.assertEqual(len(primer_df), 3)

    def test_accept_multiple_prim_lists(self):
        """
        Tests to ensure that a mp_class.create_df() can accept multiple primer lists
        as input and generate an output df.
        """
        prim_list_1 = self.mp.samp_primer_info
        prim_list_2 = MissingPrimers("test/data/primer_dump_standard.txt", 2).samp_primer_info
        primer_df = create_df([prim_list_1, prim_list_2])
        self.assertEqual(len(primer_df), 6)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
