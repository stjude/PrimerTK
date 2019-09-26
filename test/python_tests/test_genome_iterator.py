#!/usr/bin/env python3
"""
Unit tester for genome_iterator.py

Python version: Python 3.6.8 :: Anaconda, Inc.
Date: 02/28/2019
Dependencies:
 - unittest
 - pandas>=0.22.0
"""

import unittest
import pandas as pd
from primer_tk import genome_iterator as gi

class TestGenomeIterator(unittest.TestCase):
    """
    Subclass of unittest to test genome_iterator.py
    """
    def setUp(self):
        self.ref_genome = "test/data/test_standard.fa"
        self.genome = gi.genome_iterator(self.ref_genome)
        self.test_input = 'test/data/input_standard.csv'
        self.test_input2 = 'test/data/input_standard.txt'
        self.test_input3 = 'test/data/input_standard.fa'
        #self.parser = gi.get_args_iterator([])

    def test_genome_iterator(self):
        """
        Asserts object output by genome_iterator function has length of 2 (header, seq)
        """
        self.assertEqual(len(self.genome[0]), 2)

    def test_create_dataframe_csv(self):
        """
        Asserts that a pd.DataFrame is returned from function.
        """
        self.assertTrue(isinstance(gi.create_dataframe_csv(self.test_input), pd.DataFrame))

    def test_create_dataframe_txt(self):
        """
        Assets that a pd.DataFrame is returned from function.
        """
        self.assertTrue(isinstance(gi.create_dataframe_txt(self.test_input2), pd.DataFrame))

    def test_file_extension_success(self):
        """
        Asserts that a pd.DataFrame is returned only if file extension matches.
        """
        self.assertTrue(isinstance(gi.file_extension(self.test_input2), pd.DataFrame))

    def test_file_extension_failure(self):
        """
        Asserts that code exits with warning message if file extension is incorrect.
        """
        with self.assertRaises(SystemExit) as code_message:
            gi.file_extension(self.test_input3)
        self.assertEqual(code_message.exception.code,
                         "Wrong File Format, should be .txt or .csv")

    def test_match_chr_to_genome(self):
        """
        Asserts that a pd.DataFrame object is returned from function.
        """
        dataframe = gi.file_extension(self.test_input)
        self.assertTrue(isinstance(gi.match_chr_to_genome(dataframe, self.genome), pd.DataFrame))


    def test_flanking_regions_fasta(self):
        """
        Asserts that return object (list of tuples) first item is sample.
        """
        dataframe = gi.file_extension(self.test_input)
        self.assertTrue((gi.create_flanking_regions_fasta\
                             (self.genome, dataframe, 200)[0][0]).startswith('sample'))

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
