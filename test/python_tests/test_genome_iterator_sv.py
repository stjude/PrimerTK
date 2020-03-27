#!/usr/bin/env python3
"""
Unit tester for genome_iterator_sv.py

Python version: Python 3.6.8 :: Anaconda, Inc.
Date: 02/28/2019
Dependencies:
 - unittest
 - pandas >= 0.22.0
"""

import unittest
import pandas as pd
from primer_tk import genome_iterator_sv as gi

class TestGenomeIterator(unittest.TestCase):
    """
    Subclass of unittest to test genome_iterator_sv.py
    """
    def setUp(self):
        self.ref_genome = "test/data/test_sv.fa"
        self.genome = gi.genome_iterator(self.ref_genome)
        self.test_input = 'test/data/input_sv.csv'
        self.test_input2 = 'test/data/input_sv.txt'
        self.test_input3 = 'test/data/input_sv.fa'
        self.translocation_input1 = 'test/data/translocation_input.csv'
        self.translocation_input2 = 'test/data/translocation_input.txt'
        self.insertion_genome = gi.genome_iterator('test/data/insertion_test.fa')
        self.insertion_input1 = 'test/data/insertion_input.csv'
        self.insertion_input2 = 'test/data/insertion_input.txt'

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
        Asserts that a pd.DataFrame is returned from function.
        """
        self.assertTrue(isinstance(gi.create_dataframe_txt(self.test_input2), pd.DataFrame))

    def test_create_dataframe_insertion_csv(self):
        """
        Asserts that proper insertion dataframe is returned from function.
        """
        self.assertEqual(len(gi.create_dataframe_insertion_csv(self.insertion_input1)),4)

    def test_create_dataframe_insertion_txt(self):
        """
        Asserts that proper insertion dataframe is returned from function.
        """
        self.assertEqual(len(gi.create_dataframe_insertion_txt(self.insertion_input2)),4)
    def test_create_dataframe_translocation_csv(self):
        """
        Asserts that proper translocation dataframe is returned from function.
        """
        self.assertEqual(len(gi.create_dataframe_translocation_csv(self.translocation_input1)),4)

    def test_create_dataframe_translocation_txt(self):
        """
        Asserts that proper translocation dataframe is returned from function.
        """
        self.assertEqual(len(gi.create_dataframe_translocation_txt(self.translocation_input2)),4)

    def test_file_extension_insertion_success(self):
        """
        Asserts that proper df is created when sv type is insertion.
        """
        self.assertEqual(len(gi.file_extension(self.insertion_input1, 'insertion')),4)

    def test_file_extension_translocation_success(self):
        """
        Asserts that proper df is created when sv type is translocation.
        """
        self.assertEqual(len(gi.file_extension(self.translocation_input1, 'translocation')),4)

    def test_file_extension_deletion_success(self):
        """
        Asserts that a pd.DataFrame is returned only if file extension matches.
        """
        self.assertTrue(isinstance(gi.file_extension(self.test_input2, 'deletion'), pd.DataFrame))

    def test_file_extension_inversion_success(self):
        """
        Asserts that a pd.DataFrame is returned only if file extension matches.
        """
        self.assertTrue(isinstance(gi.file_extension(self.test_input2, 'inversion'), pd.DataFrame))

    def test_file_extension_failure(self):
        """
        Asserts that code exits with warning message if file extension or sv is incorrect.
        """
        with self.assertRaises(SystemExit) as code_message:
            gi.file_extension(self.test_input3, 'deletion')
        self.assertEqual(code_message.exception.code,
                         "Wrong File Format, should be .txt (tab) or .csv (comma), or check sv type.")

    def test_match_chr_to_genome(self):
        """
        Asserts that a pd.DataFrame object is returned from function.
        chr format is mismatched between genome and input file, checks for conversion.
        """
        dataframe = gi.file_extension(self.test_input2, 'deletion')
        self.assertTrue(isinstance(gi.match_chr_to_genome(dataframe, self.genome, 'deletion'), pd.DataFrame))

    def test_match_chr_to_genome_insertion(self):
        """
        Asserts that input file with format chr1 is converted to 1 to match ref genome.
        """
        dataframe = gi.file_extension(self.insertion_input1, 'insertion')
        self.assertFalse(gi.match_chr_to_genome(dataframe, self.insertion_genome, 'insertion')\
                             ["ChrNorm"].str.contains("chr").any())

    def test_flanking_regions_fasta_deletion(self):
        """
        Asserts that return object (list of tuples) first item is sample.
        """
        dataframe = gi.file_extension(self.test_input, 'deletion')
        self.assertTrue((gi.flanking_regions_fasta_deletion\
                             (self.genome, dataframe, 10)[0][0]).startswith('sample'))

    def test_flanking_regions_fasta_translocation(self):
        """
        Asserts that translocations are handled properly for both forward and reverse strand translocations.
        """
        dataframe = gi.file_extension(self.translocation_input2, 'translocation')
        flanking = gi.flanking_region_fasta_translocation(self.insertion_genome, dataframe, 5)
        self.assertTrue(len(flanking) == 4)

    def test_flanking_regions_fasta_insertion(self):
        """
        Asserts that insertions are handled properly for both forward and reverse strand insertions.
        """
        dataframe = gi.file_extension(self.insertion_input2, 'insertion')
        flanking = gi.flanking_region_fasta_insertion(self.insertion_genome, dataframe, 5)
        self.assertTrue(len(flanking) == 8)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
