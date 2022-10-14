"""
Unittest for Genome
"""

import os
import logging
import unittest

import pandas as pd

from primertk import Genome

class TestGenome(unittest.TestCase):
    """ Unittest for Genome class """
    def setUp(self) -> None:
        self.genome = Genome.Genome('./data/test_standard.fa', 'test')
        self.genome.genome_iterator()
        self.regions_file = "./data/test_input_standard.csv"
        self.header_fail = "./data/test_no_header.csv"

    def test_logger(self):
        self.assertIsInstance(self.genome.logger, logging.Logger)

    def test_genome_iterator(self):
        self.assertEqual(len(self.genome.reference), 2)
        self.assertEqual(self.genome.reference[0][0], '1')
        self.assertEqual(self.genome.reference[1][0], '2')

    def test_parse_input(self):
        self.assertIsInstance(Genome.parse_input(self.regions_file), pd.DataFrame)

    def test_no_header_fail(self):
        with self.assertRaises(AssertionError):
            Genome.parse_input(self.header_fail)
        
    def test_match_chr_to_genome(self):
        df = Genome.parse_input(self.regions_file)
        Genome.match_chr_to_genome(df, self.genome.reference)
        self.assertFalse(df['Chr'].str.contains("chr").all())

    def tearDown(self) -> None:
        os.remove('./test.log')

        