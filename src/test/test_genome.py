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
        self.genome = Genome.Fasta.from_fasta('./data/test_standard.fa')
        self.regions_file = "./data/test_input_standard.csv"
        self.header_fail = "./data/test_no_header.csv"

    def test_logger(self):
        self.assertIsInstance(self.genome.logger, logging.Logger)

    def test_from_fasta(self):
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

    def test_gets_proper_index(self):
        """ Ref genomes are 1 indexed"""
        fasta_seqs = self.genome.create_flanking_regions(self.regions_file, 0)
        self.assertEqual(fasta_seqs[0][1], 'T')
        self.assertEqual(fasta_seqs[1][1], 'G')
        self.assertEqual(fasta_seqs[2][1], 'C')

    def test_gets_proper_flanking(self):
        """ Now we want to ensure flanking regions pulled correctly"""
        fasta_seqs = self.genome.create_flanking_regions(self.regions_file, 10)
        self.assertEqual(len(fasta_seqs[0][1]), 21)
        self.assertEqual(fasta_seqs[0][1], 'ACCCTAACCCTAACCCTAACC')

    def tearDown(self) -> None:
        os.remove('./primertk.log')

        