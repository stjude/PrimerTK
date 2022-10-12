"""
Unittest for Genome
"""

import os
import logging
import unittest

from primertk import Genome

class TestGenome(unittest.TestCase):
    """ Unittest for Genome class """
    def setUp(self) -> None:
        self.genome = Genome.Genome('./data/test_standard.fa', 'test')

    def test_logger(self):
        self.assertIsInstance(self.genome.logger, logging.Logger)

    def test_genome_iterator(self):
        self.genome.genome_iterator()
        self.assertEqual(len(self.genome.reference), 2)

    def tearDown(self) -> None:
        os.remove('./test.log')

        