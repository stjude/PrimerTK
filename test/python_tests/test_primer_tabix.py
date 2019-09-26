#!/usr/bin/env python3
"""
Unit tester for primer_tabix.py

Python version: Python 3.6.8 :: Anaconda, Inc.
Dependencies:
 - unittest
 - pandas >= 0.22.0
 - pysam >= 0.15.2
"""

import unittest
import pandas as pd
from pysam import VariantFile
from primer_tk import primer_tabix as pt

class TestPrimerTabix(unittest.TestCase):
    """
    Subclass of unittest to test primer_tabix.py
    """
    def setUp(self):
        self.primer_output = "test/data/top_ranked_final_primers_standard.csv"
        self.vcf_test = 'test/data/primer_test_vcf.csv'
        self.pinfo = pt.create_tabix_df(self.primer_output)
        self.pvcf = pt.create_tabix_df(self.vcf_test)
        self.vcf = 'test/data/smallvcf.vcf.gz'
        self.vcf_in = VariantFile(self.vcf)

    def test_create_tabix_df(self):
        """
        Test dataframe creation.
        """
        self.assertTrue(isinstance(pt.create_tabix_df(self.primer_output), pd.DataFrame))

    def test_primer_range_left(self):
        """
        Test the range info is added to the dataframe.
        """
        self.assertEqual(len(pt.primer_range_left(self.pinfo['Sequence ID'], self.pinfo['Primer Rank'],
                                                  self.pinfo['Chromosome'], self.pinfo['Primer Left Seq'],
                                                  self.pinfo['Position1']).columns), 6)
    def test_primer_range_right(self):
        """
        Test the range info is added to the dataframe.
        """
        self.assertEqual(len(pt.primer_range_left(self.pinfo['Sequence ID'], self.pinfo['Primer Rank'],
                                                  self.pinfo['Chromosome'], self.pinfo['Primer Right Seq'],
                                                  self.pinfo['Position2']).columns), 6)

    def test_match_pinfo_to_vcf(self):
        """
        Normalizes the chr string to the reference vcf.
        """
        self.assertFalse(pt.match_pinfo_to_vcf(self.pinfo, self.vcf_in)['Chromosome'].str.contains("chr").any())

    def test_tabix_fetch(self):
        """
        Fetches snp info and assigns to primers.
        """
        left = pt.primer_range_left(self.pvcf['Sequence ID'], self.pvcf['Primer Rank'],
                                    self.pvcf['Chromosome'], self.pvcf['Primer Left Seq'],
                                    self.pvcf['Position1'])
        normalized = pt.match_pinfo_to_vcf(left, self.vcf_in)
        left_snps = pt.tabix_fetch(normalized['Sequence ID'], normalized['Primer Rank'],
                                   normalized['Chromosome'], normalized['Position1'],
                                   normalized['Position2'], self.vcf_in)
        self.assertEqual(len(left_snps), 20)

    def test_tabix_results_to_df(self):
        """
        Generates a pd.DataFrame from the tabix results.
        """
        left = pt.primer_range_left(self.pvcf['Sequence ID'], self.pvcf['Primer Rank'],
                                    self.pvcf['Chromosome'], self.pvcf['Primer Left Seq'],
                                    self.pvcf['Position1'])
        normalized = pt.match_pinfo_to_vcf(left, self.vcf_in)
        left_snps = pt.tabix_fetch(normalized['Sequence ID'], normalized['Primer Rank'],
                                   normalized['Chromosome'], normalized['Position1'],
                                   normalized['Position2'], self.vcf_in)
        left_df = pt.tabix_results_to_df(left_snps, "L", "Left SNP Count")
        self.assertEqual(left_df['Left SNP Count'][0], 10)

    def test_merge_left_right(self):
        """
        Merge the left and right SNP dataframes.
        """
        left = pt.primer_range_left(self.pvcf['Sequence ID'], self.pvcf['Primer Rank'],
                                    self.pvcf['Chromosome'], self.pvcf['Primer Left Seq'],
                                    self.pvcf['Position1'])
        normalized = pt.match_pinfo_to_vcf(left, self.vcf_in)
        left_snps = pt.tabix_fetch(normalized['Sequence ID'], normalized['Primer Rank'],
                                   normalized['Chromosome'], normalized['Position1'],
                                   normalized['Position2'], self.vcf_in)
        left_df = pt.tabix_results_to_df(left_snps, "L", "Left SNP Count")
        right = pt.primer_range_right(self.pvcf['Sequence ID'], self.pvcf['Primer Rank'],
                                      self.pvcf['Chromosome'], self.pvcf['Primer Right Seq'],
                                      self.pvcf['Position2'])
        normalized = pt.match_pinfo_to_vcf(right, self.vcf_in)
        right_snps = pt.tabix_fetch(normalized['Sequence ID'], normalized['Primer Rank'],
                                    normalized['Chromosome'], normalized['Position1'],
                                    normalized['Position2'], self.vcf_in)
        right_df = pt.tabix_results_to_df(left_snps, "R", "Right SNP Count")
        merged_df = pt.merge_left_right(left_df, right_df, self.pvcf)
        self.assertTrue('Left SNP Count' and 'Right SNP Count' in merged_df.columns)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
