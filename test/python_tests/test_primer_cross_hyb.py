#!/usr/bin/env python3

"""
Modules required for program.
 - python3.6+
    - pandas>=0.22.0
    - argparse
    - primer_cross_hyb
    - sequence_info
"""

import unittest
from io import StringIO
import argparse
from primer_tk.mp_class import MissingPrimers
from primer_tk.mp_class import create_df
from primer_tk import primer_cross_hyb as pch

class TestPrimerCrossHyb(unittest.TestCase):
    """
    Subclass of unittest to test primer_cross_hyb.py
    """
    def setUp(self):
        """
        Setup for tests, need to generate df as all scripts in this operate
        on a dataframe. self.primer_df is called twice to drop NA.
        """
        self.dump = "test/data/primer_dump_standard.txt"
        prim_list_0 = MissingPrimers(self.dump, 0).samp_primer_info
        prim_list_1 = MissingPrimers(self.dump, 1).samp_primer_info
        prim_list_2 = MissingPrimers(self.dump, 2).samp_primer_info
        prim_list_3 = MissingPrimers(self.dump, 3).samp_primer_info
        prim_list_4 = MissingPrimers(self.dump, 4).samp_primer_info
        self.primer_df = create_df([prim_list_0, prim_list_1, prim_list_2,
                                    prim_list_3, prim_list_4])
        self.primer_df = self.primer_df.loc[~(self.primer_df['Primer Left Seq'] == 'NA')]
        self.pa = 50
        self.fp_len = pch.get_fprimer_percent_aln(self.primer_df['Primer Left Seq'], self.pa)
        self.pd_compare_1_2 =  pch.primer_dimer_local(self.fp_len, self.primer_df['Sequence ID'],
                                                      self.primer_df['Primer Left Seq'],
                                                      self.primer_df['Primer Right Seq'])



    def test_get_fprimer_percent_aln(self):
        """
        Tests to ensure that the proper forward primer length is returned for
        the given percent alignment.
        primer lengths input = [22, 22, 22, 22, 21, 22, 24, 22, 22, 22, 22, 24]
        Should return = [11, 11, 11, 11, 10, 11, 12, 11, 11, 11, 11, 12]
        """
        self.assertEqual(self.fp_len, [11, 11, 11, 11, 10, 11, 12, 11, 11, 11, 11, 12])


    def test_primer_dimer_local(self):
        """
        For percent alignment == 50, this set of primers should contain 17 dimers in
        forward primer vs reverse primer comparison.
        """
        self.assertEqual(len(list(self.pd_compare_1_2)), 17)

    def test_list_from_gen(self):
        """
        Tests to ensure list is generated from primer_dimer_local output.
        """
        self.assertTrue(isinstance(pch.list_from_gen(self.pd_compare_1_2), list))

    def test_p_list_formatter(self):
        """
        Tests to ensure that output list is formatted properly.
        """
        pd_list = pch.list_from_gen(self.pd_compare_1_2)
        formatted = pch.p_list_formatter(pd_list)
        self.assertFalse(any('\\n' in x for x in formatted))

    def test_dimer_true(self):
        """
        Tests to ensure that exactly 3 dimers are found between forward and reverse primers.
        """
        pd_list = pch.list_from_gen(self.pd_compare_1_2)
        formatted = pch.p_list_formatter(pd_list)
        self.primer_df['Dimers1_2F'] = pch.dimer_true(self.primer_df, 2, formatted)
        self.assertEqual((self.primer_df['Dimers1_2F'][self.primer_df['Dimers1_2F'] == True].count()), 3)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
