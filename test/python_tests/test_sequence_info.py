#!/usr/bin/env python3
"""
Unit tester for sequence_info.py

Python version: Python 3.6.8 :: Anaconda, Inc.
Date: 03/22/2019
Dependencies:
 - unittest
 - biopython>=1.70
"""

import unittest
from primer_tk import sequence_info as seqinf

class TestSequence(unittest.TestCase):
    """
    Subclass of unittest to test Sequence Class for sequence_info.py
    """
    def setUp(self):
        """
        Generate a sequence object that is easy to test!
        """
        self.seq = seqinf.Sequence("ACTG")
        self.name = "seq name not seq"
        self.primer_dimer = seqinf.PrimerDimer(self.seq.sequence,
                                               seqinf.Sequence("TGCA").sequence, 2)

    def test_handles_lowercase(self):
        """
        Tests to ensure handles lowercase sequence strings.
        """
        input_seq = self.seq.sequence
        test_seq = seqinf.Sequence('actg').sequence
        self.assertEqual(input_seq, test_seq)

    def test_set_name(self):
        """
        Tests that name is properly changed when function is called.
        """
        self.name = "TestSeq1"
        self.assertEqual(self.name, "TestSeq1")

    def test_complement(self):
        """
        Tests to ensure that sequence complement is returned.
        """
        complement_seq = self.seq.complement()
        self.assertEqual(complement_seq, "TGAC")

    def test_reverse_complement(self):
        """
        Tests to ensure that reverse complement sequence is returned.
        """
        reverse_complement_seq = self.seq.reverse_complement()
        self.assertEqual(reverse_complement_seq, "CAGT")

    def test_reverse_sequence(self):
        """
        Tests to ensure that the sequence is properly reversed.
        """
        reverse_seq = self.seq.reverse_sequence()
        self.assertEqual(reverse_seq, "GTCA")


    def test_gc_percent(self):
        """
        Test to ensure that GC is properly calculated.
        """
        gc_float = self.seq.gc_percent()
        self.assertEqual(gc_float, 50.0)

    def test_0_gc_percent(self):
        """
        Test to ensure that 0 GC doesn't throw error (div0).
        """
        gc_float = seqinf.Sequence("AAAA").gc_percent()
        self.assertEqual(gc_float, 0.0)

    def test_pd_local(self):
        """
        Test to ensure pd_local is aligning exact sequences.
        ACTG & TGCA align at length 2, this is our test.
        Example:
        ACTG
        --||
          TGCA
        """
        local = self.primer_dimer.pd_local()
        self.assertEqual([item[2] for item in local], [2.0])

    def test_format_comp_align(self):
        """
        Test to ensure that complementary alignment is output.
        """
        easy_seq1 = seqinf.Sequence("AAAA").sequence
        easy_seq2 = seqinf.Sequence("TTTT").complement()
        test = seqinf.PrimerDimer(easy_seq1, easy_seq2, 4)
        local = test.pd_local()
        for item in local:
            self.assertTrue("TTTT" in seqinf.PrimerDimer.\
                                format_alignment_compl(item[0],
                                                       item[1],
                                                       item[2],
                                                       item[3],
                                                       item[4]))


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
