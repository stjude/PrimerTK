#!/usr/bin/env python3

""" Modules required to run program.

Dependencies:
 - python3.6+
    - biopython>=1.70
"""

from Bio import pairwise2

class Sequence:
    """
    Creates a sequence object from an input DNA sequence.
    """

    def __init__(self, sequence):
        """
        Initialize seq object.
        """
        self.sequence = str(sequence).upper()
        self.name = "seq_name_not_set"

    def set_name(self, name):
        """
        Sets name of sequence object.
        """
        self.name = name

    def complement(self):
        """
        Returns complement (defined in Global vars).
        Example: ATCG -> TAGC
        """
        return ''.join(map(complement_base, self.sequence))

    def reverse_complement(self):
        """
        Returns the reversed complement of the sequence.
        Example: ATCG -> TAGC -> CGAT
        """
        return ''.join(map(complement_base, self.sequence[::-1]))

    def reverse_sequence(self):
        """
        Returns the reverse of the input sequence.
        Example: ATCG -> GCTA
        """
        return self.sequence[::-1]

    def gc_percent(self):
        """
        Calculates GC as a percentage of the input sequence.
        Example: ATCG = 50
        """
        return float((self.sequence.count('G') +\
                          self.sequence.count('C')) / float(len(self.sequence))*100)


class PrimerDimer(Sequence):
    """
    Identifies dimerization between 2 sequences. Essentially,
    it identifies high degrees of complementarity between
    sequences which could lead to undesired binding behavior
    during PCR amplification. Primers bind to each other rather
    than target DNA strand if in pool.
    Example: ATCG + TAGC (high degree of complementarity).
    """
    def __init__(self, sequence, compare_sequence, pd_length):
        super(PrimerDimer, self).__init__(sequence)
        self.compare_sequence = compare_sequence
        self.pd_length = pd_length

    def pd_local(self):
        """
        Compares sequence and comparesequence for the highest number of
        identical characters within a single alignment of the two strings.
        Output is a list with:
        [sequence,
        complement(comparesequence),
        score,
        beginning-of-alignment,
        end-of-alignment]
        Use of the static method 'format_alignment_compl' to format the
        output is highly recommended to display the alignments.
        """
        for alignment in pairwise2.align.localxs(self.sequence,
                                                 self.compare_sequence,
                                                 -len(self.sequence),
                                                 -len(self.sequence)):
            if int(alignment[2]) >= self.pd_length:
                yield alignment

    @staticmethod
    def format_alignment_compl(align1, align2, score, begin, end):
        """
        A modified version of biopython.pairwise2.format_alignment.
        Instead of returning the align2 argument, the complement of align2
        is returned.
        For PrimerDimers we are trying to find complementary base pairings
        so this scores based on complementarity, whereas pairwise2 scores on
        literal matches.
        """
        seq = []
        compalign2 = Sequence(align2).complement()

        seq.append("%s\n" % align1)
        seq.append("%s%s\n" % (" "*begin, "|"*(end-begin)))
        seq.append("%s\n" % compalign2)
        seq.append("  Score=%g\n" % score)
        return ''.join(seq)

# GLOBAL VARS
_COMPLEMENTS = {'-': '-', 'N': 'N', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def is_complement(base1, base2):
    """
    Used to identify is base2 is a complement of base1
    """
    return _COMPLEMENTS[base1] == base2

def complement_base(base):
    """
    Returns the complement for a base.
    """
    if base in _COMPLEMENTS:
        return _COMPLEMENTS[base]
    print("base not in complements, something is wrong!")
