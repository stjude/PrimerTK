"""
Unittest for Genome
"""

import os
import pytest

import pandas as pd


from primertk import Genome

@pytest.fixture
def local_dir() -> str:
    return os.path.dirname(os.path.realpath(__file__))

@pytest.fixture
def fasta(local_dir: str) -> Genome.Fasta:
    ''' Returns a Fasta instance with loaded fasta '''
    return Genome.Fasta.from_fasta(os.path.join(local_dir, './data/test_standard.fa'))

@pytest.fixture
def regions_file(local_dir: str) -> str:
    local_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(local_dir, "./data/test_input_standard.csv")

@pytest.fixture
def header_fail(local_dir: str) -> str:
    return os.path.join(local_dir, "./data/test_no_header.csv")

def test_from_fasta(fasta):
    assert len(fasta.reference) == 2
    assert '1' in fasta.reference
    assert '2' in fasta.reference

def test_parse_input(regions_file):
    assert isinstance(Genome.parse_input(regions_file), pd.DataFrame)

def test_header_fail(header_fail):
    with pytest.raises(AssertionError):
        Genome.parse_input(header_fail)

def test_match_chr_to_genome(fasta, regions_file):
    df = Genome.parse_input(regions_file)
    Genome.match_chr_to_genome(df, fasta.reference)
    assert df['Chr'].str.contains("chr").any() == False

def test_gets_proper_index(fasta, regions_file):
    fasta_seqs = fasta.create_flanking_regions(regions_file, 0)
    assert fasta_seqs[0][1] == 'T'
    assert fasta_seqs[1][1] == 'G'
    assert fasta_seqs[2][1] == 'C'

def test_gets_proper_flanking(fasta, regions_file):
    """ Now we want to ensure flanking regions pulled correctly"""
    fasta_seqs = fasta.create_flanking_regions(regions_file, 10)
    assert len(fasta_seqs[0][1]) == 21
    assert fasta_seqs[0][1] == 'ACCCTAACCCTAACCCTAACC'
