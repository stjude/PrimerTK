"""
Genome object contains reference file information and manipulation
Dependencies:
 - biopython >= 1.79
"""

import logging
import sys
from typing import List

from Bio import SeqIO
import pandas as pd

class Genome:
    """
    Parses a Genome build fasta and returns a list of iterators.
    Attributes:
     - fasta: path to the reference fasta file
     - verbosity: logger verbosity, default to INFO
     - dataset: dataset name
    """
    def __init__(self, fasta: str, dataset: str, verbosity="INFO"):
        self.fasta = fasta
        self.dataset = dataset
        self.logger = self.create_logger(verbosity)
        self.logger.debug("Genome object was created.")
        self.reference = []

    def create_logger(self, verbosity: str) -> logging.Logger:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(f'{self.dataset}.log', mode='w')
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter('%(asctime)s %(process)d-%(levelname)s: %(message)s')
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)
        return logger

    def genome_iterator(self) -> None:
        """ Returns a list of iterators of a genome."""
        self.logger.info(f'Parsing reference {self.fasta}')
        try:
            with open(self.fasta, 'r') as fasta:
                for record in SeqIO.parse(fasta, 'fasta'):
                    self.reference.append((record.id, record.seq))
        except FileNotFoundError:
            self.logger.error(f'{self.fasta} does not exist or you do not have permission to access it.')
            sys.exit(1)
    
    def create_flanking_regions_fasta(self, regions_file: str, flanking_region_size: int) -> List[tuple]:
        """ Creates Fasta with flanking regions for batch primer processing 
        Args:
            regions_file: input file to design flanking seqs around
            flanking_region_size: number of bases to pad on each side of region
        Returns:
            fasta_seqs: list of tuple of (header, seq) 
        """
        fasta_seqs = []
        regions_df = parse_input(regions_file)
        self.logger.info("Attempting to create flanking seqs from input")
        match_chr_to_genome(regions_df, self.reference)
        for headers, seqs in self.reference:
            chrm = str(headers)
            seq = str(seqs)
            for sample, tag, chrom, pos in zip(regions_df['Sample'], regions_df['Tag'], regions_df['Chr'], regions_df['Pos']):
                if str(chrom) == chrm:
                    header = f"{sample}_{tag}_{chrom}:{pos}"
                    flank_seq = seq[int(pos)-int(flanking_region_size):int(pos)+int(flanking_region_size)]
                    fasta_seqs.append((header, flank_seq.upper()))
        self.logger.info("Flanking sequences created!")
        return fasta_seqs

def parse_input(regions_file: str) -> pd.DataFrame:
    """ Parses regions input file and returns dataframe """
    with open(regions_file) as f:
        header = f.readline()
        try:
            assert len(header.split(',')) ==  4
        except AssertionError:
            sys.exit("input file should be comma separated with 4 columns, Sample,Tag,Chr,Pos")
        split_header = header.rstrip().split(',')
        assert all(any(i in j for j in ["Sample", "Tag", "Chr", "Pos"]) for i in split_header), "Input file should contain header: Sample,Tag,Chr,Pos"

    regions_df = pd.read_csv(regions_file, header=0, dtype={'Chr': str})
    return regions_df

def match_chr_to_genome(dataframe: pd.DataFrame, reference: List) -> pd.DataFrame:
    """ Matches formatting between regions file and reference genome for chromosome nomenclature
    Designed to make life a bit easier on the user.
    Args:
        dataframe: dataframe created from input regions file
        reference: parsed reference list
    """
    if dataframe['Chr'].str.contains("chr").any() and "chr" not in str(reference[0][0]):
        dataframe['Chr'] = dataframe['Chr'].str.replace('chr', '')
    elif not dataframe['Chr'].str.contains("chr").any() and "chr" in str(reference[0][0]):
        dataframe['Chr'] = 'chr' + dataframe['Chr']
