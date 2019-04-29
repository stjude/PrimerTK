#!/usr/bin/env python
"""
Modules required to run program.
Dependencies:
 - python3.6+
    - pandas>=0.22.0
    - numpy>=1.14.2
    - biopython>=1.70
    - argparse
"""

import sys
import time
import argparse

from Bio import SeqIO
import pandas as pd
from primer_tk import sequence_info as seqinf

def add_iterator_subparser(subparser):
    """ Get commandline arguments
    Args:
        subparser (?): Argparse subparsers

    Returns: None
    """

    parser = subparser.add_parser("iterator_sv", help="iterator_sv subparser")

    parser.add_argument("-ref", "--ref_genome", required=True,
                        help="Reference Genome File to design primers around")
    parser.add_argument("-in", "--regions_file", required=True,
                        help="File with regions to design primers around")
    parser.add_argument("-opt_size", "--primer_opt_size",
                        dest="primer_opt_size", default="22",
                        help="The optimum primer size for output, default: 22")
    parser.add_argument("-min_size", "--primer_min_size",
                        dest="primer_min_size", default="18",
                        help="The optimum primer size for output, default: 18")
    parser.add_argument("-max_size", "--primer_max_size",
                        dest="primer_max_size", default="25",
                        help="The optimum primer size for output, default: 25")
    parser.add_argument("-opt_gc", "--primer_opt_gc",
                        dest="primer_opt_gc", default="50",
                        help="Optimum primer GC, default: 50")
    parser.add_argument("-min_gc", "--primer_min_gc",
                        dest="primer_min_gc", default="20",
                        help="Minimum primer GC, default: 20")
    parser.add_argument("-max_gc", "--primer_max_gc",
                        dest="primer_max_gc", default="80",
                        help="Maximum primer GC, default: 80")
    parser.add_argument("-opt_tm", "--primer_opt_tm",
                        dest="primer_opt_tm", default="60",
                        help="Optimum primer TM, default: 60")
    parser.add_argument("-min_tm", "--primer_min_tm",
                        dest="primer_min_tm", default="57",
                        help="minimum primer TM, default: 57")
    parser.add_argument("-max_tm", "--primer_max_tm",
                        dest="primer_max_tm", default="63",
                        help="maximum primer TM, default: 63")
    parser.add_argument("-sr", "--product_size_range",
                        dest="product_size_range", default="200-400",
                        help="Size Range for PCR Product, default=200-400")
    parser.add_argument("-flank", "--flanking_region_size",
                        dest="flanking_region_size", default="200",
                        help="This value will select how many bases up and downstream to count\
                              when flanking SNP (will do 200 up and 200 down), default: 200")
    parser.add_argument("-st", "--sequence_target",
                        dest="sequence_target", default="199,1",
                        help="default: 199,1, should be half of your flanking region size,\
                              so SNP/V will be included.")
    parser.add_argument("-mp", "--mispriming", dest="mispriming",
                        help="full path to mispriming library for primer3\
                              (/home/dkennetz/testing_p3/primers/humrep.ref")
    parser.add_argument("-tp", "--thermopath", dest="thermopath",
                        help="full path to thermo parameters for primer3 to use\
                              (/hpcf/apps/primer3/install/2.4.0/src/primer3_config/) install loc")
    parser.add_argument("-sv", "--sv-type", dest="sv",
                        choices=['deletion', 'inversion'],
                        help="currently supported SV primer generation: "
                             "deletion and inversion.")

def genome_iterator(genome):
    """
    Uses Biopython SeqIO to parse a genome Fasta and store each chr and
    sequence as a list.

    Args:
        genome (file): Full reference genome from command-line arguments.

    Returns:
        output (list): genome parsed into list of tuples of (header, seq)
        by chromosome for easy access later.
    """

    output = []
    with open(genome, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            chrm = record.id
            seq = record.seq
            output.append((chrm, seq))
        return output

def create_dataframe_csv(regions_file):
    """
    Creates a pandas DataFrame from a regions file in comma separated form
    and names the columns in the DF according to their values.

    Args:
        regions_file (file): Full path to -in input regions.csv.
        This should have coordinate positions of interest to design primers around.
        Note: the chromosome column can be of format chr1 or simply 1 (chr not necessary).
        The file should contain no headers and should be structured as follows:

        Gene1,Sample1,chr1,pos1,pos2
        Gene2,Sample2,chr2,pos1,pos2
        Gene3,Sample3,chr3,pos1,pos2
        ...
    Returns:
        regions_df (pd.DataFrame): The infile parsed to pd.DataFrame object.
    """

    regions_df = pd.read_csv(regions_file, header=None)
    regions_df.columns = ['Gene', 'Sample', 'Chr', 'PosStart', 'PosStop']
    regions_df = regions_df.astype({'Chr': str})
    return regions_df

def create_dataframe_txt(regions_file):
    """
    Creates a pandas DataFrame from a regions file in a tab separated form
    and uses the following column names in the DF according to their values.

    Args:
        regions_file (file): Full path to -in input regions.txt.
        This should have coordinate positions of interest to design primers around.
        Note: the chromosome column can be of format chr1 or simply 1 (chr not necessary).
        The file should contain no headers and should be constructed as follows:

        Gene1\tSample1\tchr1\tpos1\tpos2
        Gene2\tSample2\tchr2\tpos1\tpos2
        Gene3\tSample3\tchr3\tpos1\tpos2
        ...
    Returns:
        regions_df (pd.DataFrame): The infile parsed to pd.DataFrame object.
    """
    regions_df = pd.read_table(regions_file, header=None)
    regions_df.columns = ['Gene', 'Sample', 'Chr', 'PosStart', 'PosStop']
    regions_df = regions_df.astype({'Chr': str})
    return regions_df

def file_extension(infile):
    """
    Runs pd.read_ function to create dataframe corresponding
    to file extension, .txt and .csv accepted.
    This will call the function to create a df only if file extension is correct.
    If the file extension is wrong, it will return the error message and exit.

    Args:
        infile (file): Full path to -in input regions.txt.
    Returns:
        small_regions (pd.DataFrame): The infile parsed to pd.DataFrame object
    """
    if infile.lower().endswith('.txt'):
        small_regions = create_dataframe_txt(infile)
    elif infile.lower().endswith('.csv'):
        small_regions = create_dataframe_csv(infile)
    else:
        sys.exit("Wrong File Format, should be .txt or .csv")
    return small_regions

def match_chr_to_genome(dataframe, genome):
    """
    Used to match formatting between regionsFile input and genome for chromosome
    denotations. Some genomes contain the string "chr" and some do not.

    Args:
        dataframe (pandas object): dataframe created from create_dataframe() function.
        genome (list): genome list of tuples created from genome_iterator() function.

    Returns:
        dataframe (pandas object): dataframe formatted to match genome annotation of chromosome.
    """
    if dataframe['Chr'].str.contains("chr").any() and "chr" in str(genome[0][0]):
        pass
    elif dataframe['Chr'].str.contains("chr").any() and "chr" not in str(genome[0][0]):
        dataframe['Chr'] = dataframe['Chr'].str.replace('chr', '')
    elif not dataframe['Chr'].str.contains("chr").any() and "chr" not in str(genome[0][0]):
        pass
    elif not dataframe['Chr'].str.contains("chr").any() and "chr" in str(genome[0][0]):
        dataframe['Chr'] = 'chr' + dataframe['Chr']
    else:
        print("look at pandas error")
    return dataframe

def flanking_regions_fasta_deletion(genome, dataframe, flanking_region_size):
    """
    Makes batch processing possible, pulls down small region
    of genome for which to design primers around.
    This is based on the chromosome and position of input file.
    Each Fasta record  will contain:

    >Sample_Gene_chr:posStart-posStop
    Seq of flanking region upstream of SV + seq of flanking region downstream of SV

    Args:
        genome (list): genome list of tuples (header, seq).
        dataframe (pandas object): dataframe with sample info.
        flanking_region_size (int): length of sequence upstream and downstream of
        input coordinate position to pull as sequence to design primers around.
    """
    output = []
    for headers, seqs in genome:
        chrm = str(headers)
        seq = str(seqs)
        for gene, sample, chrom, start, stop in zip(dataframe.Gene, dataframe.Sample, dataframe.Chr,
                                                    dataframe.PosStart, dataframe.PosStop):
            if str(chrom) == chrm:
                header = str(str(sample)+"_"+str(gene)+"_"+\
                                 str(chrom)+":"+str(start)+"-"+str(stop)+"__")
                flank_seq = seq[int(start)-int(flanking_region_size):int(start)+1]\
                    +seq[int(stop):(int(stop)+int(flanking_region_size))]
                output.append((header, flank_seq.upper()))
    return output

def flanking_regions_fasta_inversion(genome, dataframe, flanking_region_size):
    """
    Makes batch processing possible, pulls down small region
    of genome for which to design primers around and generates
    flanking regions based on an inverted sequence.

    This is based on the chromosome and position of input file.
    Each Fasta record  will contain:

    >Sample_Gene_chr:posStart-posStop
    Seq of flanking region upstream of SV + seq of flanking region downstream of SV

    Args:
        genome (list): genome list of tuples (header, seq).
        dataframe (pandas object): dataframe with sample info.
        flanking_region_size (int): length of sequence upstream and downstream of
        input coordinate position to pull as sequence to design primers around.
    """
    output = []
    for headers, seqs in genome:
        chrm = str(headers)
        seq = str(seqs)
        for gene, sample, chrom, start, stop in zip(dataframe.Gene, dataframe.Sample, dataframe.Chr,
                                                    dataframe.PosStart, dataframe.PosStop):
            if str(chrom) == chrm:
                header = str(str(sample)+"_"+str(gene)+"_"+\
                                 str(chrom)+":"+str(start)+"-"+str(stop)+"__BP1")
                flank_seq = seq[int(start)-int(flanking_region_size):int(start)+1]\
                    + seqinf.Sequence(seq[int(stop):((int(stop)-int(flanking_region_size))+1):-1]).complement()
                output.append((header, flank_seq))
        for gene, sample, chrom, start, stop in zip(dataframe.Gene, dataframe.Sample, dataframe.Chr,
                                                    dataframe.PosStart, dataframe.PosStop):
            if str(chrom) == chrm:
                header = str(str(sample)+"_"+str(gene)+"_"+\
                                 str(chrom)+":"+str(start)+"-"+str(stop)+"__BP2")
                flank_seq = seqinf.Sequence(seq[int(start)+int(flanking_region_size):int(start)+1:-1]).complement()\
                    + seq[int(stop):(int(stop)+int(flanking_region_size))]
                output.append((header, flank_seq.upper()))
    return output

