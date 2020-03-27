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
    parser.add_argument("-mp", "--mispriming", dest="mispriming", required=True,
                        help="full path to mispriming library for primer3\
                              (EX: /home/dkennetz/mispriming/humrep.ref")
    parser.add_argument("-tp", "--thermopath", dest="thermopath", required=True,
                        help="full path to thermo parameters for primer3 to use\
                              (EX: /home/dkennetz/primer3/src/primer3_config/) install loc")
    parser.add_argument("-sv", "--sv-type", dest="sv",
                        choices=['deletion', 'inversion', 'insertion', 'translocation'], required=True,
                        help="currently supported SV primer generation: "
                             "deletion, inversion, insertion and translocation.")

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

def create_dataframe_insertion_csv(regions_file):
    """
    Creates a pandas DataFrame from a regions file in comma separated form
    and names the columns in the DF according to their values.

    Args:
        regions_file (file): Full path to -in input regions.csv.
        This should have coordinate positions of interest to design primers around.
        Note: the chromosome column can be of format chr1 or simply 1 (chr not necessary).
        The file should contain no headers and should be structured as follows:

        Gene1,Sample1,chrNorm1,posNorm1,posNorm2,strand,chrIns2,posIns1,posIns2,strand
        Gene2,Sample2,chrNorm2,posNorm1,posNorm2,strand,chrIns3,posIns1,posIns2,strand
        Gene3,Sample3,chrNorm3,posNorm1,posNorm2,strand,chrIns4,posIns1,posIns2,strand
        ...
    Returns:
        regions_df (pd.DataFrame): The infile parsed to pd.DataFrame object.
    """

    regions_df = pd.read_csv(regions_file, header=None)
    regions_df.columns = ['Gene', 'Sample', 'ChrNorm', 'PosNorm1', 'PosNorm2', 'StrandN',
                          'ChrIns', 'PosIns1', 'PosIns2', 'StrandI']
    regions_df = regions_df.astype({'ChrNorm': str})
    regions_df = regions_df.astype({'ChrIns': str})
    return regions_df

def create_dataframe_translocation_csv(regions_file):
    """
    Creates a pandas DataFrame from a regions file in csv form
    and uses column names in df according tot heir values.
    
    Args:
        regions_file (file): input regions.csv file.
        Note: the chromosome column can be of format chr1 or simply 1 (chr not necessary).
        The file should contain no headers and should be constructed as follows:

        Gene1,Sample1,chrNorm,posNorm,strand,chrTrans,posTrans,strand
        Gene2,Sample2,chrNorm,posNorm,strand,chrTrans,posTrans,strand
        Gene3,Sample3,chrNorm,posNorm,strand,chrTrans,posTrans,strand
        ...
    Returns:
        regions_df (pd.DataFrame): the infile parsed to pd.DataFrame object.
    """
    regions_df = pd.read_csv(regions_file, header=None)
    regions_df.columns = ['Gene', 'Sample', 'ChrNorm', 'PosNorm', 'StrandN',
                          'ChrTrans', 'PosTrans', 'StrandT']
    regions_df = regions_df.astype({'ChrNorm': str})
    regions_df = regions_df.astype({'ChrTrans': str})
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
    regions_df = pd.read_csv(regions_file, sep='\t', header=None)
    regions_df.columns = ['Gene', 'Sample', 'Chr', 'PosStart', 'PosStop']
    regions_df = regions_df.astype({'Chr': str})
    return regions_df

def create_dataframe_insertion_txt(regions_file):
    """
    Creates a pandas DataFrame from a regions file in a tab separated form
    and uses the following column names in the DF according to their values.

    Args:
        regions_file (file): Full path to -in input regions.txt.
        This should have coordinate positions of interest to design primers around.
        Note: the chromosome column can be of format chr1 or simply 1 (chr not necessary).
        The file should contain no headers and should be constructed as follows:

        Gene1\tSample1\tchrNorm1\tposNorm1\tposNorm2\tstrand\tchrIns2\tposIns1\tposIns2\tstrand
        Gene2\tSample2\tchrNorm2\tposNorm1\tposNorm2\tstrand\tchrIns3\tposIns1\tposIns2\tstrand
        Gene3\tSample3\tchrNorm3\tposNorm1\tposNorm2\tstrand\tchrIns4\tposIns1\tposIns2\tstrand
        ...
    Returns:
        regions_df (pd.DataFrame): The infile parsed to pd.DataFrame object.
    """
    regions_df = pd.read_csv(regions_file, sep='\t', header=None)
    regions_df.columns = ['Gene', 'Sample', 'ChrNorm', 'PosNorm1', 'PosNorm2', 'StrandN',
                          'ChrIns', 'PosIns1', 'PosIns2', 'StrandI']
    regions_df = regions_df.astype({'ChrNorm': str})
    regions_df = regions_df.astype({'ChrIns': str})
    return regions_df

def create_dataframe_translocation_txt(regions_file):
    """
    Creates a pandas DataFrame from a regions file in txt form
    and uses column names in df according tot heir values.
    
    Args:
        regions_file (file): input regions.csv file.
        Note: the chromosome column can be of format chr1 or simply 1 (chr not necessary).
        The file should contain no headers and should be constructed as follows:

        Gene1\tSample1\tchrNorm\tposNorm\tstrand\tchrTrans\tposTrans\tstrand
        Gene2\tSample2\tchrNorm\tposNorm\tstrand\tchrTrans\tposTrans\tstrand
        Gene3\tSample3\tchrNorm\tposNorm\tstrand\tchrTrans\tposTrans\tstrand
        ...
    Returns:
        regions_df (pd.DataFrame): the infile parsed to pd.DataFrame object.
    """
    regions_df = pd.read_csv(regions_file, sep='\t', header=None)
    regions_df.columns = ['Gene', 'Sample', 'ChrNorm', 'PosNorm', 'StrandN',
                          'ChrTrans', 'PosTrans', 'StrandT']
    regions_df = regions_df.astype({'ChrNorm': str})
    regions_df = regions_df.astype({'ChrTrans': str})
    return regions_df

def file_extension(infile, strvar):
    """
    Runs pd.read_ function to create dataframe corresponding
    to file extension, .txt and .csv accepted.
    This will call the function to create a df only if file extension is correct.
    If the file extension is wrong, it will return the error message and exit.

    Args:
        infile (file): Full path to -in input regions.txt.
        strvar (string): the structural variant type
    Returns:
        small_regions (pd.DataFrame): The infile parsed to pd.DataFrame object
    """
    if infile.lower().endswith('.txt') and strvar in ('deletion', 'inversion'):
        small_regions = create_dataframe_txt(infile)
    elif infile.lower().endswith('.csv') and strvar in ('deletion', 'inversion'):
        small_regions = create_dataframe_csv(infile)
    elif infile.lower().endswith('.txt') and strvar == 'insertion':
        small_regions = create_dataframe_insertion_txt(infile)
    elif infile.lower().endswith('.csv') and strvar == 'insertion':
        small_regions = create_dataframe_insertion_csv(infile)
    elif infile.lower().endswith('.csv') and strvar == 'translocation':
        small_regions = create_dataframe_translocation_csv(infile)
    elif infile.lower().endswith('.txt') and strvar == 'translocation':
        small_regions = create_dataframe_translocation_txt(infile)
    else:
        sys.exit("Wrong File Format, should be .txt (tab) or .csv (comma), or check sv type.")
    return small_regions

def match_chr_to_genome(dataframe, genome, strvar):
    """
    Used to match formatting between regionsFile input and genome for chromosome
    denotations. Some genomes contain the string "chr" and some do not.

    Args:
        dataframe (pandas object): dataframe created from create_dataframe() function.
        genome (list): genome list of tuples created from genome_iterator() function.

    Returns:
        dataframe (pandas object): dataframe formatted to match genome annotation of chromosome.
    """
    if strvar in ('deletion', 'inversion'):
        if dataframe['Chr'].str.contains("chr").any() and "chr" not in str(genome[0][0]):
            dataframe['Chr'] = dataframe['Chr'].str.replace('chr', '')
        elif not dataframe['Chr'].str.contains("chr").any() and "chr" in str(genome[0][0]):
            dataframe['Chr'] = 'chr' + dataframe['Chr']
    elif strvar == 'insertion':
        if dataframe['ChrNorm'].str.contains("chr").any() and "chr" not in str(genome[0][0]):
            dataframe['ChrNorm'] = dataframe['ChrNorm'].str.replace('chr', '')
            dataframe['ChrIns'] = dataframe['ChrIns'].str.replace('chr', '')
        elif not dataframe['ChrNorm'].str.contains("chr").any() and "chr" in str(genome[0][0]):
            dataframe['ChrNorm'] = 'chr' + dataframe['ChrNorm']
            dataframe['ChrIns'] = 'chr' + dataframe['ChrIns']
    elif strvar == 'translocation':
        if dataframe['ChrNorm'].str.contains("chr").any() and "chr" not in str(genome[0][0]):
            dataframe['ChrNorm'] = dataframe['ChrNorm'].str.replace('chr', '')
            dataframe['ChrTrans'] = dataframe['ChrTrans'].str.replace('chr', '')
        elif not dataframe['ChrNorm'].str.contains("chr").any() and "chr" in str(genome[0][0]):
            dataframe['ChrNorm'] = 'chr' + dataframe['ChrNorm']
            dataframe['ChrTrans'] = 'chr' + dataframe['ChrIns']
    else:
        sys.exit("Wrong SV type, please select an SV in help menu.")
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

    >Sample_Gene_chr:posStart-posStop_BP
    Seq of flanking region upstream of SV + seq of flanking region downstream of SV

    Args:
        genome (list): genome list of tuples (header, seq).
        dataframe (pandas object): dataframe with sample info.
        flanking_region_size (int): length of sequence upstream and downstream of
        input coordinate position to pull as sequence to design primers around.
    Returns:
        output (list): header + seq
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
                    + seqinf.\
                    Sequence(seq[int(stop):((int(stop)-(int(flanking_region_size)+1))):-1])\
                    .complement()
                output.append((header, flank_seq))
        for gene, sample, chrom, start, stop in zip(dataframe.Gene, dataframe.Sample, dataframe.Chr,
                                                    dataframe.PosStart, dataframe.PosStop):
            if str(chrom) == chrm:
                header = str(str(sample)+"_"+str(gene)+"_"+\
                                 str(chrom)+":"+str(start)+"-"+str(stop)+"__BP2")
                flank_seq = seqinf.\
                    Sequence(seq[int(start)+int(flanking_region_size):int(start)+1:-1])\
                    .complement()\
                    + seq[int(stop):(int(stop)+int(flanking_region_size))]
                output.append((header, flank_seq.upper()))
    return output

def flanking_region_fasta_translocation(genome, dataframe, flanking_region_size):
    """
    Pulls down small region of genome for which to design primers around and
    generates flanking regions based on strand info from input file.

    Each Fasta record will contain:
    >Sample_Gene_chrNorm:posNorm-posTrans
    Seq of flanking region upstream of posNorm + seq after posTrans based on strand
    
    Args:
        genome (list): genome list of tuples (header, seq)
        dataframe (pd.DataFrame): dataframe with sample info.
        flanking_region_size  (int): length of sequence to pad position with.
    Returns:
        result (dict): {header: seq}
    """
    output = []
    samp_norm = {}
    samp_tran = {}
    for headers, seqs in genome:
        chrm = str(headers)
        seq = str(seqs)
        for gene, sample, chrn, posn, strandn, chrt, post, strandt  in zip(dataframe.Gene, dataframe.Sample,
                                                                           dataframe.ChrNorm, dataframe.PosNorm,
                                                                           dataframe.StrandN, dataframe.ChrTrans,
                                                                           dataframe.PosTrans, dataframe.StrandT):
            if str(chrn) == chrm and strandn == '+':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+str(posn)+"-"+str(post))
                flank_seq = seq[int(posn):int(posn)+int(flanking_region_size)]
                samp_norm[header] = flank_seq
            elif str(chrn) == chrm and strandn == '-':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+str(posn)+"-"+str(post))
                flank_seq = seqinf.Sequence(seq[int(posn):int(posn)-(int(flanking_region_size)):-1])\
                            .complement()
                samp_norm[header] = flank_seq

        for gene, sample, chrn, posn, strandn, chrt, post, strandt  in zip(dataframe.Gene, dataframe.Sample,
                                                                           dataframe.ChrNorm, dataframe.PosNorm,
                                                                           dataframe.StrandN, dataframe.ChrTrans,
                                                                           dataframe.PosTrans, dataframe.StrandT):
            if str(chrt) == chrm and strandt == '+':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+str(posn)+"-"+str(post))
                flank_seq = seq[int(post)-int(flanking_region_size):int(post)]
                samp_tran[header] = flank_seq
            elif str(chrt) == chrm and strandt == '-':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+str(posn)+"-"+str(post))
                flank_seq = seqinf.Sequence(seq[int(post)+int(flanking_region_size):int(post):-1])\
                            .complement()
                samp_tran[header] = flank_seq

    result = {}
    for name1, seq1 in samp_norm.items():
        for name2, seq2 in samp_tran.items():
            if name1 == name2:
                outseq = seq2+seq1
                result[name1] = outseq

    return result
                
def flanking_region_fasta_insertion(genome, dataframe, flanking_region_size):
    """
    Makes batch processing possible, pulls down small region
    of genome for which to design primers around and generates
    flanking regions based on an inverted sequence.

    This is based on the chromosome and position of input file.
    Each Fasta record  will contain:

    Note: If strand is negative, coordinates should be in decreasing order.

    >Sample_Gene_chr:posNorm1-posNorm2_BP
    Seq of flanking region upstream of SV + seq of inserted sequence based on strand

    Args:
        genome (list): genome list of tuples (header, seq)
        dataframe (pandas object): dataframe with sample info.
        flanking_region_size (int): length of sequence upstream dna downstream of
        input coordinate position to pull as sequence to design primers around.
    Returns:
        result (dict): {header: seq}
    """
    seqsnormbp1 = {}
    seqsinsbp1 = {}
    seqsnormbp2 = {}
    seqsinsbp2 = {}
    for headers, seqs in genome:
        chrm = str(headers)
        seq = str(seqs)
        for gene, sample, chrn, startn, stopn, strandn, stopi in zip(dataframe.Gene, dataframe.Sample,
                                                                     dataframe.ChrNorm, dataframe.PosNorm1,
                                                                     dataframe.PosNorm2, dataframe.StrandN,
                                                                     dataframe.PosIns2):
            if str(chrn) == chrm and strandn == '+':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP1")
                flank_seq = seq[int(startn):int(startn)+int(flanking_region_size)]
                seqsnormbp1[header] = flank_seq
            elif str(chrn) == chrm and strandn == '-':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                             str(startn)+"-"+str(stopi)+"__BP1")
                flank_seq = seqinf.Sequence(seq[int(startn):int(startn)-(int(flanking_region_size)):-1])\
                            .complement()
                seqsnormbp1[header] = flank_seq

        for gene, sample, chrn, startn, chri, starti, stopi, strandi in zip(dataframe.Gene, dataframe.Sample,
                                                                            dataframe.ChrNorm, dataframe.PosNorm1,
                                                                            dataframe.ChrIns, dataframe.PosIns1,
                                                                            dataframe.PosIns2, dataframe.StrandI):
            if str(chri) == chrm and strandi == '+':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP1")
                flank_seq = seq[int(starti)-int(flanking_region_size):int(starti)]
                seqsinsbp1[header] = flank_seq
            elif str(chri) == chrm and strandi == '-':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP1")
                flank_seq = seqinf.Sequence(seq[int(starti)+int(flanking_region_size):int(starti):-1])\
                            .complement()
                seqsinsbp1[header] = flank_seq
        for gene, sample, chrn, startn, stopn, strandn, stopi in zip(dataframe.Gene, dataframe.Sample,
                                                                     dataframe.ChrNorm,
                                                                     dataframe.PosNorm1,
                                                                     dataframe.PosNorm2, dataframe.StrandN,
                                                                     dataframe.PosIns2):
            if str(chrn) == chrm and strandn == '+':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP2")
                flank_seq = seq[int(stopn)-int(flanking_region_size):int(stopn)]
                seqsnormbp2[header] = flank_seq
            elif str(chrn) == chrm and strandn == '-':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP2")
                flank_seq = seqinf.Sequence(seq[int(stopn)+int(flanking_region_size):int(stopn):-1])\
                            .complement()
                seqsnormbp2[header] = flank_seq

        for gene, sample, chrn, startn, chri, starti, stopi, strandi in zip(dataframe.Gene, dataframe.Sample,
                                                                            dataframe.ChrNorm, dataframe.PosNorm1,
                                                                            dataframe.ChrIns, dataframe.PosIns1,
                                                                            dataframe.PosIns2, dataframe.StrandI):
            if str(chri) == chrm and strandi == '+':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP2")
                flank_seq = seq[int(stopi):int(stopi)+(int(flanking_region_size))]
                seqsinsbp2[header] = flank_seq
            elif str(chri) == chrm and strandi == '-':
                header = str(str(sample)+"_"+str(gene)+"_"+str(chrn)+":"+\
                                 str(startn)+"-"+str(stopi)+"__BP2")
                flank_seq = seqinf.\
                    Sequence(seq[int(stopi):int(stopi)-int(flanking_region_size):-1])\
                    .complement()
                seqsinsbp2[header] = flank_seq

    result = {}
    for name1, seq1 in seqsnormbp1.items():
        for name2, seq2 in seqsinsbp1.items():
            if name1 == name2:
                outseq = seq2+seq1
                result[name1] = outseq

    for name1, seq1 in seqsnormbp2.items():
        for name2, seq2 in seqsinsbp2.items():
            if name1 == name2:
                outseq = seq1+seq2
                result[name1] = outseq
    return result
