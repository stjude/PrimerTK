Requirements:
 - gcc 6+
 - Primer3>=2.4.0
 - isPcr (compiled version with download) will run with gcc 6
 - python3.6+
   - pandas>=0.22.0
   - biopython>=1.70
   - numpy>=1.4
   - cwltool>=1.0.20190228155703
 - perl>=5.10.1

To install python packages, you will see a "requirements.txt" file next to the "Pybox" directory.

For necessary python packages:

    $pip3 install --user -r requirements.txt

You will see some loading and the packages should be installed to your user.

### The Programs

These programs may work programs less than the version specified, but they have been tested with these specific
versions, as well as newer versions of each software and are successful.

PrimerPlex has multiple built-in pipeline options, each designed to make life in the lab easier by saving the user
from manual primer design. PrimerPlex currently has 3 run modes:

1) PrimerStandard
2) PrimerPlex
3) PrimerSV

PrimerStandard (tested experimentally):

  Designed for batch primer generation around a genome.fasta file. The user will provide an input file of samples 
and positions of interest based on a certain reference build. The program will pull down flanking regions and setup
primer3 input, with flexible, user defined input parameters. It will then run primer3, generating up to 5 primers
per position. Useful information will be extracted from each primer pair and carried on to in silico pcr. In 
silico pcr will use the each primer pair as input, and attempt to find locations in the genome which will be 
amplified by the pair. After amplification, the results from PCR will be analyzed and checked for specificity. 
The top ranking primer pairs will be reported for each region. Each of the top ranking primer pairs will be written
out to a final primer output file, as well as a plate ordering format accepted by most vendors. Finally, each
primer pair will be tested for SNPs. If SNPs are found in the primer, the position in the primer will be reported 
(nucleotide) and the snp ID, frequency, and clinical relevance will be reported in a final output sheet. 

PrimerMultiplex (tested experimentally):

  The functionality of PrimerPlex is in the name: it allows for pooling of primers! PrimerPlex has the same initial
functionality of PrimerStandard. It generates regions of interest based on an input file. Prior to in silico pcr,
it compares all primers to all other primers generated and checks for high complementarity. If 2 primers meet a
percent complementarity specified by the user (default is 60%) the primer pair is dropped prior to in silico pcr.
The thought behind this is that highly complemtary primers have a good chance of binding to each other during pcr,
thus lowering their effective concentration and product expression. 
  After checking for cross-hybridization, an all-vs-all pcr input is generated. This is necessary because all
primers will be present in a pool, so any 2 primers in the pool could amplify unintended regions. in silico pcr is
run on the all-vs-all input and products all products are checked for specificty. Only primers matching the initial
input will be retained. The filtering for this is less strict, so some off-target can be expected. If a primer pair
amplifies the intended region, it will be kept. Highest ranking primers will be output, and all primers post-pcr
filtration will be analyzed for snp content.

PrimerSV (deletion tested experimentally):

  PrimerSV currently supports deletions and inversions. Translocations can be easily included; however, they have
currently been difficult for me to classify. When I have real world data to test primer generation, I will include 
translocations in the source code. 

  PrimerSV: deletion currently accepts an input file with sample, chromosome, deletion start and stop. It uses the
start and stop coordinates to parse relevant flanking sequences from the reference genome, and design primers
around the flanking region. Primers are then checked around the flanking regions and predicted product GC is
returned. This has seen high success in experiment, and I would highly recommend this tool for designing primers
around expected deletion events.

  PrimerSV: inversion currently accepts an input file with sample, chromosome, inversion start, stop. It designs a
primer pair around each inversion break point. So for a given sample and inversion, primers are designed around
inversion start and annotated BP1, and inversion stop and annotated BP2. Primers by default are designed around
flanking sequence 200bp upstream of the inversion BP1, and then 200bp of inverted sequence after BP1. Inverted
sequence is designed using BP2 as the start position (3' end of forward strand). It reads 200bp 3'-5' from BP2 and
complements the sequence (effectively making it the second strand). 

**Tutorial:**

Primer3 can be installed with gcc compiler using the following:

    $cd ~
    $git clone https://github.com/primer3-org/primer3.git primer3
    $cd primer3/src
    $make
    $make test

Then user $PATH variable can be adjusted to contain primer3 src:

    $export PATH=/home/<username>/primer3/src:$PATH

So then primer3 can be found from the command line.

Model input files are given in the "tests" directory that is included in this package for each pipeline, but a
model input file will be given for each pipeline below.
It is important that input files to the pipeline DO NOT CONTAIN HEADERS, as custom headers will be added when input
file is parsed.

PrimerStandard and PrimerMultiplex examples:

    input1.csv

        GENE1,SAMPLE1,CHR1,POS1
        GENE2,SAMPLE2,CHR2,POS2
        GENE3,SAMPLE3,CHR3,POS3

    where capitalized letters are replaced by values. A real example for a GRCh37 genome could be:

        HBA2,MyHumanSample1,16,223000
        HBA1,MyHumanSample2,16,227000
        BRCA1,MyHumanSample3,17,41196318
        ...

    The input files also support the format (the use of "chr" or not will be normalized to the reference genome
convention):

        HBA2,MyHumanSample1,chr16,223000
        HBA1,MyHumanSample2,chr16,227000
        BRCA1,MyHumanSample3,chr17,41196318
        ...

    Input files can also be tab delimited, instead of comma:

       HBA2	MyHumanSample1  chr16      223000
       HBA1	MyHumanSample2  chr16      227000
       BRCA1	MyHumanSample3  chr17      41196318
       ...

PrimerSV has a slightly different input file format, as it requires 2 positional columns (SV Start and SV stop):

    input1.csv

	GENE1,SAMPLE1,CHR1,SVSTART,SVSTOP
	GENE2,SAMPLE2,CHR2,SVSTART,SVSTOP
	GENE3,SAMPLE3,CHR3,SVSTART,SVSTOP

    Again, a real example would be:

        HBA2,MyHumanSample1,16,223000,224000
	HBA1,MyHumanSample2,16,227000,230000
	BRCA1,MyHumanSample3,17,41196318,41200000
	...
	
The same formatting conditions apply, so I will not perpetuate the input examples.

**Running the genomeStandard pipelines:**

There are two ways to run each pipeline. The pipelines can be run using the Common Workflow Language (CWL) or
script by script in a stepwise fashion.
I would recommend using CWL, as all it requires to run is the proper installs, and an input.yml file
(an example for each pipeline has been included in install).

I will first present how to run each pipeline step by step, so it is easier to understand the inputs and outputs
and test installation successes. Then I will explain how to run the CWL portion of the code
(much more straightforward).

These examples below will assume that I have reference genome Fasta in my home directory, and the input files are
in the program install directory.
Also, Primer3 is installed and you have added the install location to your path.

PrimerStandard and PrimerPlex:

First, change to the relevant directory:

    $cd ~/PrimerPlex/Pybox/genomeStandard/

**1) genome_iterator.py** (quite a lot of options, I wanted the user to have a lot of control of Primer3 setup).
 
    $python3.6 src/genome_iterator.py -h
    optional arguments:
    -h, --help            show this help message and exit
    -ref REF_GENOME, --ref_genome REF_GENOME
                          Reference Genome File to design primers around
    -in REGIONS_FILE, --regions_file REGIONS_FILE
                          File with regions to design primers around
    -opt_size PRIMER_OPT_SIZE, --primer_opt_size PRIMER_OPT_SIZE
                          The optimum primer size for output, default: 22
    -min_size PRIMER_MIN_SIZE, --primer_min_size PRIMER_MIN_SIZE
                          The optimum primer size for output, default: 18
    -max_size PRIMER_MAX_SIZE, --primer_max_size PRIMER_MAX_SIZE
                          The optimum primer size for output, default: 25
    -opt_gc PRIMER_OPT_GC, --primer_opt_gc PRIMER_OPT_GC
                          Optimum primer GC, default: 50
    -min_gc PRIMER_MIN_GC, --primer_min_gc PRIMER_MIN_GC
                          Minimum primer GC, default: 20
    -max_gc PRIMER_MAX_GC, --primer_max_gc PRIMER_MAX_GC
                          Maximum primer GC, default: 80
    -opt_tm PRIMER_OPT_TM, --primer_opt_tm PRIMER_OPT_TM
                          Optimum primer TM, default: 60
    -min_tm PRIMER_MIN_TM, --primer_min_tm PRIMER_MIN_TM
                          minimum primer TM, default: 57
    -max_tm PRIMER_MAX_TM, --primer_max_tm PRIMER_MAX_TM
                          maximum primer TM, default: 63
    -sr PRODUCT_SIZE_RANGE, --product_size_range PRODUCT_SIZE_RANGE
                          Size Range for PCR Product, default=200-400
    -flank FLANKING_REGION_SIZE, --flanking_region_size FLANKING_REGION_SIZE
                          This value will select how many bases up and
                          downstream to count when flanking SNP (will do 200 up
                          and 200 down), default: 200
    -st SEQUENCE_TARGET, --sequence_target SEQUENCE_TARGET
                          default: 199,1, should be half of your flanking region
                          size, so SNP/V will be included.
    -mp MISPRIMING, --mispriming MISPRIMING
                          full path to mispriming library for primer3
                          (/home/dkennetz/testing_p3/primers/humrep.ref
    -tp THERMOPATH, --thermopath THERMOPATH
                          full path to thermo parameters for primer3 to use
                          (/hpcf/apps/primer3/install/2.4.0/src/primer3_config/)
                          install loc

An actual command line example:

    $python3.6 src/genome_iterator.py -ref /home/dkennetz/references/GRCh37/GRCh37.fasta \
    > -in /home/dkennetz/PrimerPlex/Pybox/genomeStandard/tests/input1.txt \ 
    > -opt_size 22\
    > -min_size 18\
    > -max_size 25\
    > -opt_gc 50\
    > -min_gc 20\
    > -max_gc 80\
    > -opt_tm 60\
    > -min_tm 57\
    > -max_tm 63\
    > -sr 200-400\
    > -flank 200\
    > -st 199,1\
    > -mp /home/dkennetz/PrimerPlex/Pybox/genomeStandard/tests/humrep.ref\
    > -tp /home/dkennetz/Primer3/src/primer3_config/

Running this will produce 2 outputs, a primer3 input file and a fasta with flanking regions based on inputs:

    primer3_input.20190313-104109.txt
    flanking_regions.20190313-104109.fasta

The flanking regions are determined by -flank (the value here and the default are 200) so 200bp were pulled down
from each side of the target and will be used for primer design.

The -st (sequence target parameter) requires primer3 to include your position in the primer design, so -st should
always be (flank-1),1 to ensure your target is amplified by primers.

The GC, TM, and Size values above are defaults. None of them actually need to be included in your command line, but
the option is there incase you want to change them to different values.

The total flanking region for primer design will be (flank*2) so the default size is 400bp.

-mp is a mispriming library. Primer3 will penalize primers around these regions by default and will usually exclude
them.

-tp is the thermodynamic parameter path from primer3. This is a path to a directory which primer3 will
use to score primer pairs.

**2) Primer3:**

    $primer3_core --output=primer_dump.txt primer3_input.20190313-104109.txt

This will output a lot of intermediate files, as well as a primer_dump.txt file. The program will only use the
primer_dump.txt file going forward, which contains all relevant information.

The CWL version of this pipeline will not return the intermediate files by default, but the step by step program
will.

catting primer_dump.txt:

    $cat primer_dump.txt

SEQUENCE_ID=sample1_gene1_1:15000000__
SEQUENCE_TEMPLATE=ACACATAACCATGCCTACAACCCACAGCAAAGCTTCTTGGACACCCAGGTGGTGGTTTCAAAGTGGCCCATTGCACTGAGCCTTTGTTTCCTTCTTGGAAGCAACAGCCGTGGGCTTGTGATAATTTGATTGTTTTGAAATGGGAACAAGTGACTGCCTGAAGAAGGGCATGCATGATGTGAAGGCAGTTTGAAGAATGCCTGAAATTCAGAACGTCAGGGAGGAAGCGTCCTTCGAGATCATCAGGGCCAGCCCCTGCATTTCACAGTAGGAGCAGCTGGAGAAGGAAGACGTGTTGTTGGATCCATCAGCTGCTGTTGGATCCATTGTGCTGCCTCAAGGATGCAGGCTGTGGAAACAGGGAGCTGCAGATGCCTGAGTGACTTTGACAGAAACTTAC
SEQUENCE_TARGET=199,1
PRIMER_FIRST_BASE_INDEX=1
PRIMER_TASK=pick_detection_primers
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_PRODUCT_SIZE_RANGE=200-400
PRIMER_MAX_END_STABILITY=9.0
PRIMER_MAX_SELF_ANY_TH=45.00
PRIMER_MAX_SELF_END_TH=35.00
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00
PRIMER_PAIR_MAX_COMPL_END_TH=35.00
PRIMER_MAX_HAIRPIN_TH=24.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
PRIMER_TM_FORMULA=1
PRIMER_SALT_CORRECTIONS=1
PRIMER_SALT_MONOVALENT=50.0
PRIMER_INTERNAL_SALT_MONOVALENT=50.0
PRIMER_SALT_DIVALENT=1.5
PRIMER_INTERNAL_SALT_DIVALENT=1.5
PRIMER_DNTP_CONC=0.6
PRIMER_INTERNAL_DNTP_CONC=0.6
PRIMER_DNA_CONC=50.0
PRIMER_INTERNAL_DNA_CONC=50.0
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/hpcf/apps/primer3/install/2.4.0/src/primer3_config/
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_MAX_POLY_X=3
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_MIN_TM=57
PRIMER_OPT_TM=60
PRIMER_MAX_TM=63
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_NUM_RETURN=5
P3_FILE_FLAG=1
PRIMER_EXPLAIN_FLAG=1
PRIMER_MISPRIMING_LIBRARY=/home/dkennetz/PrimerPlex/Pybox/genomeStandard/tests/humrep.ref
PRIMER_MIN_GC=20
PRIMER_OPT_GC_PERCENT=50
PRIMER_MAX_GC=80
PRIMER_PAIR_MAX_DIFF_TM=3
PRIMER_LEFT_EXPLAIN=considered 1301, GC content failed 15, low tm 287, high tm 397, high hairpin stability 442, high repeat similarity 4, long poly-x seq 22, ok 134
PRIMER_RIGHT_EXPLAIN=considered 1242, low tm 221, high tm 399, high any compl 1, high end compl 12, high hairpin stability 259, high repeat similarity 22, long poly-x seq 22, high template mispriming score 9, ok 297
PRIMER_INTERNAL_EXPLAIN=considered 3785, GC content failed 20, low tm 720, high tm 1615, high hairpin stability 370, ok 1060
PRIMER_PAIR_EXPLAIN=considered 42, unacceptable product size 30, high mispriming library similarity 5, primer in pair overlaps a primer in a better pair 173, ok 7
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_INTERNAL_NUM_RETURNED=5
PRIMER_PAIR_NUM_RETURNED=5
PRIMER_PAIR_0_PENALTY=0.348775
PRIMER_LEFT_0_PENALTY=0.310349
PRIMER_RIGHT_0_PENALTY=0.038426
PRIMER_INTERNAL_0_PENALTY=0.034616
PRIMER_LEFT_0_SEQUENCE=ATAACCATGCCTACAACCCACA
PRIMER_RIGHT_0_SEQUENCE=TGATGGATCCAACAACACGTCT
PRIMER_INTERNAL_0_SEQUENCE=CCTTCGAGATCATCAGGGCC
PRIMER_LEFT_0=5,22
PRIMER_RIGHT_0=310,22
PRIMER_INTERNAL_0=231,20
PRIMER_LEFT_0_TM=59.690
PRIMER_RIGHT_0_TM=59.962
PRIMER_INTERNAL_0_TM=59.965
PRIMER_LEFT_0_GC_PERCENT=45.455
PRIMER_RIGHT_0_GC_PERCENT=45.455
PRIMER_INTERNAL_0_GC_PERCENT=60.000
PRIMER_INTERNAL_0_SELF_ANY_TH=0.00
PRIMER_LEFT_0_SELF_ANY_TH=0.00
PRIMER_RIGHT_0_SELF_ANY_TH=21.41
PRIMER_INTERNAL_0_SELF_END_TH=0.00
PRIMER_LEFT_0_SELF_END_TH=0.00
PRIMER_RIGHT_0_SELF_END_TH=0.00
PRIMER_LEFT_0_HAIRPIN_TH=0.00
PRIMER_RIGHT_0_HAIRPIN_TH=0.00
PRIMER_INTERNAL_0_HAIRPIN_TH=0.00
PRIMER_LEFT_0_LIBRARY_MISPRIMING=12.00, MER20 Nonautonomous DNA transposon
PRIMER_RIGHT_0_LIBRARY_MISPRIMING=11.00, L1PA7 3'-end of L1 repeat (subfamily L1PA7) - a consensus sequence
PRIMER_PAIR_0_LIBRARY_MISPRIMING=19.00, reverse Tigger2 Autonomous DNA transposon
PRIMER_LEFT_0_END_STABILITY=4.1700
PRIMER_RIGHT_0_END_STABILITY=4.1800
PRIMER_LEFT_0_TEMPLATE_MISPRIMING_TH=18.3343
PRIMER_RIGHT_0_TEMPLATE_MISPRIMING_TH=12.3592
PRIMER_PAIR_0_COMPL_ANY_TH=0.00
PRIMER_PAIR_0_COMPL_END_TH=0.00
PRIMER_PAIR_0_PRODUCT_SIZE=306
PRIMER_PAIR_0_TEMPLATE_MISPRIMING_TH=18.33
...

Will show you all the information for every primer that was designed, including the flanking sequence used, and why
primers may or may not have failed.
Good primer pairs will have a low PRIMER_PAIR_PENALTY (<=3), but I have seen success with primers with penalties up
around even 7. You don't know until you try!

**3) main_pre_pcr.py:** (formats for in silico pcr input and allows for primer multiplexing by checking all primers
vs all others for complementarity, or a simple primer standard setup)!

    $python3.6 src/main_pre_pcr.py -h

    optional arguments:
    -h, --help            show this help message and exit
    -d DUMP, --primer3_dump DUMP
                          Primer3 stdout passed into a 'dump' file to be used as
                          input
    -o OUTFILE, --outfile_name OUTFILE
                          The output filename for all primer information.
    -pa PERCENT_ALIGNMENT, --percent_alignment PERCENT_ALIGNMENT
                          Percent match between 2 primers for pair to be
                          discarded. EX: primer_len = 22, percent_aln = 60
                          dimer_len = (60/100) * 22 = 13.2 -> 13.
    -pcr {standard,multiplex}, --pcr_type {standard,multiplex}
                          perform standard or multiplex pcr on given inputs.


For standard PCR setup:

    $python3.6 src/main_pre_pcr.py -d primer_dump.txt -o total_list.csv -pcr standard

outputs 2 files, standard_pcr.txt and total_list.csv. total_list.csv has all relevant info for downstream analysis.
standard_pcr.txt will be the input to pcr.

For multiplex PCR setup:

    $python3.6 src/main_pre_pcr.py -d primer_dump.txt -o total_list.csv -pa 60 -pcr multiplex

If you examine the files "no_dimer_df.csv" and "total_list.csv", you see that they are the same length.
A quick check:

    $wc -l no_dimer_df.csv
    13
    $wc -l total_list.csv
    13

Interesting. It looks like for sample1, 5 primer pairs were generated, for sample2 only 2 were generated,
for sample3 5 were generated.

Usually, more primers can be generated for a position by adjusting input parameters, but less than
5 primers is usually the result of a difficult region (High GC, tandem repeat).

Lets drop the percent alignment to 50 and see what happens to these file:

     $wc -l no_dimer_df.csv
     7
     $wc -l total_list.csv
     13

Even more interesting! when `-pcr multiplex` is passed as input, a function is run to compare all primers against
all other primers and drop any primers that have more than 50% complementarity with any other primers, since
they will be in a pool. We see the total output by primer3 is 13, but the total output after filtering is 7.
50% complementarity is very strict for a primer pair. The primer will have > 90% complementarity to a genomic
position of interest, so it will drastically favor genomic binding to primer-primer interaction at such a
complementarity difference. However, at percent alignment > 75%, we may start to worry a bit about primer-primer
interactions. To note, we have seen success with primers using -pa at 70. I set 60 as the default to be fairly
strict in filtering.

**4) in silico pcr:**

A pre-compiled distribution of isPcr is distributed with this software package (with citations included). This has
been tested on redhat6 and ubuntu 18.04 and works well!

For isPcr, your reference genome needs to be split by chromosome. A simple tool to handle this is pyfaidx.
To install (assuming python3.6), you can simply use `pip3.6 install --user pyfaidx`.
This will make pyfaidx utilities available from your home directory.

Change directory to your reference genome builds, so for example GRCh37.fasta and run:

    $faidx -x GRCh37.fasta

This will split your fasta into all its parts!

Then, this simple (you can probably write a nicer one) shell script can be used to run isPcr against all
chromosomes:

    #!/usr/bin/env bash

    chr_dir=$1
    pcr_file=$2

    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/1.fa ${pcr_file} plex1.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/2.fa ${pcr_file} plex2.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/3.fa ${pcr_file} plex3.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/4.fa ${pcr_file} plex4.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/5.fa ${pcr_file} plex5.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/6.fa ${pcr_file} plex6.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/7.fa ${pcr_file} plex7.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/8.fa ${pcr_file} plex8.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/9.fa ${pcr_file} plex9.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/10.fa ${pcr_file} plex10.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/11.fa ${pcr_file} plex11.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/12.fa ${pcr_file} plex12.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/13.fa ${pcr_file} plex13.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/14.fa ${pcr_file} plex14.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/15.fa ${pcr_file} plex15.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/16.fa ${pcr_file} plex16.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/17.fa ${pcr_file} plex17.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/18.fa ${pcr_file} plex18.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/19.fa ${pcr_file} plex19.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/20.fa ${pcr_file} plex20.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/21.fa ${pcr_file} plex21.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/22.fa ${pcr_file} plex22.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/X.fa ${pcr_file} plexX.out
    /home/dkennetz/PrimerPlex/Pybox/isPcr/isPcr $chr_dir/Y.fa ${pcr_file} plexY.out

    cat *.out >> pcr_output.txt
    rm -rf *.out

And to run (if my reference is in my home directory and the chromosome fasta are in a subdirectory called
"chromosomes":

For standard pcr setup:

    $./ispcr_run.sh /home/dkennetz/ref_genome/GRCh37/chromosomes standard_pcr.txt

For multiplex pcr setup:

    $./ispcr_run.sh /home/dkennetz/ref_genome/GRCh37/chromosomes multiplex_pcr.txt

It takes about a minute, and as it is running, this will print to screen as each chromosome is processed:

    Loaded 249250621 letters in 1 sequences
    Loaded 243199373 letters in 1 sequences
    Loaded 198022430 letters in 1 sequences
    Loaded 191154276 letters in 1 sequences
    Loaded 180915260 letters in 1 sequences
    Loaded 171115067 letters in 1 sequences
    Loaded 159138663 letters in 1 sequences
    Loaded 146364022 letters in 1 sequences
    Loaded 141213431 letters in 1 sequences
    Loaded 135534747 letters in 1 sequences
    Loaded 135006516 letters in 1 sequences
    Loaded 133851895 letters in 1 sequences
    Loaded 115169878 letters in 1 sequences
    Loaded 107349540 letters in 1 sequences
    Loaded 102531392 letters in 1 sequences
    Loaded 90354753 letters in 1 sequences
    Loaded 81195210 letters in 1 sequences
    Loaded 78077248 letters in 1 sequences
    Loaded 59128983 letters in 1 sequences
    Loaded 63025520 letters in 1 sequences
    Loaded 48129895 letters in 1 sequences
    Loaded 51304566 letters in 1 sequences
    Loaded 155270560 letters in 1 sequences
    Loaded 59373566 letters in 1 sequences

the output will be a modified fasta of each position amplified called "pcr_output.txt".

**5) main_post_pcr.py:**

This performs filtering and analysis on pcr products, and checks for off-target amplification. It is the same
whether you used standard pcr or multiplex pcr, as it is just analyzing pcr products, and checking for off target.

    $python3.6 src/main_post_pcr.py -h

    optional arguments:
    -h, --help            show this help message and exit
    -i PCRFILE, --pcr_output PCRFILE
                          use output of isPCR
    -tp TOTAL_PRIMERS, --total_primers TOTAL_PRIMERS
                          the pre-PCR master primer file that contains all
                          sample + primer info.

This uses the file generated from `main_pre_pcr.py` "total_list.csv" and the "pcr_output.txt" file to  generate
the final post pcr result:

    $python3.6 src/main_post_pcr.py -i pcr_output.txt -tp total_list.csv

This runs fast and outputs a lot of files! "pcr_product_info.csv" shows you all primers that generated a product.
We can see that this has more products than our original set. That is because it returns all products generated
that are inside our sequencing target position. Interesting that so many combinations of primers will actually
generate product, but perhaps this can be expected since they are near each other in genomic position for
primers generated on the same sample.

We don't know the quality of unintended primer pairings, so they are filtered before the final result is
returned.

Lastly, we have 2 final files "all_final_primers.csv" and "top_final_primers.csv" remain. 
The "top_final_primers.csv" would be the primers you would order as a pool. If you are doing standard pcr,
the primers are already nicely written out to a plate format in 2 files called:

"plate_forward_primers.csv"
"plate_reverse_primers.csv"

This format contains:

Well Position, Sequence Name, Sequence (example plate_forward_primers.csv):

    Well Position,Sequence Name,Sequence
    A1,sample1_gene1_1:15000000___F,TGCATGATGTGAAGGCAGTTTG
    A2,sample2_gene2_2:30000000___F,CCCAATCTACCATCAACTAGCTCT
    A3,sample3_gene3_3:45000000___F,TCTCTTAAACACACACGTCCCA

And is acceptable by most vendors ready to order.

**6) tag_dbsnp_primer.pl:**

tag_dbsnp_primer.pl will annotate your primers with snp's if they fall on the specific primer position! 
It is a neat tool that is pretty unique and can have future use in population genetics with some background in
that area.

In order to run this, you first need a tabix indexed reference genome. I will include the ftp's or https's for
GRCh37 and GRCh38:

38:
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz (15.4 GB)
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi (2.7 MB)

To run this, perl 5.10.1 must be in your path and you must also export the prl_src directory to your path
using the following:

    $export PERL5LIB=/home/dkennetz/PrimerPlex/Pybox/prl_src:$PERL5LIB

After you have added this prl_src to your path, you can run just by simply:

    $perl5.10.1 ../prl_src/tag_dbsnp_primer.pl -file=all_final_primers.csv -dbsnp=/home/dkennetz/referece/tabix/GRCh37/All_20180418.vcf.gz -out=dbsnp.txt

This will return an output file with all primer info, annotated with dbsnp info!
