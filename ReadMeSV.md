Requirements:
 - gcc 6+
 - Primer3>=2.4.0
 - python3.6+
   - pandas>=0.22.0
   - biopython>=1.70
   - numpy>=1.4
   - cwltool>=1.0.20190228155703

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
2) PrimerMultiplex
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

       HBA2     MyHumanSample1  chr16      223000
       HBA1     MyHumanSample2  chr16      227000
       BRCA1    MyHumanSample3  chr17      41196318
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

**Running the pipelines:**

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

**Running the genomeSV pipelines:**

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

PrimerSV:

First, change to the relevant directory:

    $cd ~/PrimerPlex/Pybox/genomeSV/

**1) genome_iterator_sv.py:** (quite a lot of options, I wanted the user to have a lot of control over primer3 in).

    $python3.6 src/genome_iterator_sv.py -h
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
    -sv {deletion,inversion}, --sv-type {deletion,inversion}
                          currently supported SV primer generation: deletion and
                          inversion.

An actual command line example:

    $python3.6 src/genome_iterator_sv.py -ref /home/dkennetz/references/GRCh38/GRCh38.fasta \
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
    > -sv deletion

Running this will produce 2 outputs, a primer3 input file and a fasta with flanking regions based on inputs:

    primer3_input.20190313-104109.txt
    flanking_regions.20190313-104109.fasta

The flanking region for "deletion" type uses (flank) size upstream of pos start, and (flank) size downstream of
pos stop and designs primers around those. 1 fasta header is generated for each sample.

The flanking region for "inversion" type uses (flank) size upstream of pos start, and then uses the complement of
(flank) size upstream of pos stop. It also produces second breakpoint flanking sequence using complement 
(flank) size downstream of pos start, and (flank) size downstream of pos stop.

Inversion generates 2 fasta regions for each header, 1 at breakpoint 1 (pos start) and 1 at breakpoint 2 (pos stop).

The -st (sequence target parameter) requires primer3 to include your position in primer design, so -st should
always be (flank-1),1 to ensure your target is amplified by primers.

The GC, TM, and Size values above are defaults. If these flags are left off the command line, these values will
be used. If you'd like to change these values, simply enter a different value next to the flag.

The total flanking region for primer design will be (flank*2) so the default size is 400bp.

-mp is a mispriming library. Primer3 will penalize primers around these regions by default and will usually
exclude them because they score a high penalty.

-tp is the thermodynamic parameter path from primer3. This is the path to a directory which primer3 will use to
score primer pairs.

**2) Primer3:**

    $primer3_core --output=primer_dump.txt primer3_input.20190313-104109.txt

This will output a lot of intermediate files, as well as a primer_dump.txt file. The program will only use the
primer_dump.txt file going forward, which contains all relevant information.

The CWL version of this pipeline will not return the intermediate files by default, but the step by step program
will.

If you'd like to see an example primer_dump.txt, view the included expected_outputs file.

It will show you all the information for every primer that was designed, including the flanking sequence used, and 
why primers may or may not have failed.
Good primer pairs will have a low PRIMER_PAIR_PENALTY (<=3), but I have seen success with primers with penalties up
around even 7. You don't know until you try!

**3) main_pre_pcr.py:** (parses primer_dump.txt and uses pseudo pcr info to get "product" information).

    $python3.6 src/main_pre_pcr.py -h

    optional arguments:
    -h, --help            show this help message and exit
    -d DUMP, --primer3_dump DUMP
                          Primer3 stdout passed into a 'dump' file to be used as
                          input
    -o OUTFILE, --outfile_name OUTFILE
                          The output filename for all primer information.

to run:

    $python3.6 src/main_pre_pcr.py -d primer_dump.txt -o total_list.csv

outputs 2 files, "total_list.csv", "standard_pcr.txt"

Multiplex pcr is not currently supported with SV because I cannot amplify large deletions with in silico pcr.
It could in theory still be used to check for off-target, but I haven't found it to be too necessary for standard
primer pairs.

**4) main_post_pcr.py:** (parses exact sequence from flanking.time.fasta file that should be extracted).

    $python3.6 src/main_post_pcr.py -h

    optional arguments:
    -h, --help            show this help message and exit
    -f FLANK_FILE, --flank_file FLANK_FILE
                          use flanking_regions file from output of
                          genome_iterator_sv.py
    -tp TOTAL_PRIMERS, --total_primers TOTAL_PRIMERS
                          the pre-PCR master primer file that contains all
                          sample + primer info

to run:

    $python3.6 src/main_post_pcr.py -f flanking_regions.20190313-104109.fasta -tp total_list.csv

outputs 2 files, "total_list_gc.csv", "top_ranked_final_primers.csv".

SNP detection is also not currently detected with PrimerSV, but could be added as a feature in the future if
requested.



