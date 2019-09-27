# Standard and Multiplex Tutorial

These two pipelines have very similar workflow so I will include them together. If you followed the installation page and git cloned the repo, there is a directory called `test` in the PrimerTK repo. This is for code testing but it also has data for the user to implement the full pipeline (testing is good).

If you added primer3 and isPcr to you .bashrc in the installation guide, skip this. If not, run the following (assuming you installed primer3 and isPcr in your `~/bin` directory:

```
export PATH=~/bin/primer3:$PATH
export PATH=~/bin/isPcr33:$PATH
```

Or add the above to your `~/.bashrc` to prevent yourself from having to do this every time you log in to a new session.
So now `primer_tk`, `isPcr`, and `primer3_core` should all be accessible from the command line, anywhere in your system.
This tutorial is going to assume that you installed PrimerTK in your home directory but you may well have set it up anywhere. It is also going to assume you used the standard installation script so the path to the primer3 executable and isPcr executable will reflect that.

## 1) Setup a tutorial directory with files from test

These files are small so you can copy them over directly or just point to them as you run the program.

I will do everything in a directory called tutorial:

```
cd ~/tutorial/
cp ~/PrimerTK/test/data/humrep.ref .
cp ~/PrimerTK/test/data/input_standard.csv .
cp ~/PrimerTK/test/data/test_standard.fa .
```

## 2) Start running the program
Then let's start running the program. Step 1 called `iterator` has a lot of inputs, but that is just because I wanted to give users a lot of flexibility on thermodynamic parameters for primer3.

Let's display the help message and explain the parameters:

```
primer_tk iterator -h

usage: primer_tk iterator [-h] -ref REF_GENOME -in REGIONS_FILE
                          [-opt_size PRIMER_OPT_SIZE]
                          [-min_size PRIMER_MIN_SIZE]
                          [-max_size PRIMER_MAX_SIZE] [-opt_gc PRIMER_OPT_GC]
                          [-min_gc PRIMER_MIN_GC] [-max_gc PRIMER_MAX_GC]
                          [-opt_tm PRIMER_OPT_TM] [-min_tm PRIMER_MIN_TM]
                          [-max_tm PRIMER_MAX_TM] [-sr PRODUCT_SIZE_RANGE]
                          [-flank FLANKING_REGION_SIZE] [-st SEQUENCE_TARGET]
                          [-mp MISPRIMING] [-tp THERMOPATH]

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
                        full path to mispriming library for primer3 (EX:
                        /home/dkennetz/mispriming/humrep.ref
  -tp THERMOPATH, --thermopath THERMOPATH
                        full path to thermo parameters for primer3 to use (EX:
                        /home/dkennetz/primer3/src/primer3_config/) install
                        loc
```

As you can see, the help message is pretty detailed. The only thing that may need further explaining is the -mp (mispriming flag). This is just a file containing human repetitive elements or sequences that are common primer mispriming locations. Primer3 will use this to score primers in that location with a penalty so they are not likely to be selected.

Anything with a default value does not have to be specified on the command line (if you want to use that value as the parameter). These have been selected because I have seen good success with these values for standard pcr. In this tutorial, I will specify every value, but again you do not have to specify defaults.

*NOTE: THE -mp AND -tp FLAGS SHOULD BE THE FULL PATH TO FILE AS THIS IS WHAT PRIMER3 REQUIRES*
*NOTE: THE -st SHOULD ALWAYS BE (-flank-1,1) SO IF FLANK IS 200, -st SHOULD be 200-1,1 OR 199,1. THIS IS SO YOUR POSITION OF INTEREST IS GUARANTEED TO BE IN YOUR PRODUCT. I ALSO ALWAYS SET MY -sr UPPER LIMIT TO BE TWICE MY FLANK SIZE.*

```
primer_tk iterator \
 -ref test_standard.fa \
 -in input_standard.csv \
 -opt_size 22 \
 -min_size 18 \
 -max_size 25 \
 -opt_gc 50 \
 -min_gc 20 \
 -max_gc 80 \
 -opt_tm 60 \
 -min_tm 57 \
 -max_tm 63 \
 -sr 200-400 \
 -flank 200 \
 -st 199,1 \
 -mp /home/dkennetz/tutorial/humrep.ref \
 -tp /home/dkennetz/bin/primer3/src/primer3_config/
```

This module is just used to parse any reference genome and setup a primer3 config. The outputs will be:

```
flanking_regions.input_standard.fasta
primer3_input.input_standard.txt
```
The fasta will show us what sequences we pulled down, and the primer3_input will be the input to primer3.
The output files will maintain the same name as the input file used for the `-in` flag.

## 3) Run Primer3

This step is easy enough but can be time consuming for large datasets. In our case, it should take about 30 seconds.

```
primer3_core --output=primer_dump.txt primer3_input.input_standard.txt
```

This will output a log file named primer_dump.txt (we want) and a bunch of intermediate files (we don't want).
I like to clean the intermediates up (the intermediates also aren't kept if you use CWL).

```
rm -rf *.int *.for *.rev
```

Up to 5 primer pairs will be output per position, and all will be visible in the files. In the final step, we will have one file called "top_primers" and one called "all_primers". The top_primers file will be the top ranking primer pair that passed all filtering steps, and all_primers will output all primers that passed all filtering steps (up to 5 per position).

We can also see why primers failed at specific regions if we look at the primer dump file. For example, position 1 returned 0 primers. If we look at why:

```
PRIMER_LEFT_EXPLAIN=considered 653, low tm 476, high repeat similarity 113, long poly-x seq 64, ok 0
PRIMER_RIGHT_EXPLAIN=considered 791, GC content failed 140, low tm 304, high tm 163, high hairpin stability 47, high repeat similarity 83, long poly-x seq 54, ok 0
```
Primer3 gives us a pretty informative description. If we then take a look at the sequence it used to design primers:

```
SEQUENCE_TEMPLATE=AACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGT
```
It looks highly repetitive. This seems reasonable!


## 4) Run the pre-pcr step and the multiplex filtering

This step is used for parsing the primer3 output file and setting up the input for pcr, and more importantly for filtering for multiplex primer pools. Let's look at the help message:

```
primer_tk pre -h
usage: primer_tk pre [-h] -d DUMP -o OUTFILE [-nd NO_DIMER]
                     [-spcr STANDARD_PCR_FILE] [-mpcr MULTIPLEX_PCR_FILE]
                     [-pa PERCENT_ALIGNMENT] -pcr {standard,multiplex}

Command Line argument for total primerinput file to check if primers have a
degreeof complementarity with each other as definedby the user. Default is 60%
(fairly strict).

optional arguments:
  -h, --help            show this help message and exit
  -d DUMP, --primer3_dump DUMP
                        Primer3 stdout passed into a 'dump' file to be used as
                        input
  -o OUTFILE, --outfile_name OUTFILE
                        The output filename for all primer information.
  -nd NO_DIMER, --no_dimer NO_DIMER
                        The primers left after dimers removed.
  -spcr STANDARD_PCR_FILE, --standard_pcr_file STANDARD_PCR_FILE
                        The file to be used for standard pcr input
  -mpcr MULTIPLEX_PCR_FILE, --multiplex_pcr_file MULTIPLEX_PCR_FILE
                        The file to be used for multiplex pcr input
  -pa PERCENT_ALIGNMENT, --percent_alignment PERCENT_ALIGNMENT
                        Percent match between 2 primers for pair to be
                        discarded. EX: primer_len = 22, percent_aln = 60
                        dimer_len = (60/100) * 22 = 13.2 -> 13.
  -pcr {standard,multiplex}, --pcr_type {standard,multiplex}
                        perform standard or multiplex pcr on given inputs.
```

### For multiplex PCR:

PrimerTK does two important things for multiplexing:
1. It checks for the presence of primer dimers between any two primers in the entire dataset.
2. It sets up an all vs all PCR input.

Both of these are important because in multiplexing, all primers will be in a single pool, so they will interact with each other. If the reaction between each other is more favorable than the reaction with the DNA template, or it is competing, those primers will have less effective concentration because they will bind to each other instead of the template. This will give us lower or no representation of the expected site.

Furthermore, we do in silico PCR to check for potential off-target amplification. This has the same effect as primer dimerization to a large extent. If there are competing sites, there is less effective concentration of primer at the site of interest. 

To deal with this, a simple string matching algorithm has been implemented. This returns a score based on how well the primers align with each other. I also ask for a user input of the `-pa` or percent alignment. I use this percent alignment strategy because primers will be different lengths. If the primer length times the percent alignment is greater than the calculated alignment score, then the primer is dropped. 

This is repeated until every possible combination in the pool has been checked.

Afterwards, it returns a primer file and a pcr file that only contains the post-filter primer pairs.

to run:

```
primer_tk pre \
 -d primer_dump.txt \
 -o total_primers.csv \
 -nd no_dimers.csv \
 -mpcr multiplex_pcr.txt \
 -pa 60 \
 -pcr multiplex
```

We can see from our example if we check the `no_dimers.csv` file, that with this pa score the primer pair was not filtered. However, as the size of the dataset increased, there is much more potential for primer-primer interaction and many more will be filtered. 

Go ahead and cat the multiplex_pcr.txt output by this program as well. We can see that there are 4(!) pcr reactions that isPcr is going to simulate. This goes to show the complexity of pcr multiplexing. The reason one primer pair has 4 pcr reactions tested is because:
1. It considers the forward primers possibility to amplify using itself (this may be overkill but multiplexing results are good).
2. It obviously considers the forward with the reverse.
3. It considers the reverse with the forward.
4. It considers the reverse with the reverse.

### For standard PCR:

There is not much to consider for standard PCR input, so PrimerTK just parses information from the primer3 log and sets up a single pcr reaction for each primer pair:

primer_tk pre \
 -d primer_dump.txt \
 -o total_primers.csv \
 -spcr standard_pcr.txt \
 -pcr standard

This is fast for all datasets because it doesn't perform any all vs all comparisons.

## 5) In Silico PCR off-target check

*NOTE: IN SILICO PCR REQUIRES THE REFERENCE GENOME TO BE SPLIT INTO CHROMOSOMES FOR THE HUMAN REFERENCE (OR ANY REFERENCE ABOVE 2.1GB). I WILL POST HOW I DO THIS BELOW.*

*NOTE: THIS IS NOT REQUIRED FOR OUR TEST CASE BECAUSE THERE IS ONLY ONE CHROMOSOME.*

#### Splitting a reference genome by chromosome:

The easiest way to do this is using [pyfaidx](https://github.com/mdshw5/pyfaidx).

to install:

```
pip3 install pyfaidx --user
```

Then go to the directory with your reference genome and run:

```
faidx -x reference.fa
```

This will result in your original reference, a .fai (fasta index) file for your reference, and then each individual chromosome split out as a separate fasta file. For GRCh38 human reference this took less than 3 minutes.

### Running isPcr for the tutorial

In the PrimerTK/scripts directory, there is a script called `ispcr.sh`. This is an easy to use script that runs in silico PCR using the pcr file output from the last step. Running requires two inputs: a directory with your reference split into chromosomes, and your pcr input file (generated by PrimerTK). So from your tutorial directory (I am going to run this like I installed PrimerTK in my home directory). Since I copied the test.fa into this directory when I started, and my multiplex pcr input file is in the directory as well, I can just run like this:

```
# Usually your chromosome dir will be located somewhere else, in that case put the path.
~/PrimerTK/scripts/ispcr.sh ./ ./multiplex_pcr.txt
```

This will print:
`Running In Silico PCR on pcr input positions...`

And in our case it will take less than a second. Normally it takes a few minutes if there are a lot of targets in a real reference genome.

The script writes output to a file called `pcr_output.fa`. This will be used in the next step. If we `cat` the file, we see the following:

```
>1:1074+1381 sample3_gene3_1:1200__ 308bp GCAGAAACTCACGTCACGGT CTGCCACTACACCTTGAGCA
GCAGAAACTCACGTCACGGTggcgcggcgcagagacgggtagaacctcag
...
>1:1074-1381 sample3_gene3_1:1200__ 308bp CTGCCACTACACCTTGAGCA GCAGAAACTCACGTCACGGT
CTGCCACTACACCTTGAGCAagaggaccctgcaatgtccctagctgccag
...
```

We see the structure of the header: positions, sample name, primer1, primer2. We may notice there are two amplifications for this, but they have the same positions. The first has `1:1074+1381` a noticable (+) sign between the position coordinates, and the second `1:1074-1381` has a noticable (-) sign between the position coordinates. This indeicates that the first primer pair amplified the forward strand, while the second amplified the reverse.

In the next step, we use all positional outputs from this file to filter out off target pairs. We also use the sequence to show predicted product length and GC content. We only test amplifications up to 5kb (which may be overkill), because most products longer than even 2kb will have incomplete amplification in the lab, so I figured 5kb would cover every possible opportune amplificiations.

## 6) Post PCR Processing for multiplex samples

This is going to use the `total_primers.csv` file generated by the `pre` step, and the `pcr_output.fa` file as inputs to do some final checks. It also has some other flags, such as user specified off-target max (8 by default) and output file names:

```
primer_tk post -h
usage: primer_tk post [-h] [-i PCRFILE] [-tp TOTAL_PRIMERS]
                      [-ot OFF_TARGET_MAX] [-pcri PCR_PRODUCT_INFO]
                      [-all ALL_PRIMER_INFO] [-top TOP_FINAL_PRIMERS]
                      [-plate PLATE_BASENAME]

optional arguments:
  -h, --help            show this help message and exit
  -i PCRFILE, --pcr_output PCRFILE
                        use output of isPCR
  -tp TOTAL_PRIMERS, --total_primers TOTAL_PRIMERS
                        the pre-PCR master primer file that contains all
                        sample + primer info.
  -ot OFF_TARGET_MAX, --off_target_max OFF_TARGET_MAX
                        the maximum number of off target hits of primer.
  -pcri PCR_PRODUCT_INFO, --pcr_product_info PCR_PRODUCT_INFO
                        the information of all products generated by isPcr
  -all ALL_PRIMER_INFO, --all_primer_info ALL_PRIMER_INFO
                        all of the successful primers generated
  -top TOP_FINAL_PRIMERS, --top_final_primers TOP_FINAL_PRIMERS
                        the top primers generated for each position
  -plate PLATE_BASENAME, --plate_basename PLATE_BASENAME
                        the basename of the primers in plate format ready to
                        order.
```

The PCR product info file will include information regarding all products generated by isPcr. This will include off-target sites which are filtered out of our final primer files.

The all primers file will contain all primers that passed all the applied filters for a given position, for all positions in the file. This file will have a tag for each primer pair to indicate how many off-target sites it amplified. I have seen cases where the top ranking primer pair actually amplified more off-target than the second ranking primer pair (differences as great as 5 to 0). For this case, I go back and look at the penalty in the `primer_dump.txt` file to see if the penalty difference between the two is substantial and if they have similar tm and product gc. If the thermodynamics look good between the pair, and one has less off-target than the other, I pick the one with less off-target.

The top final primers file will display the top ranking primer pair for each position. Also, the primer positions in regards to the reference will also be output, as this can be useful to potentially trouble shoot primers that seem like they have failed.

Lastly, the plate_basename is a simple string that we want our primer plate to be named. This step will setup a ready to order file with the forward primers in a plate format and the reverse primers in a plate format. The default for these files will be:

`plated_primers_F.txt` and `plated_primers_R.txt`.
The forward and reverse primers will correspond to the same positions on the plates. So primer pair 1 will be in A1 on each plate, primer pair 2 will be in A2... etc.

To run:

primer_tk post \
 -i pcr_output.fa \
 -tp total_primers.csv \
 -ot 5 \
 -pcri product_info.csv \
 -all all_primers.csv \
 -top top_primers.csv \
 -plate plated

This will be fast and the final output files from the whole pipeline will be:

```
flanking_regions.input_standard.fasta
primer3_input.input_standard.txt
primer_dump.txt
Primer_Dimers.txt
total_primers.csv
no_dimers.csv
multiplex_pcr.txt
pcr_output.fa
product_info.csv
top_primers.csv
all_primers.csv
plated_F.csv
plated_R.csv
```

## 7) SNP Annotation of Primers <a name="tabix"></a>

This cannot be done in the tutorial, but can be done with real datasets, as long as there is a vcf file for the reference you are using. 

An example would be if I was using GRCh38 and wanted to annotate snps according to this genome. The file I may use would be (and the tabix indexed genome):

```
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20190923.vcf.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20190923.vcf.gz.tbi
```

After these are downloaded and I want to annotate my files I run:

```
primer_tk tabix -h
usage: primer_tk tabix [-h] [-vcf VCF] [-in P_INFO] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -vcf VCF, --variant-call-file VCF
                        Tabix indexed VCF.
  -in P_INFO, --primer-input-file P_INFO
                        The output of the primer pipeline.
  -o OUTPUT, --output OUTPUT
                        The name of the output file
```

The vcf file would be the tabix-indexed reference genome, the input P_INFO file would be the `all_primers.csv` file from the post-processing output, and the output file would be whatever you want!

An example of running is:

```
primer_tk tabix \
  -vcf ~/tabix_indexed_genome/GRCh38/All_20180418.vcf.gz \
  -in all_primers.csv \
  -o snp_annotated_primers.csv
```

which would return a the snp annotated primer file!

Next, proceed to the [Sructural Variants Tutorial](structural_variants.md) or the [CWL Tutorial](cwl.md) to learn how to run this pipeline as a single step workflow!
