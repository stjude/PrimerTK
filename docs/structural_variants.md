# Structural Variants Tutorial

The possible structural variants to design primers around are insertions, deletions, translocations, and inversions. These have a slightly different file format than the Standard and Multiplex inputs, and insertion + translocation are different than deletion + inversion. To view input file specifications, check out the [inputs page](inputs.md).

Although the input files are different, the pipeline inputs and outputs are the same for all SV's, so I will just perform the tutorial with a deletion.
The steps of the pipeline are also quite similar to the [Standard and Multiplex Tutorial](standard_and_multiplex.md), so if you have performed that one already, this may seem more straightforward to you.

This pipeline just parses sequences from the reference genome based on the coordinates in the inputs. If the SV type is a bit more complex than a deletion, such as an insertion of the - strand to the + strand of the normal chromosome, some more complex string manipulation is done when designing the positional flanking sequences.

*NOTE: THERE IS CURRENTLY NO OPTION TO MULTIPLEX SV'S ALTHOUGH THIS COULD BE ADDED IN A LATER RELEASE IF DESIRED.*

If you followed the installation page and git cloned the repo, there is a directory called `test` in the PrimerTK repo. This is for code testing but it also has data for the user to implement the full pipeline (testing is good).

For the SV workflow, only primer3 needs to be added to your path. If you followed standard installation and added it to your .bashrc, you can skip this. If not:

```
export PATH=~/bin/primer3:$PATH
```

Or add the above to your ~/.bashrc to prevent this in the future. Now check ot make sure `primer_tk` and `primer3_core` are accessible via the command line anywhere in your system.

## 1) Setup a tutorial directory called sv_tutorial with files from test

I will do everything in a directory called sv_tutorial:

```
mkdir -p ~/sv_tutorial/

cd ~/sv_tutorial/
cp ~/PrimerTK/test/data/humrep.ref .
cp ~/PrimerTK/test/data/test_standard.fa .
cp ~/PrimerTK/test/data/input_sv.csv .
```

## 2) Start running the program

Then let;s start running the program. Step 1 is called `iterator_sv` and it has a lot of inputs, but that is just because I wanted to give users a lot of flexibility on thermodynamic parameters for primer3.

Let's display the help message and explain the parameters:

```

primer_tk iterator_sv -h
usage: primer_tk iterator_sv [-h] -ref REF_GENOME -in REGIONS_FILE
                             [-opt_size PRIMER_OPT_SIZE]
                             [-min_size PRIMER_MIN_SIZE]
                             [-max_size PRIMER_MAX_SIZE]
                             [-opt_gc PRIMER_OPT_GC] [-min_gc PRIMER_MIN_GC]
                             [-max_gc PRIMER_MAX_GC] [-opt_tm PRIMER_OPT_TM]
                             [-min_tm PRIMER_MIN_TM] [-max_tm PRIMER_MAX_TM]
                             [-sr PRODUCT_SIZE_RANGE]
                             [-flank FLANKING_REGION_SIZE]
                             [-st SEQUENCE_TARGET] -mp MISPRIMING -tp
                             THERMOPATH -sv {deletion,inversion,insertion}

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
  -sv {deletion,inversion,insertion}, --sv-type {deletion,inversion,insertion}
                        currently supported SV primer generation: deletion,
                        inversion, and insertion.
```
The `-ref` parameter is the reference genome for which you want to design primers. The `-in` parameter is the position input file with the positions of interest. The parameters in the middle with `-opt, -min, -max` in front of them are the thermodynamic parameters. They all contain defaults that I have found to be good for most standard small fragment experiments. `-sr` is the product_size_range, `-mp` is the primer mispriming library which penalizes primers that land on sequences in the file. I have provided a copy of the file I use in the `test` directory as I think it is sufficient for most cases. `-tp` is the thermodynamics path for your primer3 program to use. This is in the primer3 install location, so if you installed primer3 in your `~/bin` like the standard installation, this path will be `~/bin/primer3/src/primer3_config/`. Lastly, `-sv` is the sv type. For this I will choose deletion. If you are doing a translocation setup, you should use insertion. This is explained in the [inputs](inputs.md) for the insertion file.

Anything with a default value does not need to be specified unless you want to change that default value to something else, but for the sake of the tutorial I will include them:

*NOTE: THE -st SHOULD ALWAYS BE (-flank-1,1) SO IF FLANK IS 200, -st SHOULD be 200-1,1 OR 199,1. THIS IS SO YOUR POSITION OF INTEREST IS GUARANTEED TO BE IN YOUR PRODUCT. I ALSO ALWAYS SET MY -sr UPPER LIMIT TO BE TWICE MY FLANK SIZE.*

```
primer_tk iterator_sv \
 -ref test_standard.fa \
 -in input_sv.csv \
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
 -mp /home/dkennetz/sv_tutorial/humrep.ref \
 -tp /home/dkennetz/bin/primer3/src/primer3_config/ \
 -sv deletion
```

This will output the following files:

```
flanking_regions.input_sv.fasta
primer3_input.input_sv.txt
```

The flanking regions file is used later in the workflow, and the primer3_input file is used in the next (primer3) step. The filename convention is always flanking_regions.<infile>.fasta and priemr3_input.<infile>.txt.

## 3) Run Primer3

This step is easy enough but can be time consuming for large datasets. In our case, it should take about 30 seconds.

```
primer3_core --output=primer_dump.txt primer3_input.input_sv.txt
```

This will output a log file named primer_dump.txt (we want) and a bunch of intermediate files (we don't want).
I like to clean the intermediates up (the intermediates also aren't kept if you use CWL).

```
rm -rf *.int *.for *.rev
```

Up to 5 primer pairs will be output per position, and all will be visible in the files. In the final step, we will have one file called "top_primers" and one called "all_primers". The top_primers file will be the top ranking primer pair that passed all filtering steps, and all_primers will output all primers that passed all filtering steps (up to 5 per position).

## 4) Run the primer3 output parsing

This is similar to the Standard Primer Design pipeline in that all we are doing is parsing the primer3 output in this step:

```
primer_tk pre_sv -h
usage: primer_tk pre_sv [-h] -d DUMP -o OUTFILE [-pcr PCRFILE]

optional arguments:
  -h, --help            show this help message and exit
  -d DUMP, --primer3_dump DUMP
                        Primer3 stdout passed into a 'dump' file to be used as
                        input
  -o OUTFILE, --outfile_name OUTFILE
                        The output filename for all primer information.
  -pcr PCRFILE, --pcrfile PCRFILE
                        The pseudopcr file
```

the -pcr flag here just sets up a pseudopcr file and could be useful in the future, so I included it. To run:

```
primer_tk pre_sv \
 -d primer_dump.txt \
 -o total_primers.csv \
 -pcr fake_pcr.txt
```

If we look at the output, we see two primer pairs were designed around 1 position. This looks like a tough region to design around! We can always go back and check out the `primer_dump.txt` file to see why a specific region was rejected. It has a line like this for each position in the file:

```
PRIMER_LEFT_EXPLAIN=considered 653, low tm 476, high repeat similarity 113, long poly-x seq 64, ok 0
PRIMER_RIGHT_EXPLAIN=considered 1308, GC content failed 557, low tm 53, high tm 521, high hairpin stability 147, high repeat similarity 21, ok 9
```

So for position 1, we see all the problems associated with the region and why primers failed there. If we actually look at the sequence it used:

```
SEQUENCE_TEMPLATE=AACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGTTGCAAAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAGGCGC
```

This is reasonable!

## 5) Pseudo PCR check

There is no actual PCR done here, as it isn't super necessary when we are not multiplexing. The likelihood of severe off-target for a single pair is pretty low, although it does happen occasionally!

To run the next step:

```
primer_tk post_sv -h
usage: primer_tk post_sv [-h] [-f FLANK_FILE] [-tp TOTAL_PRIMERS]
                         [-all ALL_FINAL_PRIMERS] [-top TOP_FINAL_PRIMERS]
                         [-plate PLATE_BASENAME]

optional arguments:
  -h, --help            show this help message and exit
  -f FLANK_FILE, --flank_file FLANK_FILE
                        use flanking_regions file from output of
                        genome_iterator_sv.py
  -tp TOTAL_PRIMERS, --total_primers TOTAL_PRIMERS
                        the pre-PCR master primer file that contains all
                        sample + primer info
  -all ALL_FINAL_PRIMERS, --all_final_primers ALL_FINAL_PRIMERS
                        all primers generated for targets
  -top TOP_FINAL_PRIMERS, --top_final_primers TOP_FINAL_PRIMERS
                        top primers generated for targets
  -plate PLATE_BASENAME, --plate_basename PLATE_BASENAME
                        the basename of the primers in plate format ready to
                        order.
```

For this step, we actually pass the flanking regions file output from step 1 back into this step because this is the sequence that was used to design the primer pair. From this we can extract the full sequence the primers are amplifying and give pseudo-pcr product information such as GC content and product length, which is why I say there is pseudo-pcr done here.

The primer positions for each pair in regards to the reference are also output, as this can be useful to potentially trouble shoot primers that seem like they have failed.

The outputs for this are all primer pairs, and top primer pair by input position, as well as an easy to order primer plate format.

To run:

primer_tk post_sv \
 -f flanking_regions.input_sv.fasta \
 -tp total_primers.csv \
 -all all_primers.csv \
 -top top_primers.csv \
 -plate plate \

And all final outputs from this pipeline are:

```
primer3_input.input_sv.txt
flanking_regions.input_sv.fasta
primer_dump.txt
total_primers.csv
fake_pcr.txt
all_primers.csv
top_primers.csv
plate_F.csv
plate_R.csv
```

## 7) SNP Annotation of Primers

To view this, please head over to [the Standard and Multiplex Tutorial Step 7](standard_and_multiplex.md#tabix).
