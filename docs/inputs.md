# Input File Formats for Various Pipelines

The following will present template files which display how the coordinate file input to the program should look. These inputs will contain the following columns which correspond to things used by the program.

*NOTE: NO INPUT FILES HAVE HEADERS. THE FIRST LINE OF THE FILE SHOULD BE THE FIRST POSITION OF INTEREST.*

*NOTE: POSITIONS ARE REFERENCE GENOME SPECIFIC*

## Multiplex / Standard Primer Generation Input File

This pipeline uses an input file with 4 columns. The file can be comma separated or tab delimited. Any other formats will cause the first step in the pipeline to error out.

The first column is "Gene" or "ID". If I am targeting a specific SNP on a specific gene, I usually use the gene name here. Otherwise, I use a random identifier. Either are okay, as this column is not matched to anything, but rather is just used in the naming convention for final output. The second column is going to be the name of the sample, or the identifier you want to give the location. The reason both column 1 and column 2 are nice is because you could look for multiple genes in the same sample (or you could design a whole set around a single sample). Again, the sample name is not checked against anything and is just used for final naming. Column 3 is the chromosome location of the target position. This can be in the format `chr1` or `1` but the `chr` should be lowercase and should have no spaces between itself and the number or X, Y. The fourth column is the position of interest. This should be the integer value only and should not contain commas.

multiplex_input.csv:

```
BRCA1,Sample1,chr17,43100000
HBA2,Sample2,chr16,172900
HBG1,SampleNice,chr11,5249000
```

The same input, only tab delimited (should end in .txt):

multiplex.txt:

```
BRCA1	Sample1	chr17	43100000
HBA2	Sample2	chr16	172900
HBG1	SampleNice	chr11	5249000
```

Personally, I think csv files are easier to use as sometimes different text editors tab things differently. If you have an issue with a txt input, try to run it again using csv input instead.

*NOTE: THE STANDARD PRIMER GEN FILE IS THE EXACT SAME AS THE MULTIPLEX INPUT.*

## Deletion and Inversion Input File

This has the same format as the standard and multiplex input, except it has a fifth column. 

Columns 1-3 are the same; however, column 4 is deletion / inversion start and column 5 is deletion / inversion stop

So if the user was looking for a 100 bp deletion / inversion in a random part of chromosome 1, the input could look like this:

deletion.csv:

```
random,mysample,1,100000,100100
...
```
This would then take flanking sequence on either end of those positions to span the breakpoint.

For an inversion, the input would be the same:

inversion.csv:

```
random,mysample,1,100000,100100
```

The program just performs different operations on the sequence when different SV types are specified.

## Insertion and Translocation Input File

This file is a bit more involved, as it includes the normal chromosome start and stop positions (where the insertion will be placed) and the inserted regions start and stop positions, as well as strand orientation of the inserted region. In theory for a translocation you can just add another entry with the positions flip-flopped, or just treat all the translocation events as individual insertions. See below:

insertion.csv:

```
gene1,sample1,chr1,1000,1100,chr2,1000,1100,+
gene2,sample2,chr1,1000,1100,chr2,1000,1100,-
```

Gene1 sample1 would have a 100 bp insertion on chr1 from chr2 using the forward strand orientation.
Gene2 sample2 would have a 100 bp insertion on chr1 from chr2 using the reverse strand orientation (the seq will be reverse complemented).

Now if this was a translocation and gene1 sample1 was a translocation (where the two regions switched places) you could just handle it in two lines by switching the positions (let's pretend the same thing happened for gene2 sample2):

translocation.csv:

```
gene1,sample1,chr1,1000,1100,chr2,1000,1100,+
gene1,sample1,chr2,1000,1100,chr1,1000,1100,+
gene2,sample2,chr1,1000,1100,chr2,1000,1100,-
gene2,sample2,chr2,1000,1100,chr1,1000,1100,-
```

Now of course you may have to think about strandedness a little bit here, but it is definitely applicable!

Before diving into the program by yourself, I would recommend heading to the tutorial of interest first:

[Standard and Multiplex Tutorial](standard_and_multiplex.md)

[Structural Variants Tutorial](structural_variants.md)
