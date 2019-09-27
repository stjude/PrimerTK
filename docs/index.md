# Welcome to PrimerTK documentation

PrimerTK is a toolkit is a python package that has the capability to design multiplex primer pools, and primer sets
around complex structural variants such as large deletions, insertions, translocation events, and inversions, and of course, standard primer sets. 

PrimerTK is used in conjunction with [primer3](https://github.com/primer3-org/primer3) and [in silico PCR](https://hgwdev.gi.ucsc.edu/~kent/src/) to extract reference sequences by coordinate position, design primers, compare them with each other to check for pooling compatibility (if multiplexing), check for off-target products, and filter for the top results. It also has the ability to detect SNP's present in primers and annotate them in results.

PrimerTK also has pipelines written in the Common Workflow Language (CWL). CWL a powerful tool to run structured pipelines. All it requires is an input yaml file, which I will provide the skeleton for in the CWL section. By using CWL, the user can simply adjust the input yaml file to meet the specific needs of the workflow and then run it and it will execute each consecutive step on its own and return a final output. This enables users to run multiple primer sets at the same time (or the same primer set under different thermodynamic conditions).

# Table of Contents
1. [Introduction](introduction.md)
2. [Install](installation.md)
3. [Input Templates](inputs.md)
4. [Standard and Multiplex Tutorial](standard_and_multiplex.md)
5. [Structural Variants Tutorial](structural_variants.md)
6. [Common Workflow Language Setup](cwl.md)
