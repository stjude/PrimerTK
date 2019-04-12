Before reading this, you should ensure that all software is installed properly and functional by doing the
walkthrough for ReadMeStandard.md, and ReadMeSV.md. The other ReadMe will walk you through installation dependencies and
will help you understand the inputs required at each step.

This will generalize a walkthrough for PrimerStandard.cwl, but the rest work in the same way. Modify the paths in the .yml file corresponding to the workflow.
Once you get to the end of this, you should understand how to generally run each CWL tool. If not feel free to reach out to me. Contact info at the bottom.

This ReadMe is applicable to both PrimerStandard.cwl and PrimerStandardSNPs.cwl.

To be able to use the cwl pipelines, python3.6 must be installed and cwl-runner must be installed.

cwl-runner would be installed if you pip install requirements.txt:

    $pip3.6 install --user cwl-runner

This will set cwl-runner up in your user environment (this will be complete if you followed the install reqs). 

From here:

1) go to the cwl pipelines directory:

    $cd genomeStandard/CWL_Pipelines/

I have provided example standard_pcr.yml and standard_pcrSNPs.yml in this directory.

See below (You will only modify what follows `path:` for each variable, each is tagged with a modify):

    ref_genome:
      class: File
      path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa #modify to your reference

    regions_file:
      class: File
      path: /home/dkennetz/PrimerPlex/Pybox/genomeStandard/tests/input1.txt # modify to your regions file location

    primer_opt_size: 22
    primer_min_size: 18
    primer_max_size: 25
    primer_opt_gc: 50
    primer_min_gc: 20
    primer_max_gc: 80
    primer_opt_tm: 60
    primer_min_tm: 57
    primer_max_tm: 63
    product_size_range: '200-400'
    flanking_region_size: 200
    sequence_target: '199,1'
    mispriming_library: '/home/dkennetz/PrimerPlex/Pybox/genomeStandard/tests/humrep.ref' # modify 
    thermodynamics_path: '/hpcf/apps/primer3/install/2.4.0/src/primer3_config/' # modify to primer3/src/primer3_config/ in your install location
    output: 'primer_dump.txt'
    outfile: 'total_list.csv'
    pcr: 'standard'
    percent_alignment: 50

    chromosome_fasta:
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr1.fa } #modify the path
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr2.fa } #modify the path
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr3.fa } #modify the path
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr4.fa } ...
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr5.fa } ...
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr6.fa } ...
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr7.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr8.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr9.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr10.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr11.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr12.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr13.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr14.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr15.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr16.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr17.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr18.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr19.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr20.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr21.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chr22.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chrX.fa }
     - {class: File, path: /research/rgs01/reference/public/genomes/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/chrY.fa }

     catted_filename: 'pcr_products.fasta'

The good part about this is that you will only have to modify the paths 1 time. I would recommend setting up an exact copy of this yml for each reference genome.

For any of the .ymls, you will have to modify the path to correspond to your file locations. Furthermore, any of the other input parameters can be adjusted to
the liking of the user.

2) To test a run and see if it works (assuming your are still in the CWL_Pipelines directory:

    $cd ../cwl_test/standard/

There will already be the expected outputs in here (if you are using GRCh37). 

Open the file "run_primer_standard.sh" and there will be 3 paths you need to change:

    export PATH=/home/<user>/PrimerPlex/Pybox/isPcr:$PATH
    export PERL5LIB=/home/<user>/PrimerPlex/Pybox/prl_src:$PERL5LIB

    cwltool --preserve-entire-environment /home/<user>/PrimerPlex/Pybox/genomeStandard/CWL_Pipelines/PrimerStandard.cwl /home/<user>/PrimerPlex/Pybox/genomeStandard/CWL_Pipelines/standard_pcr.yml

If you are unfamiliar with cwl, it is a python api that uses these "workflow description sheets" to specify what tool should run in what order, and what the expected inputs and outputs are for each.
These are known as workflows. The "standard_pcr.yml" file has all the external inputs that the user must provide to the workflow in order for successful completion.

The flag "--preserve-entire-environment" is fairly straightforward, it tells cwl to inherit your path specific variables to run the tool. The command line to run this is:

    $./run_primer_standard.sh

This will kick off the pipeline and report each step to screen.

If you have any errors running the pipeline, feel free to send them to me at:

dennis.kennetz@stjude.org

I will respond as soon as possible.
