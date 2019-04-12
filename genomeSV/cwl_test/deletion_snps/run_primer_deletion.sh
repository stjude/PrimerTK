#!/usr/bin/env bash

module load gcc/6.3.0

export PATH=$PATH:~/primer3/src/
export PATH=$PATH:~/try_install/

cwltool --preserve-entire-environment /home/dkennetz/PrimerPlex/genomeSV/CWL_Pipelines/PrimerSV_SNP.cwl /home/dkennetz/PrimerPlex/genomeSV/CWL_Pipelines/deletion_snp.yml

