#!/usr/bin/env bash

module load gcc/6.3.0

export PATH=/home/dkennetz/try_install:$PATH
export PATH=/home/dkennetz/primer3/src/:$PATH

cwltool --preserve-entire-environment /home/dkennetz/PrimerPlex/genomeStandard/CWL_Pipelines/PrimerMultiplex.cwl /home/dkennetz/PrimerPlex/genomeStandard/CWL_Pipelines/multiplex_pcr.yml

