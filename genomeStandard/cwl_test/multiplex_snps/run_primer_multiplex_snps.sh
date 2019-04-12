#!/usr/bin/env bash

module load gcc/6.3.0
module load python/3.6.1

export PATH=/home/dkennetz/try_install:$PATH
export PATH=/home/dkennetz/primer3/src:$PATH
export PATH=/home/dkennetz/.local/bin:$PATH

cwltool --preserve-entire-environment /home/dkennetz/PrimerPlex/genomeStandard/CWL_Pipelines/PrimerMultiplexSNPs.cwl /home/dkennetz/PrimerPlex/genomeStandard/CWL_Pipelines/multiplex_pcrSNPs.yml


