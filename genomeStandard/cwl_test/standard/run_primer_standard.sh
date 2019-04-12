#!/usr/bin/env bash

module load gcc/6.3.0

export PATH=$PATH:/home/dkennetz/primer3/src/
export PATH=$PATH:/home/dkennetz/try_install/

cwltool --preserve-entire-environment /home/dkennetz/PrimerPlex/genomeStandard/CWL_Pipelines/PrimerStandard.cwl /home/dkennetz/PrimerPlex/genomeStandard/CWL_Pipelines/standard_pcr.yml

