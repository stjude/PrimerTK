#!/usr/bin/env bash

module load primer3/2.4.0

cwltool --preserve-entire-environment /home/dkennetz/PrimerPlex/Pybox/genomeSV/CWL_Pipelines/PrimerSV.cwl /home/dkennetz/PrimerPlex/Pybox/genomeSV/CWL_Pipelines/inversion.yml
