#!/usr/bin/env bash

export PATH=$PATH:/home/dkennetz/primer3/src

cwltool --preserve-entire-environment /home/dkennetz/PrimerPlex/Pybox/genomeSV/CWL_Pipelines/PrimerSV.cwl /home/dkennetz/PrimerPlex/Pybox/genomeSV/CWL_Pipelines/inversion.yml
