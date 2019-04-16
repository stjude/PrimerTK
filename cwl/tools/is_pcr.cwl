#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "in silico PCR runner"

baseCommand: isPcr

inputs:
  chromosome_fasta:
    type: File
    inputBinding:
      position: 1
  pcr_input:
    type: File
    inputBinding:
      position: 2
  standard_output:
    type: string
    default: 'stdout'
    inputBinding:
      position: 3

outputs:
  chromosome_amplification:
    type: File
    outputBinding:
      glob: "*_primeramp.fa"

stdout: $(inputs.chromosome_fasta.nameroot)_primeramp.fa

