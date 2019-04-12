#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "tags all primers with SNP annotations"

baseCommand: perl

inputs:
  tag_dbsnp_primer:
    type: File
    inputBinding:
      position: 1
    default:
      class: File
      location: ../../../prl_src/tag_dbsnp_primer.pl
  total_primers_file:
    type: File
    inputBinding:
      position: 2
      prefix: -file=
      separate: false
  snp_db:
    type: File
    inputBinding:
      position: 3
      prefix: -dbsnp=
      separate: false
    secondaryFiles:
     - .tbi
  outfile_snp:
    type: string
    inputBinding:
      position: 4
      prefix: -out=
      separate: false

outputs:
  snp_indexed_file:
    type: File
    outputBinding:
      glob: $(inputs.outfile_snp)
