#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "tags all primers with SNP annotations"

baseCommand: python3.6

inputs:
  primer_tabix:
    type: File
    inputBinding:
      position: 1
    default:
      class: File
      location: ../../src/primer_tabix.py

  vcf_in:
    type: File
    inputBinding:
      position: 2
      prefix: -vcf
  primer_pipeline_output:
    type: File
    inputBinding:
      position: 3
      prefix: -in
  outfile_snp:
    type: string
    inputBinding:
      position: 4
      prefix: -o

outputs:
  snp_indexed_file:
    type: File
    outputBinding:
      glob: $(inputs.outfile_snp)
