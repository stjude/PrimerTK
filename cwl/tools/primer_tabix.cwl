#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "tags all primers with SNP annotations"

baseCommand: [primer_tk, tabix]
inputs:
  vcf_in:
    type: File
    inputBinding:
      position: 3
      prefix: -vcf
    secondaryFiles:
     - .tbi
  primer_pipeline_output:
    type: File
    inputBinding:
      position: 4
      prefix: -in
  outfile_snp:
    type: string
    inputBinding:
      position: 5
      prefix: -o

outputs:
  snp_indexed_file:
    type: File
    outputBinding:
      glob: $(inputs.outfile_snp)
