#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "main pre pcr setup"

baseCommand: [primer_tk, pre]

inputs:
  primer_dump:
    type: File
    default: primer_dump.txt
    inputBinding:
      position: 3
      prefix: -d
  outfile:
    type: string
    default: 'total_list.csv'
    inputBinding:
      position: 4
      prefix: -o
  pcr:
    type: string
    default: 'standard'
    inputBinding:
      position: 6
      prefix: -pcr
  standard_pcr_input:
    type: string
    default: 'standard_pcr_in.txt'
    inputBinding:
      position: 7
      prefix: -spcr

outputs:
  total_primers_list:
    type: File
    outputBinding:
      glob: $(inputs.outfile)
  pcr_input:
    type: File
    outputBinding:
      glob: $(inputs.standard_pcr_input)

