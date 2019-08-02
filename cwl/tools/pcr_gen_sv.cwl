#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "main pre pcr setup"

baseCommand: [primer_tk, pre_sv]
inputs:
  primer_dump:
    type: File
    default: primer_dump.txt
    inputBinding:
      position: 3
      prefix: -d
  outfile:
    type: string
    default: 'total_list_sv.csv'
    inputBinding:
      position: 4
      prefix: -o

outputs:
  total_primers_list:
    type: File
    outputBinding:
      glob: '*.csv'
  pcr_input:
    type: File
    outputBinding:
      glob: 'standard_pcr.txt'

