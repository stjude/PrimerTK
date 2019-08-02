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
  percent_alignment:
    type: int
    default: 60
    inputBinding:
      position: 5
      prefix: -pa
  pcr:
    type: string
    inputBinding:
      position: 6
      prefix: -pcr
  no_dimer:
    type: string
    default: 'no_dimers.csv'
    inputBinding:
      position: 7
      prefix: -nd
  multiplex_pcr_infile:
    type: string
    default: 'multiplex_pcr_in.txt'
    inputBinding:
      position: 8
      prefix: -mpcr

outputs:
  total_primers_list:
    type: File
    outputBinding:
      glob: $(inputs.outfile)
  pcr_multiplex_input:
    type: File
    outputBinding:
      glob: $(inputs.multiplex_pcr_infile)
  no_dimer_out:
    type: File
    outputBinding:
      glob: $(inputs.no_dimer)
