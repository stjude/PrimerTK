#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "main pre pcr setup"

baseCommand: python3.6
arguments:
 - valueFrom: primer_tk
   position: 1
   prefix: -m
 - valueFrom: pre
   position: 2

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

outputs:
  total_primers_list:
    type: File
    outputBinding:
      glob: 'total_list.csv'
  pcr_input:
    type: File
    outputBinding:
      glob: 'standard_pcr.txt'

