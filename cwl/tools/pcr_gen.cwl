#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "main pre pcr setup"

baseCommand: python3.6

inputs:
  main_pre_pcr:
    type: File
    inputBinding:
      position: 1
    default:
      class: File
      location: ../../../src/main_pre_pcr.py

  primer_dump:
    type: File
    default: primer_dump.txt
    inputBinding:
      position: 2
      prefix: -d
  outfile:
    type: string
    default: 'total_list.csv'
    inputBinding:
      position: 3
      prefix: -o
  percent_alignment:
    type: int
    default: 60
    inputBinding:
      position: 4
      prefix: -pa
  pcr:
    type: string    
    inputBinding:
      position: 5
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

