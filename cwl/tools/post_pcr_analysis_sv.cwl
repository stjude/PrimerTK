#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "gets specific region GC (no amp)"

baseCommand: python3
arguments:
 - valueFrom: primer_tk
   position: 1
   prefix: -m
 - valueFrom: post_sv
   position: 2

inputs:
  flank_file:
    type: File
    inputBinding:
      position: 3
      prefix: -f
  total_primers:
    type: File
    inputBinding:
      position: 4
      prefix: -tp
  all_primers:
    type: string
    default: 'all_primers.csv'
    inputBinding:
      position: 5
      prefix: -all
  top_primers:
    type: string
    default: 'top_primers.csv'
    inputBinding:
      position: 6
      prefix: -top
  plate_basename:
    type: string
    default: 'plated_primers'
    inputBinding:
      position: 7
      prefix: -plate

outputs:
  all_primer_info:
    type: File
    outputBinding:
      glob: $(inputs.all_primers)
  top_ranking_primers:
    type: File
    outputBinding:
      glob: $(inputs.top_primers)
  plated_primers:
    type: File[]
    outputBinding:
      glob: '$(inputs.plate_basename)*'
