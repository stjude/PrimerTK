#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "gets specific region GC (no amp)"

baseCommand: python3.6
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

outputs:
  total_output_gc:
    type: File
    outputBinding:
      glob: 'total_list_gc.csv'
  top_ranking_primers:
    type: File
    outputBinding:
      glob: 'top_ranked_final_primers.csv'
