#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "gets specific region GC (no amp)"

inputs:
  main_post_pcr:
    type: File
    inputBinding:
      position: 1
    default:
      class: File
      location: ../../src/main_post_pcr.py

  flank_file:
    type: File
    inputBinding:
      position: 2
      prefix: -f
  total_primers:
    type: File
    inputBinding:
      position: 3
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
