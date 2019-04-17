#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "post pcr analysis"

baseCommand: python3.6
arguments:
 - valueFrom: primer_tk
   position: 1
   prefix: -m
 - valueFrom: post
   position: 2
inputs:
  pcr_output:
    type: File
    inputBinding:
      position: 3
      prefix: -i
  total_primers:
    type: File
    inputBinding:
      position: 4
      prefix: -tp

outputs:
  all_product_info:
    type: File
    outputBinding:
      glob: 'pcr_product_info.csv'
  filtered_good_primers:
    type: File
    outputBinding:
      glob: 'all_final_primers.csv'
  top_ranked_primers:
    type: File
    outputBinding:
      glob: 'top_final_primers.csv'
  idt_plate_fwd:
    type: File
    outputBinding:
      glob: 'plate_forward_primers.csv'
  idt_plate_rvs:
    type: File
    outputBinding:
      glob: 'plate_reverse_primers.csv'

