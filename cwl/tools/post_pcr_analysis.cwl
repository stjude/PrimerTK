#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "post pcr analysis"

baseCommand: [primer_tk, post]
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
  pcr_product_info:
    type: string
    default: 'pcr_product_info.csv'
    inputBinding:
      position: 5
      prefix: -pcri
  all_good_primers:
    type: string
    default: 'all_primers.csv'
    inputBinding:
      position: 6
      prefix: -all
  top_primer_info:
    type: string
    default: 'top_primers.csv'
    inputBinding:
      position: 7
      prefix: -top
  plate_basename:
    type: string
    default: 'plated_primers'
    inputBinding:
      position: 8
      prefix: -plate

outputs:
  all_product_info:
    type: File
    outputBinding:
      glob: $(inputs.pcr_product_info)
  filtered_good_primers:
    type: File
    outputBinding:
      glob: $(inputs.all_good_primers)
  top_ranked_primers:
    type: File
    outputBinding:
      glob: $(inputs.top_primer_info)
  plated_primers:
    type: File[]
    outputBinding:
      glob: '$(inputs.plate_basename)*'
