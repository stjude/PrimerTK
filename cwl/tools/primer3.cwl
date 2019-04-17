#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "primer3 pipeline tool"

baseCommand: primer3_core

hints:
  SoftwareRequirement:
    packages:
      primer3:
        version: [ "2.4.0" ]

inputs:
  output:
    type: string
    default: "primer_dump.txt"
    inputBinding:
      position: 1
      prefix: --output=
      separate: false
  primer3_input_file:
    type: File
    inputBinding:
      position: 2

outputs:
  primer_dump_file:
    type: File
    outputBinding:
      glob: $(inputs.output)

