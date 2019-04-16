#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "concatenates all PCR output files"
requirements:
 - class: ShellCommandRequirement

baseCommand: cat

inputs:
  catted_filename:
    type: string

  chromosome_fastas:
    type: File[]
    streamable: true
    inputBinding:
      position: 1

outputs:
  pcr_combined:
    type: File
    streamable: true
    outputBinding:
      glob: $(inputs.catted_filename)

stdout: $(inputs.catted_filename)


