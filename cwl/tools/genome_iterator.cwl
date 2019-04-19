#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "primer3 input file generator"

baseCommand: python3.6
arguments:
 - valueFrom: primer_tk
   position: 1
   prefix: -m
 - valueFrom: iterator
   position: 2
inputs:
  ref_genome:
    type: File
    inputBinding:
      position: 3
      prefix: -ref
  regions_file:
    type: File
    inputBinding:
      position: 4
      prefix: -in
  primer_opt_size:
    type: int
    default: 22
    inputBinding:
      position: 5
      prefix: -opt_size
  primer_min_size:
    type: int
    default: 18
    inputBinding:
      position: 6
      prefix: -min_size
  primer_max_size:
    type: int
    default: 25
    inputBinding:
      position: 7
      prefix: -max_size
  primer_opt_gc:
    type: int
    default: 50
    inputBinding:
      position: 8
      prefix: -opt_gc
  primer_min_gc:
    type: int
    default: 20
    inputBinding:
      position: 9
      prefix: -min_gc
  primer_max_gc:
    type: int
    default: 80
    inputBinding:
      position: 10
      prefix: -max_gc
  primer_opt_tm:
    type: int
    default: 60
    inputBinding:
      position: 11
      prefix: -opt_tm
  primer_min_tm:
    type: int
    default: 57
    inputBinding:
      position: 12
      prefix: -min_tm
  primer_max_tm:
    type: int
    default: 63
    inputBinding:
      position: 13
      prefix: -max_tm
  product_size_range:
    type: string
    default: '200-400'
    inputBinding:
      position: 14
      prefix: -sr
  flanking_region_size:
    type: int
    default: 200
    inputBinding:
      position: 15
      prefix: -flank
  sequence_target:
    type: string
    default: '199,1'
    inputBinding:
      position: 16
      prefix: -st
  mispriming_library:
    type: string
    inputBinding:
      position: 17
      prefix: -mp
  thermodynamics_path:
    type: string
    inputBinding:
      position: 18
      prefix: -tp

outputs:
  flanking_regions_file:
    type: File
    outputBinding:
      glob: '*.fasta'
  primer3_input_file:
    type: File
    outputBinding:
      glob: '*.txt'
