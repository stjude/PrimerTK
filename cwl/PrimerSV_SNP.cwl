#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
id: "generates primers around deletion or inversion SVs"

requirements:
 - class: ScatterFeatureRequirement

inputs:
  ref_genome: File
  regions_file: File
  primer_opt_size:
    type: int
    default: 22
  primer_min_size:
    type: int
    default: 18
  primer_max_size:
    type: int
    default: 26
  primer_opt_gc:
    type: int
    default: 50
  primer_min_gc:
    type: int
    default: 20
  primer_max_gc:
    type: int
    default: 80
  primer_opt_tm:
    type: int
    default: 60
  primer_min_tm:
    type: int
    default: 57
  primer_max_tm:
    type: int
    default: 63
  product_size_range:
    type: string
    default: '200-400'
  flanking_region_size:
    type: int
    default: 200
  sequence_target:
    type: string
    default: '199,1'
  mispriming_library: string
  thermodynamics_path: string
  sv_type: string
  output:
    type: string
    default: primer3_dump.txt
  outfile: string
  all_primers: string
  top_primers: string
  plate_basename: string
  outfile_snp: string
  vcf_in:
    type: File
    secondaryFiles:
     - .tbi
outputs:
  flank_file:
    type: File
    outputSource: genome_iterator_sv/flanking_regions_file
  primer3_input_file:
    type: File
    outputSource: genome_iterator_sv/primer3_input_file
  primer_dump:
    type: File
    outputSource: primer3/primer_dump_file
  total_primer_list:
    type: File
    outputSource: pcr_setup/total_primers_list
  all_primers_output:
    type: File
    outputSource: pcr_sim/all_primer_info
  top_ranking_primer:
    type: File
    outputSource: pcr_sim/top_ranking_primers
  plate_primers:
    type: File[]
    outputSource: pcr_sim/plated_primers
  snp_file:
    type: File
    outputSource: primer_tabix/snp_indexed_file

steps:
  genome_iterator_sv:
    run: ./tools/genome_iterator_sv.cwl
    in:
      ref_genome: ref_genome
      regions_file: regions_file
      primer_opt_size: primer_opt_size
      primer_min_size: primer_min_size
      primer_max_size: primer_max_size
      primer_opt_gc: primer_opt_gc
      primer_min_gc: primer_min_gc
      primer_max_gc: primer_max_gc
      primer_opt_tm: primer_opt_tm
      primer_min_tm: primer_min_tm
      primer_max_tm: primer_max_tm
      product_size_range: product_size_range
      flanking_region_size: flanking_region_size
      sequence_target: sequence_target
      mispriming_library: mispriming_library
      thermodynamics_path: thermodynamics_path
      sv_type: sv_type
    out: [flanking_regions_file, primer3_input_file]
  primer3:
    run: ./tools/primer3.cwl
    in:
      primer3_input_file: genome_iterator_sv/primer3_input_file
      output: output
    out: [primer_dump_file]
  pcr_setup:
    run: ./tools/pcr_gen_sv.cwl
    in:
      primer_dump: primer3/primer_dump_file
      outfile: outfile
    out: [total_primers_list, pcr_input]
  pcr_sim:
    run: ./tools/post_pcr_analysis_sv.cwl
    in:
      flank_file: genome_iterator_sv/flanking_regions_file
      total_primers: pcr_setup/total_primers_list
      all_primers: all_primers
      top_primers: top_primers
      plate_basename: plate_basename
    out: [all_primer_info, top_ranking_primers, plated_primers]
  primer_tabix:
    run: ./tools/primer_tabix.cwl
    in:
      vcf_in: vcf_in
      primer_pipeline_output: pcr_sim/all_primer_info
      outfile_snp: outfile_snp
    out: [snp_indexed_file]