#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
id: "generates primers around deletion or inversion SVs"

requirements:
 - class: ScatterFeatureRequirement

inputs:
  ref_genome: File
  regions_file: File
  primer_opt_size: int
  primer_min_size: int
  primer_max_size: int
  primer_opt_gc: int
  primer_min_gc: int
  primer_max_gc: int
  primer_opt_tm: int
  primer_min_tm: int
  primer_max_tm: int
  product_size_range: string
  flanking_region_size: int
  sequence_target: string
  mispriming_library: string
  thermodynamics_path: string
  sv_type: string
  output: string
  outfile: string
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
  total_outputs_gc:
    type: File
    outputSource: pcr_sim/total_output_gc
  top_ranking_primer:
    type: File
    outputSource: pcr_sim/top_ranking_primers
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
    run: ./tools/pcr_setup.cwl
    in:
      primer_dump: primer3/primer_dump_file
      outfile: outfile
    out: [total_primers_list, pcr_standard_output]
  pcr_sim:
    run: ./tools/pcr_sim.cwl
    in:
      flank_file: genome_iterator_sv/flanking_regions_file
      total_primers: pcr_setup/total_primers_list
    out: [total_output_gc, top_ranking_primers]
  primer_tabix:
    run: ./tools/primer_tabix.cwl
    in:
      vcf_in: vcf_in
      primer_pipeline_output: pcr_sim/total_output_gc
      outfile_snp: outfile_snp
    out: [snp_indexed_file]