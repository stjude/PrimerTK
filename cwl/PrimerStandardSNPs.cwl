#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
id: "generates all standard PCR products from a given set of inputs"

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
  pcr: string
  chromosome_fasta: File[]
  catted_filename: string
  vcf_in:
    type: File
    secondaryFiles:
     - .tbi
  outfile_snp: string

outputs:
  flank_file:
    type: File
    outputSource: genome_iterator/flanking_regions_file
  primer3_input_file:
    type: File
    outputSource: genome_iterator/primer3_input_file
  primer_dump:
    type: File
    outputSource: primer3/primer_dump_file
  total_primer_list:
    type: File
    outputSource: standard_pcr_gen/total_primers_list
  catted_output:
    type: File
    outputSource: combine_pcr_outputs/pcr_combined
  all_products_info:
    type: File
    outputSource: post_pcr_analysis/all_product_info
  top_ranked_primers:
    type: File
    outputSource: post_pcr_analysis/top_ranked_primers
  snp_indexed_file:
    type: File
    outputSource: primer_tabix/snp_indexed_file

steps:
  genome_iterator:
    run: ./tools/genome_iterator.cwl
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
    out: [flanking_regions_file, primer3_input_file]
  primer3:
    run: ./tools/primer3.cwl
    in:
      primer3_input_file: genome_iterator/primer3_input_file
    out: [primer_dump_file]
  standard_pcr_gen:
    run: ./tools/standard_pcr_gen.cwl
    in:
      primer_dump: primer3/primer_dump_file
      pcr: pcr
    out: [total_primers_list, pcr_input]
  is_pcr:
    run: ./tools/is_pcr.cwl
    scatter: chromosome_fasta
    in:
      chromosome_fasta: chromosome_fasta
      pcr_input: standard_pcr_gen/pcr_input
    out: [chromosome_amplification]
  combine_pcr_outputs:
    run: ./tools/combine_pcr_output.cwl
    in:
      catted_filename: catted_filename
      chromosome_fastas: is_pcr/chromosome_amplification
    out: [pcr_combined]
  post_pcr_analysis:
    run: ./tools/post_pcr_analysis.cwl
    in:
      pcr_output: combine_pcr_outputs/pcr_combined
      total_primers: standard_pcr_gen/total_primers_list
    out: [all_product_info, filtered_good_primers, top_ranked_primers]
  primer_tabix:
    run: ./tools/primer_tabix.cwl
    in:
      vcf_in: vcf_in
      primer_pipeline_output: post_pcr_analysis/filtered_good_primers
      outfile_snp: outfile_snp
    out: [snp_indexed_file]