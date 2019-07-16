#!/usr/bin/env bats
load utilities

setup() {
    cd $BATS_TEST_DIRNAME
    GS=$(mktemp)
    AT=$(mktemp)
}
teardown() {
    cwl_teardown
    rm $GS
    rm $AT
}

@test "genome_iterator" {
    CWL_SCRIPT="../../../../cwl/tools/genome_iterator.cwl"
    CWL_INPUT="../inputs/clt/genome_iterator.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer3_input.*.txt > $AT
    cat ../data/clt/primer3_input.standard.txt > $GS
    cmp $AT $GS
    rm -rf primer3_input.*.txt
    rm -rf flanking_regions.*.fasta
}
@test "standard pcr" {
    CWL_SCRIPT="../../../../cwl/tools/standard_pcr_gen.cwl"
    CWL_INPUT="../inputs/clt/standard_pcr_gen.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat standard_pcr_in.txt > $AT
    cat ../data/clt/standard_pcr_in.txt > $GS
    cmp $AT $GS
    cat total_standard_primers.csv > $AT
    cat ../data/clt/total_standard_primers.csv > $GS
    rm -rf standard_pcr_in.txt
    rm -rf total_standard_primers.csv
}
@test "multiplex pcr" {
    CWL_SCRIPT="../../../../cwl/tools/multiplex_pcr_gen.cwl"
    CWL_INPUT="../inputs/clt/multiplex_pcr_gen.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat multiplex_pcr_in.txt > $AT
    cat ../data/clt/multiplex_pcr_input.txt > $GS
    cmp $AT $GS
    cat total_multiplex_primers.csv > $AT
    cat ../data/clt/total_multiplex_primers.csv > $GS
    cmp $AT $GS
    cat no_dimer_primers.csv > $AT
    cat ../data/clt/no_dimers_primers.csv > $GS
    rm -rf multiplex_pcr_in.txt
    rm -rf total_multiplex_primers.csv
    rm -rf Primer_Dimers.txt
    rm -rf no_dimer_primers.csv
}
@test "isPcr" {
    CWL_SCRIPT="../../../../cwl/tools/is_pcr.cwl"
    CWL_INPUT="../inputs/clt/ispcr.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat test_primeramp.fa > $AT
    cat ../data/clt/test_primeramp.fa > $GS
    cmp $AT $GS
    rm -rf test_primeramp.fa
}
@test "combine pcr output" {
    CWL_SCRIPT="../../../../cwl/tools/combine_pcr_output.cwl"
    CWL_INPUT="../inputs/clt/combine_pcr_output.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat pcr_output.fa > $AT
    cat ../data/clt/pcr_output.fa > $GS
    cmp $AT $GS
    rm -rf pcr_output.fa
}
@test "post pcr analysis" {
    CWL_SCRIPT="../../../../cwl/tools/post_pcr_analysis.cwl"
    CWL_INPUT="../inputs/clt/post_pcr_analysis.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat pcr_product_info.csv > $AT
    cat ../data/clt/pcr_product_info.csv > $GS
    cmp $AT $GS
    cat all_primers.csv > $AT
    cat ../data/clt/all_primers.csv > $GS
    cmp $AT $GS
    cat top_primers.csv > $AT
    cat ../data/clt/top_primers.csv > $GS
    cmp $AT $GS
    cat plated_primers_F.csv > $AT
    cat ../data/clt/plated_primers_F.csv > $GS
    cmp $AT $GS
    cat plated_primers_R.csv > $AT
    cat ../data/clt/plated_primers_R.csv > $GS
    cmp $AT $GS
    rm -rf pcr_product_info.csv all_primers.csv top_primers.csv plated_primers_F.csv plated_primers_R.csv
}
@test "primer tabix" {
    CWL_SCRIPT="../../../../cwl/tools/primer_tabix.cwl"
    CWL_INPUT="../inputs/clt/primer_tabix.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat db_snp.csv > $AT
    cat ../data/clt/db_snp.csv > $GS
    cmp $AT $GS
    rm -rf db_snp.csv
}