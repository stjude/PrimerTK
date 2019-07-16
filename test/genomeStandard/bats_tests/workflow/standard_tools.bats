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

@test "multiplex primers" {
    CWL_SCRIPT="../../../../cwl/PrimerMultiplex.cwl"
    CWL_INPUT="../inputs/wf/primer_multiplex.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat all_primers.csv flanking_regions.*.fasta pcr_output.fa pcr_product_info.csv plated_primers_F.csv\
    plated_primers_R.csv primer3_input.*.txt primer_multiplex_dump.txt top_primers.csv total_multiplex_primers.csv > $AT
    cat ../data/wf/all_primers.csv ../data/wf/flanking_regions.multiplex.fasta ../data/wf/pcr_output.fa ../data/wf/pcr_product_info.csv\
    ../data/wf/plated_primers_F.csv ../data/wf/plated_primers_R.csv ../data/wf/primer3_input.multiplex.txt ../data/wf/primer_multiplex_dump.txt\
    ../data/wf/top_primers.csv ../data/wf/total_multiplex_primers.csv > $GS
    cmp $AT $GS
    rm -rf all_primers.csv flanking_regions.*.fasta pcr_output.fa pcr_product_info.csv\
    plated_primers_F.csv plated_primers_R.csv primer3_input.*.txt primer_multiplex_dump.txt\
    top_primers.csv total_multiplex_primers.csv
}
@test "standard primers" {
    CWL_SCRIPT="../../../../cwl/PrimerStandard.cwl"
    CWL_INPUT="../inputs/wf/primer_standard.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat total_standard_primers.csv top_standard_primers.csv  primer_standard_dump.txt primer3_input.*.txt\
    plated_primers_standard_F.csv plated_primers_standard_R.csv flanking_regions.*.fasta pcr_standard_output.fa\
    pcr_standard_product_info.csv all_standard_primers.csv > $AT
    cat ../data/wf/total_standard_primers.csv ../data/wf/top_standard_primers.csv ../data/wf/primer_standard_dump.txt\
    ../data/wf/primer3_input.standard.txt ../data/wf/plated_primers_standard_F.csv ../data/wf/plated_primers_standard_R.csv\
    ../data/wf/flanking_regions.standard.fasta ../data/wf/pcr_standard_output.fa ../data/wf/pcr_standard_product_info.csv\
    ../data/wf/all_standard_primers.csv > $GS
    cmp $AT $GS
    rm -rf total_standard_primers.csv top_standard_primers.csv primer_standard_dump.txt primer3_input.*.txt\
    plated_primers_standard_F.csv plated_primers_standard_R.csv flanking_regions.*.fasta pcr_standard_output.fa\
    pcr_standard_product_info.csv all_standard_primers.csv
}