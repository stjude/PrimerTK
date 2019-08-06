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

@test "genome_iterator_sv_insertion" {
    CWL_SCRIPT="../../../../cwl/tools/genome_iterator_sv.cwl"
    CWL_INPUT="../inputs/sv/genome_iterator_sv_insertion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer3_input.*.txt > $AT
    cat ../data/sv/primer3_input.insertion.txt > $GS
    cmp $AT $GS
    rm -rf primer3_input.*.txt
    rm -rf flanking_regions.*.fasta
}
@test "genome_iterator_sv_deletion" {
    CWL_SCRIPT="../../../../cwl/tools/genome_iterator_sv.cwl"
    CWL_INPUT="../inputs/sv/genome_iterator_sv_deletion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer3_input.*.txt > $AT
    cat ../data/sv/primer3_input.deletion.txt > $GS
    cmp $AT $GS
    rm -rf primer3_input.*.txt
    rm -rf flanking_regions.*.fasta
}
@test "genome_iterator_sv_inversion" {
    CWL_SCRIPT="../../../../cwl/tools/genome_iterator_sv.cwl"
    CWL_INPUT="../inputs/sv/genome_iterator_sv_inversion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer3_input.*.txt > $AT
    cat ../data/sv/primer3_input.inversion.txt > $GS
    cmp $AT $GS
    rm -rf primer3_input.*.txt
    rm -rf flanking_regions.*.fasta
}
@test "primer3 insertion" {
    CWL_SCRIPT="../../../../cwl/tools/primer3.cwl"
    CWL_INPUT="../inputs/sv/primer3_insertion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer_dump.txt > $AT
    cat ../data/sv/primer_dump_insertion.txt > $GS
    cmp $AT $GS
    rm -rf primer_dump.txt
}
@test "primer3 deletion" {
    CWL_SCRIPT="../../../../cwl/tools/primer3.cwl"
    CWL_INPUT="../inputs/sv/primer3_deletion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer_dump.txt > $AT
    cat ../data/sv/primer_dump_deletion.txt > $GS
    cmp $AT $GS
    rm -rf primer_dump.txt
}
@test "primer3 inversion" {
    CWL_SCRIPT="../../../../cwl/tools/primer3.cwl"
    CWL_INPUT="../inputs/sv/primer3_inversion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer_dump.txt > $AT
    cat ../data/sv/primer_dump_inversion.txt > $GS
    cmp $AT $GS
    rm -rf primer_dump.txt
}
@test "pre pcr insertion" {
    CWL_SCRIPT="../../../../cwl/tools/pcr_gen_sv.cwl"
    CWL_INPUT="../inputs/sv/pre_pcr_sv_insertion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat total_primers_insertion.csv > $AT
    cat ../data/sv/total_primers_insertion.csv > $GS
    cmp $AT $GS
    rm -rf total_primers_insertion.csv
    rm -rf standard_pcr.txt
}
@test "pre pcr deletion" {
    CWL_SCRIPT="../../../../cwl/tools/pcr_gen_sv.cwl"
    CWL_INPUT="../inputs/sv/pre_pcr_sv_deletion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat total_primers_deletion.csv > $AT
    cat ../data/sv/total_primers_deletion.csv > $GS
    cmp $AT $GS
    rm -rf total_primers_deletion.csv
    rm -rf standard_pcr.txt
}
@test "pre pcr inversion" {
    CWL_SCRIPT="../../../../cwl/tools/pcr_gen_sv.cwl"
    CWL_INPUT="../inputs/sv/pre_pcr_sv_inversion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat total_primers_inversion.csv > $AT
    cat ../data/sv/total_primers_inversion.csv > $GS
    cmp $AT $GS
    rm -rf total_primers_inversion.csv
    rm -rf standard_pcr.txt
}
@test "post pcr insertion" {
    skip
    CWL_SCRIPT="../../../../cwl/tools/post_pcr_analysis_sv.cwl"
    CWL_INPUT="../inputs/sv/post_pcr_sv_insertion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat all_primers.csv > $AT
    cat ../data/sv/all_primers_insertion.csv > $GS
    cmp $AT $GS
    cat top_primers.csv > $AT
    cat ../data/sv/top_primers_insertion.csv > $GS
    cmp $AT $GS
    rm -rf all_primers.csv
    rm -rf top_primers.csv
    rm -rf plated_primers_F.csv
    rm -rf plated_primers_R.csv
}
@test "post pcr deletion" {
    CWL_SCRIPT="../../../../cwl/tools/post_pcr_analysis_sv.cwl"
    CWL_INPUT="../inputs/sv/post_pcr_sv_deletion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat all_primers.csv > $AT
    cat ../data/sv/all_primers_deletion.csv > $GS
    cmp $AT $GS
    cat top_primers.csv > $AT
    cat ../data/sv/top_primers_deletion.csv > $GS
    cmp $AT $GS
    rm -rf all_primers.csv
    rm -rf top_primers.csv
    rm -rf plated_primers_F.csv
    rm -rf plated_primers_R.csv
}
@test "post pcr inversion" {
    CWL_SCRIPT="../../../../cwl/tools/post_pcr_analysis_sv.cwl"
    CWL_INPUT="../inputs/sv/post_pcr_sv_inversion.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat all_primers.csv > $AT
    cat ../data/sv/all_primers_inversion.csv > $GS
    cmp $AT $GS
    cat top_primers.csv > $AT
    cat ../data/sv/top_primers_inversion.csv > $GS
    cmp $AT $GS
    rm -rf all_primers.csv
    rm -rf top_primers.csv
    rm -rf plated_primers_F.csv
    rm -rf plated_primers_R.csv
}
@test "post pcr insertion fails" {
    CWL_SCRIPT="../../../../cwl/tools/post_pcr_analysis_sv.cwl"
    CWL_INPUT="../inputs/sv/post_pcr_sv_deletion.yml"
    run cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    [ "$output" = "No primers found for any targets, try again!" ]
}
