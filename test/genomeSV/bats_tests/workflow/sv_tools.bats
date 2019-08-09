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

@test "PrimerSV Inversion Workflow" {
    CWL_SCRIPT="../../../../cwl/PrimerSV.cwl"
    CWL_INPUT="../inputs/sv/wf/inversion_sv_wf.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer3_input.*.txt > $AT
    cat ../data/sv/wf/primer3_input.inversion.txt > $GS
    cmp $AT $GS
    cat primer_dump_inversion.txt > $AT
    cat ../data/sv/wf/primer_dump_inversion.txt > $GS
    cmp $AT $GS
    cat total_inversion_list.csv > $AT
    cat ../data/sv/wf/total_inversion_list.csv > $GS
    cmp $AT $GS
    cat all_inversion_primers.csv > $AT
    cat ../data/sv/wf/all_inversion_primers.csv > $GS
    cmp $AT $GS
    cat top_inversion_primers.csv > $AT
    cat ../data/sv/wf/top_inversion_primers.csv > $GS
    cmp $AT $GS
    cat plate_inversion_primers_F.csv > $AT
    cat ../data/sv/wf/plate_inversion_primers_F.csv > $GS
    cmp $AT $GS
    cat plate_inversion_primers_R.csv > $AT
    cat ../data/sv/wf/plate_inversion_primers_R.csv > $GS
    cmp $AT $GS
    rm -rf primer3_input.*.txt
    rm -rf flanking_regions.*.fasta
    rm -rf primer_dump_inversion.txt
    rm -rf total_inversion_list.csv
    rm -rf all_inversion_primers.csv
    rm -rf top_inversion_primers.csv
    rm -rf plate_inversion_primers_F.csv
    rm -rf plate_inversion_primers_R.csv
}
@test "PrimerSV Deletion Workflow" {
    CWL_SCRIPT="../../../../cwl/PrimerSV.cwl"
    CWL_INPUT="../inputs/sv/wf/deletion_sv_wf.yml"
    cwltool --preserve-entire-environment $CWL_SCRIPT $CWL_INPUT
    cat primer3_input.*.txt > $AT
    cat ../data/sv/wf/primer3_input.deletion.txt > $GS
    cmp $AT $GS
    cat primer_dump_deletion.txt > $AT
    cat ../data/sv/wf/primer_dump_deletion.txt > $GS
    cmp $AT $GS
    cat total_deletion_list.csv > $AT
    cat ../data/sv/wf/total_deletion_list.csv > $GS
    cmp $AT $GS
    cat all_deletion_primers.csv > $AT
    cat ../data/sv/wf/all_deletion_primers.csv > $GS
    cmp $AT $GS
    cat top_deletion_primers.csv > $AT
    cat ../data/sv/wf/top_deletion_primers.csv > $GS
    cmp $AT $GS
    cat plate_deletion_primers_F.csv > $AT
    cat ../data/sv/wf/plate_deletion_primers_F.csv > $GS
    cmp $AT $GS
    cat plate_deletion_primers_R.csv > $AT
    cat ../data/sv/wf/plate_deletion_primers_R.csv > $GS
    cmp $AT $GS
    rm -rf primer3_input.*.txt
    rm -rf flanking_regions.*.fasta
    rm -rf primer_dump_deletion.txt
    rm -rf total_deletion_list.csv
    rm -rf all_deletion_primers.csv
    rm -rf top_deletion_primers.csv
    rm -rf plate_deletion_primers_F.csv
    rm -rf plate_deletion_primers_R.csv
}
