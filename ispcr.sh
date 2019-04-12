#!/usr/bin/env bash

chromosome_dir=$1
pcr_input=$2

for x in $chromosome_dir/chr*.fa
do
    isPcr $x $pcr_input stdout >> pcr_output.fa
done

