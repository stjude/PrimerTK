#!/usr/bin/env bash

#########################
# The command line help #
#########################
display_help() {
    echo "
Please provide the full path to the directory containing genome split into chromosome.fa files
as the first command line argument and the pcr input file as the second.
"
    exit 1
}

chromosome_dir=$1
pcr_input=$2

if [[ ! -d $chromosome_dir ]]; then
    >&2 echo "
ERROR: no chromosome directory was specified.
"
    display_help
    exit 1
elif [[ ! -f $pcr_input ]]; then
    >&2 echo "
ERROR: no pcr input file was specified.
"
    display_help
    exit 1
else
    echo "Running In Silico PCR on pcr input positions..."
fi

for x in $chromosome_dir/*.fa
do
    isPcr $x $pcr_input stdout >> pcr_output.fa
done

