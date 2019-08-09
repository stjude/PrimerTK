#!/usr/bin/env bash
# test developed by JRM3

# This script will ALWAYS be run from the directory in which it is locatred
# so as to ensure relative directory structures are intact.
SCRIPT_DIR=$(dirname $0)
#cd $SCRIPT_DIR

usage() {
    echo "run_tests.sh [-h]"
}

####################
### Build source ###
####################
export PYTHONPATH=$(pwd)/../../lib:$PYTHONPATH
export PATH=$(pwd)/../../scripts:$PATH

#########################
### Python unit tests ###
#########################
coverage run --source .. -m unittest discover python_tests -p "test*.py" -b
coverage report --skip-covered -m
