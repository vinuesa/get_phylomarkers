#!/usr/bin/env bash

#: version 2021-09-16 vivaMX
#: PROGRAM: run_make_tests.sh
#: AIM: runs make on the clean and test sections, the latter found in test_get_phylomarkers.t

# make sure we start tests with clean data directories
make clean

# run tests
make test

# cleanup again
make clean
