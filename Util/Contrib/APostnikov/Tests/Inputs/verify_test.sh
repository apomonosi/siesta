#!/bin/bash

program=$1
workdir=$2
output_To_Compare=$3
ref_dir=$4

cd $workdir

outname="${program}_${output_To_Compare}"

# Verify if the output reference exists
if [ -f $ref_dir/$outname ]
then
    # Get the diffs side-by-side, ignoring whitespace
    diff -y -w --suppress-common-lines $outname \
	    $ref_dir/$outname > OUT.diffs
else
    echo "$ref_dir/$outname doesn't exist"
    exit 1
fi

# Erase the diff file if it is empty
if [ ! -s OUT.diffs ]
then
	echo "Test passed for $program"
	rm -f OUT.diffs
	exit 0
else
	echo "$program: ** DIFFERENCES FOUND **"
	echo "DIFFS"
	cat OUT.diffs
	exit 1
fi