#!/bin/bash
program=$1
input=$2
workdir=$3
path_To_Inputs=$4
path_To_Executable=$5
output_To_Compare=$6

#Verify if the utility exists
if [ -f $path_To_Executable/$program ]
then
    cd $workdir
	rm *
	cp $path_To_Inputs/* .

	# Execute the utility
	$path_To_Executable/$program < $input > $program.out

	mv $output_To_Compare ${program}_${output_To_Compare}
else
	echo "Utility doesn't exist"
	exit 1
fi
