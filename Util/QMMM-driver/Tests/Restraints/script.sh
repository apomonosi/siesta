#!/bin/bash

name=restraints

mkdir work
cd work

export PS_PATH_HINT="../../../../../Tests/Pseudos"
. ../../../../../Tests/set_siesta_dir.sh "$1" $2

DRIVER=${ABS_EXEC_DIR}/siesta_qmmm

echo "Running script with SIESTA=$SIESTA"

cp ../amber.parm .

$DRIVER < ../$name.fdf > $name.out &
$SIESTA < ../$name.siestain.fdf > $name.siesta.out

