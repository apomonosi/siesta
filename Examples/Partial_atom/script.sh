#!/bin/sh

# Run "partial atom" example
#
# Argument: SIESTA execution string
# The script is passed the (probably relative) path to the siesta
# executable, and maybe with a "mpirun" prefix
SIESTA="$1"

# Extract last component of the executable, in case of mpirun-style string
REL_PATH=$(echo ${SIESTA} | awk '{print $NF}')
EXEC_PREFIX=$(echo ${SIESTA} | awk '{$NF=""; print}')
REL_PATH=$(which ${REL_PATH})
NAME=$(basename ${REL_PATH})
EXEC_DIR=$(dirname ${REL_PATH})

# Find absolute path.
pushd ${EXEC_DIR} ; ABS_EXEC_DIR=$(pwd) ; popd
ABS=${ABS_EXEC_DIR}/${NAME}

echo "Running script with SIESTA=$EXEC_PREFIX $ABS"
fractional=${ABS_EXEC_DIR}/fractional

mkdir work
cd    work
ln -s ../O.psf .

# We first run "fractional" to create the new pseudopotential file.
echo "==> Running $fractional"
$fractional O  0.5

# And then run siesta
echo "
SystemName          Oxypartial
SystemLabel         oxypartial
NumberOfAtoms       2
NumberOfSpecies     1

%block ChemicalSpeciesLabel
 1  201  O-Fraction-0.50000
%endblock ChemicalSpeciesLabel
%block SyntheticAtoms
 1  
 2 2 3 4
 1.0 2.0 0.0 0.0
%endblock SyntheticAtoms

%block PAO.basis
O-Fraction-0.50000  2
 n=2   0   2
 0.0 0.0   
   1.000      1.000   
 n=2   1   2 P   1
   4.139      2.740   
   1.000      1.000   
%endblock PAO.Basis

Spin        Polarized
MeshCutoff  200 Ry
SCF.Mixer.History 4

AtomicCoordinatesFormat  Ang
%block AtomicCoordinatesAndAtomicSpecies
 0.000  0.000  0.000  1
 0.000  0.000  1.200  1
%endblock AtomicCoordinatesAndAtomicSpecies

MD.TypeOfRun Broyden
MD.Steps     40
" > Job.fdf

echo "==> Running $SIESTA"
${ABS_EXEC_DIR}/siesta < Job.fdf > Job.out

cp Job.out ../partial-0.5.out

