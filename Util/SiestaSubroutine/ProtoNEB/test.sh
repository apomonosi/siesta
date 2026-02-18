#!/bin/sh

# NOTE: This script is of historical value only.
#       With the switch to CMake, and the support
#       for the SIESTA_PS_PATH environment variable,
#       the idioms are different (and maybe simpler)
#
# Select the appropriate run below (at the end).
#
# Make sure that you are using the right version of Siesta.
# The SIESTABIN setting below is quite naive and might not work
# in all cases. You can call this script as:
#
#            SIESTABIN=/path/to/siesta/bin test.sh
#
ROOT="../../../.."
PSEUDOS=${ROOT}/Tests/Pseudos
#
if [ -z "$SIESTABIN" ] ; then
      SIESTABIN=${ROOT}/bin/
      SIESTA=${ROOT}/bin/siesta
else
      SIESTA=${SIESTABIN}/siesta
fi
echo "Using Siesta executable: $SIESTA"
#

#rm -r work
if [ -d work ] ; then
   echo "Work directory 'work' exists. Please delete it"
   exit
fi

#
mkdir work
cd work
cp -p ../*.fdf .
cp ${PSEUDOS}/H.psf  .
cp ${PSEUDOS}/N.psf  .
ln -sf ${SIESTA} ./siesta
#

echo ""; echo "Protoneb"
mpirun -np 3 ${SIESTABIN}/SIESTA_protoNEB_exe  | tee protoneb.out
cat protoneb.out


