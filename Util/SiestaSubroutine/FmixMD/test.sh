#!/bin/sh

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

if [ -d work ] ; then
   echo "Work directory 'work' exists. Please delete it"
   exit
else
   mkdir work
fi

cp -p h2o.fast.fdf h2o.conv.fdf driver.dat work
#
cd work
cp ${PSEUDOS}/H.psf  .
cp ${PSEUDOS}/O.psf  .
#
ln -sf ${SIESTA} ./siesta

${SIESTABIN}/fmixmd-driver < driver.dat | tee driver.out

