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
cp ${PSEUDOS}/O.psf  .
ln -sf ${SIESTA} ./siesta
#

echo ""; echo "simple_pipes_serial"
${SIESTABIN}/pipes_serial    | tee simple_pipes_serial.out
mv h2o.out siesta_pipes_serial.out
echo ""; echo "simple_pipes_parallel"
${SIESTABIN}/pipes_parallel  | tee simple_pipes_parallel.out
mv h2o.out siesta_pipes_parallel.out

echo ""; echo "simple_mpi_parallel"
mpirun -np 2 ${SIESTABIN}/mpi_driver | tee simple_mpi_parallel.out

cat socket.fdf >> h2o.fdf
echo ""; echo "simple_sockets_serial"
${SIESTABIN}/sockets_serial    | tee simple_sockets_serial.out
mv h2o.out siesta_sockets_serial.out
echo ""; echo "simple_sockets_parallel"
${SIESTABIN}/sockets_parallel  | tee simple_sockets_parallel.out
mv h2o.out siesta_sockets_parallel.out

