#!/bin/sh
#
# This script is used as a wrapper for the QMMM driver and SIESTA, with the
# intention of being called from CMake/Ctest. It receives everything it needs
# as arguments, including the names for the QMMM driver and SIESTA executables.
#
# It also cleans up the output files before running the driver script and SIESTA.
#
# See the testing function in the CMakeLists.txt file for more details.

# Excecutable for the QMMM driver.
driver=$1
# FDF input for the driver.
driverfdf=$2
# Output file for the driver.
driverout=$3

# Path to the SIESTA pseudopotential files.
siestaps=$4
# Excecutable for SIESTA.
siesta=$5
# FDF input for SIESTA.
siestafdf=$6
# Output file for SIESTA.
siestaout=$7

# Preprocessing commands for SIESTA.
siestapre=$8
# Postprocessing commands for SIESTA.
siestapost=$9

echo "Running as:"
echo "${driver} < ${driverfdf} > ${driverout} &"
echo "SIESTA_PS_PATH=${siestaps} ${siestapre} ${siesta} -out ${siestaout} ${siestafdf} ${siestapost}"

# Clean up output files before running the driver script.
rm -f ${driverout} *.lst *.XV
# Run the driver script and SIESTA in parallel.
${driver} < ${driverfdf} > ${driverout} &
SIESTA_PS_PATH=${siestaps} ${siestapre} ${siesta} -out ${siestaout} ${siestafdf} ${siestapost}