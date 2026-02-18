#!/bin/sh
# This script is required to compatibilize with the driver script
# that is used in the QMMM driver tests, testscript-qmmm.sh
#
# It is a simple wrapper that runs the driver script with the
# provided input and output files. It also cleans up the
# output files before running the driver script.
#
# Since it's intended for MM-only tests, it doesn't need to
# call SIESTA as a parallel process.

# Excecutable for the QMMM driver.
driver=$1
# FDF input for the driver.
driverfdf=$2
# Output file for the driver.
driverout=$3

echo "Running as:"
echo "${driver} < ${driverfdf} > ${driverout}"

rm -f ${driverout} *.lst *.XV
${driver} < ${driverfdf} > ${driverout}