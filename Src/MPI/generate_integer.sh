#!/bin/sh
### set -x
dir=$(dirname "$0")

echo " ===> Generating integer module files from templates..."

if [ -z "$@" ] ; then
  KINDS=$(./mpi_int_explorer | sed -e 's/,/ /g')
else
  KINDS=$@
fi

rm -rf V_S_integer.uses VS_integer.uses Interfaces_integer.f90

echo "Used integer kinds: $KINDS"

for kind in ${KINDS} ; do

  echo "        USE MPI__i${kind}_V  ; USE MPI__i${kind}_S" >> V_S_integer.uses
  echo "        USE MPI__i${kind}_VS ; USE MPI__i${kind}_SV" >> VS_integer.uses

  for tag in v s sv vs ; do

    sed -e "/_type/s//_i${kind}/" -e "/type/s//integer(${kind})/" \
        "${dir}/mpi__type_${tag}.f90" >> Interfaces_integer.f90

  done

done
