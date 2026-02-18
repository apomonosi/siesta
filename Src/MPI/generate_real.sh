#!/bin/sh
### set -x
dir=$(dirname "$0")

echo " ===> Generating real, logical, character module files from templates..."

if [ -z "$@" ] ; then
  KINDS=$(./mpi_kind_explorer | sed -e 's/,/ /g')
else
  KINDS=$@
fi

rm -rf V_S_real.uses VS_real.uses Interfaces.f90

echo "Used real kinds: $KINDS"


for kind in ${KINDS} ; do

  echo "        USE MPI__r${kind}_V  ; USE MPI__r${kind}_S" >> V_S_real.uses
  echo "        USE MPI__c${kind}_V  ; USE MPI__c${kind}_S" >> V_S_real.uses
  echo "        USE MPI__r${kind}_VS ; USE MPI__r${kind}_SV" >> VS_real.uses
  echo "        USE MPI__c${kind}_VS ; USE MPI__c${kind}_SV" >> VS_real.uses

  for tag in v s sv vs ; do

    sed -e "/_type/s//_r${kind}/" -e "/type/s//real(${kind})/" \
       "${dir}/mpi__type_${tag}.f90" >> Interfaces.f90
    sed -e "/_type/s//_c${kind}/" -e "/type/s//complex(${kind})/" \
       "${dir}/mpi__type_${tag}.f90" >> Interfaces.f90

  done
done


for tag in v s sv vs ; do

  sed -e "/_type/s//_logical/" -e "/type/s//logical/" \
      "${dir}/mpi__type_${tag}.f90" >> Interfaces.f90

  sed -e "/_type/s//_character/" -e "/type/s//character(*)/" \
      "${dir}/mpi__type_${tag}.f90" >> Interfaces.f90

done
