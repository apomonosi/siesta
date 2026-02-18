#!/bin/sh
#
# The script is passed the (probably relative) path to the transiesta
# executable
TS_RAW="$1"

# Extract last component, in case of mpirun-style string
TS_REL_PATH=$(echo ${TS_RAW} | awk '{print $NF}')
TS_EXEC_PREFIX=$(echo ${TS_RAW} | awk '{$NF=""; print}')
TS_NAME=$(basename ${TS_REL_PATH})
EXEC_DIR=$(dirname ${TS_REL_PATH})

# Find absolute path
pushd ${EXEC_DIR}; ABS_EXEC_DIR=$(pwd); popd

TS_ABS=${ABS_EXEC_DIR}/${TS_NAME}
TS="$TS_EXEC_PREFIX $TS_ABS"
OBJDIR=$(basename ${ABS_EXEC_DIR})
ROOT_DIR=$(dirname ${ABS_EXEC_DIR})

echo "Running script with TranSIESTA=$TS, compiled in OBJDIR=${OBJDIR}"
if [ -z "$TBT" ] ; then
    TBT="${TS_EXEC_PREFIX} ${ABS_EXEC_DIR}/tbtrans"
else
    echo "Running script with tbtrans=$TBT"
fi

for ELEC in ts_n_terminal_elec_x ts_n_terminal_elec_z
do

    # Start with the electrode calculation
    echo "==> Electrode Calculation for $ELEC"
    mkdir Elec_$ELEC
    cd Elec_$ELEC
    ln ../C.psf .
    ln ../$ELEC.fdf .
    $TS --electrode $ELEC.fdf > $ELEC.out
    RETVAL=$?
    if [ $RETVAL -ne 0 ]; then
	echo "The electrode calculation did not go well ..."
	exit
    fi
    cp $ELEC.out ..
    # Go back to base directory
    cd ..

done

for SCAT in ts_n_terminal_3 ts_n_terminal_4
do
    echo "==> Scattering Region Calculation for $SCAT"
    mkdir Scat_$SCAT
    cd Scat_$SCAT
    ln ../C.psf .
    ln ../$SCAT.fdf .
    # Copy the electrode's .TSHS
    for ELEC in ts_n_terminal_elec_x ts_n_terminal_elec_z
    do
	ln ../Elec_$ELEC/$ELEC.TSHS .
    done	
    $TS $SCAT.fdf > $SCAT.out
    RETVAL=$?
    if [ $RETVAL -ne 0 ]; then
	echo "** The scattering region calculation for $SCAT did not go well ..."
	exit
    fi
    cp $SCAT.out ../
    rm *.TSGF*
    previous=$SCAT

    # Go back to base directory
    cd ..


    # TBTrans calculation
    echo "==> TBTrans Calculation for $SCAT"
    echo "==> Running $SCAT with tbtrans=$TBT"
    mkdir TBT_$SCAT
    cd TBT_$SCAT
    # Copy input files
    for ELEC in ts_n_terminal_elec_x ts_n_terminal_elec_z
    do
	ln ../Elec_$ELEC/$ELEC.TSHS .
    done	
    ln ../Scat_$SCAT/$SCAT.TSHS .
    ln ../$SCAT.fdf .
    $TBT $SCAT.fdf > ${SCAT}_tbt.out
    RETVAL=$?
    if [ $RETVAL -ne 0 ]; then
	echo "The scattering region calculation did not go well ..."
	exit
    fi

    cp ${SCAT}_tbt.out ../
    for f in $SCAT.TBT.AVTRANS_* $SCAT.TBT.AVTEIG_* $SCAT.TBT.AVADOS_*
    do
	if [ -e $f ]; then
	    cp $f ../
	fi
    done
    rm -f *.TBTGF*

    cd ..
done

