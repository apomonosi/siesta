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

for ELEC in ts_au_100_elec_1x1 ts_au_100_elec_3x3 ; do

    # Start with the electrode calculation
    echo "==> Electrode Calculation for $ELEC"
    mkdir Elec_$ELEC
    cd Elec_$ELEC
    ln ../Au.psf .
    ln ../$ELEC.fdf .
    $TS --electrode $ELEC.fdf > $ELEC.out
    RETVAL=$?
    if [ $RETVAL -ne 0 ]; then
	echo "The electrode calculation did not go well ..."
	exit
    fi
    cp $ELEC.out ../
    # Go back to base directory
    cd ..

    # Scattering region calculation
    if [ $ELEC == "ts_au_100_elec_1x1" ]; then
	scats="ts_au_100_3x3 ts_au_100_3x3_0.25V"
    else
	scats="ts_au_100_1x1 ts_au_100_1x1_0.25V"
    fi
    previous=
    for SCAT in $scats
    do
	echo "==> Scattering Region Calculation for $SCAT"
	mkdir Scat_$SCAT
	cd Scat_$SCAT
	ln ../Au.psf .
	ln ../$SCAT.fdf .
	# Copy the electrode's .TSHS
	ln ../Elec_$ELEC/$ELEC.TSHS .
	if [ -e ../Scat_$previous/$previous.TSDE ]; then
	    cp ../Scat_$previous/$previous.TSDE $SCAT.TSDE
	fi
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
	ln ../Elec_$ELEC/$ELEC.TSHS .
	ln ../Scat_$SCAT/$SCAT.TSHS .
	ln ../$SCAT.fdf .
	$TBT $SCAT.fdf > ${SCAT}_tbt.out
	RETVAL=$?
	if [ $RETVAL -ne 0 ]; then
	    echo "The scattering region calculation did not go well ..."
	    exit
	fi

	cp ${SCAT}_tbt.out ../
	for f in $SCAT.TBT.TRANS_* $SCAT.TBT.TEIG_* $SCAT.TBT.ADOS_*
	do
	    if [ -e $f ]; then
		cp $f ../
	    fi
	done
	rm -f *.TBTGF*

	cd ..
    done
done

