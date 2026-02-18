#!/bin/bash

# Installation script for libxc and GridXC
# This installation script has been written by:
#  Nick R. Papior, 2016-2020.
#
# The author takes no responsibility of damage done to your hardware or
# software. It is up to YOU that the script executes the correct commands.
#
# This script is released under the LGPL license.

# VERY BASIC installation script of required libraries
# for installing these packages:
#   libxc
#   libgridxc
# If you want to change your compiler version you should define the
# global variables that are used for the configure scripts to grab the
# compiler, they should be CC and FC. Also if you want to compile with
# different flags you should export those variables; CFLAGS, FFLAGS.

# If you have downloaded other versions edit these version strings
xc_v=4.3.4
gridxc_v=0.9.6

# Install path, change accordingly
# You can change this variable to control the installation path
# If you want the installation path to be a "packages" folder in
# your home directory, change to this:
# ID=$HOME/packages
if [ -z $PREFIX ]; then
    ID=$(pwd)/build
else
    ID=$PREFIX
fi
# Decide whether everything is installed in 1 directory
_single_dir=1
_make_j=2

while [ $# -gt 0 ]; do
    opt=$1 ; shift
    case $opt in
	--prefix|-p)
	    ID=$1 ; shift
	    ;;
	--libxc-version|-xc-v)
	    xc_v=$1 ; shift
	    ;;
	--libgridxc-version|-gridxc-v)
	    gridxc_v=$1 ; shift
	    ;;
	--single-directory)
	    _single_dir=1
	    ;;
	--separate-directory)
	    _single_dir=0
	    ;;
	--make-j)
	    _make_j=$1 ; shift
	    ;;
	--help|-h)
	    echo " $0 --help shows this message"
	    echo ""
	    echo "These options are available:"
	    echo ""
	    echo "  --prefix|-p <>: specify the installation directory of the libraries"
	    echo "  --libxc-version|-xc-v <>: specify the libxc version (default: $xc_v)"
	    echo "  --libgridxc-version|-gridxc-v <>: specify the libGridXC version (default: $gridxc_v)"
	    echo "  --single-directory : all libraries are installed in --prefix/{bin,lib,include} (default: YES)"
	    echo "  --separate-directory : all libraries are installed in --prefix/<package>/<version>/{bin,lib,include} (default: NO)"
	    echo "  --make-j <>: run make in parallel using <> number of cores (default: $_make_j)"
	    echo ""
	    echo "To customize compilers and flags please export these environment variables:"
	    echo "  CC"
	    echo "  FC"
	    echo "  MPICC"
	    echo "  MPIFC"
	    echo "  CFLAGS"
	    echo "  FFLAGS"
	    echo ""
	    exit 0
	    ;;
    esac
done


echo "Installing libraries in folder: $ID"
echo ""
echo "Using these software versions:"
echo "  libxc     : $xc_v"
echo "  libGridXC : $gridxc_v"
echo ""
echo "If you want other versions, please try $0 --help and check the options"
echo "  Will wait 2 sec before proceeding"
echo ""
sleep 2

mkdir -p $ID

# First we check that the user have downloaded the files
function file_exists {
    if [ ! -e $(pwd)/$1 ]; then
	echo "I could not find file $1..."
	echo "Please download the file and place it in this folder:"
	echo " $(pwd)"
	exit 1
    fi
}

# Download a file, if able and the file does not exist
which wget > /dev/null
if [ $? -eq 0 ]; then
    # success we can download using wget
    function _dwn_file {
	wget -O $1 $2
    }
else
    function _dwn_file {
	curl -o $1 $2
    }
fi

# Use download function
#  $1 is name of file
#  $2 is URL
function download_file {
    if [ ! -e $(pwd)/$1 ] ; then
	# Try and download
	_dwn_file $1 $2
    fi
}

# Check for function $?
function retval {
    local ret=$1
    local info="$2"
    shift 2
    if [ $ret -ne 0 ]; then
	echo "Error: $ret"
	echo "$info"
	exit 1
    fi
}
 
# Download files if they can
download_file libxc-${xc_v}.tar.gz http://www.tddft.org/programs/libxc/down.php?file=$xc_v/libxc-$xc_v.tar.gz
download_file libgridxc-${gridxc_v}.tar.gz https://gitlab.com/siesta-project/libraries/libgridxc/uploads/e6e4ec1e3ec648a45b3613e724c7be04/libgridxc-$gridxc_v.tar.gz

file_exists libxc-${xc_v}.tar.gz
file_exists libgridxc-${gridxc_v}.tar.gz
unset file_exists

#################
# Install libxc #
#################
if [ $_single_dir -eq 0 ]; then
    xc_dir=$ID/libxc/$xc_v
else
    xc_dir=$ID
fi
[ -d $xc_dir/lib64 ] && xc_lib=lib64 || xc_lib=lib
if [ ! -e $xc_dir/$xc_lib/libxcf90.a ]; then
    rm -rf libxc-${xc_v}
    tar xfz libxc-${xc_v}.tar.gz
    cd libxc-${xc_v}
    ./configure --enable-shared --prefix $xc_dir/
    retval $? "libxc config"
    make -j 2
    retval $? "libxc make"
    make test 2>&1 | tee libxc.test
    retval $? "libxc make test"
    make install
    retval $? "libxc  make install"
    mv libxc.test $xc_dir/
    cd ../
    rm -rf libxc-${xc_v}
    echo "Completed installing libxc"
    [ -d $xc_dir/lib64 ] && xc_lib=lib64 || xc_lib=lib
else
    echo "libxc directory already found."
fi

##################
# Install GridXC #
##################
if [ $_single_dir -eq 0 ]; then
    gridxc_dir=$ID/gridxc/$gridxc_v
else
    gridxc_dir=$ID
fi
[ -d $gridxc_dir/lib64 ] && gridxc_lib=lib64 || gridxc_lib=lib
if [ ! -e $gridxc_dir/$gridxc_lib/libgridxc_dp_mpi.a ]; then
    rm -rf libgridxc-$gridxc_v
    tar xfz libgridxc-$gridxc_v.tar.gz
    cd libgridxc-$gridxc_v

    function _build {
	rm -rf build
	mkdir build
	cd build
	../configure --enable-shared --enable-multiconfig --with-libxc=$xc_dir --prefix=$gridxc_dir $@
        retval $? "gridxc configure"
	make
        retval $? "gridxc make"
	make check > gridxc.check
        retval $? "gridxc make check"
	make install
	retval $? "gridxc make install"
    }
    [ "x$MPICC" == "x" ] && export MPICC=mpicc
    [ "x$MPIFC" == "x" ] && export MPIFC=mpifort
    [ "x$MPI_ROOT" == "x" ] && export MPI_ROOT=$(dirname $(dirname $(which $MPICC)))
    

    # Start build process
    _build --with-mpi=$MPI_ROOT CC=$MPICC FC=$MPIFC
    mv gridxc.check $gridxc_dir/gridxc_dp_mpi.check
    cd ..
    _build --with-mpi=$MPI_ROOT CC=$MPICC FC=$MPIFC --enable-single-precision
    mv gridxc.check $gridxc_dir/gridxc_sp_mpi.check
    cd ..
    _build --without-mpi
    mv gridxc.check $gridxc_dir/gridxc_dp.check
    cd ..
    _build --without-mpi --enable-single-precision
    mv gridxc.check $gridxc_dir/gridxc_sp.check
    cd ..
    
    cd ..
    rm -rf libgridxc-${gridxc_v}
    echo "Completed installing GridXC"
    [ -d $gridxc_dir/lib64 ] && gridxc_lib=lib64 || gridxc_lib=lib
else
    echo "GridXC directory already found."
fi

##########################
# Completed installation #
##########################

echo ""
echo "##########################"
echo "# Completed installation #"
echo "#   of GridXC package    #"
echo "#  and its dependencies  #"
echo "##########################"
echo ""
echo ""

echo "Please add the following to the BOTTOM of your arch.make file (for double-precision and MPI):"
echo ""
echo "GRIDXC_ROOT = $gridxc_dir"
echo "include \$(GRIDXC_ROOT)/share/org.siesta-project/gridxc_dp_mpi.mk"
if [ "$xc_dir/$xc_lib" == "$gridxc_dir/$gridxc_lib" ]; then
    echo "LDFLAGS += -Wl,-rpath,$xc_dir/$xc_lib"
else
    echo "LDFLAGS += -Wl,-rpath,$xc_dir/$xc_lib"
    echo "LDFLAGS += -Wl,-rpath,$gridxc_dir/$gridxc_lib"
fi
echo ""
