#!/bin/bash

# Installation script for xmlf90 and psml
# This installation script has been written by:
#  Nick R. Papior, 2016-2020.
#
# The author takes no responsibility of damage done to your hardware or
# software. It is up to YOU that the script executes the correct commands.
#
# This script is released under the LGPL license.

# VERY BASIC installation script of required libraries
# for installing these packages:
#   xmlf90
#   libpsml
# If you want to change your compiler version you should define the
# global variables that are used for the configure scripts to grab the
# compiler, they should be CC and FC. Also if you want to compile with
# different flags you should export those variables; CFLAGS, FFLAGS.

# If you have downloaded other versions edit these version strings
xml_v=1.5.4
psml_v=1.1.10

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
	--xmlf90-version|-xml-v)
	    xml_v=$1 ; shift
	    ;;
	--psml-version|-psml-v)
	    psml_v=$1 ; shift
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
	    echo "  --xmlf90-version|-xml-v <>: specify the xmlf90 version (default: $xml_v)"
	    echo "  --psml-version|-psml-v <>: specify the PSML version (default: $psml_v)"
	    echo "  --single-directory : all libraries are installed in --prefix/{bin,lib,include} (default: YES)"
	    echo "  --separate-directory : all libraries are installed in --prefix/<package>/<version>/{bin,lib,include} (default: NO)"
	    echo "  --make-j <>: run make in parallel using <> number of cores (default: $_make_j)"
	    echo ""
	    echo "To customize compilers and flags please export these environment variables:"
	    echo "  CC"
	    echo "  FC"
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
echo "  xmlf90  : $xml_v"
echo "  PSML    : $psml_v"
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
download_file xmlf90-${xml_v}.tar.gz https://launchpad.net/xmlf90/trunk/1.5/+download/xmlf90-$xml_v.tar.gz
download_file libpsml-${psml_v}.tar.gz https://gitlab.com/siesta-project/libraries/libpsml/uploads/95c4d83d750d1ca87659b6395c8f7ea7/libpsml-$psml_v.tar.gz

file_exists xmlf90-${xml_v}.tar.gz
file_exists libpsml-${psml_v}.tar.gz
unset file_exists

##################
# Install xmlf90 #
##################
if [ $_single_dir -eq 0 ]; then
    xml_dir=$ID/xmlf90/$xml_v
else
    xml_dir=$ID
fi
[ -d $xml_dir/lib64 ] && xml_lib=lib64 || xml_lib=lib
if [ ! -e $xml_dir/$xml_lib/libxmlf90.a ]; then
    rm -rf xmlf90-xmlf90-${xml_v} xmlf90-${xml_v}
    tar xfz xmlf90-${xml_v}.tar.gz
    if [ -d xmlf90-xmlf90-${xml_v} ]; then
	d=xmlf90-xmlf90-${xml_v}
    else
	d=xmlf90-${xml_v}
    fi
    cd $d
    ./configure --prefix $xml_dir/
    retval $? "xmlf90 config"
    make -j 2
    retval $? "xmlf90 make"
    make check 2>&1 | tee xmlf90.check
    retval $? "xmlf90 make check"
    make install
    retval $? "xmlf90  make install"
    mv xmlf90.check $xml_dir/
    cd ../
    rm -rf $d
    echo "Completed installing xmlf90"
    [ -d $xml_dir/lib64 ] && xml_lib=lib64 || xml_lib=lib
else
    echo "xmlf90 directory already found."
fi

################
# Install PSML #
################
if [ $_single_dir -eq 0 ]; then
    psml_dir=$ID/psml/$psml_v
else
    psml_dir=$ID
fi
[ -d $psml_dir/lib64 ] && psml_lib=lib64 || psml_lib=lib
if [ ! -e $psml_dir/$psml_lib/libpsml.a ]; then
    rm -rf libpsml-$psml_v
    tar xfz libpsml-$psml_v.tar.gz
    cd libpsml-$psml_v
    mkdir build ; cd build
    ../configure --prefix=$psml_dir \
		 --with-xmlf90=$xml_dir \
		 LDFLAGS="-L$xml_dir/$xml_lib -Wl,-rpath,$xml_dir/$xml_lib"
    retval $? "PSML configure"
    make -j 2
    retval $? "PSML make"
    make check 2>&1 | tee psml.test
    retval $? "PMSL make check"
    make install
    retval $? "PSML make install"
    mv psml.test $psml_dir/
    cd ../../
    rm -rf psml-${psml_v}
    echo "Completed installing PSML"
    [ -d $psml_dir/lib64 ] && psml_lib=lib64 || psml_lib=lib
else
    echo "PSML directory already found."
fi

##########################
# Completed installation #
##########################

echo ""
echo "##########################"
echo "# Completed installation #"
echo "#    of PSML package     #"
echo "#  and its dependencies  #"
echo "##########################"
echo ""
echo ""

echo "Please add the following to the BOTTOM of your arch.make file:"
echo ""
echo "XMLF90_ROOT = $xml_dir"
echo "include \$(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk"
echo "PSML_ROOT = $psml_dir"
echo "include \$(PSML_ROOT)/share/org.siesta-project/psml.mk"
if [ "$xml_dir/$xml_lib" == "$psml_dir/$psml_lib" ]; then
    echo "LDFLAGS += -Wl,-rpath,$xml_dir/$xml_lib"
else
    echo "LDFLAGS += -Wl,-rpath,$xml_dir/$xml_lib"
    echo "LDFLAGS += -Wl,-rpath,$psml_dir/$psml_lib"
fi
echo ""
