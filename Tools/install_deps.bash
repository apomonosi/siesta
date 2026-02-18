#!/bin/bash

# Installation script for (nearly) all Siesta dependencies.
# Generally the Cmake invocation will install everything, however,
# this script can be used to install stuff locally in custom locations.
#
# This installation script has been written by:
#  Nick R. Papior, 2016-2023.
#
# The author takes no responsibility of damage done to your hardware or
# software. It is up to YOU that the script executes the correct commands.
#
# This script is released under the LGPL-3.0 license.

# If you want to change your compiler version you should define the
# global variables that are used for the configure scripts to grab the
# compiler, they should be CC and FC. Also if you want to compile with
# different flags you should export those variables; CFLAGS, FFLAGS.
# Optionally you could use VENDOR=gnu|intel
# to use defaults for the specific compiler.
#
set -e


# If you have downloaded other versions edit these version strings
fdf_v=0.5.1
fdf_url=https://gitlab.com/siesta-project/libraries/libfdf/uploads/e1774d1bd63db5edfdcadbe3972bbfcc/libfdf-$fdf_v.tar.gz
xc_v=6.2.2
xc_url=https://www.tddft.org/programs/libxc/down.php?file=$xc_v/libxc-$xc_v.tar.gz
gridxc_v=2.0.2
gridxc_url=https://gitlab.com/-/project/11928876/uploads/7896893a5d2197a75d4722c665859abb/libgridxc-2.0.2.tar.gz
xml_v=1.6.3
xml_url=https://gitlab.com/siesta-project/libraries/xmlf90/uploads/8153db06dece1c0b9c38dcda31918fbf/xmlf90-${xml_v}.tar.gz
psml_v=2.0.1
psml_url=https://gitlab.com/siesta-project/libraries/libpsml/uploads/90654f5445a8d073e3c003e1735e3b5d/libpsml-$psml_v.tar.gz

flook_v=0.8.4
flook_url=https://github.com/ElectronicStructureLibrary/flook/archive/v$flook_v.tar.gz

zlib_v=1.3
zlib_url=https://zlib.net/fossils/zlib-${zlib_v}.tar.gz
hdf_v=1.14.3
hdf_url=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf_v%.*}/hdf5-${hdf_v}/src/hdf5-${hdf_v}.tar.bz2
cdfc_v=4.9.2
cdfc_url=https://github.com/Unidata/netcdf-c/archive/refs/tags/v${cdfc_v}.tar.gz
cdff_v=4.6.1
cdff_url=https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${cdff_v}.tar.gz

sdftd3_v=1.0.0
sdftd3_url=https://github.com/dftd3/simple-dftd3/releases/download/v${sdftd3_v}/s-dftd3-${sdftd3_v}-source.tar.xz

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


# Generic options for the build instructions
# Decide whether everything is installed in 1 directory
_single_dir=1
_make_j=4
_test=1
_down_only=0
_down_backend=wget
declare -a _order_packages
declare -a _packages
_order_packages=(
  libfdf
  fdf
  flook
  xmlf90
  xml
  zlib
  hdf5
  hdf
  netcdf
  cdf
  libxc
  xc
  libgridxc
  gridxc
  libpsml
  psml
  simple-dftd3
  simpledftd3
  sdftd3
  dftd3
  sd3
  d3
)
_packages=(
  libfdf
  flook
  xmlf90
  zlib
  hdf5
  netcdf
  libxc
  libgridxc
  gridxc
  libpsml
  simple-dftd3
)

_cmake_args=()
_cmake_prefix=()

function get_name {
  local in=$1
  case $1 in
    libfdf|fdf)
      echo "libfdf" ;;
    flook)
      echo "flook" ;;
    xmlf90|xml)
      echo "xmlf90" ;;
    zlib)
      echo "zlib" ;;
    hdf5|hdf)
      echo "hdf5" ;;
    netcdf|cdf)
      echo "netcdf" ;;
    libxc|xc)
      echo "libxc" ;;
    libgridxc|gridxc)
      echo "libgridxc" ;;
    libpsml|psml)
      echo "libpsml" ;;
    simple-dftd3|dftd3)
      echo "simple-dftd3" ;;
    *)
      echo "Unknown"
      exit 2
  esac
}


function get_directory {
  local path=$1
  local name=$2
  local version=$3
  shift 3
  if [ $_single_dir -eq 1 ]; then
    printf '%s' "$path"
  else
    printf '%s' "$path/$name/$version"
  fi
}


while [ $# -gt 0 ]; do
  opt=$1 ; shift
  case $opt in
    --wget)
      _down_backend=wget
      ;;
    --curl)
      _down_backend=curl
      ;;
    --prefix|-p)
      ID=$1 ; shift
      ;;
    --download-only|-do)
      _down_only=1
      ;;
    --packages|-P)
      # space separated 
      _packages=($1)
      shift
      ;;
    --skip-tests|--no-tests|--skip-test|--no-test)
      _test=0
      ;;
    --single-directory)
      _single_dir=1
      ;;
    --separate-directory|--separate-dir|-sd)
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
      echo "  --prefix|-p <>: specify the installation directory of the library"
      echo "  --packages|-P: specify which packages that should be installed, defaults to $_packages"
      echo "  --download-only|-do: only download packages, do not install anything"
      echo "  --single-directory : all libraries are installed in --prefix/{bin,lib,include} (default: YES)"
      echo "  --separate-directory : all libraries are installed in --prefix/<package>/<version>/{bin,lib,include} (default: NO)"
      echo "  --make-j <>: run make in parallel using <> number of cores (default: $_make_j)"
      echo ""
      echo "NOTE:"
      echo " Packages may adapt URLs, if needed please just download the files required and place them in this folder:"
      echo "  $(dirname $(realpath $0))"
      echo ""
      echo "To customize compilers and flags please export these environment variables:"
      echo "  CC"
      echo "  FC"
      echo "  MPICC"
      echo "  MPIFC"
      echo "  CFLAGS"
      echo "  FFLAGS"
      echo ""
      echo "For instance, to install libgridxc with libxc support (and no other packages), do:"
      echo "  $0 --packages 'libgridxc libxc'"
      echo ""
      echo "This script allows installing these packages [version]:"
      echo ""
      echo " - libfdf [$fdf_v]"
      echo " - libxc [$xc_v]"
      echo " - libgridxc [$gridxc_v]"
      echo " - xmlf90 [$xml_v]"
      echo " - libpsml [$psml_v]"
      echo " - flook [$flook_v]"
      echo " - zlib [$zlib_v]"
      echo " - hdf5 [$hdf_v]"
      echo " - netcdf-c [$cdfc_v]"
      echo " - netcdf-fortran [$cdff_v]"
      echo " - simple-dftd3 [$sdftd3_v]"
      echo ""
      exit 0
      ;;
  esac
done


echo "Installing libraries in folder: $ID"
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
has_wget=$?
which curl> /dev/null
has_curl=$?
if [ $has_wget -ne 0 ]; then
  _down_backend=curl
elif [ $has_curl -ne 0 ]; then
  _down_backend=wget
fi

case "$_down_backend" in
  wget)
    function _dwn_file {
      wget -O $1 $2
    }
    ;;
  curl)
    function _dwn_file {
      curl -o $1 $2
    }
    ;;
esac

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


function install_flook {
  flook_dir=$(get_directory $1 flook $flook_v)
  shift
  [ -d $flook_dir/lib64 ] && flook_lib=lib64 || flook_lib=lib
  if [ ! -d $flook_dir/$flook_lib ]; then
    rm -rf flook-${flook_v}
    tar xfz flook-${flook_v}.tar.gz
    cd flook-${flook_v}
    mkdir -p obj ; cd obj
    {
      echo "# Setup script creation of setup.make"
      [ "x$FC" != "x" ] && \
          echo "FC = $FC"
      [ "x$FCFLAGS" != "x" ] && \
          echo "FFLAGS = $FCFLAGS"
      [ "x$FFLAGS" != "x" ] && \
          echo "FFLAGS = $FFLAGS"
      [ "x$CC" != "x" ] && \
          echo "CC = $CC"
      [ "x$CFLAGS" != "x" ] && \
          echo "CFLAGS = $CFLAGS"
    } > setup.make
    {
      echo TOP_DIR=..
      echo include ../Makefile
    } > Makefile
    if [ "x$VENDOR" == "x" ]; then
      make liball
    else
      make VENDOR=$VENDOR liball
    fi
    retval $? "flook make liball"
    make install PREFIX=$flook_dir
    retval $? "flook make install"
    cd ../../
    rm -rf flook-${flook_v}
    echo "Completed installing flook"
    [ -d $p/lib64 ] && flook_lib=lib64 || flook_lib=lib
  else
    echo "flook directory already found."
  fi
  _cmake_args+=(
    -DWITH_FLOOK=ON
    -DFLOOK_ROOT=$flook_dir
  )
}

function install_xml {
  install_xmlf90 $@
}
function install_xmlf90 {
  xml_dir=$(get_directory $1 xmlf90 $xml_v)
  shift
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
    ./configure --prefix $xml_dir
    retval $? "xmlf90 config"
    make -j $_make_j
    retval $? "xmlf90 make"
    if [ $_test -eq 1 ]; then
      make check 2>&1 | tee xmlf90.check
      retval $? "xmlf90 make check"
    fi
    make install
    retval $? "xmlf90  make install"
    [ -e xmlf90.check ] && mv xmlf90.check $p/
    cd ../
    rm -rf $d
    echo "Completed installing xmlf90"
    [ -d $xml_dir/lib64 ] && xml_lib=lib64 || xml_lib=lib
  else
    echo "xmlf90 directory already found."
  fi
  _cmake_args+=(
    -Dxmlf90_ROOT=$xml_dir
  )
}

function install_fdf {
  install_libfdf $@
}
function install_libfdf {
  fdf_dir=$(get_directory $1 libfdf $fdf_v)
  shift
  [ -d $fdf_dir/lib64 ] && fdf_lib=lib64 || fdf_lib=lib
  if [ ! -e $fdf_dir/$fdf_lib/libfdf.a ]; then
    rm -rf libfdf-${fdf_v}
    tar xfz libfdf-${fdf_v}.tar.gz
    cd libfdf-${fdf_v}
    cmake -B _build -G Ninja -DCMAKE_INSTALL_PREFIX=${fdf_dir}
    retval $? "libfdf cmake setup"
    cmake --build _build/
    retval $? "libfdf cmake build"
    if [ $_test -eq 1 ]; then
      pushd _build && ctest && popd
      retval $? "libfdf tests"
    fi
    cmake --install _build/
    retval $? "libfdf make install"
    cd ../
    rm -rf libfdf-$fdf_v
    echo "Completed installing limfdf"
    [ -d $fdf_dir/lib64 ] && fdf_lib=lib64 || fdf_lib=lib
  else
    echo "libfdf directory already found."
  fi
  _cmake_args+=(
    -Dlibfdf_ROOT=$fdf_dir
  )
  _cmake_prefix+=($fdf_dir)
}

function install_dftd3 {
  install_simple_dftd3 $@
}
function install_sdftd3 {
  install_simple_dftd3 $@
}
function install_d3 {
  install_simple_dftd3 $@
}
function install_sd3 {
  install_simple_dftd3 $@
}
function install_simpledftd3 {
  install_simple_dftd3 $@
}
function install_simple_dftd3 {
  sdftd3_dir=$(get_directory $1 s-dftd3 $sdftd3_v)
  shift
  [ -d $sdftd3_dir/lib64 ] && sdftd3_lib=lib64 || sdftd3_lib=lib
  if [ ! -e $sdftd3_dir/$sdftd3_lib/libs-dftd3.a ]; then
    rm -rf s-dftd3-${sdftd3_v}
    tar xfJ s-dftd3-${sdftd3_v}.tar.xz
    cd s-dftd3-${sdftd3_v}
    cmake -B _build -G Ninja -DWITH_OpenMP=0 -DCMAKE_INSTALL_PREFIX=${sdftd3_dir}
    retval $? "simple-dftd3 cmake setup"
    cmake --build _build/
    retval $? "simple-dftd3 cmake build"
    if [ $_test -eq 1 ]; then
      pushd _build && ctest && popd
      retval $? "simple-dftd3 tests"
    fi
    cmake --install _build/
    retval $? "simple-dftd3 make install"
    cd ../
    rm -rf s-dftd3-$sdftd3_v
    echo "Completed installing simple-dftd3"
    [ -d $sdftd3_dir/lib64 ] && sdftd3_lib=lib64 || sdftd3_lib=lib
  else
    echo "simple-dftd3 directory already found."
  fi
  _cmake_args+=(
    -DWITH_DFTD3=ON
    -Ds-dftd3_ROOT=$sdftd3_dir
  )
  _cmake_prefix+=($sdftd3_dir)
}


function install_psml {
  install_libpsml $@
}
function install_libpsml {
  install_xmlf90 $1
  psml_dir=$(get_directory $1 libpsml $psml_v)
  shift
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
    make -j $_make_j
    retval $? "PSML make"
    if [ $_test -eq 1 ]; then
      make check 2>&1 | tee psml.test
      retval $? "PMSL make check"
    fi
    make install
    retval $? "PSML make install"
    [ -e psml.test ] && mv psml.test $psml_dir/
    cd ../../
    rm -rf psml-${psml_v}
    echo "Completed installing PSML"
    [ -d $psml_dir/lib64 ] && psml_lib=lib64 || psml_lib=lib
  else
    echo "PSML directory already found."
  fi
  _cmake_args+=(
    -Dlibpsml_ROOT=$psml_dir
  )
}

function install_xc {
  install_libxc $@
}
function install_libxc {
  xc_dir=$(get_directory $1 libxc $xc_v)
  shift
  [ -d $xc_dir/lib64 ] && xc_lib=lib64 || xc_lib=lib
  if [ ! -e $xc_dir/$xc_lib/libxcf90.a ]; then
    rm -rf libxc-${xc_v}
    tar xfz libxc-${xc_v}.tar.gz
    cd libxc-${xc_v}
    ./configure --enable-shared --prefix $xc_dir
    retval $? "libxc config"
    make -j $_make_j
    retval $? "libxc make"
    if [ $_test -eq 1 ]; then
      make test 2>&1 | tee libxc.test
      retval $? "libxc make test"
    fi
    make install
    retval $? "libxc  make install"
    [ -e libxc.test ] && mv libxc.test $xc_dir/
    cd ../
    rm -rf libxc-${xc_v}
    echo "Completed installing libxc"
    [ -d $xc_dir/lib64 ] && xc_lib=lib64 || xc_lib=lib
  else
    echo "libxc directory already found."
  fi
  _cmake_prefix+=($xc_dir)
}

function install_gridxc {
  install_libgridxc $@
}
function install_libgridxc {
  gridxc_dir=$(get_directory $1 libgridxc $gridxc_v)
  shift
  [ -d $gridxc_dir/lib64 ] && gridxc_lib=lib64 || gridxc_lib=lib
  if [ ! -e $gridxc_dir/$gridxc_lib/libgridxc_dp.a ]; then
    rm -rf libgridxc-$gridxc_v
    tar xfz libgridxc-$gridxc_v.tar.gz
    cd libgridxc-$gridxc_v

    function _build {
      rm -rf build
      mkdir build
      cd build
      if [ -d $xc_dir ]; then
        ../configure --enable-shared --enable-multiconfig --with-libxc=$xc_dir --prefix=$gridxc_dir $@
      else
        ../configure --enable-shared --enable-multiconfig --prefix=$gridxc_dir $@
      fi
      retval $? "gridxc configure"
      make -j $_make_j
      retval $? "gridxc make"
      if [ $_test -eq 1 ]; then
        make check > gridxc.check
        retval $? "gridxc make check"
      fi
      make install
      retval $? "gridxc make install"
    }
    [ "x$MPICC" == "x" ] && export MPICC=mpicc
    [ "x$MPIFC" == "x" ] && export MPIFC=mpifort
    [ "x$MPI_ROOT" == "x" ] && export MPI_ROOT=$(dirname $(dirname $(which $MPICC)))

    # Start build process
    _build --with-mpi=$MPI_ROOT CC=$MPICC FC=$MPIFC
    [ -e gridxc.check ] && mv gridxc.check $gridxc_dir/gridxc_dp_mpi.check
    cd ..
    _build --with-mpi=$MPI_ROOT CC=$MPICC FC=$MPIFC --enable-single-precision
    [ -e gridxc.check ] && mv gridxc.check $gridxc_dir/gridxc_sp_mpi.check
    cd ..
    _build --without-mpi
    [ -e gridxc.check ] && mv gridxc.check $gridxc_dir/gridxc_dp.check
    cd ..
    _build --without-mpi --enable-single-precision
    [ -e gridxc.check ] && mv gridxc.check $gridxc_dir/gridxc_sp.check
    cd ..
    
    cd ..
    rm -rf libgridxc-${gridxc_v}
    echo "Completed installing GridXC"
    [ -d $gridxc_dir/lib64 ] && gridxc_lib=lib64 || gridxc_lib=lib
  else
    echo "GridXC directory already found."
  fi
  _cmake_args+=(
    -Dlibgridxc_ROOT=$gridxc_dir
  )
  _cmake_prefix+=($gridxc_dir)
}


function install_zlib {
  zlib_dir=$(get_directory $1 zlib $zlib_v)
  [ -d $zlib_dir/lib64 ] && zlib_lib=lib64 || zlib_lib=lib
  if [ ! -e $zlib_dir/$zlib_lib/libz.a ]; then
    rm -rf zlib-${zlib_v}
    tar xfz zlib-${zlib_v}.tar.gz
    cd zlib-${zlib_v}
    ./configure --prefix $zlib_dir/
    retval $? "zlib config"
    make -j $_make_j
    retval $? "zlib make"
    if [ $_test -eq 1 ]; then
      make test 2>&1 | tee zlib.test
      retval $? "zlib make test"
    fi
    make install
    retval $? "zlib make install"
    [ -e zlib.test ] && mv zlib.test $zlib_dir/
    cd ../
    rm -rf zlib-${zlib_v}
    echo "Completed installing zlib"
    [ -d $zlib_dir/lib64 ] && zlib_lib=lib64 || zlib_lib=lib
  else
    echo "zlib directory already found."
  fi
  _cmake_prefix+=($zlib_dir)
}


function install_hdf {
  install_hdf5 $@
}
function install_hdf5 {
  install_zlib $1
  hdf_dir=$(get_directory $1 hdf5 $hdf_v)
  shift
  [ -d $hdf_dir/lib64 ] && hdf_lib=lib64 || hdf_lib=lib
  if [ ! -e $hdf_dir/$hdf_lib/libhdf5.a ]; then
    rm -rf hdf5-$hdf_v
    tar xfj hdf5-$hdf_v.tar.bz2
    cd hdf5-$hdf_v
    mkdir build ; cd build
    ../configure --prefix=$hdf_dir \
      --enable-shared --enable-static \
      --enable-fortran --with-zlib=$zlib_dir \
      LDFLAGS="-L$zlib_dir/$zlib_lib -Wl,-rpath,$zlib_dir/$zlib_lib"
    retval $? "hdf5 configure"
    make -j $_make_j
    retval $? "hdf5 make"
    if [ $_test -eq 1 ]; then
      make check-s 2>&1 | tee hdf5.test
      retval $? "hdf5 make check-s"
    fi
    make install
    retval $? "hdf5 make install"
    [ -e hdf5.test ] && mv hdf5.test $hdf_dir/
    cd ../../
    rm -rf hdf5-${hdf_v}
    echo "Completed installing hdf5"
    [ -d $hdf_dir/lib64 ] && hdf_lib=lib64 || hdf_lib=lib
  else
    echo "hdf5 directory already found."
  fi
  _cmake_prefix+=($hdf_dir)
}


function install_netcdf_c {
  install_hdf5 $1
  cdf_dir=$(get_directory $1 netcdf $cdfc_v)
  shift
  [ -d $cdf_dir/lib64 ] && cdf_lib=lib64 || cdf_lib=lib
  if [ ! -e $cdf_dir/$cdf_lib/libnetcdf.a ]; then
    rm -rf netcdf-c-$cdfc_v
    tar xfz netcdf-c-$cdfc_v.tar.gz
    cd netcdf-c-$cdfc_v
    mkdir build ; cd build
    ../configure --prefix=$cdf_dir \
      --enable-shared --enable-static \
      --enable-netcdf-4 --disable-dap \
      CPPFLAGS="-I$hdf_dir/include -I$zlib_dir/include" \
      LDFLAGS="-L$hdf_dir/$hdf_lib -Wl,-rpath,$hdf_dir/$hdf_lib \
      -L$zlib_dir/$zlib_lib -Wl,-rpath,$zlib_dir/$zlib_lib"
    retval $? "netcdf configure"
    make -j $_make_j
    retval $? "netcdf make"
    if [ $_test -eq 1 ]; then
      make check 2>&1 | tee netcdf-c.test
      retval $? "netcdf make check"
    fi
    make install
    retval $? "netcdf make install"
    [ -e netcdf-c.test ] && mv netcdf-c.test $cdf_dir/
    cd ../../
    rm -rf netcdf-c-$cdfc_v
    echo "Completed installing C NetCDF library"
    [ -d $cdf_dir/lib64 ] && cdfc_lib=lib64 || cdfc_lib=lib
  else
    echo "netcdf directory already found."
  fi
  _cmake_args+=(
    -DWITH_NETCDF=ON
    -DNetCDF_ROOT=$cdf_dir
  )
  _cmake_prefix+=($cdf_dir)
}

function install_netcdf_f {
  # this will install directly into the C-library location
  # This makes sense since they are tightly coupled.
  if [ ! -e $cdf_dir/$cdf_lib/libnetcdff.a ]; then
    rm -rf netcdf-fortran-$cdff_v
    tar xfz netcdf-fortran-$cdff_v.tar.gz
    cd netcdf-fortran-$cdff_v
    mkdir build ; cd build
    ../configure CPPFLAGS="-DgFortran -I$zlib_dir/include \
      -I$hdf_dir/include -I$cdf_dir/include" \
      LIBS="-L$zlib_dir/$zlib_lib -Wl,-rpath,$zlib_dir/$zlib_lib \
      -L$hdf_dir/$hdf_lib -Wl,-rpath,$hdf_dir/$hdf_lib \
      -L$cdf_dir/$cdf_lib -Wl,-rpath,$cdf_dir/$cdf_lib \
      -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz" \
      --prefix=$cdf_dir --enable-static --enable-shared
    retval $? "netcdf-fortran configure"
    make -j $_make_j
    retval $? "netcdf-fortran make"
    if [ $_test -eq 1 ]; then
      make check 2>&1 | tee netcdf-fortran.serial
      retval $? "netcdf-fortran make check"
    fi
    make install
    retval $? "netcdf-fortran make install"
    [ -e netcdf-fortran.test ] && mv netcdf-fortran.test $cdf_dir/
    cd ../../
    rm -rf netcdf-fortran-$cdff_v
    echo "Completed installing Fortran NetCDF library"
  else
    echo "netcdf-fortran library already found."
  fi
}


# First the script will download the files so that we can error
# out fast
function has_pack {
  local s=$1
  local pack
  for pack in "${_packages[@]}"
  do
    if [ "$pack" == "$s" ]; then
      return 0
    fi
  done
  return 1
}

for pack in "${_order_packages[@]}"
do
  # Check if it is in the packages
  if $(has_pack $pack) ; then
    case $pack in
      libfdf|fdf)
        download_file libfdf-$fdf_v.tar.gz $fdf_url ;;
      flook)
        download_file flook-$flook_v.tar.gz $flook_url ;;
      libxc|xc)
        download_file libxc-$xc_v.tar.gz $xc_url ;;
      simple-dftd3|sdftd3|sd3)
        download_file s-dftd3-${sdftd3_v}.tar.xz $sdftd3_url ;;
      libgridxc|gridxc)
        download_file libgridxc-$gridxc_v.tar.gz $gridxc_url ;;
      xmlf90|xml)
        download_file xmlf90-$xml_v.tar.gz $xml_url ;;
      libpsml|psml)
        download_file libpsml-$psml_v.tar.gz $psml_url ;;
      zlib)
        download_file zlib-$zlib_v.tar.gz $zlib_url ;;
      hdf5|hdf)
        download_file hdf5-$hdf_v.tar.bz2 $hdf_url ;;
      netcdf|cdf)
        download_file netcdf-c-$cdfc_v.tar.gz $cdfc_url
        download_file netcdf-fortran-$cdff_v.tar.gz $cdff_url ;;
      *)
        echo "Error in program... for $pack"
        ;;
    esac
  fi
done

if [ $_down_only -eq 1 ]; then
  echo "Asked to only download, quitting..."
  exit 0
fi

for pack in "${_order_packages[@]}"
do
  # Check if it is in the packages
  if $(has_pack $pack) ; then
    case $pack in
      netcdf|cdf)
        install_netcdf_c $ID
        install_netcdf_f
        ;;
      *)
        install_${pack//-/_} $ID
        ;;
    esac
  fi
done


echo ""
echo "##############################"
echo "# Completed installation     #"
echo "# of the following packages: #"
for pack in "${_packages[@]}"
do
  printf '%-29s%s\n' "# - ${pack}" "#"
done
echo "##############################"
echo ""
  
if [ ${#_cmake_args[@]} -gt 1 ]; then
  _cmake_args=($(echo "${_cmake_args[@]}" | tr ' ' '\n' | sort -u))
fi

# Remove duplicates in the cmake_prefix path
if [ ${#_cmake_prefix[@]} -eq 0 ]; then
  _cmake_prefix=\$CMAKE_PREFIX_PATH
else
  _cmake_prefix=$(echo "${_cmake_prefix[@]}" | tr ' ' '\n' | sort -u | tr '\n' ';')
  _cmake_prefix="$_cmake_prefix\$CMAKE_PREFIX_PATH"
fi

echo "Please add this to the cmake invocation (please append additional paths if needed):"
echo ""
echo "  cmake ${_cmake_args[@]} -DCMAKE_PREFIX_PATH=\"$_cmake_prefix\""
echo ""
