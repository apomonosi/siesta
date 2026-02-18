# Building Siesta with CMake

Siesta requires CMake >= 3.20, and (if used) the ninja (>=1.10) backend.
Both cmake and ninja can be installed easily in most systems, in
particular with conda or pip.

This document provides the basic instructions for building, but
the details depend on the way the software stack is organized in
the target machine. We provide more information about what to do
in various cases in the [build-tools][_build-tools] repository.

Siesta's CMake scripts have borrowed ideas from the DFTB+, SIRIUS, and
other open-source projects, and we fully acknowledge their
contributions and efforts towards streamlining the build-system experience.


## Quick and go

The most basic compilation of Siesta can be done simply by:

```shell
   cmake -S. -B_build -DCMAKE_INSTALL_PREFIX=/path/to/installation  #<MORE_OPTIONS>
   cmake --build _build -j 4
   cmake --install _build
```

If your system is relatively standard and all required dependencies
are found this will result in the `/path/to/installation/bin`
directory being populated with executables for `siesta`, `tbtrans`,
and all the Siesta utilties.

Most likely, you will have to work a bit more to provide the
`<MORE_OPTIONS>` part of the above command.  Follow the instructions
below, and note the sections below on toolchains for some examples.


## HPC sites/admins

Here are some clear instructions that should benefit the performance
of Siesta when running on HPC infrastructure.

- Build your executable for the hardware it is meant to run on.
  If you are cross-compiling, ensure it's being built with the CPU extensions
  that is available on the running machine.
  When cross-compiling, it will automatically remove any host-specific
  default flags in Siesta.

- Use interprocedural optimizations, this will gain ~5% of performance.
  (`-DSIESTA_WITH_IPO=ON`)

- Use any optimization flags you know works good on your architecture.
  The Siesta built-system has flag corrections for certain source files
  which are known in having difficulties with optimizations. Hence,
  flags are generally *safe to use*.

- Avoid using `-ffast-math` (or the compiler equivalents) as it tends to
  mess up the SCF cycles (at least in some cases...).

- One can install both `siesta` with and without OpenMP in the same binary
  location to enable users the choice of using OpenMP or not.

  This is mainly beneficial for TranSiesta+TBtrans runs.

  It is recommended to add this to the build flags to suffix
  the executables with `_omp`:

  ```shell
  cmake ... -DSIESTA_WITH_OPENMP=ON \
            -DSIESTA_EXECUTABLE_SUFFIX=_omp \
            -DSIESTA_SHARED_LIBRARY_SUFFIX=_omp \
            -DSIESTA_STATIC_LIBRARY_SUFFIX=_omp
  ```
  which will create the executable `siesta_omp`.

  The OpenMP enabled executable should, for the moment, only
  be used when it's known to be beneficial (tests are required!).

- Lastly, Siesta is mainly relying on external libraries, BLAS, LAPACK,
  ScaLAPACK/ELPA/ELSI.

  Please ensure these are built for optimum performance!


## Dependencies

Siesta relies on a number of libraries, some required and some
optional. The CMake approach facilitates the handling of these dependencies,
and several modes of operation are available:

* You can pre-install (some of) the libraries and
  provide CMake with hints about where to find them. This can be done
  by setting the CMAKE_PREFIX_PATH environment variable to the paths
  of the installed packages. For example:

   ```
      export CMAKE_PREFIX_PATH=/path/libxc/share/cmake:/path/libgridxc/share/cmake
      cmake ...
   ```

   CMake to look for the required packages in the directories specified in CMAKE_PREFIX_PATH.

   Alternatively, you can set CMAKE_PREFIX_PATH as a CMake variable when running the cmake command. For example:

   ```
      cmake ... -DCMAKE_PREFIX_PATH=/path/libxc/share/cmake;/path/libgridxc/share/cmake
   ```

   Please note the different delimiters used: ":" for Unix OS and ";" for
   CMake list separator, and remember to replace `/path/...` with the
   actual path to the installed libraries on your system.

* The source for (some of) the external libraries can be made
  available in the External/<package> subdirectories by activating git
  submodules within the Siesta distribution:

   ```shell
      git submodule update --init --recursive
   ```

   This option requires versions of git 2.13 and above, and is only
   available for versions of Siesta that have been accessed via `git
   clone`. Release and source code packages obtained through direct
   downloads from Gitlab do not feature this capability. If CMake
   detects source code is those directories, it can compile it
   on-the-fly.

* The source for (some of) the external libraries can be downloaded
  automatically and compiled on-the-fly

   If users do not have internet access on the compiling machine one
   must send the sources by other means. To aid this procedure one may
   use the `stage_submodules.sh` script to gather all sources for
   later uploading. (More on this in the [build-tools][_build-tools]  repo)

Here is the list of dependencies together with their options, with notes on whether they can
be compiled on-the-fly or not.

### System software

These are functionalities that should be already installed in most
systems devoted to scientific computing. Please check your site's
documentation about available versions.

Note that this "system software stack" should be consistent, and the
libraries (e.g. LAPACK, ScaLAPACK, etc) should have versions that
match the compiler and MPI system used. Software-stack support can get
complicated pretty quickly, and most sites use package managers such
as Spack or EasyBuild, offering users the possibility to load
environment modules to select the appropriate software
stack. Depending on the installation, either the libraries are installed
in "well-known" places or environment modules set their
variables in a CMake-friendly way, so CMake will detect the
appropriate libraries and no further action is needed on the part of
the user. If not, the CMake search can be helped by setting "hint"
variables, as explained below. We suggest that CMake is run first
without any extra libraries, and the log inspected to see if the
appropriate libraries have been found.

#### BLAS (required)

- `BLAS_LIBRARY=<name of library>|NONE` specifies the library name
  for linking.

   The `NONE` setting should can be used when BLAS is implicitly linked through other libraries (for example, when OpenBLAS provides LAPACK),
   or by the compiler itself (e.g. the Cray compiler )
   
- `BLAS_LIBRARY_DIR=<path to library>` place where to find the library  `BLAS_LIBRARY`

- `BLAS_LINKER_FLAG` flags to use when linking

Example:
```shell
   cmake ... -DBLAS_LIBRARY=blis \
             -DBLAS_LIBRARY_DIR=/opt/blis/lib
```


#### LAPACK (required)

- `LAPACK_LIBRARY=<name of library>|NONE` specifies the library name
  for linking. If `NONE` LAPACK is implicitly linked through other
  libraries/flags or the compiler itself (e.g. Cray).
  
- `LAPACK_LIBRARY_DIR=<path to library>` place where to find the library
  `LAPACK_LIBRARY`
- `LAPACK_LINKER_FLAG` flags to use when linking

Example:
```shell
   cmake ... -DLAPACK_LIBRARY=openblas \
             -DLAPACK_LIBRARY_DIR=/opt/openblas/lib \
             -DBLAS_LIBRARY=NONE
```

#### MPI (highly recommended)

- Use `SIESTA_WITH_MPI=ON|OFF` to enable, disable support respectively.

This MPI option defaults to be `ON` when an MPI installation is found. If one
is found but MPI is not wanted, the option must be explicitly set to OFF in
the command line:

```shell
   cmake  -DSIESTA_WITH_MPI=OFF ...
```

If the MPI library supports it, the `mpi_f08` module will be used for
the MPI interfaces. If not supported, the legacy interfaces
implemented in Src/MPI will be used. Full user control can be achieved
with the variable

```
SIESTA_WITH_MPI_INTERFACES [f08|legacy|none|auto]
```

The `none` value must be set if the compiler/MPI library combination
cannot use the f08 interfaces nor the legacy interfaces. The 'auto'
choice (that can be omitted) indicates no preference, and exists due
to a technical issue with older versions of CMake.


#### ScaLAPACK (required for MPI support)

- `SCALAPACK_LIBRARY=<name of library>|NONE` specifies the library name
  for linking. If `NONE` ScaLAPACK is implicitly linked through other
  libraries/flags or the compiler itself.
- `SCALAPACK_LIBRARY_DIR=<path to library>` place where to find the library
  `SCALAPACK_LIBRARY`
- `SCALAPACK_LINKER_FLAG` flags to use when linking

Example:
```shell
   cmake ... -DSCALAPACK_LIBRARY="-lmkl=cluster" \
             -DBLAS_LIBRARY=NONE -DLAPACK_LIBRARY=NONE
```
For certain versions of unix-based Scalapack installations, this flag needs to
be set explicitly to `-DSCALAPACK_LIBRARY="-lscalapack-openmpi"` (or the variant
you want to use). If not, CMake will fail to find it.


#### NetCDF (highly recommended)

Enable writing NetCDF files for faster (parallel) IO and
also for easier post-processing utilities.


- `SIESTA_WITH_NETCDF=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if NetCDF can be found (i.e. specifying `NetCDF_PATH`
  should be enough).
- `NetCDF_ROOT|PATH=<path to installation>` to specificy the location
  of the NetCDF installation. Generally should be enough with this
  flag.
- `NetCDF_INCLUDE_DIR` to manually specify include directories
  for modules etc.

There are some extra support libraries based on NetCDF which are shipped
with Siesta and are required for TBtrans:

- `SIESTA_WITH_NCDF=ON|OFF` enable NCDF support (a wrapper around NetCDF
  that makes it easier to work with).
  This is automatically detected, and the default is sufficient.

- `SIESTA_WITH_NCDF_PARALLEL=ON|OFF` allow parallel IO through NCDF library.
  This is automatically detected, and the default is sufficient.


#### OpenMP

Enable threading support using OpenMP.

Both TranSiesta/TBtrans may benefict from OpenMP when running large systems
which may result in performance gains plus memory reductions. The rest of
Siesta will not benefit as much.

- `SIESTA_WITH_OPENMP=OFF|ON` to disable, enable support respectively.
  By default it is `OFF`.

Users are recommended to test whether it makes sense for them to
use the threading support.

Be aware of `OMP_NUM_THREADS` and `OMP_PROC_BIND` variables
which may highly influence the performance gains.


### Required domain-specific libraries

These libraries were formerly part of the Siesta distribution, but now
they are packaged separately and made available to third parties. They
are part of the offerings of the Electronic Structure Library (ESL),
and they feature as modules in the MaX ("Materials at the eXascale")
Center of Excellence.

There might be pre-installed versions of (some of) these libraries. If
not, they can be fetched and compiled on-the-fly if there is an
internet connection.  Alternatively, the source code can be made
available by instantiation of a git submodule (it will be placed in
the External/<package> folder, or by manual deposit (in that one or in a
custom source directory). All these options are controlled by
the setting of the `<PACKAGE>_FIND_METHOD` CMake variable, where
`<PACKAGE>` stands for the uppercase form of any of the package names
below. By default, its value is:

```
   <PACKAGE>_FIND_METHOD="cmake;pkgconf;source;fetch"
```

allowing all available kinds of library-search strategies, which
will be tried in this order:

- `cmake` will search using CMake's `find_package`
- `pkgconf` will search using the pkg-config wrappers within CMake.

     `cmake` and `pkgconf` are generically implemented using package
     finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH` and
     `PKG_CONFIG_PATH` are (prepended)/appended with the directories
     that should be searched.

- `source` will compile the sources if found in `<PACKAGE>_SOURCE_DIR`,
  which typically defaults to `External/<package>`. The sources can be
  checked out via git submodule instantiation, manually cloned,
  unpacked from a release archive, etc.

- `fetch` will fetch the source , using the variables:

   * `<PACKAGE>_GIT_REPOSITORY`: the URL of the Git repository to
      clone the sources.  It defaults to the original development
      site. Changing it may be useful for testing forks with
      different functionaties (mostly for developers or advanced users).
       
   * `<PACKAGE>_GIT_TAG`: points to the revision to be checked out
     (tag, branch, or commit hash).


NOTE: The value of the `<PACKAGE>_FIND_METHOD` variables defaults to that of `SIESTA_FIND_METHOD`,
which in turn defaults to `cmake;pkgconf;source;fetch`. 

#### libfdf (required)

(Submodule) sources in the External/libfdf folder. 

#### libpsml (required)

Enables the reading of pseudopotential files in the `PSML` file
format, such as those downloaded from www.pseudo-dojo.org.

(Submodule) sources in the External/libpsml folder.

#### xmlf90 (required)

It is used in Siesta as a dependency of libpsml, and also to produce XML output.

(Submodule) sources in the External/xmlf90 folder.

#### libgridxc (required)

Evaluates the XC energy and potentials on the grids
where the density is calculated. It can leverage the libxc library, which
is optional (see below) but highly recommended.

The libgridxc library used must be compatible with the MPI and
grid-precision settings in SIESTA (e.g. those controlled by the
`SIESTA_WITH_MPI` and `SIESTA_WITH_GRID_SP` variables). This is automatic
if it is compiled on-the-fly. If if it pre-installed, CMake will verify
the compatibility.

Extra variables:

- `LIBGRIDXC_MIN_VERSION` can be used to allow the use
  of a (pre-compiled)  lower minimum version (the default is "2.0.1").
  Note that the use of older versions of libgridxc might
  cause some malfunctions. For example, the new default of
  reparametrizing the pseudopotential grids might not
  work well with vdw functionals. The absolute minimum
  version is 0.10.0.

(Submodule) sources in the External/libgridxc folder.

### Optional libraries

#### simple-DFTD3 (optional)

Add support for DFTD3 dispersion corrections as suggested by Grimme et.al.

(Submodule) sources in the External/DFTD3/s-dftd3 folder.

- `SIESTA_WITH_DFTD3=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if the `External/DFTD3/` directory contains
  directories with the appropriate sources.

The variables controlling this package's handling are prefixed by `S-DFTD3`
(e.g. `S-DFTD3_FIND_METHOD`). The package depends on a number of sub-packages
that are handled automatically. Their names and their associated prefixes are:

- mctc-lib (MCTC-LIB)
- mstore (MSTORE)
- test-drive (TEST-DRIVE)
- toml-f (TOML-F)


#### libxc (highly recommended)

libxc is an XC functional library implementing a very large variety of functionals.

libgridxc can leverage libxc and use its functionals.

- `SIESTA_WITH_LIBXC=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if the library can be found.
  
This library must be pre-installed, with a valid set of pkg-config `.pc`
files. Its search is controlled by the variables:

   * `CMAKE_PREFIX_PATH`: One of its components must point to the root
     of the libxc installation.
   
   * `LIBXC_Fortran_INTERFACE=f03;f90` to search for a specific
     interface.  The default is to accept both, with preference for
     the f03 interface. To only search for f90, set
     `-DLIBXC_Fortran_INTERFACE=f90`.


#### ELPA (native interface) (recommended)

[ELPA](https://elpa.mpcdf.mpg.de/) provides efficient and scalable
direct eigensolvers for symmetric (hermitian) matrices.

This refers to the native interface to ELPA in Siesta.  ELPA is also
offered by the ELSI library of solvers, both as built-in and as
'external library'.  If ELSI is compiled with an external ELPA, the
same ELPA library must be used.

- `SIESTA_WITH_ELPA=ON|OFF` to enable, disable support respectively.
  A search for an ELPA library is carried out, and the option set
  to `ON` if found.

The ELPA library must be pre-compiled, and it will be found if its root
installation directory is in `CMAKE_PREFIX_PATH`.  See
`Config/cmake/Modules/FindCustomElpa.cmake` for further details if problems
arise, particularly with older versions of ELPA.

The CMake system will try to determine automatically the features of
the ELPA library found, in particular whether it has GPU support, and
the names to be used in the interface for GPU kernels and options.
This is achieved by analyzing the output of the 'elpa2_print_kernels'
program that is installed together with the ELPA library. If for some
reason this program cannot be used or trusted, there is the
possibility of setting the relevant variables automatically. For
example:

```
-DELPA_SET_GPU_STRING="nvidia-gpu"
-DELPA_SET_REAL_GPU_KERNEL="ELPA_2STAGE_REAL_NVIDIA_GPU"
-DELPA_SET_COMPLEX_GPU_KERNEL="ELPA_2STAGE_COMPLEX_NVIDIA_GPU"
```

Note that if an external ELPA library is not available, Siesta can
still use its native ELPA interface through the built-in ELPA library
in ELSI (see below). This is automatically set by the build system if
ELSI support is configured (which is the default).

#### ELSI (recommended, used by default)

[ELSI](https://wordpress.elsi-interchange.org/) is a library of
electronic-structure solvers, offering a unified interface. 

- `SIESTA_WITH_ELSI=ON|OFF` to enable, disable support respectively.
  Defaults to `ON`.

The ELSI library interface is now by default compiled on-the-fly. The
source can be downloaded automatically, or it can reside in a
submodule (instantiated in the External/ELSI-project/elsi_interface
directory).

The compilation options are controlled by the following variables
(default values are in parenthesis):

   - `SIESTA_WITH_ELSI_PEXSI`: Enable PEXSI support in ELSI (ON by default)
     It might be deactivated at configuration time if proper versions of
     bison/flex are not found. Specifically, the compilation of the
     scotch library dependency needs flex 2.6.4.
   
   - `ELSI_WITH_GPU`: Enable GPU support in ELSI (OFF)
   
   - `ELSI_WITH_EXTERNAL_ELPA`: Use an external ELPA library in ELSI.
      (ON if an ELPA library is available. Highly recommended)
      
   - `ELSI_WITH_ELPA2021`: Use ELPA 2021 (NOT as well tested as the default ELPA2020, OFF).
   
   - `ELSI_WITH_TESTS`: Compile (Fortran) tests in ELSI (ON)
   
   - `ELSI_Fortran_FLAGS`,`ELSI_C_FLAGS`, `ELSI_CXX_FLAGS`: Compiler
      flags for ELSI compilation. Their values are passed directly to the
      `CMAKE_<LANG>_FLAGS_<CONFIG>` set of variables for the duration
      of the compìlation of the ELSI subproject. If the top-level CMAKE_BUILD_TYPE
      is not set for any reason, a custom ELSI config type is created.

Note that if an external ELPA library is not available, Siesta can
still use its native ELPA interface (usually requested with SIESTA_WITH_ELPA)
through the built-in ELPA library in ELSI. This is automatically set by
the build system.

GPU support in ELSI can come from an external ELPA library (which can
be compiled for an appropriate GPU architecture, i.e. nvidia, amd, etc)
or from the built-in ELPA library (in this case only Nvidia cards with CUDA
are supported). This is controlled by the `ELSI_WITH_GPU`
variable. Additionally, the following variables are relevant:

   - `CMAKE_CUDA_ARCHITECTURES`: Default `80` (A100 class)
   - `CMAKE_CUDA_FLAGS`  (default `-arch=sm_${CMAKE_CUDA_ARCHITECTURES} -O3 -std=c++11`)
   

For advanced users only: A precompiled ELSI library can be used if the
`ELSI_ALLOW_PRECOMPILED` variable is set to `ON`.  The library will be
found if its root installation directory is in `CMAKE_PREFIX_PATH`.
Note that pre-compiled libraries do not in general provide mechanisms
for the easy discovery of information about their capabilities
(i.e. external ELPA library, GPU support). Extra variables might have
to be set from the command line in this case:

- ELSI_HAS_EXTERNAL_ELPA
- ELSI_ELPA_HAS_GPU_SUPPORT
- ELSI_ELPA_GPU_STRING
- ELSI_ELPA_REAL_GPU_KERNEL
- ELSI_ELPA_COMPLEX_GPU_KERNEL
- SIESTA_ELPA_HAS_GPU_SUPPORT

The use of pre-compiled libraries is discouraged for this reason.
On-the-fly compilation is preferred for its flexibility.

#### PEXSI (native interface)

Note: This might be mildly outdated regarding build options. In general,
the PEXSI subsystem of ELSI is to be preferred.

The Pole EXpansion and Selected Inversion (PEXSI) method evaluates
certain selected elements of matrix functions, e.g., the Fermi-Dirac
function of the KS Hamiltonian, yielding the density matrix. It can be
used as a highly scalable and efficient alternative to diagonalization
methods.

- `SIESTA_WITH_PEXSI=ON|OFF` to enable, disable support respectively.
  Defaults to `OFF`.

This refers to the native interface to PEXSI in Siesta.  PEXSI is also
offered by the ELSI library of solvers, both as built-in and as
'external library'.  If ELSI is compiled with an external PEXSI, the
same PEXSI library must be used.


The PEXSI library must be pre-compiled, and
it can be found in several ways:

- If an ELSI library (with PEXSI) is available, it will be used if
CMAKE_PREFIX_PATH includes the path to the ELSI installation.

- For (pre-release) PEXSI 2.1 (which uses a new scheme with subordinate
libraries), the environment variable PEXSI_ROOT must point to the
PEXSI installation directory.

- For PEXSI v2.0, CMAKE_PREFIX_PATH must include the path to the
PEXSI installation, and it will use the CMake config package.

The PEXSI libraries can be obtained from

- V2.1 (pre-release): https://bitbucket.org/berkeleylab/pexsi/src/for2.1/
- V2.0: https://bitbucket.org/berkeleylab/pexsi/src/master/

See the file `Config/cmake/search_for_native_PEXSI_provider.cmake`
for more details.

#### CheSS (experimental)

CheSS implements a linear-scaling Fermi-Operator-Expansion method with
Chebyshev polynomials. It was developed by the BigDFT project.

- `SIESTA_WITH_CHESS=ON|OFF` to enable, disable support respectively.
  Defaults to `OFF`.

See `Config/cmake/Modules/FindCustomCHESS.cmake` for details on
how to link against CheSS. You might have to do some work to compile
the needed dependencies from the BigDFT project, and tweak the finder module.

#### FFTW

The FFTW library is currently only used in the Util/STM/ol-stm utility.

- `SIESTA_WITH_FFTW=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if found.

#### FLOOK

[Flook](https://github.com/ElectronicStructureLibrary/flook)
implements an interface to access internal Siesta variables from a
built-in Lua interpreter, enabling the creation of custom scripted
functionality (e.g. control of molecular dynamics trajectories, novel
geometry-relaxation methods) without recompiling the code.

- `SIESTA_WITH_FLOOK=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if found.

(Submodule) sources in the External/Lua-Engine folder. 

An installed version might be found if its root directory is contained
in `CMAKE_PREFIX_PATH` or in the *environment* variable `FLOOK_ROOT`.

If not found in the default directory, the source can be fetched
automatically from the flook development site, or from the
URL contained in the *environment* variable `FLOOK_PACKAGE`.

If you already have an installation of lua that is compatible with flook,
you can set the `LUA_DIR` *environment* variable to the directory containing
`lib/liblua.*`. Then, flook will not compile its own version of lua and
will use the provided one.


Note: Currently FLOOK cannot be compiled in the submodule to be usable
      as a shared library. Hence, if you need flook, then:

        SIESTA_SHARED_LIBS=OFF

      must be used (default).

#### Wannier90 "wrapper-library" interface

A wrapper interface between Siesta and wannier90 (version 3.1.0) so
that wannier90 can be called from Siesta on-the-fly.

- `SIESTA_WITH_WANNIER90=ON|OFF` to enable, disable support respectively.
  Defaults to `OFF`.

Turn on the `SIESTA_WITH_WANNIER90` option and use an *environment*
variable `WANNIER90_PACKAGE` that points to a pristine
`wannier90-3.1.0.tar.gz` file (a remote URL is also allowed).  This
file will be unpacked, patched, and compiled into a wrapper library
that Siesta can use directly.


### Other Siesta Options

Siesta provides other options to control some internal features. The
generic Siesta executable should be sufficient for most, but some
users may have different needs.

- `SIESTA_WITH_GRID_SP=OFF|ON` use single-precision grid operations (Default `OFF`).
  It can reduce the memory requirements for large mesh-cutoffs
  and/or large unit-cells, at the expense of some precision.
  The default is to use double precision `-DSIESTA_WITH_GRID_SP=OFF`

If you enable this feature you might consider adding the variable
`SIESTA_SUFFIX=_grid_sp` so that the Siesta executable and libraries
are properly tagged in this case.


### Shared libraries

It is possible to build Siesta in *shared* library mode.
Simply pass:

```
  cmake ... -DSIESTA_SHARED_LIBS=ON
```
and everything will be compiled in shared mode.
This will not set the runtime paths (RPATH) of the executables/libraries.
If you want the installed targets to have the runtime path specified
simply add `-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON` as well.

Note: FLOOK does *not* work with shared library mode.


### Profiling

Support for profiling is overloaded in the "timer" interface. The currently enabled
profiling options are:

- `SIESTA_WITH_PROFILE_NVTX=ON|OFF`  (Defaults to `OFF`)
  If enabled, the variable `SIESTA_PROFILE_NVTX_LIBRARY` must
  contain the path to the `libnvToolsExt` library by nVidia.

The EXTRAE library by BSC can be used to generate MPI traces, but its configuration is not
yet supported by the CMake build system. See the next section.

### Use of options not supported directly by CMake

Some options are not yet supported directly by the CMake framework, but
they can still be supported indirectly.

For example, if it is desired to use the MUMPS library in TranSiesta, the library
needs to be pre-compiled, and the following CMake variable specified:

```
   -DFortran_FLAGS="-DSIESTA__MUMPS <other flags>"
   -DSIESTA_LINKER_FLAGS="-L/path/to/mumps/lib -lzmumps -lmumps_common <other libraries>"
```
where '<>' are any libraries that MUMPS depends on, and '/path/to/mumps' is the root
of the MUMPS installation.

This mechanism is general, even if a bit coarse. It is expected that native CMake support
for all options will be offered soon.

### Cross compiling

Cross compilation is generally automatically detected and will remove
any host-specific flags for the compilation. If one wishes to disable
host-specific flags, one can use the flag:
```
   # When cmake does not detect cross-compilation
   cmake ... -DSIESTA_WITH_HOST_OPTIMIZATION=OFF
   # When cmake detects cross-compilation (i.e. -DCMAKE_SYSTEM_NAME=something)
   cmake ... -DCMAKE_SYSTEM_NAME=alternate-system
```

When cross-compiling Siesta one needs to pass the appropriate kinds
for the compiler used. This is necessary because the automatic
generation of kind information would require the execution of
binary auxiliary programs, which cannot be done when cross-compiling.

Simply define these flags
```
 -DSIESTA_INTEGER_KINDS="4;8" -DSIESTA_REAL_KINDS="4;8"
```
where one might change 4 and 8 with the appropriate integer32, integer64, real32 and
real64 kinds. Then the generation of binary executables will be omitted.


## Settings for compiler flags and build types

### Compilation flags

Compilation flags are generally managed through the environment
variables (NOT CMake variables).

- `FC` for specifying the fortran compiler
- `FFLAGS` for specifying the compilation flags

An invocation might be:
```shell
   FC=gfortran FFLAGS='-O3 -march=native' cmake ...
```

Alternatively, the flags can be supplied on the command line
```shell
   cmake -DFortran_FLAGS=-Os -DC_FLAGS=-Os
```
This enables fine tuning of the compiler flags.

Customarily, CMake uses the `CMAKE_<LANG>_FLAGS`.
These may also be used

### Toolchains for compiler flags

Siesta's build system supports the use of `SIESTA_TOOLCHAIN` files to
gather compiler flags. (This is different from, although related to
the functionality provided by the `CMAKE_TOOLCHAIN_FILE` variable; see below)

This feature can be used in several ways:

```shell
   # Use default toolchain files located in Config/cmake/toolchains
   
     cmake ... -DSIESTA_TOOLCHAIN=gnu
     
   # or a full path (for local edited files)
   
     cmake ... -DSIESTA_TOOLCHAIN=/path/to/toolchain/file/gnu.cmake
```

When using `SIESTA_TOOLCHAIN` one can use multiple toolchains.
This can be valuable for overwriting or adding variables from various
toolchains (mainly useful for developers). For example:

```shell
   cmake -DSIESTA_TOOLCHAIN=gnu;local ...
```

to use `./Config/cmake/toolchains/gnu.cmake` and `./local.cmake`.

Currently the default toolchain to use will be decided internally as follows:

- GNU compilers will use the `Config/cmake/toolchains/gnu.cmake`
  toolchain file.
  
- Intel (and the newer Intel LLVM backend) compilers will use the
  `Config/cmake/toolchains/intel.cmake` toolchain file.
  
- Otherwise a _generic_ toolchain file will be used, which uses the
  default CMake variables.

To gain complete control of the compiler flags (without adding the toolchain
ones) you will have to select the `none` toolchain and set the flags by hand:

```shell
   cmake -DSIESTA_TOOLCHAIN=none -DFortran_FLAGS="-Os -DSOME__VAR"
```

### Custom toolchain files and pre-loading of cache variables

A custom toolchain may contain any setting of variables. They can
be thought of as an `arch.make` file with default parameters.

Toolchains for a few sample computer systems can be found in
`Config/cmake/toolchains`, but note that it is impractical to cover
all possibilities. Further examples and documentation for building can
be found in [this repository](https://gitlab.com/siesta-project/ecosystem/build-tools).

There are two ways to load a custom toolchain file:

```
   # Use of the -C option to pre-load the (cache) variables contained in the file:
   
     cmake ... -C Config/cmake/toolchains/leonardo.cmake
     
   # Use the form:
   
     cmake ... -DCMAKE_TOOLCHAIN_FILE=Config/cmake/toolchains/leonardo.cmake
```

The second form is somewhat of an abuse of the CMake functionality intended for cross-compilation. The first form loads
the variables at the beginning of the CMake configuration run, before the `project()` call is executed, so custom
toolchain files loaded with `-C` can contain settings for the build type (see next section).

Parameters that exists in a toolchain file can be overridden on
the command-line with `cmake -D<VAR>=<VALUE>`.

Developers are suggested to create custom toolchain files with the
appropriate variables. CMake's `presets` functionality could also be
used, although it is not officially supported by Siesta yet.

### Build type

CMake uses a build-type to help determine the flags to be used. The default
build type (`Release`) should be sufficient for most (if not all users), but
other types can be enabled with, for example:

```shell
   cmake -DCMAKE_BUILD_TYPE=Debug
```

Currently the default Siesta toolchain files support these build types:

- `Release`: the default and recommended build type, it uses a high optimization
  level without sacrificing accuracy.
  
- `Debug`: used for debugging Siesta, or if there are runs that show problems.
  *Bug reports* should use this build
  
- `Check`: used for debug + checking code execution runs, primarily
  useful for developers; equally good for bug-reports.
  
- `RelWithDebInfo`: a release mode with debug symbols.

- `MinSizeRel`: optimizes the executables for minimum size (`-Os`)

Compiler flags associated to the different build types
can also be specified in the command line. For example:

```shell
   cmake -DFortran_FLAGS="-Wall" -DFortran_FLAGS_DEBUG=-g -DCMAKE_BUILD_TYPE=Debug
```

will result in the internal variable `CMAKE_Fortran_FLAGS` evaluating to the string `-Wall -g`, because
`Fortran_FLAGS` is sort of a global setting, while `Fortran_FLAGS_DEBUG` applies only to the `Debug` build type.

Custom build types can also be specified if needed for particular uses:

```shell
   cmake  -DFortran_FLAGS_EXOTIC="-fcrazy-unroll -ffind-bugs"  -DCMAKE_BUILD_TYPE=Exotic
```

The currently supported build-types in the shipped toolchain files are:

- `Fortran_FLAGS`
- `Fortran_FLAGS_RELEASE`
- `Fortran_FLAGS_DEBUG`
- `Fortran_FLAGS_CHECK`
- `Fortran_FLAGS_RELWITHDEBINFO`
- `Fortran_FLAGS_MINSIZEREL`


### Interprocedural optimizations

Siesta enables compilation of its executables using *interprocedural optimization*.
This can give some performance by optimizing code across libraries/modules/objects.
For production systems you are encouraged to utilize this option for best performance:

```shell
   cmake ... -DSIESTA_WITH_IPO=true
```

Note, your compiler has to support this, CMake will stop if it has
problems with it. Expect a *very long* and memory demanding linker step!


### Custom suffixes for executables and/or libraries

By default, Siesta executables and libraries are installed with plain
names such as `siesta`, `tbtrans`, `libsiesta.a`, etc.  A custom
suffix can be added by using the variables

- SIESTA_SUFFIX (for executables and libraries)
- SIESTA_EXECUTABLE_SUFFIX
- SIESTA_SHARED_LIBRARY_SUFFIX
- SIESTA_STATIC_LIBRARY_SUFFIX

For instance, to change the suffixes for the executables:

```shell
   cmake ... -DSIESTA_EXECUTABLE_SUFFIX=.x
```
will create `siesta.x`, and similarly for all other executables.

This can be used to have multiple executables in the same `PATH` using
different backends. For example, a GPU enabled executable can be build
with `-DSIESTA_EXECUTABLE_SUFFIX=_gpu`.


### Inspecting the configuration log file

The first part of the CMake invocation:

```shell
   cmake -S. -B_build -DCMAKE_INSTALL_PREFIX=/path/to/installation  #<MORE_OPTIONS>
```

will produce lenghty output related to the configuration steps, the variables used, packages
found (or not found), etc. A final summary of all the options applied is also shown.

This output can appear overwhelming at first, but it pays to read it with some care, particularly if errors appear. It
is also important to keep it in case a building issue has to be brought to the developers.
In fact, in some cases *even more output* might be useful.

Using

```shell
   cmake --log-level=debug ...
```

will show more internal steps. 


### Building in parallel (recommended)

To build the project in parallel simply add the `-j N` option. For example:

```shell
cmake --build _build -j 4
```
will build the project using 4 processes.


## Tests

Siesta's build system integrates a testing framework, which is exercised
by executing the ``ctest`` program under the build directory.
This will first run tests, and then subsequently use a Python script
to compare the reference output with the test output.
```shell
cd _build
ctest <options>
```

If the required external libraries have been compiled as part of
the current CMake invocation, installation tests for them will
also be executed.

If you want to run only a minimal set of Siesta tests, you can select
the simplest representative cases by choosing the label "simple":
```shell
ctest -L simple <options>
```

One can test separately the execution of the examples and the
verification of the numerical precision of the results (when compared
to the provided reference outputs). The latter test involves the
setting of appropriate tolerances, and its results might sometimes
depend on the libraries used or on whether MPI is used or not. Further
work is being done on the architecture of the verification process, to
make it more robust and complete.

For output verification you will need a local Python
installation with the `ruamel.yaml` package installed.
As the verification is *on* by default, you should only disable this
verification step if you don't have access to Python3 with the above package.
```shell
ctest -LE verify # Runs and skips verification for all tests.
```

You can also run a single test using the -R option with the name
of the test itself; for example:

```shell
ctest -R default_basis
```

## Spack support

[Spack](https://spack.io) is a package manager that can ease the task of maintaining
and configuring the software stack.

We make available spack recipes and sample spack environment files for Siesta and
relevant needed libraries in [this auxiliary repo][_build-tools].

## Deprecated options

### Legacy form of the options

Most of the options above are prefixed by `SIESTA_`. 
The simpler forms without the prefix (e.g. `WITH_MPI` instead of
`SIESTA_WITH_MPI`) are **no longer accepted** by default.

### Units definitions

Since Siesta 5.0 the internal units in Siesta are based on
CODATA-2018, so some numerical differences from older versions can be
expected. In general these differences will be quite small, but
sometimes the results might be appreciably different from older
calculations, particularly for molecular dynamics runs.

One can select the units employed by using this cmake invocation:

```shell
# this will use the same units as in the 4 series and older
cmake ... -DSIESTA_WITH_UNIT_CONVENTION=legacy
# this will use the new units (the default)
cmake ... -DSIESTA_WITH_UNIT_CONVENTION=codata2018
```

It is **NOT** recommended to use the `legacy` scheme for anything but
tests. The new units should *always* be preferred.


[_build-tools]: https://gitlab.com/siesta-project/ecosystem/build-tools

