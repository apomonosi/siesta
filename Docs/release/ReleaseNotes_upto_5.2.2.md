# Release Notes for Siesta

This file contains, in an expanded ChangeLog style, notes
regarding the evolution of the Siesta code. The changes are
grouped under headings representing past (and upcoming) releases.

# 5.2.X versions

### 5.2.2 (2025-02-04) Bugfix release

* Clarify the conditions under which PEXSI support in ELSI is compiled.

* Refactor citation support. Added more references.

* Add CI support for automatic generation of .rst files for siesta-docs.

* Enabled GPU activation in ELSI build DM/EDM operations with external ELPAs.

* Bumped default electronic temperature from 299.9863 K to a
  more standard 300 K (changed from default in Ry to default in K).

* Make Fortran required in the find_package() call for MPI.

* Added documentation for TBT.k.File for completeness.

* Bugfix: Fixed memory leaks for molecular dynamics simulations that
  came from the new mesh algorithm.

* Bugfix: fixed ELSI not respecting "Number-of-Eigenstates" option.

* Bugfix: Now the N of a polarization orbital is bumped if it's lower than the
  one of the orbital generating it.

* Bugfix: IntDos2Eig util expecting integers instead of reals in intdos file.


### 5.2.1 (2024-12-11) Bugfix release

* Fixed bug for parallel runs where a segmentation fault would occur.

* Added the vdw_D3 energy to TranSiesta energy.

* Reworked D3 periodicity guessing to rely on cell_periodic.
  Now it should cover most cases.

* Further removal of backspace calls in utilities.

* Removed legacy build instructions for tbtrans.

* Improved documentation in many aspects.


### 5.2.0 (2024-11-05) Feature release

* Improve the performance of grid-related communications by using
  asynchronous techniques (by Rogeli Grima, BSC). The BSC
  modifications have been made compatible with the variable-nsm
  redistribution scheme by Federico Pedron, who also added documentation.
  Also, replace the algorithm for pre-computation of the communication
  pattern for the redistribution of orbital-based matrices from
  block-cyclic to grid-based for better scalability.

* New on-the-fly compilation of ELSI with automatic handling of its
  dependencies.

* Interface to the ELSI library of solvers, in "DM mode" (i.e.,
  calling directly the integrated solvers that take H and S and return
  DM and the EDM). Support for collinear spin and k-point sampling is
  available (in fully parallel mode). The non-collinear spin case is
  not yet supported.

* Streamlined atomic charge output: now Mulliken, Voronoi and Hirshfeld
  charges exhibit the same behaviour across different parts of the code,
  with the possibility of printing any of the three during, before or after the
  SCF cycle. In this regard, new options are available for atomic charge
  printing.

* Prettified atomic charge output, while also making it easier to post-process.

* Added new toolchain to fix Intel LLVM-based compilations (ifx).

* Added -Monline option for NVHPC compilations.

* Added profiling support for NVTX.

* Bumped DFT-D3 version to avoid issues with Cray compilers.

* Added workaround on flook interface for Cray compilers.

* Fixed total energy output in MD.nc output.

* Fixed ioeig being called from all nodes, thus producing an I/O bottleneck.

* Fixed a bug introduced in 5.0 that prevented the printing of the TDEIG file.

* Removed references to k-cutoff in test suite.

* Removed m_hs_matrix, which was never used.

* Modularized manual to use in tandem with the Siesta Docs webpage.

* Added thread-count to times file.

* Expanded information offered in the build info.

* Denchar now prints dummy atom for compatibility with VMD and XCrysden.


# 5.0.X versions

### 5.0.2 (2024-10-29) Bugfix Release

* Improvements to the handling of the GPU set-string and kernel names for ELPA.

* Fixed bug in pressure annealing scaling factor which resulted in incorrect
  behaviour of the method, as presented in
  https://gitlab.com/siesta-project/siesta/-/issues/30.

* Fixed the version generation script to print proper version numbers.

* Fixed missing rcore initialization in atom.F.

* Improved documentation for GF file options.

* Bumped libgridxc version to 2.0.2 to prevent issues with VdW functionals.

### 5.0.1 (2024-06-25) Bugfix Release

* Further improvements to the cmake infrastructure:
  * Enabled compatibility with shared projects (Note that building SIESTA as a
    shared library is not compatible with FLOOK yet).
  * Enabled manual integer and real kind especification for MPI interfaces.
  * Fixed libhsx name (previously "liblibhsx").
  * Fixed Scalapack/MPI priorities, now CMake will not look for Scalapack if a
    serial compilation was asked for.

* Added detection of missing valence PAO:
  * When users fail to provide a valence orbital in the PAO.Basis block that
    should be occupied, SIESTA will print an error.
  * If the shell is empty, SIESTA will print a warning and populate the existing
    orbitals accordingly.
  * This also includes checks for the presence of semicore states that might be
    required by the pseudopotential.

* Fixed the fallback mechanism for polarization orbital generation when extra
  semicore states are present ("pol_ptr not associated" error).

* Further fixes for CODATA standardization in Wannier-related routines.

* Added clarifications on Spin magnetic moment definition and units.

* Added further clarifications on backwards compatibility for basis set
  generation defaults.

* Fixed TSHS version determination errors with certain compilers.

* Increased robustness for ION file read, without relying on iostat to see
  if some piece of information is present in the file.

* Fixed SiestaAsSubroutine examples and updated test scripts.

* Removed E_bs debug printing.

* Fixed involuntary printing of the output EPSIMG file.

### 5.0.0 (2024-05-18) Feature Release

* Use new, prefixed, names for the CMake options and variables, e.g.
  SIESTA_WITH_MPI instead of WITH_MPI.  Old-style variables
  **are not accepted** by default.

* Enabled custom linker flags:

    SIESTA_LINKER_FLAGS_PRE
    SIESTA_LINKER_FLAGS_POST
    SIESTA_LINKER_FLAGS (same as the above flag)

  These may be used to pass additional flags to the linker stage of the
  binaries (currently required for Metis and MUMPS support).

  They can be used together with a new scheme that wraps the definition
  of targets and their dependencies.

  Some target names are now prefixed with SIESTA_, and some other variables
  are also prefixed.

  Improvements in CMake messages and summary information.

* Corrected an issue where the XV files were not read when only UseSaveData
  was set, #361

* Correct an issue that resulted in slightly inaccurate Mulliken
  overlap populations.

* Bring more up-to-date the interface to the CheSS library (experimental).

* Add support for profiling with the NVTX tools by nVidia. The calls
  are overloaded on the timing interface.

* Updated the native PEXSI interface to use new (post 1.2.0) routine names.
  A custom interface file allows using both standalone PEXSI (>= 2.0) and ELSI
  libraries.

* fixed serial writing/reading of NetCDF files, basically all .grid.nc files
  in spin calculations was stored wrong, only affecting serial compilations.
  Fixes #310

* fixed denchar bug for non-allocated arrays, #321

## 5.0.0-beta1 (2023-09-29) Beta Feature Release

### Backward compatibility issues

* Unit precisions has changed. From version 5 and onwards the units used will be based on the
  CODATA values. This change results in print-out changes and possibly some changes in many-step
  MD simulations. Printouts of energies in eV and lengths in Ang results in a change.
   v5 ENERGY [eV]  = 0.9999921447466303 * v4 ENERGY [eV]
   v5 LENGTH [Ang] = 1.0000003985490675 * v4 LENGTH [Ang]
  More details can be found in Src/libunits directory.

* The default FFT ProcessorY value is determined with a different scheme. Executions that are
  sensitive to this value can set the FFT.Processory.traditional fdf variable.

* The labels in the Mulliken analysis CML blocks have been changed to use "population" instead of "charge".

* The HSX file format has changed to reduce disk-space and increase precision.

* Mixtures of exchange-correlation functionals should now be specified in the XC.mix block.

* Defaults for basis sets have been changed:
   PAO.EnergyShift has been halved to 0.01 Ry (changed form 0.02 Ry)
   PAO.SoftDefault is now true (changed from false)
   PAO.OldStylePolOrbs is now false (changed from true)
   PAO.FixSplitTable is now true (changed from false)
   PAO.SplitNormH is now 0.45 (changed from PAO.SplitNorm)
   ReparametrizePseudos is now true (changed from false)

### Changes

* Bring the stand-alone PEXSI interface up to date with versions 2.X of the PEXSI library.

* Removed legacy makefile-based building system. CMake (>=3.20) can be used for configuration and building.

* New "wrapper-library" interface to Wannier90, allowing Siesta to drive Wannier90 on-the-fly
  instead of using the customary workflow of preparing the inputs, running Wannier90, and re-running
  Siesta. A separate download of the pristine Wannier90 sources is required.

* Added a new lua checkpoint, AFTER_MOVE, right after the new atomic coordinates are
  determined (Pol Febrer).

* New default of the units, precision increase and following CODATA-2018 values, see Backwards
  compatibility issues

* Implement a new unit-handling paradigm in fdf: the units table and associated 'inquire_unit'
  function are implemented by the client code (Siesta). The libfdf library processes units by
  calling an appropriate procedure pointer.

* Use an external libfdf library. This is a new required dependency, which can be satisfied via
  an installed package, a submodule, or a direct fetch from the libfdf repository.

* Extend the options for pseudopotential handling with ps-file-spec in the ChemicalSpecies block
  and SIESTA_PS_PATH environment variable (Alberto Garcia).

* Add XC information to CML file (@william-chem)

* Added support for D3 dispersion corrections (Federico Pedron)

   * Most, if not all, of the GGA functionals are supported by default (PBE, LYP and their
     variants), but not LDA or VDW functionals. Note that D3 corrections should not be used
     in conjuction with VDW functionals in any case. Mixed "cocktail" functionals
     are not supported either.
   * Note, however, that inputing custom D3 parameters is a feature of this implementation.
     This can be useful for those rare cases that are not supported.
   * A higher degree of control of the model is also posible, by changing all relevant cut-off
     values (2-body and 3-body interactions, and coordination number calculations). These
     already have reasonable values by default.
   * The implementation relies on the simple-DFTD3 library (https://dftd3.readthedocs.io).
     The relevant code can be incorporated as git submodules in External/DFTD3, and compiled
     on the fly (see INSTALL and INSTALL.CMake).

* New building framework, including CMake support and streamlining of makefile-based scheme

   * The fdf code in Src/fdf is now serial-only, with clients taking
     care of the broadcast of the fdf database object. This removes
     the need for multiple compilations of the library to suit each
     client.

   * The rest of the internal libraries (in directories under Src)
     have been made to follow the same scheme of single compilation:
     MPI, psoplib, ncps, MatrixSwitch, and the linear-algebra
     fallbacks in Src/Libs.

   * New versions of the ncdf and fdict libraries (renamed as
     'easy-fdict' and 'easy-ncdf') have been prepared and used to
     follow the single-compilation scheme.

   * Programs in Utils and Pseudo are now compiled under the same
     building directory as Siesta itself. The main user-facing change
     is the location of the setup script for the build directory. The
     typical step is now 'sh ../Config/obj_setup.sh' from the build
     directory. A new target 'utils' builds the whole suite, and new
     'install' and 'install_utils' targets copy everything to an
     installation directory.  This approach removes the need to handle
     OBJDIR when dealing with different builds, and simplifies the
     installation of all the programs.

   * Automatic compilation of required library dependencies can be
     achieved by setting a flag in arch.make, using library code in
     git submodules in ExtLibs/. CMake is used for the building of the
     libraries.

   * In addition to the above upgrades of the makefile-based building
     scheme, there is now an independent CMake framework, which enables
     mostly seamless compilation of the suite of programs. A variety of
     mechanisms for discovery and use of dependencies are supported.

   * New experimental spack recipes (intended to be used for now in a
     separately configured spack repo) are available.

   * CMake framework also supports the automated running of tests via ctest.


* Made the DFT+U scheme compatible with spin-orbit coupling.

* Added support for PSML pseudopotential files:

   * Changes to offer more information about PAO generation, and to implement
     new features such as automatic handling of semicore states.
   * New psoplib and ncps libraries in Libs handle KB generation and 'Froyen' interface
   * New stand-alone program Pseudo/vnl-operator/psop to generate PSML files with Siesta-style KB  projectors
   * New program Pseudo/converters/psml2psf
   * An external LibGridXC replaces the built-in SiestaXC
   * Libxc support provided through libGridXC (with xc specs in XC.mix block)
   * New mandatory dependencies: xmlf90, libpsml, libgridxc (with optional libxc)
   * Alloc module revamped to make logging work with external libraries
   * Memory logging revamped
   * Streamlined 'Util' makefiles with much smaller dependency lists
   * New, more modular, building specification with slimmer arch.make file.

* Allow parallelization over orbitals instead of over k-points in the unfold utility

* Allow mesh parallel distributions with different numbers of subpoint multiplicity.

* Minor changes for f2003 conformance

* Provide more options for the shift of the origin of coordinates.

* Some fixes for library operation (avoid stopping when 'onlyS' is set; re-opening of unit 6)

* Cleanup of vibra/fcbuild code: use dynamic arrays, and LAPACK solver by default

* Extended Hirshfeld and Voronoi partition analysis to spin. Update of output, including new CML blocks.

* Performance increase for tall+skinny matrices in TS

* Added band-unfolding utility (J.M. Soler)

* Overlap gradient now saved in *.nc file

* Allowing finer grained k-point input from users

* Off-site spin-orbit (full) coupling (Ramon Cuadrado)

* Cleaned data-structures for atom code

* Interface to the CheSS linear-scaling library (Stephan Mohr)

* Update internal siesta_forces logic wrt TDDFT

* Implementation of Real-Time TDDFT (Rafi Ullah)

* Added OpenMP nesting information (only print-out)

* Removed unused python scripts across all utilities.

* Fixed a lot of issues with Cray compilers.

* Now a default value of 1.0 Ang is used for the lattice constant.
  This ensures that a proper unit cell is recognized and Siesta
  does not default to the isolated molecule treatment.

* Changed the default value for BlockSize, to enable both better
  runtimes and improve the default compatibility with GPU-enabled
  versions of ELPA.

* Extended allocation functions to support 5-dimensional arrays.

* Huge rework of Siesta Tests directory:

  * All of the tests that were intended as physics-relevant examples can
    now be found under the /Examples folder.
  * Tests now are categorized according to the functionality they are
    testing within Siesta. Total runtimes have been shortened by a lot,
    and new reference outputs are provided.
  * Tests support ctest automation to run during after CMake builds.
  * Tests have a run script within each folder that provides an easy way
    to run with an already compiled version of Siesta.

* New defaults available for the following variables, which now reflect better
  their intended usage:
  * DM.UseSaveDM now defaults to true, unless UseSaveData (deprecated option)
    is specified. In that case, it will default to the value of UseSaveData.
  * SaveHS now defaults to true.
  * ChangeKgridInMD now defaults to true.

# 4.1.X versions

### Backward compatibility issues

* TranSiesta/TBtrans eta values for the device region are now defaulted
  to the minimal electrode eta value / 100

### Changes

* Sort shells in PAO.Basis by n quantum number

* Fixed reading Target.Pressure, reported on mailing list 24 March 2021.


## 4.1.5 (2021-01-27)  Feature release

### Backward compatibility issues

* Siesta is now developed on the GitLab platform: www.gitlab.com/siesta-project/siesta
  A number of Siesta-related packages are developed here: www.gitlab.com/siesta-project

* TranSiesta/TBtrans eta values now default to 1 meV

* TranSiesta equilibrium contours now default to using the "right" scheme

* Maximum l for KB projectors is now automatically set to the highest
  l in the PS file. This may result in slight changes.

* Bug in Voigt representation output, !25
  Users with scripts that parses stress-tensor-Voigt will
  need to adapt.

### Changes

* Removed memory leaks for k-point samplings in MD runs, #58

* If compiling with Gfortran 10, there will be lots of errors, please see manual

* Changed LDA+U to DFT+U to clarify its usage is not restricted to LDA

* Fixed regression in OMM routines due to resetting of qtots
  (introduced when fixing the NC/SOC occupation problem)

* Enabled the spin-spiral code, !20
  Feedback is requested from expert users of spin-spiral systems.
  _Not_ production stable.

* Added final break-point for LUA, LUA_FINALIZE just before Siesta stops.

* Added dipole calculation from vacuum region, !22
  For gated calculations this is the preferred method.

* Bug, fixed nnzs == nspin parallelization calculations which
  could revert order of DM reads.

* Removed some minor memory leaks in mesh-subs

* Bug, fixed Ha/Bohr unit in fdf

* Updated flook installation script to 0.8.1

* Update Docs/REPORTING_BUGS

* Document the setting of 'neigwanted' and print them if the diag solver allows it.

* Fix computation of NC/SOC occupations when the (optional) number of
  eigenstates handled ('neigwanted') is less than the number of orbitals, and
  enabled parallel k calculations !9

* Fix reading of wave-functions in Util/COOP/fat.f90

* Added Obj/ARCH-EXPERIMENTAL for suggested more modular building scheme

* Added interface code to use GPUs with the ELPA library.

* Removed a small memory leak (from siesta_init)

* Enabled parallel-over-k for NC/SOC calculations

* Fixed bug for writing eigenvalues when NumberOfOrbitals was used (together with NC/SOC)

* Fixed printing Voronoi and Hirshfeld charges if user requested LDOS calculations

* Fixed cell transpose when using socket calculations (only important for skewed cells)

  * Added citation information to output (end of run)

* Allowed Mesh.Sizes as list input so users can specify their own Mesh size

* Fixed *.nc writes, fixes lp:1810279.

* Fixes polarization issue with Bessel orbitals

* Made the optional BSC_CELLXC code accessible at run time, instead of
  through pre-processing at compile time.

* Lots of updates for documentation

* Precision problem in the vibra utility is fixed, lp:1816719

* Added wavefunction tools for spin-orbit calculations.
  Now denchar, COOP/COHP, fatbands, and stm utilities can process spin-orbit output.

* Removed *ALL* OMP collapse statements, Intel 2019 is buggy.

#### TranSiesta / TBtrans

* Added charge tolerance after SCF to ensure TS didn't go out of the basin
  This is basically to catch heavy losses of electrons/protons which results
  in incorrect results.

* TS.Elecs.DM.Init now defaults to diagon, while this yields worse results
  it has little influence for correctly setup systems and it will generally
  make SCF a little easier.

* Much more efficient dq implementation for fixing charge fluctuations

* Fixed TBtrans calculations with spin-flags (i.e. TB-only)

* Allow E-field for 1-electrode calculations, this allows 1-electrode
   capacitor setups

* Fixed bug in voltage potential which was severe for capacitor like setups.
   It also changes regular TS runs.

* Added TSFA[C] files which show forces on atoms in the device region

* Added TS-only energies, that is energies calculated from DM/H are
   now only using the updated elements of the respective regions.

* Fixed a bug in TBtrans when piping input

* Much faster Bloch expansions in transiesta/tbtrans calculations
   Using tiling will now result in extremely fast self-energy calculations

* Allowed complex contours in tbtrans calculations via negative eta values
   and user defined energy contours, the precision of outputted contours
   is also increased to the maximum field width.

* Changed Eta defaults in electrodes

* Enabled usage of external truly infinite electrodes (real-space self-energies)

* Fixed contour file output which is now usable for various post-processing
   utilities.

* Fixed tbtrans command line arguments, lp:1829974

* Added delta-Ef to the electrode block to specify off-set in electrodes.
   Mainly useful for semi-conducting electrodes.

* Added separate energies for NEGF calculation, now the TranSiesta energies
   are more divided and should be more comperable since they are calculated
   on the updated sub-set.


## 4.1-b4 (2018-11-07)  Bugfix beta release

### Changes to default operational parameters

* MeshCutoff has been increased to 300 Ry (from 100)
* MaxSCFIterations has been increased to 1000 (from 50)
* SCFMustConverge is now default true (from false)

### Changes

* Added developer documentation found in Docs/developer
    Both ford (preferred) and Doxygen may be used

* Generally increased precision in many output files

* Lots of fixes and updates for the Lua/flook interaction

* Auxiliary supercell handling when reading DM matrices:
    Siesta can now read and convert nearly *any* DM matrix and make it
    match the used sparse pattern.

* Fixed minor inconsistencies when handling Bessel basis

* Updated all diagonalization routines
    - ELPA and MRRR for k-point sampling.
    - Less memory usage

* Fixed bug on reading *.ion* files (lp:1751723)

* Updated internal integration table sizes (slightly increased precision)

* PDOS files now also contain the fermi-level such that tools may easily
    align the energy to 0.

* Added more digits to dDmax which may be relevant when performing
    Spin-Orbit/Non-Collinear calculations.

* Fixed bug related to writing out non-collinear spin eigenvalues,
    and also for spin-orbit. (lp:1708634)

* Fixed parallel PDOS calculations of non-colinear and spin-orbit.
    (lp:1718162)

* Added calculated charges to the Lua interface (check the charges
    while running).

* Fixed lots of compilation issues related to the utilities
    (lp:1712317, lp:1712319, lp:1711850)

* Fix for reading a ghost basis (lp:1736455, lp:1738425)

* Fix when fdf-input lines are too long. Instead of discarding the
    remaining line, fdf now "dies" to inform users of possible erroneous
    input. (lp:1728281)

* Fixed Monkhorst-Pack displacements when the displacement was larger
    than 1 (lp:1721479)

* Fix for possible heap allocated arrays (Intel compilers) (lp:1704370)

* Ensured many files to be closed properly at the end of the runs.

* Added basic compiler information to the siesta/transiesta/tbtrans
    header (compiler output)

* Performing SOC calculations does not not require all species
    to have SOC contributions.

#### TranSiesta / TBtrans

* Disk-space reduction when mixing non-periodic and periodic electrodes

* Now tiling is also enabled for Bloch expansions. This is actually faster
   than repetitions, so users should prefer tiling the electrodes

* TranSiesta is now intrinsic to the Siesta executable. An
   electrode should now be calculated using 'siesta --electrode'
   The TranSiesta executable still exists but is nothing but 'siesta --electrode'

* Many bug-fixes related to pivoting tables; this should only
   change the effective BTD matrices, and has no relevance to the
   accuracy of the calculations

* Huge performance increase in TBtrans in many different parts of the code

* Bug-fix for out-of-core calculations for spin-polarized TBT.Spin 2 calculations

* Fixed the default 2-terminal Poisson ramp. The ramp is now
   defaulted to be positioned in the central region.
     TS.Poisson ramp-central

* Small memory reduction by de-allocating unused siesta memory when
   entering transiesta.

* Fixed the box Poisson for N-electrode calculations when using
   skewed electrodes. Thanks to Azar Ostovan and Thomas Frederiksen.

* Fixed tbtrans setup for bias-window-only calculations. Now the contours
   are correctly interpreted.

* Fixed tbtrans AVCEIG output.

* Change TBtrans DOS output such that there is no normalization

* Enabled tbtrans 1-orbital calculations in the BTD matrices.

* Fixed sign-convention changes in orbital-currents. Now they are
   checked and works together with sisl (>0.9).

* Allowed external GF files for the self-energies. This is mainly beneficial
   for TBtrans as we can add external electrodes *anywhere* in the device.
   Say Buttiker-probes.

* Bugfix when the left electrode was set to -|V|/2 (the default |V|/2 is
   unaffected).

* Added much more output to the TBT*.nc files; electrode information is now
   complete, and also the BTD matrices are written.

* Enabled tbtrans -fdf TBT.Analyze which runs all pivoting schemes, this
   may be very beneficial to run with tbtrans before performing calculations.
   Choosing the correct pivoting scheme can be really important!

* Enabled output file on tbtrans command line:
     tbtrans --out TBT.out RUN.fdf
   is (more or less) equivalent to:
     tbtrans RUN.fdf > TBT.out

* Made Fermi charge correction more aggressive for faster convergence.

* TBtrans can now calculate DM, COOP and COHP curves. They are calculated
   in the supercell picture and can thus be analyzed cross-boundary as well.
   They are calculated both from the Green function and the spectral function.
   The coming >0.9.3 release of sisl will enable this analysis.

* Fixed TBtrans DOS (Green) calculations when performing k-point calculations. There
   can be small differences when comparing Green function DOS between this version
   and prior versions. The bug is only present when time-reversal-symmetry is applied.

## 4.1-b3  (2017-07-03) Bugfix beta release

* Manual greatly overhauled and updated in various parts

* Fixed DOS and PDOS for non-colinear and spin-orbit

* Fixed bug when printing initial spin-configuration

* Enabled restarting calculations with different spin-configurations,
    i.e. one may go from an unpolarized calculation to a polarized, or
    from a polarized to an unpolarized (also non-colinear and spin-orbit).

* Lots of bug-fixes for transiesta and tbtrans

* Bug-fix for spin-orbit coupling normalization

* Fixed minor memory leaks

* Many improvements for Lua enabled runs

* Added installation scripts of
    netcdf/hdf5/zlib/flook

* Fixes to the <>.nc file for high spin configuration >= non-colinear

## 4.1-b2  (2016-11-28) Bugfix beta release

* The configure script has been removed.  Its use was discouraged and
  would often yield erroneous arch.make files.  To circumvent any
  confusions it has been obsoleted until further notice.

* Instead of the configure script two default arch.make files now
  exist in the Obj directory (gfortran.make, intel.make) which should
  be guidelines for creating one's own arch.make file.

* Several fixes for bugs reported for the b1 release. See Docs/CHANGES

## 4.1-b1  (2016-08-31) Beta release

Please see the Manual for full details

### Backward-compatibility issues

* The mixing routines have completely changed, hence the same
    convergence path cannot be expected. This, unfortunately, makes
    comparison difficult with prior versions. However, the final
    converged system should be transferable.

* SIESTA now defaults to mixing the Hamiltonian instead of the
    density matrix. To revert to density-matrix mixing, use
    "SCF.Mix DM". The option to mix after the initial scf step is now
    'on' by default.

* SIESTA now defaults to monitoring convergence for both the
    density matrix AND the Hamiltonian. To revert to only density
    matrix convergence, use: "SCF.Converge.H false"

* A major number of fdf-flags concerning mixing
    parameters have changed to a more consistent naming scheme.
    However, all previous flags are still in effect but the newer
    flags have precedence. The previous flags are the default values
    for the newer flag-names.

* Two additional files are created (H_DMGEN and H_MIXED), which
      contain the Hamiltonian at various stages through the SCF.
      Currently they are intended for developers and may be removed in
      the final 4.1 release.  You may delete these without problems.

### New features

* LDA+U (Javier Junquera)

    * Full incorporation of the LDA+U implementation in SIESTA
    * Two different LDA+U projectors available
    * Estimate the best U according to: Cococcioni and Gironcoli in PRB, 2005

* Spin-Orbit coupling (Ramon Cuadrado)

    * On-site approximation for spin-orbit-coupling matrix elements.

* MRRR method for diagonalization (Alberto Garcia)

    * This will typically be faster than divide-and-conquer algorithm
    and may be the future default. For Gamma-point calculations.

* ELPA method for diagonalization (Alberto Garcia)

    * This provides better scalability compared to ScaLAPACK for large
    # of processors. For Gamma-point calculations.

* Added interface to the PEXSI library for calculating the density
   matrix, DOS, and LDOS (Alberto Garcia)

    * This library provides massive parallelism and better
    scalability, but should only be used for very large systems.

* SIESTA is now hybrid-parallelised (Nick R. Papior)

    * One may compile Siesta/Transiesta in serial, OpenMP, MPI, or
    MPI+OpenMP mode.

* Re-write of non-equilibrium Green function code (Nick R. Papior)

    * N>=1 terminal support in transiesta
    * Improved convergence
    * Different ways of handling charge-reductions in SCF
    * All electrodes have settings _per electrode_ for full customization
    * Greatly reduced memory usage
    * Skewed axis are enabled for further reduction of complex systems
    * Implemented MUMPS method for inversion
    * Implemented LAPACK for inversion
    * Implemented BTD method for inversion (extremely fast)
    * Fully OpenMP threaded
    * Start directly from transiesta enabled
    * Temperature gradients as non-equilibrium a possibility

* Complete rewrite of tbtrans utility (Nick R. Papior)

    * Made tbtrans a stand-alone utility for user-defined tight-binding method
    * EXTREME SCALE version (BTD-only)
      - Memory consumption _only_ dependent on "device" region
    * N>=1 electrode support
    * Region projections for transmission through "eigenstates"
    * Custom change of the self-energies
    * k -> k' transmissions
    * Interpolation of bias Hamiltonians
    * Bond-currents
    * Fully OpenMP threaded and/or MPI parallelized
    * DOS and bulk transmission for electrodes
    * Gf-DOS/spectral-DOS for device region

* Mixing routines rewritten (Nick R. Papior)

    * New mixing schemes Pulay (Guarenteed Reduction)
    * Custom mixing restart options (full user customizability)

* Added more constraints (Nick R. Papior)

    * Constraints are verbose and many new ways of using constraints exists

* NetCDF4 file format for siesta -> parallel IO (Nick R. Papior)

    * Provides a standard intrinsically grouped file for retaining
    nearly all siesta related information in one file.

* Enabled convergence control of density-, energy density matrices,
   Hamiltonian and energy.

* LUA scripting in siesta (Nick R. Papior)

    * This is an experimental feature

* Gate calculations (Nick R. Papior)

    * Charge and Hartree gate

* Utilities (Nick R. Papior)

    * All make-files are now prepared to enable parallel builds
      - this makes compilation *MUCH* faster. For example:
         make -j4
        will compile using 4 cores.
    * Grimme utility
      - easy creation of FDF block for Grimme parameters
    * SpPivot, pivot sparsity pattern using several different routines
    * TS/** different utilities related to transiesta


# 4.0.X versions


## 4.0.3 (2021- ) FUTURE Bug fix release

### Changes

* Remove some more memory leaks when using k-point sampling.

* Fix for honoring the 'gen_zval' field in psf files

* Removed some minor memory leaks in mesh-subs

* Update and add .md extension to main README

* Update Docs/REPORTING_BUGS

* Add ChangeLog.md since 4.0

* Fixed cell transpose for socket communication, fixes lp:1835196

* Make matel registry pool re-sizeable

* Fixed grid output when using cell-sampling, lp:1799991

* Fixed vibra utility precision issue, fixed lp:1816719

* Increase flexibility in the handling of pseudopotentials

* Fix some issues with polarization orbitals. Bessel orbs improvements.

* Fix C10 computation in molecularmechanics


## 4.0.2 (2018-07-19) Bug fix release

### Backward compatibility issues

* This release increases the size of the internal tables for
  two-center integrals used in some matrix element calculations. This
  means that calculations are slightly more heavy, but the accuracy is
  also superior. One can regain the *old* less accurate behaviour by
  setting Compat.Matel.NRTAB to true in the fdf input file (this is
  ONLY recommended for testing purposes).

### Changes

* Enabled ion.nc files for ghost atoms (#1738425)

* Enabled ghost atoms to read ion.nc files (#1736455)

* Forced libfdf to "die" when too long input strings are passed (#1728281)

* Monkhorst-Pack grids not properly shifted to [0;1[ when user specified large displacements (#1721479)

* Updated README content in Util directory (#1712319)

* Fixed building all utilities (#1712317)

* Fixed EIG file format for non-collinear spin (#1708634)

* Fixed possible segfault when using too large unit-cells (#1704370)

* Fixed ghost orbital energies (#1695130)

* Removed vpsa2bin and vpsb2asc codes (part of ATOM)

* Added new units in fdf (Hartree [Ha] and milli-Hartree [mHa])

* Added more tests

* Fixed a missing close when using Siesta in "master" mode

* Fixed several memory accounting errors and a couple of missing close
  statements. Also removed all deprecated PASTE calls.

* Added compiler version information to the version.F90 file

* Fixed writing the Bessel ghost atom to ion.xml/nc files

* Extended information in PDOS xml files. Now the atomic number
  as well as whether the orbital is a polarization shell is present

* Increased the default size of two-center integrals tabulation arrays from 128 to 1024.
  This is a change that results in more accurate values (of basically everything). Set
  'Compat.Matel.NRTAB true' to use previous value

* Added spin-monitoring in the SCF: for spin-polarized and non-collinear calculations
  the total spin-moment is written

* Fixed non-collinear Mulliken populations in parallel

* Updated gnubands to the new code that has been present since the 4.0 release.
  The new updated gnubands provides more functionalities

* Extended precision output in EIG, KP and PDOS files

* Fixed possible non-optimal DM initialization when unit-cell folding of orbital
  connections is performed.

* Bugfix for writing out too many wavefunction coefficients,
  see https://www.mail-archive.com/siesta-l@uam.es/msg10291.html


## 4.0.1 (2017-07-04) Bug fix release

### Changes

* Better standard compliance in code structure

* Fix bug related to SlabDipoleCorrection which couldn't be turned off (lp:1630827)

* Fix non-collinear bandstructure calculations (lp:1636100)

* Fix for 'nodes' basis generation option (lp:1625725)

* Fixes VCA mixing of pseudos (lp:1633039)

* Fixes integer energy specifications in ProjectedDensityOfStates block (lp:1657584)

* Added print-outs when GridCellSampling is used

* TBtrans_rep changed the written DOS units to 1/eV (they were in 1/Ry)

* Fix for Bader charge analysis (lp:1656273)

* Fix memory problem when memory usage close to limit (lp:1665294)

* Forced Diag.ParallelOverK to false if non-collinear spin configuration (lp:1666428)

* Updated Eig2DOS to be more like gnubands (options are the same)

* Added Geometry.Must.Converge flag to ensure the relaxation has converged

* Enabled internal walltime check to forcefully stop SIESTA after a
  certain limit.

* Updates and fixes for Util/STM/ol-STM
   - Fix wrong fftw call
   - Enabled direct reading of WFSX files

* Added new interpolation option to Util/Macroave

## 4.0.0 (2016-06-23) Feature release

This version includes the van der Waals functionals, the new
load-balancing code for real-space grid operations, a Wannier90
interface, a new Orbital-Minimization-Method solver, and other
improvements and bug fixes that have been part of the development
version for some time.

### Backward compatibility issues

Please take into account the following changes in behavior (more
details in the TECHNICAL NOTES section in Docs/release_notes.4.0)

* The grid functions (charge densities, potentials, etc) were in
single precision by default in the 3.X versions, but are in double
precision by default for post-3.X versions. The 'phi' array that holds
the values of the basis orbitals on the real-space grid is kept in
single precision. Please take this into account if you compare the
results with those of siesta-3.X runs. See the manual in both versions
for more information.

* Changes in the geometry used for the analysis of the electronic
structure, as well as in the handling of the density matrix (DM) and
hamiltonian (H). This will _slightly_ change the output of most
calculations and the detailed results of any post-processing. Keep
this in mind if you need to maintain coherency within a project.

* The default 'dynamics' option has been changed from 'verlet' to 'CG'.
  There should really be a new 'single-point' default which completely
  avoids 'siesta_move'. The old behavior can be recovered by using the
  'compat-pre-v4-dynamics' switch.

* Single-point calculations do not write .STRUCT_NEXT_ITER files, and
  the coordinates in the XV file are the current ones, unmoved.
  Extra output in siesta_options is avoided for this case.

* Electric field and dipole correction for slab calculations

  - Older versions applied an incorrect dipole correction when also
    using an electric field (old behavior may be recovered by forcing
    SlabDipoleCorrection to .false.)

  - Older versions over-estimated the energy contribution from the
    dipole correction by a factor of 2 (old behavior cannot be
    recovered).

### New Features

Please see the relevant section of the manual for more information.

* New SiestaXC library for exchange and correlation, implementing several
  van der Waals functionals and some newer GGA ones.

* New code to improve the load-balancing of the operations in the
  real-space grid when running in parallel.

* A new interface to the Wannier90 code for the generation of
  maximally localized wannier functions.

* New electronic-structure solver implementing the the Orbital-Minimization-Method

* New mixing options, including hamiltonian and charge-density mixing

* Charge-confinement and "filteret" basis-set generation options.

* Improved MPI version of the Siesta-as-subroutine code for executing
  independent calculations.

* New code for 'server' operation via sockets and i-PI interface

* New JobList utility to organize and run multiple jobs.

* Timer with call-tree awareness.

* Performance and usability enhancements in TranSiesta.

* Enhancements to the restart capabilities in molecular-dynamics runs.

* New options for wave-function output. New WFSX format.

* 'Fatbands' analysis.

* Enhancements to the COOP/COHP analyzer

* Hirshfeld and Voronoi charges. Bader analysis output

* Force-convergence diagnostics.

* New HSX file format for Hamiltonian/Overlap files

* Calculation of the vacuum level for non-bulk systems.

* Updates to other analysis tools in Util/

* Replacement of license-encumbered routines by new ones.

* Enhancements to the manual and build system.

### Other changes

* A number of bugs have been fixed, and there have been
numerous cosmetic changes to the output and the code itself. See Docs/CHANGES
for a full list.

## 3.X 2.X 1.X 0.X  Older releases

See the files in `Docs/older_release_notes`

<!--
Local Variables:
mode: fundamental
fill-column: 70
indent-tabs-mode: nil
coding: utf-8
End:
-->
