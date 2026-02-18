This library can be used as an external component for MM engines
in order to use SIESTA as a QM workhorse.

Compile and link your code of choice, adding this folder location
to the linking paths and using -lsiestaqmmm.

## The interface

The interface is done via a single subroutine:

```fortran
   subroutine get_siesta_forces_and_stress( args )
```

### INPUTS (positions and distances in Bohr, charges in e-):
  n_qm (integer): Number of atoms in the QM region.
  n_mm (integer): Number of atoms in the MM region.
  r_qm (double, size = (3 x n_qm) ): Positions of atoms in the QM region.
  r_mm (double, size = (3 x n_mm) ): Positions of atoms in the MM region.
  pc   (double, size = (n_mm) ): Classical partial charges of MM atoms.
  cell (double, size = (3 x 3) ): Unit cell vectors.
  last (logical): Indicates whether this is the last MD step.

### OUTPUTS (forces in Ry/Bohr, energy in Ry ):
  energy (double): Siesta total QM/MM energy.
  f_qm   (double, size = (3 x n_qm) ): Siesta forces over QM atoms.
  f_mm   (double, size = (3 x n_mm) ): Siesta forces contribution over MM atoms.
  stress (double, size = (3 x 3) ): Cell stress contribution from Siesta.

In the first call this library sets up a SIESTA call according to
the following options that MUST BE INCLUDED in a siestaqmmm.fdf file
within the working directory. This is in addition to the fdf file
REQUIRED for SIESTA (i.e. this library won't set up an fdf file for
the QM region, it must be provided by the user).

```bash
 ParallelCommand 'mpirun' # Can be others, such as mpiexec or srun.
 nCpus 2                  # Amount of CPUs to dedicate to SIESTA.
 runSerial F              # Run SIESTA in serial mode, ingoring the above options.
 SiestaInput  NAME        # Name of the fdf file for SIESTA.
```
The name of this specific input file, siestaqmmm.fdf, can be changed using
the set_qmmm_input_file( fname ) method included in this same library. You can find 
a sample file under the sample_driver folder.
