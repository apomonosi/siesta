## Simple "Siesta as subroutine/server examples"

The files here (and in Src) implement three different kinds of functionality:

* Use of a MPI-based API  (using fsiesta_mpi)
* Communication with an external Siesta process using pipes
* Communication with an external Siesta process using sockets

The descriptions and documentation are somewhat outdated, mostly for the pipes and sockets cases.
A proper re-design of the examples hierarchy will be done shortly.

### MPI-based API

See the example in Src/phonons.F90

### Pipes and sockets-based API

Examples in Src using this mode of operation are compiled by CMake. The `test.sh` script in this directory
contains some hints to run them properly.


