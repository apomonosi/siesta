! ---
! Copyright (C) 1996-2021       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module version_info

implicit none

character(len=*), parameter :: version_str =  &
"TBTRANS_VERSION"
character(len=*), parameter :: compiler_version = &
"COMPILER_VERSION"
character(len=*), parameter :: siesta_arch= &
"SIESTA_ARCH"
character(len=*), parameter :: fflags= &
"FFLAGS"
character(len=*), parameter :: fppflags= &
"FPPFLAGS"
character(len=*), parameter :: libs= &
"LIBS"

private
public :: version_str
public :: siesta_arch, fflags, fppflags, libs
public :: compiler_version

public :: prversion

!================================================================

CONTAINS

subroutine prversion

!$ use omp_lib, only: openmp_version

implicit none

logical :: has_parallel(2)

#ifdef TBT_PHONON
write(6,'(2a)') "PHtrans Version : ", trim(adjustl(version_str))
#else
write(6,'(2a)') "TBtrans Version : ", trim(adjustl(version_str))
#endif
write(6,'(2a)') 'Architecture    : ', trim(adjustl(siesta_arch))
write(6,'(2a)') 'Compiler version: ', trim(adjustl(compiler_version))
write(6,'(2a)') 'Compiler flags  : ', trim(adjustl(fflags))
write(6,'(2a)') 'PP flags        : ', trim(adjustl(fppflags))
write(6,'(2a)') 'Libraries       : ', trim(adjustl(libs))

has_parallel(:) = .false.
write(6,'(a)',ADVANCE='NO') 'Parallelisations: '

! Check for MPI
#ifdef MPI
has_parallel(1) = .true.
write(6,'(a)',ADVANCE='NO') 'MPI'
#endif

! Check for OpenMP
!$ if (has_parallel(1)) write(6,'(a)', ADVANCE='NO') ', '
!$ write(6,'(a)',ADVANCE='NO') 'OpenMP'
!$ has_parallel(2) = .true.

! Complete parallel line
if ( any(has_parallel) ) then
  write(6,'(a)') ''
else
  write(6,'(a)') 'none'
end if

! Simply write out the version as given by the library
!$ write(6,'(a,i0)') '* OpenMP version: ', openmp_version

#ifdef USE_GEMM3M
write(6,'(a)') 'GEMM3M support'
#endif
#ifdef CDF
write(6,'(a)') 'NetCDF support'
#endif
#ifdef NCDF_4
write(6,'(a)') 'NetCDF-4 support'
#ifdef NCDF_PARALLEL
write(6,'(a)') 'NetCDF-4 MPI-IO support'
#endif
#endif
#if defined(ON_DOMAIN_DECOMP) || defined(SIESTA__METIS)
write(6,'(a)') 'METIS ordering support'
#endif

end subroutine prversion

end module version_info
