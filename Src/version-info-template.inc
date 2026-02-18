!
! Copyright (C) 1996-2021       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module siesta_version_info

implicit none

character(len=*), parameter :: siesta_version_str = &
"@SIESTA_VERSION@"
character(len=*), parameter :: siesta_compiler_version = &
"@CMAKE_Fortran_COMPILER_ID@-@CMAKE_Fortran_COMPILER_VERSION@"

character(len=*), parameter :: siesta_arch= &
"@CMAKE_SYSTEM_PROCESSOR@"
character(len=*), parameter :: siesta_fflags= &
"@Fortran_FLAGS_CURRENT_multiline@"
character(len=*), parameter :: siesta_fppflags= &
"<not-set>"
character(len=*), parameter :: siesta_libs= "<not-set>"

private
public :: siesta_version_str
public :: siesta_arch, siesta_fflags, siesta_fppflags, siesta_libs
public :: siesta_compiler_version

public :: siesta_print_version

!================================================================

CONTAINS

subroutine siesta_print_version

!$ use omp_lib, only: openmp_version

! Simple routine to print the version string. Could be extended to
! provide more information, if needed.

! Use free format in file to make more room for long option strings...
implicit none

character(len=256) :: cmd
logical :: has_parallel(2), header_printed

cmd = " "
call get_command_argument(0, cmd)
write(6,'(2a)') 'Executable      : ', trim(cmd)
write(6,'(2a)') 'Version         : ', trim(adjustl(siesta_version_str))
write(6,'(2a)') 'Architecture    : ', trim(adjustl(siesta_arch))
write(6,'(2a)') 'Compiler version: ', trim(adjustl(siesta_compiler_version))
write(6,'(2a)') 'Compiler flags  : ', trim(adjustl(siesta_fflags))
!write(6,'(2a)') 'PP flags        : ', trim(adjustl(siesta_fppflags))
!write(6,'(2a)') 'Libraries       : ', trim(adjustl(siesta_libs))

header_printed = .false.
call process_env_var("SIESTA_PS_PATH",header_printed)
call process_env_var("SIESTA_ELPA_GPU_STRING",header_printed)

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
#ifdef SIESTA__FLOOK
write(6,'(a)') 'Lua support'
#endif
#ifdef SIESTA__PEXSI
write(6,'(a)') 'Native PEXSI support'
#endif
#ifdef SIESTA__ELPA

#ifdef ELPA_HAS_GPU_SUPPORT
   write(6,'(a,a)') " --- ELPA GPU support: ", trim("@ELPA_GPU_STRING@")
#endif
#endif
#ifdef SIESTA__ELSI
write(6,'(a)') 'ELSI support. Solvers:'

#ifdef ELSI_HAS_EXTERNAL_ELPA
#include "elpa/elpa_version.h"

#ifdef ELPA_HAS_GPU_SUPPORT
write(6,'(a)') '   ELPA (external, GPU support: ' // trim("@ELPA_GPU_STRING@") // ')'
#else
write(6,'(a)') '   ELPA (external) '
#endif
write(6,'(a,i0)') '     ELPA API version (from .h file):', ELPA_API_VERSION
write(6,'(a,i0)') '     ELPA API version (from CMake var):', @ELPA_API_VERSION@

#else

#ifdef ELSI_HAS_CUDA
write(6,'(a)') '   ELPA (internal, GPU support: ' // trim("@ELSI_ELPA_GPU_STRING@") // ')'
#else
write(6,'(a)') '   ELPA (internal) '
#endif

#endif

write(6,'(a)') '   NTPoly'
write(6,'(a)') '   OMM'
#ifdef ELSI_HAS_PEXSI
write(6,'(a)') '   PEXSI'
#endif
#ifdef ELSI_HAS_MAGMA
write(6,'(a)') '   MAGMA'
#endif
#ifdef ELSI_HAS_EIGENEXA
write(6,'(a)') '   EIGENEXA'
#endif
#ifdef ELSI_HAS_SIPS
write(6,'(a)') '   SIPS'
#endif
#endif
#ifdef SIESTA__DFTD3
write(6,'(a)') 'DFT-D3 support'
#endif
#ifdef SIESTA__WANNIER90
write(6,'(a)') 'Wannier90 wrapper support'
#endif

end subroutine siesta_print_version

subroutine process_env_var(name,header_printed)
  character(len=*), intent(in) :: name  ! env var name
  logical, intent(inout)        :: header_printed
  
  integer :: status, len
  character(len=:), allocatable :: ps_path

  call get_environment_variable(trim(name),status=status,length=len)
  if ((status == 0) .and. (len > 0)) then
    allocate(character(len=len) :: ps_path)
    if (.not. header_printed) then
      write(6,'(a)') 'Environment variables: '
      header_printed = .true.
    endif
    call get_environment_variable(trim(name),value=ps_path)
    write(*,'(2x,a)') trim(name) // ': '//trim(ps_path)
    deallocate(ps_path)
  endif

end subroutine process_env_var

end module
