
module qmmm_version_info

implicit none


integer, dimension(3), save  :: num_version = (/0,0,0/)
character(len=*), parameter :: version_str =  &
"SIESTA_VERSION"
character(len=*), parameter :: qmmm_version_str =  &
"QMMM_DRIVER_VERSION"
character(len=*), parameter :: siesta_arch= &
"SIESTA_ARCH"
character(len=*), parameter :: fflags= &
"FFLAGS"

private
public :: num_version, version_str, qmmm_version_str
public :: siesta_arch, fflags

end module qmmm_version_info
!================================================================

subroutine qmmm_prversion

! Simple routine to print the version string. Could be extended to
! provide more information, if needed.

! Use free format in file to make more room for long option strings...

use qmmm_version_info
implicit none

write(6,'(2a)') "QMMM-driver Version: ", trim(qmmm_version_str)
write(6,'(2a)') "Siesta Library Version: ", trim(version_str)
write(6,'(2a)') 'Architecture  : ', siesta_arch
write(6,'(2a)') 'Compiler flags: ', fflags

#ifdef MPI
write(6,'(a)') 'PARALLEL version'
#else
write(6,'(a)') 'SERIAL version'
#endif

#ifdef TRANSIESTA
write(6,'(a)') 'TRANSIESTA support'
#endif
#ifdef CDF
write(6,'(a)') 'NetCDF support'
#endif

end subroutine qmmm_prversion
!----------------------------------------------------------

subroutine get_qmmm_version(v)
  use qmmm_version_info
  implicit none
  integer, intent(out)  :: v(3)
  v = num_version
end subroutine get_qmmm_version

