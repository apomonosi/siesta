!     
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!    

! Fully created by Nick Papior Andersen to conform with the io_s
! library.
#include "mpi_macros.f"

module m_ts_iodm
  
  use precision, only: dp
  use parallel, only : Node
  
  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData2D

  use m_os, only : file_exist
  use io_sparse_m, only: io_read, io_write

  implicit none
  
  private
  public :: write_ts_dm, read_ts_dm

contains
  
  subroutine read_ts_dm( file, dit, nsc, DM, EDM, Ef, found, Bcast)

#ifdef MPI
  use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(inout) :: nsc(3)
    type(dSpData2D), intent(inout) :: DM, EDM
    real(dp), intent(inout) :: Ef
    logical, intent(out) :: found
    logical, intent(in), optional :: Bcast

! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity) :: sp
    character(len=256) :: fn
    logical :: lBcast
    integer :: iu, five(5), no_u, nspin, ierr
    integer, allocatable, target :: gncol(:)
#ifdef MPI
    MPI_COMM_TYPE :: comm
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    found = file_exist(file, Bcast = Bcast )

    if ( .not. found ) then

       ! Clean-output
       call delete(DM)
       call delete(EDM)

       return

    end if

    ! Make name for the readed TSDE
    fn = 'IO-TSDE: '//trim(file)

    if ( Node == 0 ) then
       call io_assign(iu)
       open( iu, file=file, form='unformatted', status='old' )
       rewind(iu)
       read(iu,iostat=ierr) five
       if ( ierr /= 0 ) then
         rewind(iu)
         read(iu) five(1), five(2)
         five(3:5) = 0
       end if
    end if
    
#ifdef MPI
    if ( lBcast ) then
      call MPI_Bcast(five,5,MPI_integer,0,MPI_Comm_World,MPIerror)
    else
      comm = dist_comm(dit)
      call MPI_Bcast(five,5,MPI_integer,0,comm,MPIerror)
    end if
#endif

    no_u = five(1)
    nspin = five(2)
    nsc(1) = five(3)
    nsc(2) = five(4)
    nsc(3) = five(5)

    allocate(gncol(no_u))
    ! When reading the sparsity pattern we will forcefull read the gncol from Node == 0.
    ! And a negative value will force it to be distributed upon the read_d2D call.
    ! So ensure it is not negative.
    gncol(1) = 0

    ! Read in the sparsity pattern (distributed)
    call io_read(iu, no_u, sp, trim(fn), gncol=gncol, dit=dit, Bcast=Bcast)
    call io_read(iu, sp, DM, nspin, trim(fn), gncol=gncol, dit=dit, Bcast=Bcast)

    ! Read EDM
    call io_read(iu, sp, EDM, nspin, trim(fn), gncol=gncol, dit=dit, Bcast=Bcast)

    ! Clean-up
    call delete(sp)

    deallocate(gncol)

    ! Read Ef and close
    if ( Node == 0 ) then

       read(iu) Ef

       call io_close(iu)

    end if

#ifdef MPI
    if ( lBcast ) then
      call MPI_BCast(Ef,1,MPI_Double_Precision, &
          0,MPI_Comm_World, MPIerror)
    else
      comm = dist_comm(dit)
      call MPI_Bcast(Ef,1,MPI_Double_Precision, &
          0,comm,MPIerror)
    end if
#endif

  end subroutine read_ts_dm

  subroutine write_ts_dm(file, nsc, DM, EDM, Ef )
    
! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    integer, intent(in) :: nsc(3)
    type(dSpData2D), intent(inout) :: DM, EDM
    real(dp), intent(in) :: Ef
    
! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    integer, allocatable, target :: gncol(:)
    integer :: no_u, nspin
    integer :: iu

    external :: io_assign, io_close

    ! Gather sparse pattern
    dit => dist(DM)
    sp  => spar(DM)
    call attach(sp,nrows_g=no_u)
    ! Retrieve number of spin-components
    nspin = size(DM, 2)
    
    if ( Node == 0 ) then

       ! Open file
       call io_assign( iu )
       open( iu, file=file, form='unformatted', status='unknown' )
       rewind(iu)
       
       write(iu) no_u, nspin, nsc

    end if

    allocate(gncol(no_u))
    ! Signal that the first call will b-cast the gncol (speeds up
    ! subsequent writes)
    gncol(1) = -1

    ! Write sparsity pattern...
    call io_write(iu,sp, dit=dit, gncol=gncol)

    ! Write density matrix
    call io_write(iu, DM , gncol=gncol)

    ! Write energy density matrix
    call io_write(iu, EDM, gncol=gncol)

    deallocate(gncol)

    ! Write Ef and close
    if ( Node == 0 ) then
       
       write(iu) Ef
       
       call io_close(iu)
       
    end if
    
  end subroutine write_ts_dm
  
end module m_ts_iodm
