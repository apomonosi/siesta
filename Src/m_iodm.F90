!     
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!     
     
! This module has been rewritten to conform to Nick Papior Andersens IO routines
! It has also been coded by Nick Papior Andersen
#include "mpi_macros.f"

module m_iodm

  use parallel, only : Node

  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData2D

  use m_os, only: file_exist
  use io_sparse_m, only: io_read, io_write 

  implicit none
  
  private
  public :: write_dm, read_dm

  interface write_dm
   module procedure write_dm_1
   module procedure write_dm_2
  end interface write_dm

contains
        
  subroutine read_dm( file, dit, nsc, DM, found, Bcast )

#ifdef MPI
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    ! The orbital distribution that should be attached to
    ! DM
    type(OrbitalDistribution), intent(in) :: dit
    ! The supercell information (if present, otherwise return 0s)
    integer, intent(inout) :: nsc(3)
    ! The density matrix
    type(dSpData2D), intent(inout) :: DM
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
    MPI_COMM_TYPE :: dit_comm
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    found = file_exist(file, Bcast=.true.)

    ! If the file has not been found
    if ( .not. found ) then

       ! Clean-output
       call delete(DM)
       
       return
       
    end if

    ! Make name for the readed DM
    fn = 'IO-DM: '//trim(file)

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
      dit_comm = dist_comm(dit)
      call MPI_Bcast(five,5,MPI_integer,0,dit_comm,MPIerror)
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

    ! Read DM
    call io_read(iu,sp,DM ,nspin, trim(fn), gncol=gncol, dit=dit, Bcast=Bcast)

    ! Clean-up (sp is not fully deleted, it just only resides in DM)
    call delete(sp)

    ! All have this allocated (Node == 0 have just a larger one...)
    deallocate(gncol)

    ! Close
    if ( Node == 0 ) then

       call io_close(iu)

    end if

  end subroutine read_dm
  
  subroutine write_dm_1( file, nsc, DM )
    
! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    integer, intent(in) :: nsc(3)
    type(dSpData2D), intent(inout) :: DM
    
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
    sp => spar(DM)
    call attach(sp, nrows_g=no_u)
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
    ! signal we need to globalize and distribute it
    gncol(1) = -1

    ! Write sparsity pattern...
    call io_write(iu, sp, dit=dit, gncol=gncol)

    ! Write density matrix
    call io_write(iu, DM, gncol=gncol)

    deallocate(gncol)

    ! Close
    if ( Node == 0 ) then
       call io_close(iu)
    end if
    
  end subroutine write_dm_1
  
  
  subroutine write_dm_2( file, nsc, DM1, DM2 )
    
   ! **********************
   ! * INPUT variables    *
   ! **********************
       character(len=*), intent(in) :: file
       integer, intent(in) :: nsc(3)
       type(dSpData2D), intent(inout) :: DM1
       type(dSpData2D), intent(inout) :: DM2
       
   ! ************************
   ! * LOCAL variables      *
   ! ************************
       type(Sparsity), pointer :: sp
       type(OrbitalDistribution), pointer :: dit
       integer, allocatable, target :: gncol(:)
       integer :: no_u, nspin1, nspin2
       integer :: iu
   
       external :: io_assign, io_close
       
       ! Gather sparse pattern
       dit => dist(DM1)
       sp => spar(DM1)
       call attach(sp, nrows_g=no_u)
       ! Retrieve number of spin-components
       nspin1 = size(DM1, 2)
       nspin2 = size(DM2, 2)
       
       if ( Node == 0 ) then
   
          ! Open file
          call io_assign( iu )
          open( iu, file=file, form='unformatted', status='unknown' )
          rewind(iu)
          
          write(iu) no_u, nspin1+nspin2, nsc
   
       end if
   
       allocate(gncol(no_u))
       ! signal we need to globalize and distribute it
       gncol(1) = -1
   
       ! Write sparsity pattern...
       call io_write(iu, sp, dit=dit, gncol=gncol)
   
       ! Write density matrix
       call io_write(iu, DM1, gncol=gncol)
       call io_write(iu, DM2, gncol=gncol)
   
       deallocate(gncol)
   
       ! Close
       if ( Node == 0 ) then
          call io_close(iu)
       end if
       
     end subroutine write_dm_2

end module m_iodm
