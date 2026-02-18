! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------------------------------------
module siesta_handlers_m

  public :: siesta_set_handlers
  
CONTAINS

  subroutine siesta_set_handlers()

#ifdef MPI
    use mpi_siesta, only: set_mpi_timer_handler
#endif
    use sys, only: set_die_handler, set_bye_handler
    use sys, only: set_message_handler, set_reset_message_handler

    use alloc_handlers_m, only: alloc_memory_event
    use alloc, only: set_alloc_event_handler

#ifdef SIESTA__GRIDXC_HAS_PP
    ! libgridxc handlers
    use gridxc, only: gridxc_set_timer_start_handler
    use gridxc, only: gridxc_set_timer_stop_handler
    use gridxc, only: gridxc_set_error_handler
    use gridxc, only: gridxc_set_alloc_event_handler
    use gridxc, only: gridxc_set_alloc_error_handler
#endif

#ifdef SIESTA__PSML_HAS_PP
    ! psml handlers
    use m_psml, only: ps_set_error_handler
#endif
    
    call set_die_handler(die)
    call set_bye_handler(bye)
    call set_message_handler(message)
    call set_reset_message_handler(reset_messages_file)
    call set_alloc_event_handler(alloc_memory_event)

#ifdef MPI
    call set_mpi_timer_handler(timer_mpi)
#endif

#ifdef SIESTA__GRIDXC_HAS_PP
    call gridxc_set_error_handler(die)
    call gridxc_set_alloc_error_handler(alloc_die)
    call gridxc_set_alloc_event_handler(alloc_memory_event)
    call gridxc_set_timer_start_handler(gridxc_timer_start)
    call gridxc_set_timer_stop_handler(gridxc_timer_stop)
#endif
#ifdef SIESTA__PSML_HAS_PP
    call ps_set_error_handler(die)
#endif
  
end subroutine siesta_set_handlers


! Stand-alone 'die' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The use of a Siesta-specific 'mpi_siesta' module.
!   - The need to have the pxf functionality.
!   - The use of 'unit 6' as output.
!
!------------------------------------------------------

      subroutine die(str)

      use siesta_cml
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in)  :: str
      
      integer Node
      
      external ::  io_assign, io_close
      integer  ::  lun
#ifdef MPI
      integer MPIerror
#endif

      external :: pxfabort
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif
      
! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)
      write(6,'(a,i4)') 'Stopping Program from Node: ', Node
      write(0,'(a,i4)') 'Stopping Program from Node: ', Node
!!                                       endif
      if (Node .eq. 0) then
         call io_assign( lun )
         open(lun,file="MESSAGES",status="unknown",  &
             position="append",action="write")
         write(lun,"(a)") 'FATAL: ' // trim(str)
         call io_close(lun)
         flush(6)
         flush(0)
         If (cml_p) Then
            Call cmlFinishFile(mainXML)
         Endif                  !cml_p
      endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine die

!------------------------------------------------
      ! Note extra argument here
      subroutine alloc_die(str,code)
        character(len=*), intent(in) :: str
        integer, intent(in)          :: code
        write(*,*) "gridxc alloc error: "// trim(str)
        if (code == 0) call die(str)
      end subroutine alloc_die
!------------------------------------------------
      
      subroutine message(level,str)

      use parallel, only : Node

      ! One of INFO, WARNING, FATAL
      character(len=*), intent(in)  :: level

      character(len=*), intent(in)  :: str

      external ::  io_assign, io_close
      integer  ::  lun
      
      if (Node .eq. 0) then
         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)
         call io_assign(lun)
         open(lun,file="MESSAGES",status="unknown",   &
              position="append",action="write")
         write(lun,"(a)") trim(level) // ": " // trim(str)
         call io_close(lun)
         flush(6)
         flush(0)
      endif

      end subroutine message

      subroutine reset_messages_file()
      use parallel, only : Node

      integer :: lun
      external ::  io_assign, io_close
      
      if (Node .eq. 0) then
         call io_assign(lun)
         ! Open with 'replace' to clear content
         open(lun,file="MESSAGES",status="replace",  &
              position="rewind",action="write")
         call io_close(lun)
      endif
      end subroutine reset_messages_file
      
!---------------------------------------------------------
      subroutine bye(str)

      use siesta_cml
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in)  :: str

      integer  :: Node
      external :: pxfflush
#ifdef MPI
      integer rc, MPIerror
#endif

#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if (Node.eq.0) then
         write(6,'(a)') trim(str)
         write(6,'(a)') 'Requested End of Run. Bye!!'
         flush(6)
         If (cml_p) Then
            Call cmlFinishFile(mainXML)
         Endif                  !cml_p
      endif

#ifdef MPI
      call MPI_Finalize(rc)
#endif
      stop
      
    end subroutine bye

!------------------------
  SUBROUTINE timer_mpi( name, opt )
    character(len=*), intent(in):: name
    integer,          intent(in):: opt

#ifdef MPI_TIMING
    external timer
    call timer( name, opt )
#endif

  END SUBROUTINE timer_mpi

#ifdef SIESTA__GRIDXC_HAS_PP
  !----------
  subroutine gridxc_timer_start(str)
    character(len=*), intent(in)  :: str
    call timer("gridxc@"//trim(str),1)
  end subroutine gridxc_timer_start
  !
  subroutine gridxc_timer_stop(str)
    character(len=*), intent(in)  :: str
    call timer("gridxc@"//trim(str),2)
  end subroutine gridxc_timer_stop
#endif

  end module siesta_handlers_m
