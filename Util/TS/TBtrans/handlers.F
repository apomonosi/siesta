!--------------------------------------------------
! Stand-alone routine to capture error messages from
! the alloc module
!
! This functionality could be made more general,
! and use a uniform interface for all the utility
! modules developed in-house. (Let's say, call it
! 'error_report' with severity arguments, etc).
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The tagging of the message by Node.
!   - The use of 'unit 6' as output and '0' as error.
!
! Routines using this functionality should include
! the following
!
!   subroutine alloc_error_report(str,code)
!     character(len=*), intent(in) :: str
!     integer, intent(in)          :: code
!   end subroutine alloc_error_report
!
!------------------------------------------------------

      subroutine alloc_error_report(str,code)

      character(len=*), intent(in)  :: str
      integer, intent(in)  :: code

      if (code == 0) then
        call die(str)
      else
         write(0,*) trim(str)
         write(6,*) trim(str)
      endif

      end subroutine alloc_error_report
!--------------------------------------------------
! Stand-alone routine to capture memory events from
! the alloc module
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta
!
!   - The use of the memory_log module for reports
!
! Routines using this functionality should include
! the following
!
!   subroutine alloc_memory_event(str,code)
!     character(len=*), intent(in) :: str
!     integer, intent(in)          :: code
!   end subroutine alloc_memory_event
!
!------------------------------------------------------

      subroutine alloc_memory_event(bytes,name)

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      ! Do nothing

      end subroutine alloc_memory_event

! The PSML library calls a "die" routine when it encounters an
! error. This routine should take care of carrying out any needed
! cleaning and terminating the program.  As the details would vary with
! the client program, each program has to provide its own.
! 

      subroutine psml_die(str)
      character(len=*), intent(in) :: str

      write(0,"(a)") str
      STOP

      end subroutine psml_die

      subroutine die(str)

#ifdef MPI
      use mpi_siesta
#endif
      character(len=*), intent(in)  :: str

      integer :: Node
#ifdef MPI
      integer MPIerror
#endif
      
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
         call pxfflush(6)
         call pxfflush(0)
      endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine die

!--------------
      subroutine message(level,str)

#ifdef MPI
      use mpi_siesta
#endif

      ! One of INFO, WARNING, FATAL
      character(len=*), intent(in)  :: level
      character(len=*), intent(in)  :: str

      integer :: Node
#ifdef MPI
      integer MPIerror
#endif
      
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif
      if (Node .eq. 0) then
         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)
         call pxfflush(6)
         call pxfflush(0)
      endif

      end subroutine message

!--------------

      subroutine bye(str)

#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in)  :: str

#ifdef MPI
      integer rc, MPIerror
#endif
      integer :: Node
      
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if (Node.eq.0) then
         write(6,'(a)') trim(str)
         write(6,'(a)') 'Requested End of Run. Bye!!'
         call pxfflush(6)
      endif

#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye

      
