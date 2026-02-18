module unfold_handlers_m

  public :: unfold_set_handlers
  
CONTAINS

subroutine unfold_set_handlers()

  use sys, only: set_die_handler, set_bye_handler
  use sys, only: set_message_handler

!  use alloc_handlers_m, only: alloc_memory_event
!  use alloc, only: set_alloc_event_handler

  call set_die_handler(die)
  call set_bye_handler(bye)
  call set_message_handler(message)

!  call set_alloc_event_handler(alloc_memory_event)
  
end subroutine unfold_set_handlers

! -- Handlers

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
     flush(6)
     flush(0)
  endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call abort()
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
         flush(6)
         flush(0)
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
         flush(6)
      endif

#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
    end subroutine bye

! auxiliary routine to provide a traceback
  subroutine abort()
    interface
       subroutine c_abort() bind(C,name="abort")
       end subroutine c_abort
    end interface

    call c_abort()
  end subroutine abort

!----------
end module unfold_handlers_m
