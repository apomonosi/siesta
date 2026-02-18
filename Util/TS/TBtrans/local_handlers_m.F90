module local_handlers_m

  public :: die
  public :: message
  public :: bye

CONTAINS
  
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

end module local_handlers_m
