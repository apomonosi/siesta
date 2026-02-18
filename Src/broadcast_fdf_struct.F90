!
!   Broadcast complete fdf structure
!
#include "mpi_macros.f"
    SUBROUTINE broadcast_fdf_struct(reading_node,comm)
#ifndef MPI
      implicit none

      integer, intent(in)       :: reading_node         ! Node which contains the struct
      MPI_COMM_TYPE, intent(in)       :: comm
#else
      USE_MPI

      use fdf, only: fdf_serialize_struct, fdf_recreate_struct
      use fdf, only: fdf_set_started
      
      implicit none

      integer, intent(in)       :: reading_node         ! Node which contains the struct
      MPI_COMM_TYPE, intent(in)       :: comm

      character(len=1), allocatable    :: bufferFDF(:)
      integer                   :: ierr, nchars, rank

      call MPI_Comm_Rank( Comm, rank, ierr )
      
      if (rank == reading_node) then
         call fdf_serialize_struct(bufferFDF)
         nchars = size(bufferFDF)
      endif

      call MPI_Bcast(nchars, 1,                                 &
                     MPI_INTEGER, reading_node, comm, ierr)
      if (ierr .ne. MPI_SUCCESS) then
        call die("Error broadcasting size of fdf struct")
      endif

      if (rank /= reading_node) then
         ALLOCATE(bufferFDF(nchars), stat=ierr)
         if (ierr .ne. 0) then
            call die("Error allocating buffer in non-root node")
         endif
      endif

      call MPI_Bcast(bufferFDF, nchars,              &
                     MPI_CHARACTER, reading_node, comm, ierr)
      if (ierr .ne. MPI_SUCCESS) then
        call die("Error Broadcasting bufferFDF")
      endif

      if (rank /= reading_node) then
         call fdf_recreate_struct(bufferFDF)
         call fdf_set_started(.true.)
      endif

      DEALLOCATE(bufferFDF)
#endif
    END SUBROUTINE broadcast_fdf_struct
