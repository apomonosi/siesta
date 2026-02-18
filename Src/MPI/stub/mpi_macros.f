! All these are *fake* variables and should never do anything.
! They are thus considered *fault on use* and should result in
! some kind of compiler error upon usage. We are not
! there yet.

#define MPI_COMM_TYPE           INTEGER
#define MPI_REQUEST_TYPE        INTEGER
#define MPI_STATUS_TYPE         INTEGER, DIMENSION(1)
#define MPI_STATUS_TYPE_2D      INTEGER, DIMENSION(:,:), allocatable
#define MPI_GROUP_TYPE          INTEGER
#define MPI_DATA_TYPE           INTEGER
#define MPI_COMM_ID(X)          X
#define MPI_COMM_RESET          -1

#define USE_MPI
#define USE_MPI_ONLY_STATUS     use mpi_siesta, only: MPI_STATUS_SIZE
#define USE_MPI_ONLY_COMM
#define USE_MPI_ONLY_REQUEST
#define USE_MPI_ONLY_NULL
#define USE_MPI_ONLY_GROUP

#define ALLOC_MPI_STATUS_2D(x,n) allocate(x(1,n))
#define MPI_STATUS_SOURCE(X)    X(1)

