#if defined(MPI_INTERFACE_F08)
#define MPI_COMM_TYPE           TYPE(MPI_Comm)
#define MPI_REQUEST_TYPE        TYPE(MPI_Request)
#define MPI_STATUS_TYPE         TYPE(MPI_Status)
#define MPI_STATUS_TYPE_2D      TYPE(MPI_Status), dimension(:), allocatable
#define MPI_GROUP_TYPE          TYPE(MPI_Group)
#define MPI_DATA_TYPE           TYPE(MPI_Datatype)
#define MPI_COMM_ID(X)          X%mpi_val
#define MPI_COMM_RESET          MPI_Comm_Null

#define USE_MPI                 use mpi_f08
#define USE_MPI_ONLY_STATUS     use mpi_siesta, only: MPI_Status
#define USE_MPI_ONLY_COMM       use mpi_siesta, only: MPI_Comm
#define USE_MPI_ONLY_REQUEST    use mpi_siesta, only: MPI_Request
#define USE_MPI_ONLY_NULL       use mpi_siesta, only: MPI_Comm_Null
#define USE_MPI_ONLY_GROUP      use mpi_siesta, only: MPI_Group

#define ALLOC_MPI_STATUS_2D(x,n) allocate(x(n))
#define MPI_STATUS_SOURCE(X)    X%mpi_source

#else
#define MPI_COMM_TYPE           INTEGER
#define MPI_REQUEST_TYPE        INTEGER
#define MPI_STATUS_TYPE         INTEGER, DIMENSION(MPI_STATUS_SIZE)
#define MPI_STATUS_TYPE_2D      INTEGER, DIMENSION(:,:), allocatable
#define MPI_GROUP_TYPE          INTEGER
#define MPI_DATA_TYPE           INTEGER
#define MPI_COMM_ID(X)          X
#define MPI_COMM_RESET          -1

#define USE_MPI                 use mpi
#define USE_MPI_ONLY_STATUS     use mpi_siesta, only: MPI_STATUS_SIZE
#define USE_MPI_ONLY_COMM
#define USE_MPI_ONLY_REQUEST
#define USE_MPI_ONLY_NULL
#define USE_MPI_ONLY_GROUP

#define ALLOC_MPI_STATUS_2D(x,n) allocate(x(MPI_STATUS_SIZE,n))
#define MPI_STATUS_SOURCE(X)    X(MPI_SOURCE)
#endif



