! Module for easy reading of sparsity patterns and 
! data distributed in sparsity patterns.

! Ideally all IO routines should utilize these routines

! Fully implemented by Nick Papior Andersen

! SIESTA io-module
module io_sparse_m

  use class_OrbitalDistribution
  use class_Sparsity
  use precision, only : sp, dp
#ifdef MPI
  use mpi_siesta, only : MPI_AllReduce, MPI_Sum, MPI_Max
  use mpi_siesta, only : MPI_Bcast, MPI_IBcast
  use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self
  use mpi_siesta, only : MPI_Send, MPI_Recv
  use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
  use mpi_siesta, only : MPI_Success, MPI_Status_Size
  use mpi_siesta, only : MPI_STATUSES_IGNORE
#endif

  implicit none 

  private

  public :: Node_Sp_gncol
  public :: io_read_Sp
  public :: io_write_Sp


  public :: io_read_d1D, io_read_d2D
  public :: io_write_d1D, io_write_d2D
  public :: io_write_r1D, io_write_r2D

  public :: io_read
  interface io_read
    module procedure io_read_Sp
    module procedure io_read_d1D
    module procedure io_read_d2D
  end interface io_read

  public :: io_write
  interface io_write
    module procedure io_write_Sp
    module procedure io_write_d1D
    module procedure io_write_d2D
  end interface io_write

  public :: io_write_r
  interface io_write_r
    module procedure io_write_Sp
    module procedure io_write_r1D
    module procedure io_write_r2D
  end interface io_write_r

  ! The counting functions
  public :: count_blocks
  public :: count_consecutive, count_consecutive_sum
  public :: max_consecutive, max_consecutive_sum

contains

  ! Returns a consecutive number of contributions
  ! starting from the specified index
  function count_blocks(dit,n_t) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t
    integer :: n
    ! Local variables
    integer :: cur_node, i

    n = 1
    cur_node = node_handling_element(dit,1)
    do i = 2 , n_t
      ! if the idx is not present, just return
      if ( cur_node /= node_handling_element(dit,i) ) then
        n = n + 1
        cur_node = node_handling_element(dit,i)
      end if
    end do

  end function count_blocks

  ! Returns a consecutive number of contributions
  ! starting from the specified index
  function count_consecutive(dit,n_t,idx) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t, idx
    integer :: n
    ! Local variables
    integer :: cur_node, i

    n = 1
    cur_node = node_handling_element(dit,idx)
    do i = idx + 1 , n_t
      ! if the idx is not present, just return
      if ( cur_node /= node_handling_element(dit,i) ) return
      n = n + 1
    end do

  end function count_consecutive

  function count_consecutive_sum(dit,n_t,nidx,idx) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t, nidx(n_t), idx
    integer :: n

    n = count_consecutive(dit,n_t,idx)
    n = sum(nidx(idx:idx-1+n))

  end function count_consecutive_sum

  function max_consecutive(dit,n_t) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t
    integer :: n
    ! Local variables
    integer :: i, cur_n

    i = 1
    n = 0
    do while ( i <= n_t )

      ! Count number of consecutive numbers
      cur_n = count_consecutive(dit,n_t,i)
      if ( cur_n > n ) then
        n = cur_n
      end if

      ! Step counter
      i = i + cur_n

    end do

  end function max_consecutive

  function max_consecutive_sum(dit,n_t,nidx) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t, nidx(n_t)
    integer :: n
    ! Local variables
    integer :: i, cur_n, tmp

    i = 1
    n = 0
    do while ( i <= n_t )

      ! Count number of consecutive numbers
      cur_n = count_consecutive(dit,n_t,i)
      tmp = sum(nidx(i:i-1+cur_n))
      if ( tmp > n ) then
        n = tmp
      end if

      ! Step counter
      i = i + cur_n

    end do

  end function max_consecutive_sum

  !> Distribute a global n-col variable
  !!
  !! This routine enables a local-node collection of the n-col
  !! variable, or a B-cast.
  !!
  !! if lNode < 0, it means a broad cast, otherwise the corresponding
  !! node will retain the elements upon return.
  subroutine Node_Sp_gncol(toNode,sp,dit,no,gncol)
    integer, intent(in) :: toNode
    type(Sparsity), intent(inout) :: sp
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: no
    integer, intent(inout) :: gncol(no)

    integer, pointer :: ncol(:)
    integer :: lNode, comm

    ! grab ncol
    call attach(sp,n_col=ncol)

    ! Get current node in distribution
    lNode = dist_node(dit)
    comm = dist_comm(dit)

#ifdef MPI
    if ( toNode < 0 ) then
      ! Here everybody collects, so no need to handle it
      call broadcast()
    else
      ! Signal that no further action is needed
      gncol(1) = 0
      call collect()
    end if
#else
    gncol(:) = ncol(:)
#endif

#ifdef MPI
  contains

    subroutine collect()
      integer :: gio, n, io, nb
      integer, allocatable :: ibuf(:)
      integer :: BNode, MPIerror
    
      nb = count_blocks(dit,no)
      allocate(ibuf(nb))

      gio = 1
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then
          io = index_global_to_local(dit,gio,lNode)
          if ( lNode == toNode ) then
            gncol(gio:gio-1+n) = ncol(io:io-1+n)
          else
            nb = nb + 1
            call MPI_ISend(ncol(io), n, MPI_Integer, &
                toNode, gio, comm, ibuf(nb), MPIerror)
          end if
        else if ( lNode == toNode ) then
          nb = nb + 1
          call MPI_IRecv(gncol(gio), n, MPI_Integer, &
              BNode, gio, comm, ibuf(nb), MPIerror)
        end if
        gio = gio + n
      end do

      if ( nb > 0 ) then
        call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
      end if
      deallocate(ibuf)

    end subroutine collect

    subroutine broadcast()
      integer :: gio, n, io, nb
      integer, allocatable :: ibuf(:)
      integer :: BNode, MPIerror

      nb = count_blocks(dit,no)
      allocate(ibuf(nb))

      gio = 1
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)
      
        if ( BNode == lNode ) then
          ! copy over to global array
          io = index_global_to_local(dit,gio,lNode)
          gncol(gio:gio-1+n) = ncol(io:io-1+n)
        end if

        nb = nb + 1
        ! The MPI standard says that as long as collective
        ! are called in the same order, there is no ambiguity
        ! in the recieving end (i.e. tags are not necessary).
        call MPI_IBcast(gncol(gio), n, MPI_Integer, &
            BNode, comm, ibuf(nb), MPIerror)
        gio = gio + n
        
      end do

      if ( nb > 0 ) then
        call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
      end if
      deallocate(ibuf)

    end subroutine broadcast
#endif
    
  end subroutine Node_Sp_gncol


  ! Parse options as provided to super routines
  ! This converts the dit and Bcast arguments into:
  !   - lNode
  !     dist_node(dit), otherwise Node
  !   - comm
  !     dist_comm(dit), otherwise MPI_COMM_WORLD
  subroutine parallel_default_args(comm, lNode, dit, ldit, Bcast, lBcast)
    integer, intent(out) :: comm, lNode
    type(OrbitalDistribution), intent(in), optional :: dit
    logical, intent(out), optional :: ldit
    logical, intent(in), optional :: Bcast
    logical, intent(out), optional :: lBcast

    integer :: ierr

    ! Default node
    lNode = 0

#ifdef MPI
    ! Default communicator
    comm = MPI_COMM_WORLD

    ! Default to not use a distribution, nor b-casting
    ldit = present(dit)
    
    if ( present(lBcast) ) then
      lBcast = .false.
      if ( present(Bcast) ) lBcast = Bcast

      ! Do not signal a distribution
      if ( lBcast ) then
        ldit = .false.
      end if
      
    end if

    ! Now parse options
    if ( ldit ) then
      comm = dist_comm(dit)
      lNode = dist_node(dit)

      ! Only use a distributed read in case there are
      ! more than 1 node available, otherwise not needed
      ! Also bcasting is not needed
      if ( dist_nodes(dit) == 1 ) then
        ldit = .false.
      end if

    else if ( lBcast ) then

      ! get node from the communicator
      call MPI_Comm_rank(comm, lNode, ierr)
      
    end if

#else
    ! Not used, or needed in serial compilations
    comm = -1
    if ( present(ldit) ) ldit = .false.
    if ( present(lBcast) ) lBcast = .false.
#endif

  end subroutine parallel_default_args


  ! Reads in a sparsity pattern at the
  ! current position in the file (iu)
  ! The sparsity pattern "sp" will be returned
  ! as populated.
  ! If dist is supplied it will distribute
  ! the sparsity pattern as supplied (this implies Bcast = .true.)
  ! Else if Bcast is true it will b-cast the sparsity 
  ! pattern fully.
  subroutine io_read_Sp(iu, no, sp, tag, dit, Bcast, gncol)

    ! File handle
    integer, intent(in) :: iu
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(no)

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, ind, nl, n, i, nb
    integer :: lNode, comm
    logical :: ldit, lBcast
    integer, pointer :: lgncol(:) => null()
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: MPIerror, BNode
#endif

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit, Bcast, lBcast)

    if ( lNode == 0 ) then
      ! Local node

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
      end if

      ! First read in number of non-zero 
      ! entries per orbital
      read(iu) lgncol

    end if

    ! local number of elements
    nl = no

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit ) then

      ! Number of local elements
      nl = num_local_elements(dit,no,lNode)
      allocate(ncol(nl))

      ! allocate all requests
      nb = count_blocks(dit,no)
      allocate(ibuf(nb))

      ! distribute ncol
      gio = 1
      nb = 0
      do while ( gio <= no ) 

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then

          io = index_global_to_local(dit,gio,lNode)

          if ( lNode == 0 ) then
            ncol(io:io-1+n) = lgncol(gio:gio-1+n)
          else
            nb = nb + 1
            call MPI_IRecv(ncol(io), n, MPI_Integer, &
                0, gio, comm, ibuf(nb), MPIerror)
          end if

        else if ( lNode == 0 ) then
          nb = nb + 1
          call MPI_ISend(lgncol(gio), n, MPI_Integer, &
              BNode, gio, comm, ibuf(nb), MPIerror)

        end if

        gio = gio + n

      end do

      if ( nb > 0 ) then
        call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
      end if
      deallocate(ibuf)


    else if ( lBcast ) then

      ! Everything should be b-casted
      if ( lNode == 0 ) then
        ncol => lgncol
      else
        allocate(ncol(nl))
      end if

      ! Bcast everything
      call MPI_Bcast(ncol(1), nl, MPI_Integer, 0, comm, MPIError)


    else if ( lNode == 0 ) then

      ncol => lgncol

    end if
#else
    ! Point to the buffer
    ncol => lgncol
#endif


    !>
    !> At this point we have distributed ncol
    !> Now distribute the index pointer
    !> 
    
    ! Allocate (local) pointer
    allocate(l_ptr(nl))

    l_ptr(1) = 0
    do io = 2 , nl
      l_ptr(io) = l_ptr(io-1) + ncol(io-1)
    end do

    ! Number of local non-zero elements
    ! (also works for any bcast methods)
    n_nzs = l_ptr(nl) + ncol(nl)

    ! Allocate space for column indices
    allocate(l_col(n_nzs))

#ifdef MPI
    if ( ldit ) then

      ! We have a distributed read
      if ( lNode == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(ibuf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Read in columns and distribute
      gio = 1
      ind = 0
      nb = 0
      do while ( gio <= no ) 

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,lNode)

          if ( lNode == 0 ) then

            do i = io , io - 1 + n
              read(iu) l_col(ind+1:ind+ncol(i))
              ind = ind + ncol(i)
            end do

          else

            ! count the number of received entities
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_IRecv(l_col(ind+1), i, MPI_Integer, &
                0, gio, comm, ibuf(nb), MPIerror)
            ind = ind + i

          end if

        else if ( lNode == 0 ) then

          i = 0
          do io = gio , gio + n - 1
            read(iu) ibuf(i+1:i+lgncol(io))
            i = i + lgncol(io)
          end do

          call MPI_Send(ibuf(1), i, MPI_Integer, &
              BNode, gio, comm, MPIerror)

        end if

        gio = gio + n

      end do

      if ( lNode == 0 ) then
        if ( .not. present(gncol) ) deallocate(lgncol)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
        end if
      end if

      ! Clean up memory
      deallocate(ibuf)

    else if ( lBcast ) then

      if ( lNode == 0 ) then

        ind = 0
        do gio = 1 , no
          read(iu) l_col(ind+1:ind+ncol(gio))
          ind = ind + ncol(gio)
        end do

      end if

      ! Bcast
      call MPI_Bcast(l_col(1), n_nzs, MPI_Integer, 0, comm, MPIError)

    else if ( lNode == 0 ) then
#endif

      ind = 0
      do io = 1 , no
        read(iu) l_col(ind+1:ind+ncol(io))
        ind = ind + ncol(io)
      end do

#ifdef MPI       
    end if
#endif

    ! Create the sparsity pattern
    call newSparsity(sp,nl,no, n_nzs, ncol, l_ptr, l_col, trim(tag))

    ! de-allocate
    deallocate(l_ptr,l_col)
    if ( ldit ) deallocate(ncol)
    if ( lBcast .and. lNode /= 0 ) deallocate(ncol)
    if ( lNode == 0 .and. .not. present(gncol) ) deallocate(lgncol)

  end subroutine io_read_Sp

  ! Writes a sparsity pattern at the
  ! current position in the file (iu)
  ! If dist is supplied it will write a distributed sparsity pattern
  subroutine io_write_Sp(iu, sp, dit, gncol)

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! distribution
    type(OrbitalDistribution), intent(in), optional :: dit
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    integer, pointer :: lgncol(:) => null()
    integer, pointer :: ncol(:), l_col(:) => null()

    integer :: lno, no, io, max_n, ind, n, i, nb
    integer :: comm, lNode
    logical :: ldit
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: gio
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    ! Get the sparsity sizes
    call attach(sp,n_col=ncol, list_col=l_col, nrows=lno,nrows_g=no)

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit)

    if ( ldit ) then

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
        if ( size(gncol) < no ) then
          call die("io_write_Sp: wrong size argument gncol")
        end if
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

#else
      call die('Error in code, non-full contained sp')
#endif

    else
      lgncol => ncol
    end if

    if ( lNode == 0 ) then

      write(iu) lgncol

    end if

#ifdef MPI
    ! Write the list_col array
    if ( ldit ) then

      nb = count_blocks(dit,no)

      ! The ionode now has the maximum retrieved array
      if ( lNode == 0 ) then
        ! Retrive the maximum number of non-zero
        ! elements in each row
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(ibuf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then

          io = index_global_to_local(dit,gio,lNode)

          if ( lNode == 0 ) then
            do i = io , io - 1 + n
              write(iu) l_col(ind+1:ind+ncol(i))
              ind = ind + ncol(i)
            end do
          else
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_ISend(l_col(ind+1), i, MPI_Integer, &
                0, gio, comm, ibuf(nb), MPIerror)
            ind = ind + i
          end if
        else if ( lNode == 0 ) then
          call MPI_Recv(ibuf(1), max_n, MPI_Integer, &
              BNode, gio, comm, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: io_write_Sp')
          i = 0
          do io = gio , gio - 1 + n
            write(iu) ibuf(i+1:i+lgncol(io))
            i = i + lgncol(io)
          end do

        end if
        gio = gio + n
      end do

      if ( .not. present(gncol) ) deallocate(lgncol)
      if ( lNode /= 0 .and. nb > 0 ) then
        call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
      end if
      deallocate(ibuf)

    else if ( lNode == 0 ) then

      ind = 0
      do io = 1 , no
        write(iu) l_col(ind+1:ind+lgncol(io))
        ind = ind + lgncol(io)
      end do

    end if

#else

    ind = 0
    do io = 1 , no
      write(iu) l_col(ind+1:ind+lgncol(io))
      ind = ind + lgncol(io)
    end do

#endif

  end subroutine io_write_Sp

  subroutine io_read_d1D(iu, sp, dSp1D, tag, dit, Bcast, gncol)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData1D), intent(inout) :: dSp1D
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: io, lno, no, ind, n_nzs, n, i, nb
    integer :: comm, lNode
    logical :: ldit, lBcast
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: max_n, gio
    integer :: MPIerror, BNode
#endif

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit, Bcast, lBcast)

    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)

    if ( ldit ) then
      call newdSpData1D(sp,dit,dSp1D,name=trim(tag))

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if
#else
      call die('Error in distribution, io_read_d1D')
#endif

    else

      ! Create the Fake distribution
#ifdef MPI
      call newDistribution(no, MPI_Comm_Self, fdit, name='Fake dist')
#else
      call newDistribution(no, -1, fdit, name='Fake dist')
#endif
      call newdSpData1D(sp, fdit, dSp1D, name=trim(tag))
      ! Clean up the distribution again
      call delete(fdit)
    end if

    ! retrieve data placement
    a => val(dSp1D)

    ! Only read distributed if we are not b-casting
    if ( ldit ) then
#ifdef MPI

      nb = count_blocks(dit,no)

      ! Allocate the maximum number of entries
      if ( lNode == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1 
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then

          io = index_global_to_local(dit,gio,lNode)

          if ( lNode == 0 ) then

            do i = io , io - 1 + n
              read(iu) a(ind+1:ind+ncol(i))
              ind = ind + ncol(i)
            end do

          else

            ! count the number of received entities
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_IRecv(a(ind+1), i, MPI_Double_Precision, &
                0, gio, comm, ibuf(nb), MPIerror)
            ind = ind + i

          end if

        else if ( lNode == 0 ) then

          i = 0
          do io = gio , gio - 1 + n
            read(iu) buf(i+1:i+lgncol(io))
            i = i + lgncol(io)
          end do

          call MPI_Send(buf(1), i, MPI_Double_Precision, &
              BNode, gio, comm, MPIerror)

        end if

        gio = gio + n

      end do

      if ( .not. present(gncol) ) deallocate(lgncol)
      if ( lNode /= 0 ) then
        if ( nb > 0 ) then
          call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
        end if
      end if
      deallocate(ibuf)
#else
      call die('Error in distribution for, io_read_d1D')
#endif
    else if ( lNode == 0 ) then

      ind = 0
      do io = 1 , no
        read(iu) a(ind+1:ind+ncol(io))
        ind = ind + ncol(io)
      end do

    end if

#ifdef MPI
    if ( lBcast ) then

      call MPI_Bcast(a(1), n_nzs, MPI_Double_Precision, &
          0, comm, MPIError)

    end if
#endif

  end subroutine io_read_d1D

  subroutine io_write_d1D(iu, dSp1D, gncol)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData1D), intent(inout) :: dSp1D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: io, lno, no, ind, n, i, nb
    integer :: comm, lNode
    logical :: ldit
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    call attach(sp,nrows=lno,nrows_g=no,n_col=ncol)
    
    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit)

    ! Retrieve data
    a => val(dSp1D)

    if ( ldit ) then

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if
#else
      call die('Error in distribution, io_write_d1D')
#endif

#ifdef MPI

      nb = count_blocks(dit,no)

      ! The ionode now has the maximum retrieved array
      if ( lNode == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Write the data...
      gio = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then

          io = index_global_to_local(dit,gio,lNode)

          if ( lNode == 0 ) then
            do i = io , io - 1 + n
              write(iu) a(ind+1:ind+ncol(i))
              ind = ind + ncol(i)
            end do
          else
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_ISend(a(ind+1), i, MPI_Double_Precision, &
                0, gio, comm, ibuf(nb), MPIerror)
            ind = ind + i
          end if
        else if ( lNode == 0 ) then
          call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
              BNode, gio, comm, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: io_write_d1D')
          i = 0
          do io = gio , gio - 1 + n
            write(iu) buf(i+1:i+lgncol(io))
            i = i + lgncol(io)
          end do

        end if
        gio = gio + n
      end do

      if ( .not. present(gncol) ) deallocate(lgncol)
      if ( lNode == 0 ) then
        deallocate(buf)
      else 
        ! Wait for the last one to not send
        ! two messages with the same tag...
        if ( nb > 0 ) then
          call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
        end if
        deallocate(ibuf)
      end if
#endif
    else if ( lNode == 0 ) then

      ind = 0
      do io = 1 , no
        write(iu) a(ind+1:ind+ncol(io))
        ind = ind + ncol(io)
      end do

    end if

  end subroutine io_write_d1D

  subroutine io_write_r1D(iu, dSp1D, gncol)

    use precision, only: psp => sp
    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData1D), intent(inout) :: dSp1D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: io, lno, no, ind, n, i, nb
    integer :: comm, lNode
    logical :: ldit
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    call attach(sp,nrows=lno,nrows_g=no,n_col=ncol)

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit)

    ! Retrieve data
    a => val(dSp1D)

    if ( ldit ) then

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if
#else
      call die('Error in distribution, io_write_r1D')
#endif

#ifdef MPI

      nb = count_blocks(dit,no)

      ! The ionode now has the maximum retrieved array
      if ( lNode == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Write the data...
      gio = 1
      ind = 0
      nb = 0 
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == lNode ) then

          io = index_global_to_local(dit,gio,lNode)

          if ( lNode == 0 ) then
            do i = io , io - 1 + n
              write(iu) real(a(ind+1:ind+ncol(i)), psp)
              ind = ind + ncol(i)
            end do
          else
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_ISend(a(ind+1), i, MPI_Double_Precision, &
                0, gio, comm, ibuf(nb), MPIerror)
            ind = ind + i
          end if
        else if ( lNode == 0 ) then
          call MPI_Recv( buf(1), max_n, MPI_Double_Precision, &
              BNode, gio, comm, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: io_write_r1D')
          i = 0
          do io = gio , gio - 1 + n
            write(iu) real(buf(i+1:i+lgncol(io)), psp)
            i = i + lgncol(io)
          end do

        end if
        gio = gio + n
      end do

      if ( .not. present(gncol) ) deallocate(lgncol)
      if ( lNode == 0 ) then
        deallocate(buf)
      else 
        ! Wait for the last one to not send
        ! two messages with the same tag...
        if ( nb > 0 ) then
          call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
        end if
        deallocate(ibuf)
      end if
#endif
    else if ( lNode == 0 ) then

      ind = 0
      do io = 1 , no
        write(iu) real(a(ind+1:ind+ncol(io)), psp)
        ind = ind + ncol(io)
      end do

    end if

  end subroutine io_write_r1D

  subroutine io_read_d2D(iu, sp, dSp2D, dim2, tag, &
      sparsity_dim, dit, Bcast, gncol)

    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData2D), intent(inout) :: dSp2D
    ! The non-sparse dimension
    integer, intent(in) :: dim2
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! This denotes the sparsity dimension (either 1 or 2) 1=default
    integer, intent(in), optional :: sparsity_dim
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: lno, no, n_nzs, sp_dim
    integer :: lNode, comm
    logical :: ldit, lBcast
#ifdef MPI
    integer :: MPIerror
#endif

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit, Bcast, lBcast)

    sp_dim = 1
    if ( present(sparsity_dim) ) sp_dim = sparsity_dim

    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,nnzs=n_nzs)

    if ( ldit ) then
      call newdSpData2D(sp,dim2,dit,dSp2D,name=trim(tag), &
          sparsity_dim=sp_dim)

#ifdef MPI
      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
      if ( lgncol(1) < 0 ) then
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if
#else
      call die('Error in distribution, io_read_d2D')
#endif

    else
      ! Create the Fake distribution
#ifdef MPI
      call newDistribution(no, MPI_Comm_Self, fdit, name='Fake dist')
#else
      call newDistribution(no, -1, fdit, name='Fake dist')
#endif
      call newdSpData2D(sp,dim2,fdit,dSp2D,name=trim(tag), &
          sparsity_dim=sp_dim)
      call delete(fdit)
    end if

    a => val(dSp2D)

    if ( sp_dim == 1 ) then
      call read_sp_dim1()
    else
      call read_sp_dim2()
    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then

      call MPI_Bcast(a(1,1), dim2*n_nzs, MPI_Double_Precision, &
          0, comm, MPIError)

    end if
#endif

  contains

    subroutine read_sp_dim1()
      integer :: io, ind, n, i, nb, s
#ifdef MPI
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: gio, max_n
      integer :: BNode
#endif

      if ( ldit ) then
#ifdef MPI

        nb = count_blocks(dit,no)

        ! Allocate maximum number of entries
        if ( lNode == 0 ) then
          max_n = max_consecutive_sum(dit,no,lgncol)
          allocate(buf(max_n))
        else
          allocate(ibuf(nb))
        end if

        do s = 1 , dim2

          gio = 1
          ind = 0
          nb = 0
          do while ( gio <= no )

            BNode = node_handling_element(dit,gio)

            ! Get number of consecutive orbitals
            n = count_consecutive(dit,no,gio)

            if ( BNode == lNode ) then

              io = index_global_to_local(dit,gio,lNode)

              if ( lNode == 0 ) then
                do i = io , io - 1 + n
                  read(iu) a(ind+1:ind+ncol(i),s)
                  ind = ind + ncol(i)
                end do
              else
                i = sum(ncol(io:io-1+n))
                nb = nb + 1
                call MPI_IRecv(a(ind+1,s), i, MPI_Double_Precision, &
                    0, gio, comm, ibuf(nb), MPIerror)
                ind = ind + i
              end if
            else if ( lNode == 0 ) then
              i = 0
              do io = gio , gio - 1 + n
                read(iu) buf(i+1:i+lgncol(io))
                i = i + lgncol(io)
              end do
              call MPI_Send(buf(1), i, MPI_Double_Precision, &
                  BNode, gio, comm, MPIerror)
            end if
            gio = gio + n
          end do

          if ( lNode /= 0 .and. nb > 0 ) then
            call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
          end if

        end do

        if ( lNode == 0 ) then
          deallocate(buf)
        else
          deallocate(ibuf)
        end if

        if ( .not. present(gncol) ) deallocate(lgncol)

#else
        call die('Error in distribution for, io_read_d2D[sp=1]')
#endif
      else if ( lNode == 0 ) then

        do s = 1 , dim2 
          ind = 0
          do io = 1 , no
            read(iu) a(ind+1:ind+ncol(io),s)
            ind = ind + ncol(io)
          end do
        end do

      end if

    end subroutine read_sp_dim1
  
    subroutine read_sp_dim2()
      integer :: io, ind, n, i, nb
#ifdef MPI
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: gio, max_n
      integer :: BNode
#endif

      if ( ldit ) then
#ifdef MPI

        nb = count_blocks(dit,no)

        ! Allocate maximum number of entries
        if ( lNode == 0 ) then
          max_n = max_consecutive_sum(dit,no,lgncol) * dim2
          allocate(buf(max_n))
        else
          allocate(ibuf(nb))
        end if

        ! Read sparse blocks and distribute
        gio = 1 
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( BNode == lNode ) then

            io = index_global_to_local(dit,gio,lNode)

            if ( lNode == 0 ) then
              do i = io , io - 1 + n
                read(iu) a(1:dim2,ind+1:ind+ncol(i))
                ind = ind + ncol(i)
              end do
            else
              ! count the number of received entities
              i = sum(ncol(io:io-1+n))
              nb = nb + 1
              call MPI_IRecv(a(1,ind+1), dim2*i, MPI_Double_Precision, &
                  0, gio, comm, ibuf(nb), MPIerror)
              ind = ind + i
            end if
          else if ( lNode == 0 ) then
            i = 0
            do io = gio , gio - 1 + n
              read(iu) buf(i+1:i+dim2*lgncol(io))
              i = i + dim2*lgncol(io)
            end do
            call MPI_Send(buf(1), i, MPI_Double_Precision, &
                BNode, gio, comm, MPIerror)
          end if
          gio = gio + n
        end do

        if ( lNode == 0 ) then
          deallocate(buf)
        else
          if ( nb > 0 ) then
            call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
          end if
          deallocate(ibuf)
        end if

        if ( .not. present(gncol) ) deallocate(lgncol)

#else
        call die('Error in distribution for, io_read_d2D[sp=2]')
#endif
      else if ( lNode == 0 ) then

        ind = 0
        do io = 1 , no
          read(iu) a(1:dim2,ind+1:ind+ncol(io))
          ind = ind + ncol(io)
        end do

      end if

    end subroutine read_sp_dim2
    
  end subroutine io_read_d2D

  
  subroutine io_write_d2D(iu, dSp2D, gncol)

    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData2D), intent(inout) :: dSp2D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: lno, no, n_nzs, sp_dim, dim2
    integer :: comm, lNode
    logical :: ldit

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit)

    ! Retrieve data
    a => val(dSp2D)
    sp_dim = spar_dim(dSp2D)
    if ( sp_dim == 1 ) then
      dim2 = size(a, dim=2)
    else
      dim2 = size(a, dim=1)
    end if

    if ( sp_dim == 1 ) then
      call write_sp_dim1()
    else
      call write_sp_dim2()
    end if

  contains

    subroutine write_sp_dim1()
      integer :: io, s, ind, n, i, nb
#ifdef MPI
      integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: gio, max_n
#endif

      if ( ldit ) then

#ifdef MPI
        if ( present(gncol) ) then
          lgncol => gncol
        else
          allocate(lgncol(no))
          lgncol(1) = -1
        end if
        if ( lgncol(1) < 0 ) then
          call Node_Sp_gncol(0, sp, dit, no, lgncol)
        end if

#else
        call die('Error in distribution, io_write_d2D[sp=1]')
#endif

#ifdef MPI
        nb = count_blocks(dit,no)

        if ( lNode == 0 ) then
          max_n = max_consecutive_sum(dit,no,lgncol)
          allocate(buf(max_n))
        else
          allocate(ibuf(nb))
        end if

        do s = 1 , dim2
          gio = 1
          ind = 0
          nb = 0
          do while ( gio <= no )

            BNode = node_handling_element(dit,gio)

            ! Get number of consecutive orbitals
            n = count_consecutive(dit,no,gio)

            if ( BNode == lNode ) then

              io = index_global_to_local(dit,gio,lNode)

              if ( lNode == 0 ) then
                do i = io , io - 1 + n
                  write(iu) a(ind+1:ind+ncol(i),s)
                  ind = ind + ncol(i)
                end do
              else
                i = sum(ncol(io:io-1+n))
                nb = nb + 1
                call MPI_ISend(a(ind+1,s), i, MPI_Double_Precision, &
                    0, gio, comm, ibuf(nb), MPIerror)
                ind = ind + i
              end if
            else if ( lNode == 0 ) then
              call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
                  BNode, gio, comm, MPIstatus, MPIerror)
              if ( MPIerror /= MPI_Success ) &
                  call die('Error in code (1): io_write_d2D[sp=1]')
              i = 0
              do io = gio , gio - 1 + n
                write(iu) buf(i+1:i+lgncol(io))
                i = i + lgncol(io)
              end do

            end if

            gio = gio + n
          end do

          if ( lNode /= 0 .and. nb > 0 ) then
            call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
          end if

        end do

        if ( .not. present(gncol) ) deallocate(lgncol)
        if ( lNode == 0 ) then
          deallocate(buf)
        else
          deallocate(ibuf)
        end if
#else
        call die('Error in io_write_d2D[sp=1]')
#endif
      else if ( lNode == 0 ) then

        do s = 1 , dim2
          ind = 0
          do io = 1 , no
            write(iu) a(ind+1:ind+ncol(io),s)
            ind = ind + ncol(io)
          end do
        end do

      end if

    end subroutine write_sp_dim1

    subroutine write_sp_dim2()
      integer :: io, ind, n, i, nb
#ifdef MPI
      integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: gio, max_n
#endif

      if ( ldit ) then

#ifdef MPI
        if ( present(gncol) ) then
          lgncol => gncol
        else
          allocate(lgncol(no))
          lgncol(1) = -1
        end if
        if ( lgncol(1) < 0 ) then
          call Node_Sp_gncol(0, sp, dit, no, lgncol)
        end if

#else
        call die('Error in distribution, io_write_d2D[sp=2]')
#endif

#ifdef MPI
        nb = count_blocks(dit,no)

        if ( lNode == 0 ) then
          max_n = max_consecutive_sum(dit,no,lgncol) * dim2
          allocate(buf(max_n))
        else
          allocate(ibuf(nb))
        end if

        gio = 1
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( BNode == lNode ) then

            io = index_global_to_local(dit,gio,lNode)

            if ( lNode == 0 ) then
              do i = io , io - 1 + n
                write(iu) a(1:dim2,ind+1:ind+ncol(i))
                ind = ind + ncol(i)
              end do
            else
              i = sum(ncol(io:io-1+n))
              nb = nb + 1
              call MPI_ISend(a(1,ind+1), dim2*i, MPI_Double_Precision, &
                  0, gio, comm, ibuf(nb), MPIerror)
              ind = ind + i
            end if
          else if ( lNode == 0 ) then
            call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
                BNode, gio, comm, MPIstatus, MPIerror)
            if ( MPIerror /= MPI_Success ) &
                call die('Error in code (2): io_write_d2D[sp=2]')
            i = 0
            do io = gio , gio - 1 + n
              write(iu) buf(i+1:i+dim2*lgncol(io))
              i = i + dim2*lgncol(io)
            end do
          end if
          gio = gio + n
        end do

        if ( .not. present(gncol) ) deallocate(lgncol)
        if ( lNode == 0 ) then
          deallocate(buf)
        else
          if ( nb > 0 ) then
            call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
          end if
          deallocate(ibuf)
        end if
#else
        call die('Error in io_write_d2D[sp=2]')
#endif
      else if ( lNode == 0 ) then

        ind = 0
        do io = 1 , no
          write(iu) a(1:dim2,ind+1:ind+ncol(io))
          ind = ind + ncol(io)
        end do

      end if

    end subroutine write_sp_dim2

  end subroutine io_write_d2D

  
  subroutine io_write_r2D(iu, dSp2D, gncol)

    use precision, only: psp => sp
    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData2D), intent(inout) :: dSp2D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: lgncol(:) => null(), ncol(:)

    integer :: lno, no, n_nzs, sp_dim, dim2
    integer :: comm, lNode
    logical :: ldit

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)

    ! Get default parameters
    call parallel_default_args(comm, lNode, dit, ldit)

    ! Retrieve data
    a => val(dSp2D)
    sp_dim = spar_dim(dSp2D)
    if ( sp_dim == 1 ) then
      dim2 = size(a, dim=2)
    else
      dim2 = size(a, dim=1)
    end if

    if ( sp_dim == 1 ) then
      call write_sp_dim1()
    else
      call write_sp_dim2()
    end if

  contains

    subroutine write_sp_dim1()
      integer :: io, s, ind, n, i, nb
#ifdef MPI
      integer :: MPIstatus(MPI_STATUS_SIZE)
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: gio, max_n
      integer :: MPIerror, BNode
#endif

      if ( ldit ) then

#ifdef MPI
        if ( present(gncol) ) then
          lgncol => gncol
        else
          allocate(lgncol(no))
          lgncol(1) = -1
        end if
        if ( lgncol(1) < 0 ) then
          call Node_Sp_gncol(0, sp, dit, no, lgncol)
        end if

#else
        call die('Error in distribution, io_write_r2D[sp=1]')
#endif

#ifdef MPI
        nb = count_blocks(dit,no)

        if ( lNode == 0 ) then
          max_n = max_consecutive_sum(dit,no,lgncol)
          allocate(buf(max_n))
        else
          allocate(ibuf(nb))
        end if

        do s = 1 , dim2
          gio = 1
          ind = 0
          nb = 0
          do while ( gio <= no )

            BNode = node_handling_element(dit,gio)

            ! Get number of consecutive orbitals
            n = count_consecutive(dit,no,gio)

            if ( BNode == lNode ) then

              io = index_global_to_local(dit,gio,lNode)

              if ( lNode == 0 ) then
                do i = io , io - 1 + n
                  write(iu) real(a(ind+1:ind+ncol(i),s), psp)
                  ind = ind + ncol(i)
                end do
              else
                i = sum(ncol(io:io-1+n))
                nb = nb + 1
                call MPI_ISend(a(ind+1,s), i, MPI_Double_Precision, &
                    0, gio, comm, ibuf(nb), MPIerror)
                ind = ind + i
              end if
            else if ( lNode == 0 ) then
              call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
                  BNode, gio, comm, MPIstatus, MPIerror)
              if ( MPIerror /= MPI_Success ) &
                  call die('Error in code (1): io_write_r2D[sp=1]')
              i = 0
              do io = gio , gio - 1 + n
                write(iu) real(buf(i+1:i+lgncol(io)), psp)
                i = i + lgncol(io)
              end do
            end if
            gio = gio + n
          end do

          if ( lNode /= 0 .and. nb > 0 ) then
            call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
          end if

        end do

        if ( .not. present(gncol) ) deallocate(lgncol)
        if ( lNode == 0 ) then
          deallocate(buf)
        else
          deallocate(ibuf)
        end if
#else
        call die('Error in io_write_r2D[sp=1]')
#endif
      else if ( lNode == 0 ) then

        do s = 1 , dim2
          ind = 0
          do io = 1 , no
            write(iu) real(a(ind+1:ind+ncol(io),s), psp)
            ind = ind + ncol(io)
          end do
        end do

      end if

    end subroutine write_sp_dim1

    subroutine write_sp_dim2()
      integer :: io, ind, n, i, nb
#ifdef MPI
      integer :: MPIstatus(MPI_STATUS_SIZE)
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: gio, max_n
      integer :: MPIerror, BNode
#endif

      if ( ldit ) then

#ifdef MPI
        if ( present(gncol) ) then
          lgncol => gncol
        else
          allocate(lgncol(no))
          lgncol(1) = -1
        end if
        if ( lgncol(1) < 0 ) then
          call Node_Sp_gncol(0, sp, dit, no, lgncol)
        end if

#else
        call die('Error in distribution, io_write_r2D[sp=2]')
#endif

#ifdef MPI
        nb = count_blocks(dit,no)

        if ( lNode == 0 ) then
          max_n = max_consecutive_sum(dit,no,lgncol) * dim2
          allocate(buf(max_n))
        else
          allocate(ibuf(nb))
        end if

        gio = 1
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( BNode == lNode ) then

            io = index_global_to_local(dit,gio,lNode)

            if ( lNode == 0 ) then
              do i = io , io - 1 + n
                write(iu) real(a(1:dim2,ind+1:ind+ncol(i)), psp)
                ind = ind + ncol(i)
              end do
            else
              i = sum(ncol(io:io-1+n))
              nb = nb + 1
              call MPI_ISend(a(1,ind+1), dim2*i, MPI_Double_Precision, &
                  0, gio, comm, ibuf(nb), MPIerror)
              ind = ind + i
            end if
          else if ( lNode == 0 ) then
            call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
                BNode, gio, comm, MPIstatus, MPIerror)
            if ( MPIerror /= MPI_Success ) &
                call die('Error in code (2): io_write_r2D[sp=2]')
            i = 0
            do io = gio , gio - 1 + n
              write(iu) real(buf(i+1:i+dim2*lgncol(io)), psp)
              i = i + dim2*lgncol(io)
            end do
          end if
          gio = gio + n
        end do

        if ( .not. present(gncol) ) deallocate(lgncol)
        if ( lNode == 0 ) then
          deallocate(buf)
        else
          if ( nb > 0 ) then
            call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
          end if
          deallocate(ibuf)
        end if
#else
        call die('Error in io_write_r2D[sp=2]')
#endif
      else if ( lNode == 0 ) then

        ind = 0
        do io = 1 , no
          write(iu) real(a(1:dim2,ind+1:ind+ncol(io)), psp)
          ind = ind + ncol(io)
        end do

      end if

    end subroutine write_sp_dim2

  end subroutine io_write_r2D

end module io_sparse_m

