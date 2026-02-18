! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module ncdf_io_m

  use precision, only : sp, dp, grid_p
  use parallel, only : Node, Nodes

  use class_OrbitalDistribution
  use class_Sparsity
  use io_sparse_m, only: max_consecutive_sum, Node_Sp_gncol
  use io_sparse_m, only: count_consecutive, count_blocks

#ifdef NCDF_4
  use netcdf_ncdf, ncdf_parallel => parallel
#endif
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  private

#ifdef NCDF_4
  public :: cdf_init_mesh
  public :: cdf_r_grid
  public :: cdf_w_grid

  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

  public :: cdf_r_sp
  public :: cdf_r_d1D, cdf_r_d2D

#endif

  ! Create the type to hold the mesh distribution data
  type :: tMeshDist
    private
    ! number of mesh divisions in each axis
    integer :: m(3)
    ! mesh box bounds of each node in each direction
    ! box(1,iAxis,iNode)=lower bounds
    ! box(2,iAxis,iNode)=upper bounds
    integer, allocatable :: box(:,:,:)
  end type tMeshDist
  type(tMeshDist), save :: distr

contains

#ifdef NCDF_4

  subroutine cdf_init_mesh(mesh,nsm)

    use parallel, ProcY => ProcessorY

    ! mesh divisions (fine points), number of fine points per big points
    integer, intent(in) :: mesh(3), nsm

    ! Local quantities
    integer :: nm(3), nmesh, ntot
    integer :: ProcZ, blocY, blocZ, nremY, nremZ
    integer :: dimX, dimY, dimZ
    integer :: PP, iniY, iniZ, PY, PZ

    if ( .not. allocated(distr%box) ) then
      allocate(distr%box(2,3,0:Nodes-1))
    end if

    nm(1:3) = mesh(1:3) / nsm
    nmesh = product(mesh)
    distr%m(1:3) = mesh(1:3)

    ProcZ = Nodes/ProcY

    blocY = nm(2)/ProcY
    nremY = nm(2) - blocY*ProcY
    blocZ = nm(3)/ProcZ
    nremZ = nm(3) - blocZ*ProcZ

    dimX = nm(1) * nsm

    ntot = 0

    PP   = 0
    iniY = 1
    do PY = 1, ProcY

      dimY = blocY
      if ( PY <= nremY ) dimY = dimY + 1  ! Add extra points starting from the first nodes
      dimY = dimY * nsm                 ! For fine points

      iniZ = 1
      do PZ = 1, ProcZ
        dimZ = blocZ
        if ( PZ <= nremZ ) dimZ = dimZ + 1
        dimZ = dimZ*nsm                 ! For fine points

        distr%box(1,1,PP) = 1
        distr%box(2,1,PP) = dimX
        distr%box(1,2,PP) = iniY
        distr%box(2,2,PP) = iniY + dimY - 1
        distr%box(1,3,PP) = iniZ
        distr%box(2,3,PP) = iniZ + dimZ - 1

        ntot = ntot + dimX * dimY * dimZ

        iniZ = iniZ + dimZ
        PP   = PP + 1

      end do

      iniY = iniY + dimY

    end do

    if (ntot /= nmesh) then
      if (Node == 0) then
        write(6,*) "Nominal npt: ", nmesh, " /= assigned npt:", ntot
      end if
      call die()
    end if

  end subroutine cdf_init_mesh

  ! Reads in a sparsity pattern at the
  ! current position in the file (iu)
  ! The sparsity pattern "sp" will be returned
  ! as populated.
  ! If dist is supplied it will distribute
  ! the sparsity pattern as supplied (this implies Bcast = .true.)
  ! Else if Bcast is true it will b-cast the sparsity 
  ! pattern fully.
  subroutine cdf_r_Sp(ncdf, no, sp, tag, dit, Bcast, gncol)

    ! File handle
    type(hNCDF), intent(inout) :: ncdf
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(no)

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, gind, ind, nl, n, i, nb
    logical :: ldit, lBcast, lIO
    integer, pointer :: lgncol(:) => null()
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: MPIerror, BNode
#endif

#ifdef MPI
    ldit = present(dit)
#else
    ldit = .false.
#endif
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    if ( lIO ) then

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
      end if

      ! First read in number of non-zero 
      ! entries per orbital
      call ncdf_get_var(ncdf,'n_col',lgncol)

    end if

    nl = no

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit ) then

      ! Number of local elements
      nl = num_local_elements(dit,no,Node)
      allocate(ncol(nl))

      ! allocate all requests
      nb = count_blocks(dit,no)
      allocate(ibuf(nb))

      ! Distribute it
      gio = 1
      nb = 0
      do while ( gio <= no ) 

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == Node ) then

          io = index_global_to_local(dit,gio,Node)

          if ( Node == 0 ) then
            ncol(io:io-1+n) = lgncol(gio:gio-1+n)
          else
            nb = nb + 1
            call MPI_IRecv(ncol(io), n, MPI_Integer, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if

        else if ( Node == 0 ) then
          nb = nb + 1
          call MPI_ISend(lgncol(gio), n, MPI_Integer, &
              BNode, gio, MPI_Comm_World, ibuf(nb), MPIerror)

        end if

        gio = gio + n

      end do

      if ( nb > 0 ) then
        call MPI_WaitAll(nb,ibuf,MPI_STATUSES_IGNORE,MPIerror)
      end if
      deallocate(ibuf)

    else if ( lBcast ) then

      ! Everything should be b-casted
      if ( Node == 0 ) then
        ncol => lgncol
      else
        allocate(ncol(nl))
      end if

      ! Bcast everything
      call MPI_Bcast(ncol(1), nl, MPI_Integer, &
          0, MPI_Comm_World, MPIError)

    else if ( lIO ) then

      ncol => lgncol

    end if
#else
    ! Point to the buffer
    ncol => lgncol
#endif

    ! Allocate pointer
    allocate(l_ptr(nl))

    l_ptr(1) = 0
    do io = 2 , nl
      l_ptr(io) = l_ptr(io-1) + ncol(io-1)
    end do

    ! Number of local non-zero elements
    ! (also works for any bcast methods)
    n_nzs = l_ptr(nl) + ncol(nl)

    ! Allocate space
    allocate(l_col(n_nzs))

#ifdef MPI
    if ( ldit ) then

      ! We have a distributed read
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(ibuf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Read in columns
      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no ) 

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == Node ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,Node)

          if ( Node == 0 ) then

            i = sum(ncol(io:io-1+n))
            call ncdf_get_var(ncdf, 'list_col', l_col(ind+1:ind+i), &
                start=(/gind/), count=(/i/))
            ind = ind + i
            gind = gind + i

          else

            ! count the number of received entities
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_IRecv(l_col(ind+1), i, MPI_Integer, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
            ind = ind + i

          end if

        else if ( Node == 0 ) then

          i = sum(lgncol(gio:gio-1+n))
          call ncdf_get_var(ncdf, 'list_col', ibuf(1:i), &
              start=(/gind/), count=(/i/))

          call MPI_Send(ibuf(1), i, MPI_Integer, &
              BNode, gio, MPI_Comm_World, MPIerror)

          gind = gind + i

        end if

        gio = gio + n

      end do

      if ( Node == 0 ) then
        if ( .not. present(gncol) ) deallocate(lgncol)
      else
        if ( nb > 0 ) then
          call MPI_Waitall(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
      end if
      deallocate(ibuf)

    else if ( lBcast ) then

      if ( Node == 0 ) then

        ind = sum(ncol(1:no))
        call ncdf_get_var(ncdf,'list_col',l_col(1:ind), &
            count=(/ind/))
        if ( ind /= n_nzs ) then
          call die('Error in reading sparsity pattern, &
              &size not equivalent.')
        end if
      end if

      ! Bcast
      call MPI_Bcast(l_col(1), n_nzs, MPI_Integer, &
          0, MPI_Comm_World, MPIError)

    else
#endif

      ind = sum(ncol(1:no))
      call ncdf_get_var(ncdf, 'list_col', l_col(1:ind), &
          count=(/ind/) )

#ifdef MPI       
    end if
#endif

    ! Create the sparsity pattern
    call newSparsity(sp,nl,no, n_nzs, ncol, l_ptr, l_col, trim(tag))

    ! de-allocate
    deallocate(l_ptr,l_col)
    if ( ldit ) deallocate(ncol)
    if ( lBcast .and. Node /= 0 ) deallocate(ncol)
    if ( lIO .and. .not. present(gncol) ) deallocate(lgncol)

  end subroutine cdf_r_Sp

  subroutine cdf_w_Sp(ncdf, sp, dit, gncol)

    type(hNCDF), intent(inout) :: ncdf
    type(Sparsity), intent(inout) :: sp
    type(OrbitalDistribution), intent(in), optional :: dit
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    integer, pointer :: lgncol(:) => null()
    integer, pointer :: ncol(:), l_col(:) => null()

    logical :: ldit
    integer :: lno, no, n_nzs
#ifdef MPI
    integer :: MPIerror
#endif
    integer :: n_nnzs, g_nzs

    ! Write the sparsity to the file...
    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,list_col=l_col,nnzs=n_nzs)

    ldit = present(dit)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then
#ifdef MPI
      call MPI_Reduce(n_nzs,g_nzs,1,MPI_Integer, &
          MPI_Sum, 0, MPI_Comm_World, MPIerror)

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if
#else
      call die('Error in code, non-full contained sp')
#endif
    else
      ! serial (no distribution)
      lgncol => ncol
      g_nzs = n_nzs
    end if

    ! Check dimension nnzs
    call ncdf_inq_dim(ncdf,'nnzs',len=n_nnzs)
    if ( Node == 0 .and. n_nnzs /= g_nzs ) then
      call die('Number of non-zero elements is not equivalent.')
    end if

#ifdef MPI
    if ( ldit ) then
      if ( parallel_io(ncdf) ) then
        ! Write out in parallel
        call write_parallel()
      else
        call write_parallel_single()
      end if

      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      call write_sequential()
    end if
#else
    call write_sequential()
#endif

  contains

#ifdef MPI
    subroutine write_parallel()

      integer :: BNode
      integer :: gio, io, ind, gind, n, i, nb

      if ( lgncol(1) < 0 ) then
        ! this means broadcast lgncol on all nodes
        call Node_Sp_gncol(-1, sp, dit, no, lgncol)
      end if

      ! we write it using MPI
      call ncdf_par_access(ncdf,name='n_col',access=NF90_COLLECTIVE)

      ! Calculate indices for this node
      io = no / Nodes
      ind = io * Node + 1
      if ( Node == Nodes - 1 ) then
        ! correct the last node to write the rest, probably not
        ! balanced, but shouldn't be too problematic for
        ! very large matrices
        io = io + mod(no, Nodes)
      end if

      ! Write in parallel (faster than in below loop)
      call ncdf_put_var(ncdf,'n_col',lgncol(ind:ind+io-1), &
          start=(/ind/), count=(/io/))

      call ncdf_par_access(ncdf,name='list_col',access=NF90_COLLECTIVE)

      ! Determine the number of blocks
      gio = 1
      gind = 1
      ind = 0
      nb = 1
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        ! Number of consecutive elements
        i = sum(lgncol(gio:gio-1+n))

        if ( BNode == Node ) then
          ! Count number of actual writes
          nb = nb + 1
          call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+i), &
              count=(/i/),start=(/gind/))
          ind = ind + i
        end if

        gind = gind + i
        gio = gio + n

      end do

      ind = count_blocks(dit, no)
      do i = nb, ind ! only runned if some iterations are missed
        ! Fake access
        call ncdf_put_var(ncdf,'list_col',l_col(1:1), &
            start=(/1/), count=(/0/) )
      end do

    end subroutine write_parallel
    
    subroutine write_parallel_single()

      integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
      integer :: gio, io, max_n, ind, gind, n, i, nb
      integer, allocatable :: ibuf(:)

      if ( lgncol(1) < 0 ) then
        ! this means gather lgncol on IO-node
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

      ! Write gncol
      call ncdf_put_var(ncdf,'n_col',lgncol)

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(ibuf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( Node == BNode ) then
          
          io = index_global_to_local(dit,gio,Node)
          i = sum(ncol(io:io-1+n))

          if ( Node == 0 ) then
            call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+i), &
                count=(/i/),start=(/gind/))
            gind = gind + i
          else
            nb = nb + 1
            call MPI_ISend(l_col(ind+1), i, MPI_Integer, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i
          
        else if ( Node == 0 ) then
          call MPI_Recv(ibuf(1), max_n, MPI_Integer, &
              BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: cdf_w_Sp')
          
          call MPI_Get_Count(MPIstatus, MPI_Integer, i, MPIerror)
          call ncdf_put_var(ncdf,'list_col',ibuf(1:i), &
              count=(/i/),start=(/gind/))
          gind = gind + i
        end if

        gio = gio + n
      end do

      if ( Node /= 0 .and. nb > 0 ) then
        call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
      end if

      deallocate(ibuf)

    end subroutine write_parallel_single
#endif

    subroutine write_sequential()
      call ncdf_put_var(ncdf,'n_col',ncol)
      call ncdf_put_var(ncdf,'list_col',l_col)
    end subroutine write_sequential
    
  end subroutine cdf_w_Sp


  subroutine cdf_r_d1D(ncdf,vname,sp, dSp1D, tag, dit, Bcast, gncol)

    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(Sparsity), intent(inout) :: sp
    type(dSpData1D), intent(inout) :: dSp1D
    ! The tag of the data format
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

    integer :: lno, no, n_nzs
    logical :: ldit, lBcast, lIO
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: MPIerror, BNode
    integer :: max_n
#endif

    integer :: gind, gio, ind, n, nb, i, io

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)
    if ( ldit ) ldit = lno /= no

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
      call die('Error in distribution, cdf_r_d1D')
#endif

    else
      ! Create the Fake distribution
#ifdef MPI
      call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
      call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
      call newdSpData1D(sp,fdit,dSp1D,name=trim(tag))
      ! Clean up the distribution again
      call delete(fdit)
    end if

    a => val(dSp1D)

    if ( ldit ) then
#ifdef MPI

      nb = count_blocks(dit,no)

      ! Allocate the maximum number of entries
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1 
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( BNode == Node ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,Node)

          if ( Node == 0 ) then

            ! Count number of elements
            i = sum(ncol(io:io-1+n))
            call ncdf_get_var(ncdf, vname, a(ind+1:ind+i), &
                start=(/gind/), count=(/i/))
            ind = ind + i
            gind = gind + i

          else

            ! count the number of received entities
            i = sum(ncol(io:io-1+n))
            nb = nb + 1
            call MPI_IRecv(a(ind+1), i, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
            ind = ind + i

          end if

        else if ( Node == 0 ) then

          i = sum(lgncol(gio:gio-1+n))
          call ncdf_get_var(ncdf, vname, buf(1:i), &
              start=(/gind/), count=(/i/))

          call MPI_Send(buf(1), i, MPI_Double_Precision, &
              BNode, gio, MPI_Comm_World, MPIerror)

          gind = gind + i

        end if

        gio = gio + n

      end do

      if ( .not. present(gncol) ) deallocate(lgncol)
      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if
#else
      call die('Error in distribution for, cdf_r_d1D')
#endif
    else if ( lIO ) then

      i = sum(ncol(1:no))
      call ncdf_get_var(ncdf, vname, a(:), count=(/i/))

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then

      call MPI_Bcast(a(1),n_nzs,MPI_Double_Precision, &
          0, MPI_Comm_World, MPIError)

    end if
#endif

  end subroutine cdf_r_d1D
  
  subroutine cdf_w_d1D(ncdf, vname, dSp1D, gncol)

    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData1D), intent(inout) :: dSp1D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:), lgncol(:)
    integer :: lno, no, n_nzs
    real(dp), pointer :: a(:)

    logical :: ldit

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,list_col=l_col,nnzs=n_nzs)

    ! If they are different we should 
    ! use the distribution setting
    ldit = lno /= no

    if ( ldit ) then

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if

    else
      ! serial (no distribution)
      lgncol => ncol
    end if
    
    a => val(dSp1D)

#ifdef MPI
    if ( ldit ) then
      if ( parallel_io(ncdf) ) then
        ! Write out in parallel
        call write_parallel()
      else
        call write_parallel_single()
      end if
      
      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      call write_sequential()
    end if
#else
    call write_sequential()
#endif

  contains

#ifdef MPI
    subroutine write_parallel()

      integer :: BNode
      integer :: gio, ind, gind, n, i, nb

      if ( lgncol(1) < 0 ) then
        ! this means broadcast lgncol on all nodes
        call Node_Sp_gncol(-1, sp, dit, no, lgncol)
      end if

      call ncdf_par_access(ncdf, name=vname, access=NF90_COLLECTIVE)

      ! Loop blocks
      gio = 1
      gind = 1
      ind = 0
      nb = 1
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        ! Number of consecutive elements
        i = sum(lgncol(gio:gio-1+n))

        if ( BNode == Node ) then
          ! Count number of actual writes
          nb = nb + 1
          call ncdf_put_var(ncdf, vname, a(ind+1:ind+i), &
              count=(/i/), start=(/gind/))
          ind = ind + i
        end if

        gind = gind + i
        gio = gio + n

      end do

      ind = count_blocks(dit, no)
      do i = nb, ind ! only runned if some iterations are missed
        ! Fake access
        call ncdf_put_var(ncdf, vname, a(1:1), start=(/1/), count=(/0/))
      end do

    end subroutine write_parallel
    
    subroutine write_parallel_single()

      integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
      integer :: gio, io, max_n, ind, gind, n, i, nb
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)

      if ( lgncol(1) < 0 ) then
        ! this means gather lgncol on IO-node
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        if ( Node == BNode ) then
          
          io = index_global_to_local(dit,gio,Node)
          i = sum(ncol(io:io-1+n))

          if ( Node == 0 ) then
            call ncdf_put_var(ncdf, vname, a(ind+1:ind+i), &
                count=(/i/), start=(/gind/))
            gind = gind + i
          else
            nb = nb + 1
            call MPI_ISend(a(ind+1), i, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i
          
        else if ( Node == 0 ) then
          call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
              BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: cdf_w_d1D')
          
          call MPI_Get_Count(MPIstatus, MPI_Double_Precision, i, MPIerror)
          call ncdf_put_var(ncdf, vname, buf(1:i), &
              count=(/i/), start=(/gind/))
          gind = gind + i
        end if

        gio = gio + n
      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if

    end subroutine write_parallel_single
#endif

    subroutine write_sequential()
      call ncdf_put_var(ncdf, vname, a)
    end subroutine write_sequential
    
  end subroutine cdf_w_d1D


  subroutine cdf_r_d2D(ncdf, vname, sp, dSp2D, dim2, tag, &
      sparsity_dim, dit, Bcast, gncol)

    use class_dSpData2D

    ! File handle
    type(hNCDF), intent(inout) :: ncdf
    ! Variable name
    character(len=*), intent(in) :: vname
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

    integer :: lno, no, n_nzs
    integer :: sp_dim
    logical :: ldit, lBcast, lIO
#ifdef MPI
    integer :: MPIerror
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    sp_dim = 1
    if ( present(sparsity_dim) ) sp_dim = sparsity_dim

    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,nnzs=n_nzs)
    if ( ldit ) ldit = lno /= no

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
      call die('Error in distribution, cdf_r_d2D')
#endif

    else
      ! Create the Fake distribution
#ifdef MPI
      call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
      call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
      call newdSpData2D(sp,dim2,fdit,dSp2D,name=trim(tag), &
          sparsity_dim=sp_dim)
    end if

    a => val(dSp2D)

#ifdef MPI
    if ( ldit ) then
      if ( sp_dim == 1 ) then
        call read_sp_dim1()
      else
        call read_sp_dim2()
      end if

      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      ! this is both lio and serial read
      call ncdf_get_var(ncdf, vname, a)
    end if
#else
    call ncdf_get_var(ncdf, vname, a)
#endif

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then
      
      call MPI_Bcast(a(1,1), dim2*n_nzs, MPI_Double_Precision, &
          0, MPI_Comm_World, MPIError)

    end if
#endif

#ifdef MPI
  contains

    subroutine read_sp_dim1()

      integer :: io, ind, n, i, nb, s, gio, gind, max_n
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: MPIerror, BNode

      nb = count_blocks(dit,no)

      ! Allocate maximum number of entries
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      do s = 1 , dim2

        gio = 1 
        gind = 1
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( BNode == Node ) then

            ! Get the local orbital
            io = index_global_to_local(dit,gio,Node)

            i = sum(ncol(io:io-1+n))
            if ( Node == 0 ) then

              call ncdf_get_var(ncdf, vname, a(ind+1:ind+i,s), &
                  start=(/gind,s/), count=(/i/))
              gind = gind + i

            else

              ! count the number of received entities
              nb = nb + 1
              call MPI_IRecv(a(ind+1,s), i, MPI_Double_Precision, &
                  0, gio, MPI_Comm_World, ibuf(nb), MPIerror)

            end if
            ind = ind + i

          else if ( Node == 0 ) then

            i = sum(lgncol(gio:gio-1+n))
            call ncdf_get_var(ncdf, vname, buf(1:i), &
                start=(/gind,s/), count=(/i/))
            call MPI_Send(buf(1), i, MPI_Double_Precision, &
                BNode, gio, MPI_Comm_World, MPIerror)
            gind = gind + i
          end if
          
          gio = gio + n
        end do

        if ( Node /= 0 .and. nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if

      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        deallocate(ibuf)
      end if

    end subroutine read_sp_dim1

    subroutine read_sp_dim2()
      
      integer :: io, ind, n, i, nb, gio, gind, max_n
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)
      integer :: MPIerror, BNode

      nb = count_blocks(dit,no)

      ! Allocate maximum number of entries
      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol) * dim2
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Read sparse blocks and distribute
      gind = 1
      ind = 0
      gio = 1 
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)
        if ( BNode == Node ) then

          ! Get the local orbital
          io = index_global_to_local(dit,gio,Node)

          i = sum(ncol(io:io-1+n))
          if ( Node == 0 ) then
            call ncdf_get_var(ncdf, vname, a(1:dim2,ind+1:ind+i), &
                start=(/1,gind/), count=(/dim2,i/) )
            gind = gind + i

          else

            ! count the number of received entities
            nb = nb + 1
            call MPI_IRecv(a(1,ind+1), dim2*i, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i

        else if ( Node == 0 ) then

          i = sum(lgncol(gio:gio-1+n))
          call ncdf_get_var(ncdf, vname, buf(1:i), &
              start=(/1,gind/), count=(/dim2,i/))
          call MPI_Send(buf(1), i, MPI_Double_Precision, &
                BNode, gio, MPI_Comm_World, MPIerror)
          gind = gind + i

        end if

        gio = gio + n
      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if

    end subroutine read_sp_dim2

#endif
    
  end subroutine cdf_r_d2D

  subroutine cdf_w_d2D(ncdf, vname, dSp2D, gncol)

    use class_dSpData2D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData2D), intent(inout) :: dSp2D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:), lgncol(:)
    integer :: lno, no, n_nzs
    integer :: dim2, sp_dim
    real(dp), pointer :: a(:,:)

    logical :: ldit

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,list_col=l_col,nnzs=n_nzs)

    ! If they are different we should 
    ! use the distribution setting
    ldit = lno /= no

    if ( ldit ) then

      if ( present(gncol) ) then
        lgncol => gncol
      else
        allocate(lgncol(no))
        lgncol(1) = -1
      end if

    else
      ! serial (no distribution)
      lgncol => ncol
    end if

    a => val(dSp2D)
    sp_dim = spar_dim(dSp2D)
    if ( sp_dim == 1 ) then
      dim2 = size(a, dim=2)
    else
      dim2 = size(a, dim=1)
    end if

#ifdef MPI
    if ( ldit ) then
      if ( parallel_io(ncdf) ) then
        ! Write out in parallel
        call write_parallel()
      else
        call write_parallel_single()
      end if

      if ( .not. present(gncol) ) deallocate(lgncol)
    else
      call write_sequential()
    end if
#else
    call write_sequential()
#endif

  contains

#ifdef MPI
    subroutine write_parallel()
      integer :: BNode
      integer :: gio, ind, gind, n, i, nb

      if ( lgncol(1) < 0 ) then
        ! this means broadcast lgncol on all nodes
        call Node_Sp_gncol(-1, sp, dit, no, lgncol)
      end if

      call ncdf_par_access(ncdf, name=vname, access=NF90_COLLECTIVE)

      gio = 1
      gind = 1
      ind = 0
      nb = 1
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)

        ! Number of consecutive elements
        i = sum(lgncol(gio:gio-1+n))

        if ( BNode == Node ) then
          ! Count number of actual writes
          nb = nb + 1
          if ( sp_dim == 1 ) then
            ! this will make a temporary copy... :(
            call ncdf_put_var(ncdf, vname, a(ind+1:ind+i,:), &
                count=(/i,dim2/), start=(/gind,1/))
          else
            call ncdf_put_var(ncdf, vname, a(:,ind+1:ind+i), &
                count=(/dim2,i/), start=(/1,gind/))
          end if
          ind = ind + i
        end if

        gind = gind + i
        gio = gio + n

      end do

      ind = count_blocks(dit, no)
      do i = nb, ind ! only runned if some iterations are missed
        ! Fake access
        call ncdf_put_var(ncdf, vname, a(1:1,1:1), start=(/1,1/), count=(/0,0/))
      end do

    end subroutine write_parallel

    subroutine write_parallel_single()
      if ( lgncol(1) < 0 ) then
        ! this means gather lgncol on IO-node
        call Node_Sp_gncol(0, sp, dit, no, lgncol)
      end if

      if ( sp_dim == 1 ) then
        call write_parallel_single_sp_dim1()
      else
        call write_parallel_single_sp_dim2()
      end if

    end subroutine write_parallel_single

    subroutine write_parallel_single_sp_dim1()
      integer :: BNode, MPIstatus(MPI_STATUS_SIZE), MPIerror
      integer :: gio, io, max_n, ind, gind, n, i, nb, is
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol)
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      ! Loop size
      do is = 1, dim2
        
        gio = 1
        gind = 1
        ind = 0
        nb = 0
        do while ( gio <= no )

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)

          if ( Node == BNode ) then

            io = index_global_to_local(dit,gio,Node)
            i = sum(ncol(io:io-1+n))

            if ( Node == 0 ) then
              call ncdf_put_var(ncdf, vname, a(ind+1:ind+i,is), &
                  count=(/i,1/), start=(/gind,is/))
              gind = gind + i
            else
              nb = nb + 1
              call MPI_ISend(a(ind+1,is), i, MPI_Double_Precision, &
                  0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
            end if
            ind = ind + i
            
          else if ( Node == 0 ) then
            call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
                BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
            if ( MPIerror /= MPI_Success ) &
                call die('Error in code: cdf_w_d2D')
          
            call MPI_Get_Count(MPIstatus, MPI_Double_Precision, i, MPIerror)
            call ncdf_put_var(ncdf, vname, buf(1:i), &
                count=(/i,1/), start=(/gind,is/))
            gind = gind + i
          end if
          
          gio = gio + n
        end do

        if ( Node /= 0 .and. nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if

      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        deallocate(ibuf)
      end if

    end subroutine write_parallel_single_sp_dim1

    subroutine write_parallel_single_sp_dim2()
      integer :: BNode, MPIstatus(MPI_STATUS_SIZE), MPIerror
      integer :: gio, io, max_n, ind, gind, n, i, nb
      real(dp), allocatable :: buf(:)
      integer, allocatable :: ibuf(:)

      nb = count_blocks(dit,no)

      if ( Node == 0 ) then
        max_n = max_consecutive_sum(dit,no,lgncol) * dim2
        allocate(buf(max_n))
      else
        allocate(ibuf(nb))
      end if

      gio = 1
      gind = 1
      ind = 0
      nb = 0
      do while ( gio <= no )

        BNode = node_handling_element(dit,gio)

        ! Get number of consecutive orbitals
        n = count_consecutive(dit,no,gio)
        
        if ( Node == BNode ) then

          io = index_global_to_local(dit,gio,Node)
          i = sum(ncol(io:io-1+n))

          if ( Node == 0 ) then
            call ncdf_put_var(ncdf, vname, a(:,ind+1:ind+i), &
                count=(/dim2,i/), start=(/1,gind/))
            gind = gind + i
          else
            nb = nb + 1
            call MPI_ISend(a(1,ind+1), i*dim2, MPI_Double_Precision, &
                0, gio, MPI_Comm_World, ibuf(nb), MPIerror)
          end if
          ind = ind + i

        else if ( Node == 0 ) then
          call MPI_Recv(buf(1), max_n, MPI_Double_Precision, &
              BNode, gio, MPI_Comm_World, MPIstatus, MPIerror)
          if ( MPIerror /= MPI_Success ) &
              call die('Error in code: cdf_w_d2D')

          call MPI_Get_Count(MPIstatus, MPI_Double_Precision, i, MPIerror)
          call ncdf_put_var(ncdf, vname, buf(1:i), &
              count=(/dim2,i/), start=(/1,gind/))
          gind = gind + i
        end if
        
        gio = gio + n
      end do

      if ( Node == 0 ) then
        deallocate(buf)
      else
        if ( nb > 0 ) then
          call MPI_WaitAll(nb, ibuf, MPI_STATUSES_IGNORE, MPIerror)
        end if
        deallocate(ibuf)
      end if

    end subroutine write_parallel_single_sp_dim2
#endif
    
    subroutine write_sequential()
      call ncdf_put_var(ncdf, vname, a)
    end subroutine write_sequential

  end subroutine cdf_w_d2D
  
  subroutine cdf_w_grid(ncdf,name,nmeshl,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: nmeshl(3)
    real(grid_p), intent(in) :: grid(:)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), allocatable :: gb(:)
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

      ! Ensure collective writing
      call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

      lb(:)  = distr%box(1,:,Node)
      nel(:) = distr%box(2,:,Node) - lb(:) + 1
      if ( .not. all(nel == nmeshl) ) then
        call die('cdf_r_grid: cannot assert the grid-size from the &
            &stored grid and the denoted size')
      end if

      if ( present(idx) ) then
        call ncdf_put_var(ncdf,name,grid, &
            start=(/lb(1),lb(2),lb(3),idx/), &
            count=(/nel(1),nel(2),nel(3),1/) )          
      else
        call ncdf_put_var(ncdf,name,grid, start=lb, count=nel )
      end if

    else

      mnpt = 0
      do iN = 0 , Nodes - 1
        nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
        mnpt = max(mnpt,product(nel))
      end do

      ! The main node can safely write the data...
      if ( Node == 0 ) then

        allocate(gb(mnpt))

        ! First save it's own data
        lb(:) = distr%box(1,:,0) 
        nel(:) = distr%box(2,:,0) - lb(:) + 1

        if ( present(idx) ) then
          call ncdf_put_var(ncdf,name,grid, &
              start=(/lb(1),lb(2),lb(3),idx/), &
              count=(/nel(1),nel(2),nel(3),1/) )
        else
          call ncdf_put_var(ncdf,name,grid, &
              start=lb, count=nel )
        end if

        ! Loop on the remaining nodes
        do iN = 1 , Nodes - 1

          ! we retrieve data from the iN'th node.
          call MPI_Recv(gb,mnpt,MPI_grid_real,iN,iN, &
              MPI_Comm_World,MPIstat, MPIerror)

          ! Just make sure we only pass the correct size
          call MPI_Get_Count(MPIstat, MPI_Grid_Real, inpt, MPIerror)

          lb(:) = distr%box(1,:,iN) 
          nel(:) = distr%box(2,:,iN) - lb(:) + 1
          if ( inpt /= product(nel) ) then
            call die('Error when receiving the distributed &
                &grid data for writing to the NetCDF file.')
          end if

          if ( present(idx) ) then
            call ncdf_put_var(ncdf,name,gb(1:inpt), &
                start=(/lb(1),lb(2),lb(3),idx/), &
                count=(/nel(1),nel(2),nel(3),1/) )
          else
            call ncdf_put_var(ncdf,name,gb(1:inpt), &
                start=lb, count=nel )
          end if

        end do

        deallocate(gb)
      else
        mnpt = product(nmeshl)
        call MPI_Send(grid,mnpt,MPI_grid_real,0, &
            Node,MPI_Comm_World,MPIerror)
      end if

    end if

#else
    if ( present(idx) ) then
      call ncdf_put_var(ncdf,name,grid,start=(/1,1,1,idx/), &
          count=nmeshl)
    else
      call ncdf_put_var(ncdf,name,grid, &
          count=nmeshl)
    end if
#endif

  end subroutine cdf_w_grid

  subroutine cdf_r_grid(ncdf,name,nmeshl,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: nmeshl(3)
    real(grid_p), intent(inout), target :: grid(:)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), pointer :: gb(:) => null()
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

      ! Ensure collective writing
      call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

      lb(:)  = distr%box(1,:,Node)
      nel(:) = distr%box(2,:,Node) - lb(:) + 1
      if ( .not. all(nel == nmeshl) ) then
        call die('cdf_r_grid: cannot assert the grid-size from the &
            &stored grid and the denoted size')
      end if

      if ( present(idx) ) then
        call ncdf_get_var(ncdf,name,grid, &
            start=(/lb(1),lb(2),lb(3),idx/), &
            count=(/nel(1),nel(2),nel(3),1/) )          
      else
        call ncdf_get_var(ncdf,name,grid, start=lb, count=nel )
      end if

    else

      mnpt = 0
      do iN = 0 , Nodes - 1
        nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
        mnpt = max(mnpt,product(nel))
      end do

      ! The main node can safely read the data...
      if ( Node == 0 ) then

        if ( mnpt > product(nmeshl) ) then
          ! allocate, the IO node have a smaller
          ! space than the others
          allocate(gb(mnpt))
        else
          gb => grid
        end if

        ! Loop on the remaining nodes (we may possibly
        ! reuse the current grid array, and hence
        ! we first read the other nodes...)
        do iN = 1 , Nodes - 1

          lb(:) = distr%box(1,:,iN) 
          nel(:) = distr%box(2,:,iN) - lb(:) + 1

          ! get total number of points
          inpt = product(nel)

          if ( present(idx) ) then
            call ncdf_get_var(ncdf,name,gb(1:inpt), &
                start=(/lb(1),lb(2),lb(3),idx/), &
                count=(/nel(1),nel(2),nel(3),1/) )
          else
            call ncdf_get_var(ncdf,name,gb(1:inpt), &
                start=lb, count=nel )
          end if

          ! we send data to the iN'th node.
          call MPI_Send(gb,inpt,MPI_grid_real,iN,iN, &
              MPI_Comm_World,MPIerror)

        end do

        if ( mnpt > product(nmeshl) ) then
          deallocate(gb)
        end if

        ! Retrieve it's own data
        lb(:) = distr%box(1,:,0) 
        nel(:) = distr%box(2,:,0) - lb(:) + 1

        if ( present(idx) ) then
          call ncdf_get_var(ncdf,name,grid, &
              start=(/lb(1),lb(2),lb(3),idx/), &
              count=(/nel(1),nel(2),nel(3),1/) )
        else
          call ncdf_get_var(ncdf,name,grid, &
              start=lb, count=nel )
        end if

      else
        mnpt = product(nmeshl)
        call MPI_Recv(grid,mnpt,MPI_grid_real,0, &
            Node,MPI_Comm_World,MPIstat,MPIerror)
      end if

    end if

#else
    if ( present(idx) ) then
      call ncdf_get_var(ncdf,name,grid,start=(/1,1,1,idx/), &
          count=nmeshl)
    else
      call ncdf_get_var(ncdf,name,grid,start=(/1,1,1/), &
          count=nmeshl)
    end if
#endif

  end subroutine cdf_r_grid

#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module ncdf_io_m
