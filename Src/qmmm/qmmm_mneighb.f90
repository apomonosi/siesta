!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
!
! Some modifications were introduced by c. Sanz-navarro (2009) in order
! to adapt it to the Hybrid code
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
module QMMM_neighbour
  !! This module deals with QM-MM neighbour lists.
  use precision, only : dp

  public :: update_QMMM_mneighb
  public :: QMMM_mneighb

  private
  integer , pointer, public :: QMMM_jan(:)
    !! Atom-index of QM-MM neighbours.
  real(dp), pointer, public :: QMMM_xij(:,:)
    !! Squared distances from atom ia to MM neighbours.

  logical :: QMMM_pointers_allocated = .false.
    !! Whether QMMM_* arrays are allocated.

  integer, public :: num_mmvecs
    !! Total number of MM neighbours.
  integer, pointer, public :: grid_veclist(:)
    !! Neighbour atom index in grid.
  integer, pointer, public :: grid_nr(:,:)
    !! Grid number of neighbour (i.e. the amount of cells we have to move).

  integer , public :: MAXNNA = 200
    !! The maximum number of neighbours for a given atom.
  real(dp) :: rcoor  = 0.0_dp
    !! Maximum distance for neighbours.
contains

  subroutine update_QMMM_mneighb( na_qm, na_mm, natot, r, qmcell, pc, rcut_QMMM )
    use alloc      , only : re_alloc, de_alloc
    use precision  , only : dp
    use QMMM_helper, only : get_pbc_vectors, vec_norm
    use sys        , only : die

    !! Updates QM-MM neighbour lists. All outputs are stored in this module.
    implicit none
    integer         , intent(in)    :: na_qm
      !! Total number of QM atoms.
    integer         , intent(in)    :: na_mm
      !! Total number of MM atoms.
    integer         , intent(in)    :: natot
      !! Total number of atoms (QM+MM).
    real(dp)        , intent(in)    :: rcut_QMMM
      !! Cut-off radius for QM-MM interactions.
    real(dp)        , intent(in)    :: r(3,natot)
      !! Atomic coordinates.
    real(dp)        , intent(in)    :: pc(na_mm)
      !! Partial charges of MM atoms.
    real(dp)        , intent(inout) :: qmcell(3,3)
      !! Unit cell vectors.

    real(dp) :: kcell(3,3)
    integer :: grid_nrx_max, grid_nry_max, grid_nrz_max, nx, ny, nz
    integer :: i, j, n, in
    integer, pointer :: pbc_condition(:,:,:,:)

    ! Neighbour list variables.
    integer  :: dimvec, grid_nr_indx(3), nna
    real(dp) :: drij(3)
    integer, pointer :: grid_veclistemp(:), grid_nrtemp(:,:)
    logical, pointer :: mm_present(:)

    dimvec = na_qm * 3000
    if ( dimvec > 100000000 ) &
      call die( 'Solvent Energy and Forces: "dimvec" too large!' )

    call reclat( qmcell, kcell, 0 )
    rcoor = rcut_QMMM

    ! Initialize routine for neighbour search. This call creates the cells.

    call QMMM_mneighb( qmcell, rcoor, natot, r, 0, 1, nna )

    call de_alloc( QMMM_jan , 'QMMM_jan' , 'update_QMMM_mneighb' )
    call de_alloc( QMMM_xij , 'QMMM_xij' , 'update_QMMM_mneighb' )
    nullify( QMMM_jan, QMMM_xij )
    call re_alloc( QMMM_jan , 1, MAXNNA, 'QMMM_jan' , 'update_QMMM_mneighb' )
    call re_alloc( QMMM_xij , 1, 3, 1, MAXNNA, 'QMMM_xij' , &
                   'update_QMMM_mneighb' )

    call de_alloc( grid_veclist   , 'grid_veclist'   , 'update_QMMM_mneighb' )
    call de_alloc( grid_nr        , 'grid_nr'        , 'update_QMMM_mneighb' )
    nullify( grid_veclist, grid_nr )

    nullify( grid_veclistemp, grid_nrtemp )
    call re_alloc( grid_veclistemp, 1, dimvec, 'grid_veclistemp', &
                   'update_QMMM_mneighb' )
    call re_alloc( grid_nrtemp    , 1, 3, 1, dimvec, 'grid_nrtemp', &
                   'update_QMMM_mneighb' )

    grid_veclistemp(:) = 0
    grid_nrtemp(:,:)   = 0

    grid_nrx_max = 1 + INT( 2.0_dp * rcoor / vec_norm( qmcell(:,1), 3 ) )
    grid_nry_max = 1 + INT( 2.0_dp * rcoor / vec_norm( qmcell(:,2), 3 ) )
    grid_nrz_max = 1 + INT( 2.0_dp * rcoor / vec_norm( qmcell(:,3), 3 ) )

    nullify( pbc_condition )
    allocate( pbc_condition( na_mm, -grid_nrx_max:grid_nrx_max, &
                                    -grid_nry_max:grid_nry_max, &
                                    -grid_nrz_max:grid_nrz_max) )

    pbc_condition(:,:,:,:) = 0

    n = 1
    do i = 1, na_qm
      ! Look for neighbours of atom ia
      call QMMM_mneighb( qmcell, rcoor, natot, r, i, 0, nna )
      if ( nna > MAXNNA ) call die( 'QMMM_mneighb: MAXNNA too small.' )

      do in = 1, nna
        j = QMMM_jan(in)
        if ( j <= na_qm ) cycle

        drij(1:3) = QMMM_xij(1:3,in) - r(1:3,j) + r(1:3,i)

        call get_pbc_vectors( qmcell, kcell, drij, grid_nr_indx )

        n = n +1
        if ( abs(pc(j-na_qm)) > 0.0_dp ) &
          pbc_condition( j-na_qm, grid_nr_indx(1), grid_nr_indx(2), &
                                  grid_nr_indx(3) ) = 1
      enddo ! End of loop over neighbour list for each atom.

      ! Verifies whether veclist dimension is sufficient.
      if ( (n-1) > dimvec ) then
        write(6,*)'Dimension Neighbour list (required, used)=', n-1, dimvec
        call die( 'update_QMMM_mneighb: Stopping Program' )
      endif
    enddo ! Loop over QM atoms.

    nullify( mm_present )
    call re_alloc( mm_present, 1, na_mm, 'mm_present' , 'update_QMMM_mneighb' )
    mm_present(:) = .false.

    n = 1
    do j = 1, na_mm
    do nx = -grid_nrx_max, grid_nrx_max
    do ny = -grid_nry_max, grid_nry_max
    do nz = -grid_nrz_max, grid_nrz_max
      if ( pbc_condition(j, nx, ny, nz) /= 1 ) cycle
      if ( mm_present(j) ) cycle

      grid_veclistemp(n) = j
      grid_nrtemp(1,n)   = nx
      grid_nrtemp(2,n)   = ny
      grid_nrtemp(3,n)   = nz
      n = n +1

      mm_present(j) = .true.
    enddo
    enddo
    enddo
    enddo
    deallocate( pbc_condition )
    nullify( pbc_condition )

    call de_alloc( mm_present, 'mm_present' , 'update_QMMM_mneighb' )

    ! Allocates the true veclist.
    call re_alloc( grid_veclist , 1, n-1, 'grid_veclist', &
                   'update_QMMM_mneighb' )
    call re_alloc( grid_nr, 1, 3, 1, n-1, 'grid_nr'     , &
                   'update_QMMM_mneighb' )
    grid_veclist(1:n-1) = grid_veclistemp(1:n-1)
    grid_nr(1:3,1:n-1)  = grid_nrtemp(1:3,1:n-1)
    call de_alloc( grid_veclistemp, 'grid_veclistemp', 'update_QMMM_mneighb' )
    call de_alloc( grid_nrtemp    , 'grid_nrtemp'    , 'update_QMMM_mneighb' )

    num_mmvecs = n -1
  end subroutine update_QMMM_mneighb

  subroutine QMMM_mneighb( cell, range, na, xa, ia, isc, nna )
    !! Finds the neighbours of an atom in a cell with periodic boundary
    !! conditions. This subroutine behaves exactly as mm_mneigh, but
    !! setting up the QM-MM neighbour arrays. See the other routine
    !! for more details.
    use alloc    , only : re_alloc, de_alloc
    use precision, only : dp

    implicit none
    integer , parameter :: nx = 3
      !! Space dimentionality.

    integer , intent(in)    :: ia
      !! Atom whose neighbours are needed. A routine initialization must be
      !! done by a first call with ia = 0.
    integer , intent(in)    :: isc
      !! Single-counting switch (0=No, 1=Yes). if isc=1, only neighbours
      !! with ja <= ia are included in mm_jan.
    integer , intent(in)    ::  na
      !! Total number of atoms.
    real(dp), intent(in)    :: range
      !! Maximum distance of neighbours required.
    real(dp), intent(in)    :: xa(nx,na)
      !! Atomic positions in cartesian coordinates.
    integer , intent(out)   :: nna
        !! Number of neighbour atoms within range of ia
    real(dp), intent(inout) :: cell(nx,nx)
      !! Unit cell vectors.

    logical , save :: frstme        = .true.
    integer , save :: iamove(1)     = 0
    real(dp), save :: celast(nx, nx) = 0.0_dp, &
                      rglast         = 0.0_dp, &
                      x0(nx)         = 0.0_dp

    logical :: samcel
    integer :: ix, jx

    if ( .not. QMMM_pointers_allocated ) then
      nullify( QMMM_jan, QMMM_xij )
      !  Dimension arrays to initial size MAXNNA
      call re_alloc( QMMM_jan ,       1, MAXNNA, 'QMMM_jan' , 'mm_mneighb' )
      call re_alloc( QMMM_xij , 1, 3, 1, MAXNNA, 'QMMM_xij' , 'mm_mneighb' )
      QMMM_pointers_allocated = .true.
    endif

    ! Initialization section
    if ( frstme .or. (ia <= 0) .or. (range > RGLAST) ) then
      ! Find if cell or range have changed
      samcel = .true.
      if ( frstme ) then
        samcel = .false.
      else
        do ix = 1, nx
        do jx = 1, nx
          if ( abs(cell(jx,ix) - celast(jx,ix)) > 0.0_dp ) then
            samcel = .false.
            exit
          endif
        enddo
        enddo
        if ( abs(range - rglast) > 0.0_dp ) samcel = .false.
      endif

      if ( .not. samcel ) then ! Cell initializations
        ! Store cell and range for comparison in subsequent calls
        celast(:,:) = cell(:,:)
        rglast      = range
        frstme      = .false.

        ! Notify to ranger that cell has changed
        call QMMM_mranger( 'CELL', nx, cell, range, na, xa, na, iamove, &
                           ia, isc, x0, nna, MAXNNA )
      endif

      ! Notify to ranger that atoms have moved
      call QMMM_mranger( 'MOVE', nx, cell, range, na, xa, na, iamove, &
                         ia, isc, x0, nna, MAXNNA )
    endif

    ! Find neighbours of atom ia
    if ( ia > 0 ) call QMMM_mranger( 'FIND', nx, cell, range, na, xa, na, &
                                     iamove, ia, isc, x0, nna, MAXNNA )
  end subroutine QMMM_mneighb

  subroutine QMMM_mranger( mode, nx, cell, range, na, xa, namove, iamove, &
                           ia0, isc, x0, nna, MAX_NNA )
    !! Finds the neighbours of an atom in a cell with periodic boundary
    !! conditions. Alternatively, it finds the atoms within a sphere
    !! centered at an arbitrary point. It also allows to update the atomic
    !! positions one at a time, what is useful in Montecarlo simulations.
    !! Works exactly the same as mm_ranger, so see the routine above for more
    !! details.
    !!
    !! Outputs stored in the module are: QMMM_jan(nna) and QMMM_xij(3,nna)
    !! which are only calculated if mode = 'FIND'.
    use alloc    , only : de_alloc, re_alloc
    use precision, only : dp

    implicit none
    character(len=4), intent(in)    :: mode
      !! 'CELL' => Initialize or reshape cell; 'MOVE' => Move atom(s);
      !! 'FIND' => Find neighbours.
    integer         , intent(in)    :: nx
      !! Space dimension.
    real(dp)        , intent(in)    :: range
      !! Maximum distance of neighbours required.
    integer         , intent(in)    :: na
      !! Total number of atoms.
    real(dp)        , intent(in)    :: xa(nx, na)
      !! Atomic positions in cartesian coordinates.
    integer         , intent(in)    :: namove
      !! Number of atoms to be moved when mode='MOVE'.
    integer         , intent(in)    :: iamove(namove)
      !! Indices of atoms to be moved when mode='MOVE' and namove < na.
    integer         , intent(in)    :: ia0
      !! Atom whose neighbours are needed. if ia0 = 0, X0 is used as origin.
      !! This last case only happens when mode='FIND'.
    integer         , intent(in)    :: isc
      !! Single-counting switch. When isc = 1, only neighbours with ja <= ia0
      !! are included; it is not used unless mode='FIND' and ia0 /= 0.
    real(dp)        , intent(in)    :: x0(nx)
      !! Origin from which atoms are to be found, when mode='FIND' and ia0 = 0.
    integer         , intent(out)   :: nna
      !! Number of 'neighbour' atoms within range of ia0/x0, when mode='FIND'.
    real(dp)        , intent(inout) :: cell(*)
      !! Unit cell vectors, in truth (nx, nVec).
    integer         , intent(inout) :: MAX_NNA
      !! Size of arrays QMMM_jan, _xij and _r2ij, used when mode='FIND'.

    ! NCR is the ratio between range radius and mesh-planes distance. It fixes
    ! the size (and number) of mesh cells. Recommended values are between 1-3.
    integer, parameter :: ncr = 2

    ! dxmarg and dxrang are used for automatic cell generation. dxmarg is the
    ! minimum margin relative to coordinate range; dxrang is the minimum margin
    ! relative to range and EPS is a small number to be subtracted from 1.
    real(dp), parameter :: dxmarg = 0.1_dp
    real(dp), parameter :: dxrang = 1.0_dp
    real(dp), parameter :: EPS    = 1.0e-14_dp

    ! Internals
    integer  :: ia , ja,       & ! Atom index.
                jm ,           & ! Mesh index
                jem,           & ! Extended-mesh index
                in ,           & ! Neighbour-mesh-cell index
                ix , jx , ixx, & ! Cartesian / double-cartesian coordinate index
                nam,           & ! Number of atoms to move
                j_aux,         & ! Auxiliary to avoid compiler warnings
                nnmmax           ! Maximum number of neighbour mesh cells

    logical  :: movall, & ! Move all atoms?
                nulcel    ! Null cell?

    real(dp) :: dplane, &    ! Distance between lattice or mesh planes
                r2    , &    ! Squared distance between atoms
                xdiff , &    ! Range of atom coordinates
                xmarg , &    ! Distance to margin
                xmax  , xmin ! Min/Max atom coordinates

    ! Arrays.
    integer , pointer ::  &
      inx(:)   ,          & ! (nx) Neighbour-cell-coordinate indices
      i1nx(:)  , i2nx(:), & ! (nx) Min/Max neighbour-cell-coordinate indices
      j1nx(:)  , j2nx(:)    ! (nx) Min/Max vertex-cell-coordinate indices

    real(dp), pointer ::  &
      dmx(:)   , & ! (nx) In-cell atomic position in mesh coordinates
      dx(:)    , & ! (nx) Vector between two atoms
      dx0m(:)      ! (nx) Origin position within mesh cell

    ! Saved values.
    integer , save :: maxna = 0, & ! Maximum number of atoms.
                      iam,       & ! Atom-to-move index
                      iem,       & ! Extended-mesh index
                      im ,       & ! Mesh index
                      nm ,       & ! Number of mesh cells
                      nem,       & ! Number of extended-mesh cells
                      nnm          ! Number of neighbour mesh cells

    logical , save :: frstme = .true. ! First time calling

    real(dp), save :: range2,    & ! Squared range
                      RNGMAX,    & ! Maximum range
                      rrange       ! Slightly reduced range

    ! Saved arrays.
    integer, pointer, save :: &
      ianext(:), iaprev(:), & ! (na) Pointers to next/previous atoms in cell
      iema(:),              & ! (na) Extended-mesh index of atoms
      i1emx(:) , i2emx(:) , & ! (nx) Min/Max value of extended-mesh-crd indices
      imx(:),               & ! (nx) Mesh-cell index for each mesh vector
      i1mx(:)  , i2mx(:)  , & ! (nx)  Min/Max value of mesh-coordinate indices
      nemx(:)  ,            & ! (nx) Extended-mesh cells in each mesh direction
      nmx(:)   ,            & ! (nx) Mesh cells in each mesh direction
      nnx(:)   ,            & ! (nx) Neighbour-cell ranges
      ia1m(:)  ,            & ! (nm) Pointer to first atom in mesh cell
      imesh(:) ,            & ! (nem) Extended to normal mesh correspondence
      idnm(:)                 ! (maxnm) Index-distance between neighbour
                              !         mesh points

    real(dp), pointer, save :: &
      celmsh(:), & ! (nx*nx) Mesh-cell vectors
      rcell(:) , & ! (nx*nx) Reciprocal cell vectors
      rmcell(:), & ! (nx*nx) Reciprocal mesh-cell vectors
      dxam(:,:), & ! (nx,na) Atom position within mesh cell
      dxnm(:,:)    ! (nx,maxnm) Cartesian vector between neighbour mesh points

    ! Externals.
    real(dp), external :: DISMin, DDOT

    ! Allocate local memory - check for change in number of atoms and if
    ! there has been one then re-initialise.
    if ( frstme ) then
      nullify( i1emx, i2emx, imx, i1mx, i2mx, nemx, nmx, nnx )
      call re_alloc( i1emx , 1, nx   , 'i1emx' , 'QMMM_mranger' )
      call re_alloc( i2emx , 1, nx   , 'i2emx' , 'QMMM_mranger' )
      call re_alloc( imx   , 1, nx   , 'imx'   , 'QMMM_mranger' )
      call re_alloc( i1mx  , 1, nx   , 'i1mx'  , 'QMMM_mranger' )
      call re_alloc( i2mx  , 1, nx   , 'i2mx'  , 'QMMM_mranger' )
      call re_alloc( nemx  , 1, nx   , 'nemx'  , 'QMMM_mranger' )
      call re_alloc( nmx   , 1, nx   , 'nmx'   , 'QMMM_mranger' )
      call re_alloc( nnx   , 1, nx   , 'nnx'   , 'QMMM_mranger' )

      nullify( celmsh, rcell, rmcell )
      call re_alloc( celmsh, 1, nx*nx, 'celmsh', 'QMMM_mranger' )
      call re_alloc( rcell , 1, nx*nx, 'rcell' , 'QMMM_mranger' )
      call re_alloc( rmcell, 1, nx*nx, 'rmcell', 'QMMM_mranger' )
    endif

    if ( na > maxna ) then
      call de_alloc( ianext, 'ianext', 'QMMM_mranger' )
      call de_alloc( iaprev, 'iaprev', 'QMMM_mranger' )
      call de_alloc( iema  , 'iema'  , 'QMMM_mranger' )
      call de_alloc( dxam  , 'dxam'  , 'QMMM_mranger' )
      nullify( ianext, iaprev, iema, dxam )

      call re_alloc( ianext,        1, na, 'ianext', 'QMMM_mranger' )
      call re_alloc( iaprev,        1, na, 'iaprev', 'QMMM_mranger' )
      call re_alloc( iema  ,        1, na, 'iema'  , 'QMMM_mranger' )
      call re_alloc( dxam  , 1, nx, 1, na, 'dxam'  , 'QMMM_mranger' )

      maxna  = na
      frstme = .false.
    endif

    ! The follwing are only used locally.
    nullify( inx, i1nx, i2nx, j1nx, j2nx, dmx, dx, dx0m )
    call re_alloc( inx , 1, nx, 'inx' , 'QMMM_mranger' )
    call re_alloc( i1nx, 1, nx, 'i1nx', 'QMMM_mranger' )
    call re_alloc( i2nx, 1, nx, 'i2nx', 'QMMM_mranger' )
    call re_alloc( j1nx, 1, nx, 'j1nx', 'QMMM_mranger' )
    call re_alloc( j2nx, 1, nx, 'j2nx', 'QMMM_mranger' )
    call re_alloc( dmx , 1, nx, 'dmx' , 'QMMM_mranger' )
    call re_alloc( dx  , 1, nx, 'dx'  , 'QMMM_mranger' )
    call re_alloc( dx0m, 1, nx, 'dx0m', 'QMMM_mranger' )

    ! Cell-mesh initialization section
    if ( (mode == 'CELL') .or. (mode == 'cell') .or. frstme .or. &
         (range > rngmax) ) then
      rngmax = range ! Store range for comparison in subsequent calls

      ! Reduce the range slitghtly to avoid numerical-roundoff ambiguities
      rrange = range * (1.0_dp - EPS)
      range2 = rrange * rrange

      ! Check if cell must be generated automatically
      nulcel = .true.
      do ixx = 1, nx*nx
        if ( abs( cell(ixx) ) > 0.0_dp ) then
          nulcel = .false.
          exit
        endif
      enddo

      if ( nulcel ) then
        do ix = 1, nx
          ! Find atom position bounds
          xmin =  1.e30_dp
          xmax = -1.e30_dp
          do ia = 1, na
            xmin = Min( xmin, xa(ix,ia) )
            xmax = MAX( xmax, xa(ix,ia) )
          enddo

          ! Determine 'cell margins' to prevent intercell interactions
          xdiff = xmax - xmin
          xmarg = MAX( range * dxrang, xdiff * dxmarg )

          ! Define orthorrombic cell
          ixx       = ix + nx * ( ix -1 )
          cell(ixx) = xdiff + 2.0_dp * xmarg
        enddo
      endif

      ! Find reciprocal cell vectors (not multiplied by 2*pi)
      call reclat( cell, rcell, 0 )

      ! Find number of mesh divisions
      nm = 1
      do ix = 1, nx
        ixx     = 1 + nx * ( ix -1 )
        dplane  = 1.0_dp / SQRT( DDOT( nx, rcell(ixx:ixx+nx-1), 1, &
                                           rcell(ixx:ixx+nx-1), 1 ) )
        nmx(ix) = NCR * INT ( 0.999_dp * dplane / rrange )
        if ( nmx(ix) <= 0 ) nmx(ix) = 1
        nm = nm * nmx(ix)
      enddo

      ! Find mesh-cell vectors
      ixx = 0
      do ix = 1, nx
      do jx = 1, nx
        ixx = ixx + 1
        celmsh(ixx) = cell(ixx)  / nmx(ix)
        rmcell(ixx) = rcell(ixx) * nmx(ix)
      enddo
      enddo

      ! Find index-range of neighbour mesh cells and of extended mesh
      nnm = 1
      nem = 1
      do ix = 1, nx
        ixx       = 1 + nx * ( ix -1 )
        dplane    = 1.0_dp / SQRT( DDOT(nx, rmcell(ixx:ixx+nx-1), 1, &
                                            rmcell(ixx:ixx+nx-1), 1) )
        nnx(ix)   = INT( rrange / dplane ) +1
        j1nx(ix)  = 0
        j2nx(ix)  = 1
        i1nx(ix)  = -nnx(ix)
        i2nx(ix)  =  nnx(ix)
        i1mx(ix)  = 0
        i2mx(ix)  = nmx(ix) - 1
        i1emx(ix) = -nnx(ix)
        i2emx(ix) =  nmx(ix) + nnx(ix) - 1
        nemx(ix)  =  nmx(ix) + 2 * nnx(ix)

        nnm = nnm * ( 1 + 2 * nnx(ix) )
        nem = nem * nemx(ix)
      enddo

      ! Allocate arrays whose dimensions are now known
      call de_alloc( ia1m , 'ia1m' , 'QMMM_mranger' )
      call de_alloc( idnm , 'idnm' , 'QMMM_mranger' )
      call de_alloc( dxnm , 'dxnm' , 'QMMM_mranger' )
      call de_alloc( imesh, 'imesh', 'QMMM_mranger' )
      nullify( ia1m, idnm, dxnm, imesh )
      call re_alloc( ia1m ,        1, nm , 'ia1m' , 'QMMM_mranger' )
      call re_alloc( idnm ,        1, nnm, 'idnm' , 'QMMM_mranger' )
      call re_alloc( dxnm , 1, nx, 1, nnm, 'dxnm' , 'QMMM_mranger' )
      call re_alloc( imesh,        1, nem, 'imesh', 'QMMM_mranger' )

      ! Find which mesh cells are actually within range
      nnmmax = nnm
      nnm    = 0

      do in = 1, nnmmax
        j_aux = in
        call indarr( -1, nx, i1nx, i2nx, inx, 1, j_aux )
        nnm = nnm + 1

        ! idnm is the extended-mesh-index distance between
        ! neighbour mesh cells
        idnm(nnm) = inx(nx)
        do ix = nx-1, 1, -1
          idnm(nnm) = inx(ix) + nemx(ix) * idnm(nnm)
        enddo

        ! dxnm is the vector distance between neighbour mesh cells
        do ix = 1, nx
          dxnm(ix,nnm) = 0.0_dp

          do jx = 1, nx
            ixx          = ix + nx * ( jx -1 )
            dxnm(ix,nnm) = dxnm(ix,nnm) + celmsh(ixx) * inx(jx)
          enddo
        enddo
      enddo

      ! Find correspondence between extended and reduced (normal) meshes
      do iem = 1,nem
        j_aux = iem
        call indarr( -1, nx, i1emx, i2emx, imx, 1, j_aux )
        call indarr( +1, nx, i1mx,  i2mx,  imx, 1, im  )
        imesh(iem) = im
      enddo

      movall = .true. ! We should move atoms.
    else
      movall = .false.
    endif ! mode = cell
    ! End of cell initialization section

    ! Atom-positions (relative to mesh) initialization section
    if ( (mode == 'MOVE') .or. (mode == 'move') .or. movall) then
      if ( namove == na ) movall = .true.

      if ( movall ) then ! Initialize 'atoms in mesh-cell' lists
        nam = na
        ianext(:) = 0
        iaprev(:) = 0
        ia1m(:)   = 0
      else
        nam = namove
      endif

      do iam = 1, nam ! Loop on moved atoms
        ! Select atoms to move
        if (movall) then
          ia = iam
        else
          ia = iamove(iam)
          ja = iaprev(ia) ! Supress atom from its previous mesh-cell
          if ( ja /= 0 ) ianext(ja) = ianext(ia)
          ja = ianext(ia)
          if ( ja /= 0 ) iaprev(ja) = iaprev(ia)

          iem = iema(ia)
          im  = imesh(iem)
          if ( ia1m(im) == ia ) ia1m(im) = ja
        endif

        ! Find mesh-cell in which atom is.
        do ix = 1, nx
          ixx = 1 + nx * (ix-1)
          dmx(ix) = DDOT( nx, rmcell(ixx:ixx+nx-1), 1, xa(1:nx,ia), 1)
          imx(ix) = INT( dmx(ix) + 1000.0_dp ) - 1000
          dmx(ix) = dmx(ix) - imx(ix)
          imx(ix) = MOD( imx(ix) + 1000 * nmx(ix), nmx(ix) )
        enddo
        call indarr( +1, nx, i1emx, i2emx, imx, 1, iem )
        call indarr( +1, nx, i1mx , i2mx , imx, 1, im  )
        iema(ia) = iem

        ! Put atom first in its new mesh-cell
        ja = ia1m(im)
        if ( ja /= 0 ) iaprev(ja) = ia
        ianext(ia) = ja
        ia1m(im)   = ia

        ! Find atomic position relative to mesh
        do ix = 1, nx
          dxam(ix,ia) = 0.0_dp
          do jx = 1, nx
            ixx         = ix + nx * ( jx -1 )
            dxam(ix,ia) = dxam(ix,ia) + celmsh(ixx) * dmx(jx)
          enddo
        enddo
      enddo
    endif ! mode = move
    ! End of atom-positions initialization section

    ! Proper search section
    if ( (mode == 'FIND') .or. (mode == 'find') ) then
      rrange = range * ( 1.0_dp - EPS )
      range2 = rrange * rrange

      ! Find the mesh cell of the center of the sphere
      if ( ia0 <= 0 ) then
        ! Find mesh cell of position X0
        do ix = 1, nx
          ixx = 1 + nx * ( ix -1 )
          dmx(ix) = DDOT( nx, rmcell(ixx:ixx+nx-1), 1, X0, 1 )
          imx(ix) = INT( dmx(ix) + 1000.0_dp ) - 1000
          dmx(ix) = dmx(ix) - imx(ix)
          imx(ix) = MOD( imx(ix) + 1000 * nmx(ix), nmx(ix) )
        enddo
        call indarr( +1, nx, i1emx, i2emx, imx, 1, iem )

        do ix = 1, nx
          dx0m(ix) = 0.0_dp
          do jx = 1, nx
            ixx = ix + nx * ( jx -1 )
            dx0m(ix) = dx0m(ix) + celmsh(ixx) * dmx(jx)
          enddo
        enddo
      else
        ! Find mesh cell of atom ia0
        iem = iema(ia0)
        do ix = 1, nx
          dx0m(ix) = dxam(ix,ia0)
        enddo
      endif

      ! Loop on neighbour mesh cells and on the atoms within them
      ! This is usually the only time-consuming loop
      nna = 0
      do in = 1,nnm
        jem = iem + idnm(in)
        jm  = imesh(jem)

        ! Loop on atoms of neighbour cell. Try first atom in this mesh-cell.
        ja = ia1m(jm)

        do while ( ja /= 0 )
          ! Check that single-counting exclusion does not apply
          if ( (ia0 <= 0) .or. (isc == 0) .or. (ja <= ia0) ) then
            ! Find vector and distance to atom ja
            r2 = 0.0_dp
            do ix = 1, nx
              dx(ix) = dxnm(ix,in) + dxam(ix,ja) - dx0m(ix)
              r2     = r2 + dx(ix) * dx(ix)
            enddo

            ! Check if atom ja is within range
            if ( r2 <= range2 ) then
              nna = nna + 1
              ! Check that array arguments are not overflooded
              if ( nna > MAX_NNA ) then
                MAX_NNA = MAX_NNA + na
                call re_alloc( QMMM_jan ,        1, MAX_NNA, 'QMMM_jan', &
                              'QMMM_mranger' )
                call re_alloc( QMMM_xij , 1, nx, 1, MAX_NNA, 'QMMM_jan', &
                              'QMMM_mranger')
              endif

              QMMM_jan(nna)    = ja
              QMMM_xij(:, nna) = dx(:)
            endif
          endif
          ! Take next atom in this mesh-cell and go to begining of loop
          ja = ianext(ja)
        enddo ! ja /= 0
      enddo   ! in
    endif     ! mode = FIND
    ! End of search section

    call de_alloc( inx , 'inx' , 'QMMM_mranger' )
    call de_alloc( i1nx, 'i1nx', 'QMMM_mranger' )
    call de_alloc( i2nx, 'i2nx', 'QMMM_mranger' )
    call de_alloc( j1nx, 'j1nx', 'QMMM_mranger' )
    call de_alloc( j2nx, 'j2nx', 'QMMM_mranger' )
    call de_alloc( dmx , 'dmx' , 'QMMM_mranger' )
    call de_alloc( dx  , 'dx'  , 'QMMM_mranger' )
    call de_alloc( dx0m, 'dx0m', 'QMMM_mranger' )

    frstme = .false.
  end subroutine QMMM_mranger

  subroutine indarr( IOPT, ND, Imin, Imax, Iind, Jmin, Jind )
    !! Finds the global index in a multidimensional array from the indexes
    !! in each dimension, or viceversa (the first is an explicit solution
    !! of the standard index-resolution problem that the compiler solves
    !! each time an array element is referenced).
    ! Written by J.M.Soler. 1996.
    !
    ! Indexes I() are taken as periodic, i.e. their modulus Imax(ID)-Imin(ID)+1
    ! is taken before using them. This simplifies its use as indexes of a
    ! mesh with periodic boundary conditions. This modulus operation is
    ! also done with J, so that the output I() are always within range.
    ! if IOPT=0, nothing is done.
    !
    ! Sample usage to find the Laplacian of a function defined in a mesh
    ! with periodic boundary conditions in a space of variable dimension.
    !
    ! subroutine laplacian( N, dx, F, FLAPL )
    !   real(dp) :: F(*), FLAPL(*), dx(3)
    !   integer  :: N(3), I1(3), I2(3), I(3)
    !
    !   nmesh = N(1) * N(2) * N(3)
    !   I1(:) = 1; I2(:) = N(:)
    !   do imesh = 1, nmesh
    !     call indarr( -1, 3, I1, I2, I, 1, imesh )
    !     FLAPL(imesh) = 0.0
    !     do ID = 1,3
    !       do K = -1,1,2
    !         I(ID) = I(ID) + K
    !         call indarr( +1, 3, I1, I2, I, 1, jmesh )
    !         FLAPL(imesh) = FLAPL(imesh) + F(jmesh) / dx(ID)**2
    !         I(ID) = I(ID) - K
    !       enddo
    !       FLAPL(imesh) = FLAPL(imesh) - 2 * F(imesh) / dx(ID)**2
    !     enddo
    !   enddo
    ! end subroutine laplacian

    implicit none
    integer, intent(in) :: iopt
      !! (IOPT > 0) => from I to J. (IOPT < 0) => from J to I.
    integer, intent(in) :: nd
      !! Number of array dimensions.
    integer, intent(in) :: Imin(nd)
      !! Minimum value of array indices.
    integer, intent(in) :: Imax(nd)
      !! Maximum value of array indices.
    integer, intent(in) :: Jmin
      !! Minimum value of the global index (usually 1).
    integer, intent(inout) :: Iind(nd)
      !! Array indices in each dimension.
    integer, intent(inout) :: Jind
      !! Global index.

    integer :: ID, Ktmp, Ntmp

    if ( IOPT > 0 ) then
      Jind = 0
      do ID = ND, 1, -1
        Ntmp = Imax(ID) - Imin(ID) + 1
        Ktmp = Iind(ID)  - Imin(ID)
        Ktmp = MOD( Ktmp + 1000 * Ntmp, Ntmp )
        Jind = Ktmp + Ntmp * Jind
      enddo
      Jind = Jind + Jmin
    elseif (IOPT < 0) then
      Ktmp = Jind - Jmin
      do ID = 1,ND
        Ntmp     = Imax(ID) - Imin(ID) + 1
        Iind(ID) = Imin(ID) + MOD( Ktmp, Ntmp )
        Ktmp     = Ktmp / Ntmp
      enddo
    endif
  end subroutine indarr
end module QMMM_neighbour

