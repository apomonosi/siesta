program lindhard_function
  !! This program computes the Lindhard function for a given set of k-points and
  !! eigenvalues.
  !!
  !! Written by P. Ordejón, 2014. Updated by B. Guster.
  !! See DOI: 10.1088/1361-648X/ab8522

  use files            , only : slabel, label_length
  use lindhard_kgrid_m , only : setup_kgrid, evaluate, kpoint_crd
  use lindhard_reinit_m, only : lindhard_reinit
  use lindhard_tools_m , only : readlindhard, inver3, readeig, fermif
  use m_io             , only : io_assign, io_close
  use precision        , only : dp
  use sys              , only : die
#ifdef MPI
  use mpi
#endif

  implicit none
  character(len=150)                  :: sname
    !! System name, used to initialise read
  character(len=label_length+4), save :: fname

  real(dp), parameter :: tiny = 1.0e-8_dp
  integer  :: mscell(3,3), iu, nkps, no, ns, nk, nkp, nkpoints, myid, numprocs,&
              ierr, nperproc, ix, iy, iz, ivx, ivy, ivz, ikx, iky, ikz, &
              firstband, lastband, nblind, ikqx, ikqy, ikqz, io, jo, nkx, nky, &
              nkz, nq1, nq2, ngridx, ngridy, ngridz, icrd, bsize, ng(3), igmin(3)
  logical  :: lstop
  real(dp) :: alat, ucell(3,3), scell(3,3), kv(3,3), kvi(3,3), kvmin(3), &
              k(3), kshift(3), ef, x(3), ek, ekq, q0(3), q1(3), kdist,&
              q2(3), q(3), susc, mysusc, ocup,dif, temp, gridk(3,3), dkg(3)

  real(dp), allocatable :: ks(:,:), es(:,:,:), e(:,:,:,:), &
                           ei(:,:,:,:), fe(:,:,:,:)

#ifdef MPI
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
#else
  myid = 0
  numprocs = 1
#endif

  if ( myid == 0 ) write(6,*) 'Total number of MPI processes = ',numprocs

  ! Find out SIESTA grid
  if ( myid == 0 ) then
    call lindhard_reinit( sname )
    call redcel( alat, ucell, scell, mscell )
    call setup_kgrid( ucell, gridk, dkg, igmin, ng )

    ! gridk(i,j) = Vectors defining the k-grid parallelepiped (i coordinate of j vector)
    ! ng(i) = number of grid divisions of the i reciprocal lattice vector (RLV)
    ! kv(i,j) = i cartessian component of j RLV
    do icrd = 1, 3
      kv(:,icrd) = gridk(:,icrd) * ng(icrd)
    enddo

    ! kvmin = coordinates in k-space of corner of the RLV grid
    do icrd = 1,3
      kvmin(icrd) = gridk(icrd,1) * (igmin(1) + dkg(1)) + &
                    gridk(icrd,2) * (igmin(2) + dkg(2)) + &
                    gridk(icrd,3) * (igmin(3) + dkg(3))
    enddo

    nkx = ng(1)
    nky = ng(2)
    nkz = ng(3)
    nkpoints = nkx*nky*nkz
    call inver3( kv, kvi )

    ! Read SIESTA KP file (to compare to grid computed here).  Only main node.
    fname = slabel
    fname = trim(fname) // '.KP'

    call io_assign( iu )
    open( iu, file=fname, form='formatted', status='unknown' )

    read( iu,* ) nkps
    if ( nkps /= nkpoints ) then
      write(6,*) 'Wrong number of kpoints in KP file', nkps, nkpoints
      call die('Error: incorrect number of k-points in KP file.')
    endif

    allocate ( ks(nkps,3) )
    do ikx = 1, nkps
      ! The first component is just a dummy.
      read(iu,*) iky, (ks(ikx,icrd),icrd=1,3)
    enddo
    call io_close( iu )
  endif

  ! Communicate number of grid points to all the nodes
#ifdef MPI
  call MPI_BCAST(nkx, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST(nky, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST(nkz, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
  nkpoints = nkx*nky*nkz
#endif

  ! Read eigenvalues from EIG file. Only main node.
  if ( myid == 0 ) then
    allocate( es(1,1,1) )
    no=1
    ns=1
    nk=1

    call readeig( es, ef, no, ns, nk, .true. )
    deallocate( es )
    allocate( es(no,ns,nk) )
    call readeig( es, ef, no, ns, nk, .false. )

    if ( nk /= nkpoints ) then
      write(6,*) 'Wrong number of kpoints', nk, nkpoints
      call die('Error: incorrect number of k-points in EIG file.')
    endif
  endif

  ! Communicate number of orbitals and fermi energy to all the nodes.
#ifdef MPI
  call MPI_BCAST( no, 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ef, 1, MPI_double_precision, 0, MPI_COMM_WORLD, ierr )
#endif

  if ( myid == 0 ) then
    call readlindhard( temp, q0, q1, q2, ngridx, ngridy, ngridz, nq1, nq2, &
                       firstband, lastband )
    if ( ngridx == 0 ) ngridx = nkx
    if ( ngridy == 0 ) ngridy = nky
    if ( ngridz == 0 ) ngridz = nkz
    if ( firstband == 0 ) firstband = 1
    if ( lastband  == 0 ) lastband  = no

    write(6,*) 'Lindhard calculation parameters:'
    write(6,*) '  Temperature = ', temp, " K"
    write(6,*) '  firstband   = ', firstband
    write(6,*) '  lastband    = ', lastband
    write(6,*) '  Grid Size   = ', ngridx, ngridy, ngridz
    write(6,*) '  nq1, nq2    = ', nq1, nq2

    if ( (ngridx / 2) * 2 /= ngridx ) call die('ERROR: ngridx must be even.')
    if ( (ngridy / 2) * 2 /= ngridy ) call die('ERROR: ngridy must be even.')
    if ( (ngridz / 2) * 2 /= ngridz ) call die('ERROR: ngridz must be even.')

    if ( (ngridx / nq1) * nq1 /= ngridx ) &
      call die('ERROR: ngridx must be multiple of nq1')
    if ( (ngridy / nq2) * nq2 /= ngridy ) &
      call die('ERROR: ngridy must be multiple of nq2')
  endif

  ! We paralelize splitting k-points over nodes in the x-direction.
  ! Could be revised for further performance optimization.

  lstop = .false.
  if ( myid == 0 ) then
    nperproc = ngridx / numprocs
    write(6,*) "Number of grid points per processor:", nperproc
    if ( nperproc*numprocs /= ngridx ) then
      lstop = .true.
      write(6,*) 'Number of MPI processes must be a multiple of the number'//&
                 ' of grid points in the x direction'
    endif
  endif

#ifdef MPI
  call MPI_BCAST( lstop, 1, MPI_logical, 0, MPI_COMM_WORLD, ierr )
#endif
  if ( lstop ) &
    call die('Error: Mismatch between number of MPI processes and grid points')

#ifdef MPI
  call MPI_BCAST( temp     , 1, MPI_double_precision, 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( firstband, 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( lastband , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ngridx   , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ngridy   , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ngridz   , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nq1      , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nq2      , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nperproc , 1, MPI_integer         , 0, MPI_COMM_WORLD, ierr )
#endif

  nblind = lastband - firstband +1

  if ( myid == 0 ) then
    allocate ( e(nkx,nky,nkz,nblind) )

    ! Save energies and check that generated grid coincides with Siesta grid
    nkp = 0
    do iz = 1, nkz
    do iy = 1, nky
    do ix = 1, nkx

      nkp = nkp +1
      jo  = 0
      do io = firstband, lastband
        jo = jo +1

        ! Beware: we work on a single spin channel.
        e(ix,iy,iz,jo) = es(io,1,nkp) - ef
      enddo

      ! Compute coordinates of k point to check if k-point grid is computed
      ! correctly.
      call kpoint_crd( kv, kvmin, nkx, nky, nkz, ix, iy, iz, k )

      lstop = .true.
      kdist = (k(1) - ks(nkp,1))**2 + (k(2) - ks(nkp,2))**2 + &
              (k(3) - ks(nkp,3))**2
      if ( kdist > tiny ) then
        do ivx = 0,1
        do ivy = 0,1
        do ivz = 0,1
          do icrd = 1, 3
            kshift(icrd) = ks(nkp,icrd) + ivx*ng(1)*gridk(icrd,1) + &
                           ivy*ng(2)*gridk(icrd,2) + ivz*ng(3)*gridk(icrd,3)
          enddo
          kdist = (k(1) - kshift(1))**2 + (k(2) - kshift(2))**2 + &
                  (k(3) - kshift(3))**2
          if ( kdist < tiny ) lstop = .false.
          if ( .not. lstop ) exit
        enddo
          if ( .not. lstop ) exit
        enddo
          if ( .not. lstop ) exit
        enddo

        if ( lstop ) then
          write(6,*) 'Mismatch in kpoint grid.'
          write(6,*) 'nkp = ', nkp
          do icrd = 1, 3
            write(6,'(i4,3f10.5)') icrd, k(icrd), ks(nkp,icrd)
          enddo
          call die( 'Error: k-point grid does not match k-points.' )
        endif
      endif
    enddo
    enddo
    enddo

    deallocate(ks)
    deallocate(es)
  endif

  ! Interpolate energies on grid points and save in array ei
  ! Compute and save also Fermi occupation, to avoid multiple computing
  ! The origin of the resulting interpolation grid is at the first k-point
  ! of the input grid.

  allocate ( ei(nblind, 0:ngridx-1, 0:ngridy-1, 0:ngridz-1) )
  allocate ( fe(nblind, 0:ngridx-1, 0:ngridy-1, 0:ngridz-1) )

  if ( myid == 0 ) then
    do ikz = 0, ngridz-1
      x(3) = ikz / real(ngridz,kind=dp)
    do iky = 0, ngridy-1
      x(2) = iky  /real(ngridy,kind=dp)
    do ikx = 0, ngridx-1
      x(1) = ikx / real(ngridx,kind=dp)

      do io=1,nblind
        call evaluate( e(1,1,1,io), nkx, nky, nkz, x, ek )
        ei(io,ikx,iky,ikz) = ek
        fe(io,ikx,iky,ikz) = fermif( ek, temp )
      enddo
    enddo
    enddo
    enddo

    deallocate( e )
  endif

  ! Communicate ei and fe to all the nodes
#ifdef MPI
  bsize = nblind*ngridx*ngridy*ngridz
  call MPI_BCAST( ei(1,0,0,0), bsize, MPI_double_precision, 0, MPI_COMM_WORLD, &
                  ierr )
  call MPI_BCAST( fe(1,0,0,0), bsize, MPI_double_precision, 0, MPI_COMM_WORLD, &
                  ierr )
#endif

  ! Open output file.  Only main node.
  if ( myid == 0 ) then
    fname = slabel
    fname = trim(fname)//'.LINDHARD'

    call io_assign( iu )
    open( iu, file=fname, form='formatted', status='unknown' )
  endif

  ! Loop over q points (only qz so far)
  do ix = -ngridx / 2, ngridx / 2, nq1
  do iy = -ngridy / 2, ngridy / 2, nq2

    q(1) = ix / real(ngridx,kind=dp)
    q(2) = iy / real(ngridy,kind=dp)
    q(3) = 0.0_dp

    mysusc = 0.0_dp
    susc   = 0.0_dp

    do ikx = myid*nperproc, (myid+1)*nperproc-1
      ikqx = modulo( ix+ikx, ngridx )
    do iky = 0, ngridy-1
      ikqy = modulo( iy+iky, ngridy )
    do ikz = 0, ngridz-1
      ikqz = modulo( ikz, ngridz )

      do io = 1, nblind
        ek = ei(io,ikx,iky,ikz)
        do jo = 1, nblind
          ekq = ei(jo,ikqx,ikqy,ikqz)
          dif = ekq-ek

          if ( abs(dif) > 1.0e-30_dp ) then
            ocup   = ( fe(io,ikx,iky,ikz) - fe(jo,ikqx,ikqy,ikqz) )
            mysusc = mysusc + ocup / dif
          endif
        enddo
      enddo

    enddo ! kz
    enddo ! ky
    enddo ! kx

#ifdef MPI
    call MPI_REDUCE( mysusc, susc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                     MPI_COMM_WORLD, ierr )
#else
    susc = mysusc
#endif
    if( myid == 0 ) then
      write(6 ,'(3f20.10)') q(1), q(2), 0.5_dp * susc / (ngridx*ngridy*ngridz)
      write(iu,'(3f20.10)') q(1), q(2), 0.5_dp * susc / (ngridx*ngridy*ngridz)
    endif
  enddo
  enddo

  deallocate( ei )
  deallocate( fe )

#ifdef MPI
  call MPI_FINALIZE( ierr )
#endif

end program lindhard_function