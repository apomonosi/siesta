module QMMM_mmforces_m
   !! This module contains the subroutines needed to calculate
   !! electron and nuclei forces over point charges (i.e. atoms).
   !
   ! There are two steps covered in this module:
   !   1) Forces over MM atoms due to the SCF density-point charge interaction.
   !   2) Forces over MM atoms AND QM atoms due to the neutral atom potential-
   !      point charge interaction.
   !
   ! If we split the QM charge density into nuclei + atomic(initial) + SCF, we
   ! see that the last item covers the first two contributions.
   !
   ! NOTE: Derivations here coincide with the gradient, but the minus sign
   !       is missing (i.e., forces should be "- gradients"). This is probably
   !       due to the fact that the electron density is positive, while in
   !       reality the electron charge is negative. Thus, a -1 appears in front
   !       of everything and the forces are actual equal to the gradients.
   !
   !
   ! FORCES DUE TO SCF DENSITY (DRHO)
   !
   ! The forces due to SCF are calculated with Ewald sums in the same fashion
   ! as the electrostatic potential; in fact, the equations stem from the fact
   ! that the total energy is the grid integral of rho * Vpc.
   ! As usual, we have two contributions to the forces, the reciprocal and the
   ! real (r_ij is the distance vector, while d_ij = |r_ij|):
   !
   ! Freal(i) = q_i * sum_j dVol * rho(j) / d_ij^2 * ( A - B ) * r_ij
   !      ->> A = erfc( sqrt(a) * d_ij) / d_ij
   !      ->> B = 2 * sqrt(a/pi) * exp( -a * d_ij^2)
   !
   ! Frecip(i) = ( q_i / V ) * sum_k (4pi / k^2) exp(-k^2 / 4a) * S_kr
   !   ->> S_kr = cos(k * r_i) * S_imag - sin(k * r_i) * S_real
   !         ->> S_real = sum_j dVol * rho(j) * cos(k * r_ij)
   !         ->> S_imag = sum_j dVol * rho(j) * sin(k * r_ij)
   !
   ! Where a is the alpha gaussian term in Ewald, and j loops over grid
   ! points. It's important to note that since the density is normalized
   ! over volume, we need to take into account the volume element dVol
   ! in the summation (with dVol * rho being the actual "charge").
   ! In the actual implementation, dVol is post-multiplied after all
   ! summations are done, since actually adding it within the summations
   ! can result in severe precision and numerical errors.
   !
   ! The stress tensor follows a similar derivation, starting from
   ! - 1/V * dEqmmm / de_uv, with e being the strain.
   !
   ! Sreal_uv = sum_i q_i * sum_j dVol * rho(j) * ( A - B ) * C
   !      ->> A = erfc( sqrt(a) * d_ij) / d_ij
   !      ->> B = 2 * sqrt(a/pi) * exp( -a * d_ij^2)
   !      ->> C = r_ij(u) * r_ij(v) / d_ij^2
   !
   ! Srecip_uv = ( 1 / V^2 ) * sum_k (4pi / k^2) exp(-k^2 / 4a) * A * S_kr
   !   ->>    A = delta(u,v) - 2 * k(u) * k(v) (1 + k^2 / 4a) / k^2
   !   ->> S_kr = Smm_real * Sqm_imag + Smm_imag * Sqm_real
   !         ->> Sqm_real = sum_j dVol * rho(j) * cos(k * r_ij)
   !         ->> Sqm_imag = sum_j dVol * rho(j) * sin(k * r_ij)
   !         ->> Smm_real = sum_i q_i * cos(k * r_ij)
   !         ->> Smm_imag = sum_i q_i * sin(k * r_ij)
   !
   !
   ! FORCES DUE TO NEUTRAL ATOM POTENTIAL (Vna)
   !
   ! These are actually much simpler since SIESTA already provides Vna
   ! and its gradients as a function of distance with respect to the nuclei,
   ! via the atmfuncs module.
   !
   ! We essentially have:
   !   1) Fmm_i = q_i * sum_j( gradVa_j(r_i) )
   !   2) Fqm_j = sum_i ( q_i * gradVa_j(r_i) )
   !
   ! And the stress follows a similar derivation.
   !
   implicit none
   public :: QMMM_mmforces
   public :: QMMM_mmforce_vna

   private
contains
  subroutine QMMM_mmforces( ucell, Rho, ntpl, ntml, ntm, dvol, lattice_volume )
    !! Computes forces exerted by electrons and nuclei over point charges.
    !! This subroutine is an interface for existing methods.
    !
    ! Dev Note: An MPI barrier is called at the end of the routine to
    ! properly compute timings in cases where grid distribution is unbalanced.
    use precision     , only : dp, grid_p
    use mesh          , only : meshLim, nsm
    use QMMM_core     , only : doing_QMMM, QMMM_coulomb_type, &
                               COULOMB_EWALD, COULOMB_REAL_CUTOFF
    use QMMM_structure, only : mm_charges
#ifdef MPI
    use mpi_siesta    , only : MPI_Barrier, MPI_Comm_World
#endif
    implicit none
    integer     , intent(in)    :: ntpl
      !! Total number of grid points for this node.
    integer     , intent(in)    :: ntml(3)
      !! Number of mesh points along each axis for this node.
    integer     , intent(in)    :: ntm(3)
      !! Number of mesh points along each axis.
    real(dp)    , intent(in)    :: ucell(3,3)
      !! Unit cell vectors.
    real(dp)    , intent(in)    :: dvol
      !! Volume of a point of the unit cell.
    real(dp)    , intent(in)    :: lattice_volume
      !! Total lattice volume.
    real(grid_p), intent(in)    :: Rho(ntpl)
      !! Electron density in the grid.

    integer :: box_offset(3)
#ifdef MPI
    integer :: MPIerr
#endif

    if ( mm_charges%n < 1 ) return

    box_offset(:) = ( meshLim(1,:) -1 ) * nsm

    mm_charges%f(:,:)      = 0.0_dp
    mm_charges%stress(:,:) = 0.0_dp
    select case ( QMMM_coulomb_type() )
    case ( COULOMB_EWALD )
      call mmforce_ewald( mm_charges%n, mm_charges%r, mm_charges%f,      &
                          mm_charges%pc, mm_charges%stress, ucell, Rho,  &
                          ntpl, ntml, ntm, box_offset, lattice_volume, dvol )
    case ( COULOMB_REAL_CUTOFF )
      call mmforce_direct( mm_charges%n, mm_charges%r, mm_charges%f,      &
                           mm_charges%pc, mm_charges%stress, ucell, Rho,  &
                           ntpl, ntml, ntm, box_offset, lattice_volume, dvol )
    case default
      call die( "ERROR: QMMM_mmforces - Wrong type of Coulomb interaction.")
    end select

    ! Fix timer issues.
#ifdef MPI
    call MPI_Barrier( MPI_COMM_WORLD, MPIerr)
#endif
  end subroutine QMMM_mmforces

  subroutine mmforce_direct( na_mm, rmm, fmm, pc, stress, ucell, Rho, &
                             ntpl, ntml, ntm, box0, lattice_volume, dvol )
    !! Computes forces exerted by electrons and nuclei over point charges using
    !! a cut-off scheme.
    use precision     , only : dp, grid_p
    use QMMM_core     , only : QMMM_density_rcut
    use QMMM_helper   , only : get_lattice_type, pbc_displ_vector
    use QMMM_neighbour, only : num_mmvecs, grid_veclist, grid_nr
#ifdef MPI
    use alloc         , only : re_alloc, de_alloc
    use mpi_siesta    , only : MPI_AllReduce, MPI_double_precision, &
                               MPI_sum, MPI_Comm_World
#endif

    implicit none
    integer     , intent(in) :: na_mm
      !! Number of MM atoms.
    integer     , intent(in) :: ntpl
      !! Total number of grid points.
    integer     , intent(in) :: ntml(3)
      !! Number of mesh points along each axis for this node.
    integer     , intent(in) :: ntm(3)
      !! Number of mesh points along each axis.
    integer     , intent(in) :: box0(3)
      !! The starting mesh cell index for this node.
    real(dp)    , intent(in) :: rmm(3,na_mm)
      !! Atomic positions.
    real(dp)    , intent(in) :: pc(na_mm)
      !! Atomic (classical) partial charges.
    real(dp)    , intent(in) :: ucell(3,3)
      !! Unit cell vectors.
    real(dp)    , intent(in) :: dvol
      !! Volume of a point of the unit cell.
    real(dp)    , intent(in) :: lattice_volume
      !! Total lattice volume.
    real(grid_p), intent(in) :: Rho(ntpl)
      !! Electron density in the grid.
    real(dp)    , intent(inout) :: fmm(3,na_mm)
      !! Atomic forces.
    real(dp)    , intent(inout) :: stress(3,3)
      !! Cell stress.

    character :: lattice_type
    integer   :: ix, iy, iz, imesh, iat, icrd, ivec, js, ix0, iy0, iz0
    real(dp)  :: drij(3), rcut_rho2, stress_fact, dist, xm, ym, zm, dE, &
                 kcell(3,3)
#ifdef MPI
    integer :: MPIerr
    real(dp), pointer :: flocal(:,:), stresslocal(:,:)
#endif

    if ( num_mmvecs == 0 ) return
    rcut_rho2   = ( QMMM_density_rcut() ) ** 2
    stress_fact = 1.0_dp / lattice_volume
    call reclat( ucell, kcell, 0 )

    lattice_type = get_lattice_type( ucell )
    ! For now, we put a range of orbital equals to rmax0 (bohrs).
    ! Loop over the mesh points
    do iz0 = 0, ntml(3) -1
    do iy0 = 0, ntml(2) -1
    do ix0 = 0, ntml(1) -1
      imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

      ! We use only grid points where there is density of charge.
      if ( abs( Rho(imesh) ) > 0.0_dp ) then
        ix = ix0 + box0(1)
        iy = iy0 + box0(2)
        iz = iz0 + box0(3)

        xm = ucell(1,1) * ix / ntm(1) + ucell(1,2) * iy / ntm(2) &
           + ucell(1,3) * iz / ntm(3)
        ym = ucell(2,1) * ix / ntm(1) + ucell(2,2) * iy / ntm(2) &
           + ucell(2,3) * iz / ntm(3)
        zm = ucell(3,1) * ix / ntm(1) + ucell(3,2) * iy / ntm(2) &
           + ucell(3,3) * iz / ntm(3)

        do iat = 1, num_mmvecs ! Loop over MM atoms
          js = grid_veclist(iat)

          if ( lattice_type == 'D' ) then
            drij(1) = xm - rmm(1,js) + grid_nr(1,iat) * ucell(1,1)
            drij(2) = ym - rmm(2,js) + grid_nr(2,iat) * ucell(2,2)
            drij(3) = zm - rmm(3,js) + grid_nr(3,iat) * ucell(3,3)
          else
            drij(1) = xm - rmm(1,js)
            drij(2) = ym - rmm(2,js)
            drij(3) = zm - rmm(3,js)

            do icrd = 1, 3
            do ivec = 1, 3
              drij(icrd) = drij(icrd) + grid_nr(ivec,iat) * ucell(icrd,ivec)
            enddo
            enddo
          endif
          call pbc_displ_vector( ucell, kcell, drij )
          dist = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)

          ! Forces exerted on point charges by electrons.
          De = 0.0_dp
          if ( dist < rcut_rho2 ) cycle

          ! If our density point and point charge are too close,
          ! this might diverge. Thus, we avoid it.
          dist = ( 1.0_dp / sqrt( dist ) ) ** 3
          De   = -2.0_dp * dvol * Rho(imesh) * pc(js) * dist

          do ivec = 1, 3
            fmm(ivec,js) = fmm(ivec,js) - De * drij(ivec)

            do icrd = 1, 3
              stress(icrd,ivec) = stress(icrd,ivec) &
                               - stress_fact * rmm(icrd,js) * De * drij(ivec)
            enddo
          enddo
        enddo ! MM atoms
      endif   ! abs(Rho(imesh)) > 0
    enddo     ! grid X
    enddo     ! grid Y
    enddo     ! grid Z

#ifdef MPI
    nullify( flocal, stresslocal )

    call re_alloc( flocal, 1, 3, 1, na_mm, 'flocal', 'mmforce_direct' )
    call re_alloc( stresslocal, 1, 3, 1, 3,  'stresslocal', 'mmforce_direct' )

    flocal(:,:)      = fmm(:,:)
    stresslocal(:,:) = stress(:,:)

    fmm(:,:)    = 0.0_dp
    stress(:,:) = 0.0_dp
    call MPI_AllReduce( flocal(1,1), fmm(1,1), 3*na_mm, MPI_double_precision, &
                        MPI_sum, MPI_Comm_World, MPIerr )
    call MPI_AllReduce( stresslocal(1,1), stress(1,1), 9, MPI_double_precision,&
                        MPI_sum, MPI_Comm_World, MPIerr )

    call de_alloc( flocal     , 'flocal'     , 'mmforce_direct' )
    call de_alloc( stresslocal, 'stresslocal', 'mmforce_direct' )
#endif

  end subroutine mmforce_direct

  subroutine mmforce_ewald( na_mm, rmm, fmm, pc, stress, ucell, Rho, ntpl, &
                            ntml, ntm, box0, lattice_volume, dvol )
    !! Computes forces exerted by electrons and nuclei over point charges using
    !! the Ewald method.
    use alloc         , only : re_alloc, de_alloc
    use precision     , only : dp, grid_p
    use QMMM_core     , only : QMMM_density_rcut
    use QMMM_ewald_m  , only : QMMM_ewald
    use QMMM_helper   , only : pbc_displ_vector
    use units         , only : pi2, pi
    use sys           , only : die
#ifdef MPI
    use alloc         , only : re_alloc, de_alloc
    use mpi_siesta    , only : MPI_AllReduce, MPI_double_precision, &
                               MPI_sum, MPI_Comm_World
#endif

    implicit none
    integer     , intent(in) :: na_mm
      !! Number of MM atoms.
    integer     , intent(in) :: ntpl
      !! Total number of grid points.
    integer     , intent(in) :: ntml(3)
      !! Number of mesh points along each axis for this node.
    integer     , intent(in) :: ntm(3)
      !! Number of mesh points along each axis.
    integer     , intent(in) :: box0(3)
      !! The starting mesh cell index for this node.
    real(dp)    , intent(in) :: rmm(3,na_mm)
      !! Atomic positions.
    real(dp)    , intent(in) :: pc(na_mm)
      !! Atomic (classical) partial charges.
    real(dp)    , intent(in) :: ucell(3,3)
      !! Unit cell vectors.
    real(dp)    , intent(in) :: dvol
      !! Volume of a point of the unit cell.
    real(dp)    , intent(in) :: lattice_volume
      !! Total lattice volume.
    real(grid_p), intent(in) :: Rho(ntpl)
      !! Electron density in the grid.
    real(dp)    , intent(inout) :: fmm(3,na_mm)
      !! Atomic forces.
    real(dp)    , intent(inout) :: stress(3,3)
      !! Cell stress.

    integer   :: iat, icrd, ivec, imesh, ix, iy, iz, ikew, ix0, iy0, iz0, &
                 nmesh
    real(dp)  :: rcut_rho3, kcell(3,3), dist, dist2, krecip(3), kronij, &
                 kr, drij(3), dE, xm, ym, zm, sqAlPi, pi4vol, &
                 expterm, volInv, inv4a, kpref, rcut2, f_at(3), ew_pre, &
                 Sqm_r, Sqm_i, str_at(3,3)

    real(dp), pointer :: S_qm_real(:), S_qm_imag(:), S_mm_real(:), &
                         S_mm_imag(:), stress_atom(:,:,:), meshcrd(:,:)
    integer , pointer :: meshidx(:)

#ifdef MPI
    integer :: MPIerr, buff_size
    real(dp), pointer :: frc_buf(:,:), str_buf(:,:,:), S_buffer(:)
#endif
    volInv = 1.0_dp / lattice_volume
    pi4vol = 2.0_dp * pi2 * volInv

    ! We first pre-filter the mesh points to simplify further loops.
    !
    ! The total cost of the current operation is:
    !   nmesh * 2 + nmesh_filtered * na_mm + nmesh_filtered * n_kewald
    !
    ! Compare to the standard apporach, where it is:
    !   nmesh * na_mm + nmesh * n_kewald
    !
    ! So essentially, so long as n_kewald + na_mm is larger than 2,
    ! it is safe to assume that this approach is much more efficient.
    nmesh = 0
    do iz0 = 0, ntml(3) -1
    do iy0 = 0, ntml(2) -1
    do ix0 = 0, ntml(1) -1
      imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

      if ( .not. (abs(Rho(imesh)) > 0.0_dp) ) cycle
      nmesh = nmesh + 1
    enddo
    enddo
    enddo

    nullify( meshcrd, meshidx )
    call re_alloc( meshcrd, 1, 3, 1, nmesh, 'meshcrd', 'mmforce_ewald' )
    call re_alloc( meshidx, 1, nmesh, 'meshidx', 'mmforce_ewald' )

    ! We pre-compute the coordinates for each mesh point.
    nmesh = 0
    do iz0 = 0, ntml(3) -1
    do iy0 = 0, ntml(2) -1
    do ix0 = 0, ntml(1) -1
      imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

      if ( .not. (abs(Rho(imesh)) > 0.0_dp) ) cycle
      nmesh = nmesh + 1

      ix = ix0 + box0(1)
      iy = iy0 + box0(2)
      iz = iz0 + box0(3)

      xm = ucell(1,1) * ix / ntm(1) + ucell(1,2) * iy / ntm(2) &
         + ucell(1,3) * iz / ntm(3)
      ym = ucell(2,1) * ix / ntm(1) + ucell(2,2) * iy / ntm(2) &
         + ucell(2,3) * iz / ntm(3)
      zm = ucell(3,1) * ix / ntm(1) + ucell(3,2) * iy / ntm(2) &
         + ucell(3,3) * iz / ntm(3)

      meshidx(nmesh)   = imesh
      meshcrd(1,nmesh) = xm
      meshcrd(2,nmesh) = ym
      meshcrd(3,nmesh) = zm
    enddo
    enddo
    enddo

    ! We calculate the reciprocal lattice vectors
    call reclat( ucell, kcell, 0 )

    nullify( S_qm_real, S_qm_imag )
    call re_alloc( S_qm_real, 1, QMMM_ewald%n_kpoints, 'S_qm_real', &
                   'mmforce_ewald' )
    call re_alloc( S_qm_imag, 1, QMMM_ewald%n_kpoints, 'S_qm_imag', &
                   'mmforce_ewald' )

    ! We use this to store each MM atom's contribution to stress, to then
    ! add it up properly in a reduction.
    nullify( stress_atom )
    call re_alloc( stress_atom, 1, 3, 1, 3, 1, na_mm, 'stress_atom', &
                   'mmforce_ewald' )
    stress_atom(:,:,:) = 0.0_dp

    ! Real part of the Ewald summation
    ! Loop over the mesh points
    rcut_rho3 = QMMM_density_rcut()
    rcut_rho3 = rcut_rho3 * rcut_rho3 * rcut_rho3
    rcut2     = QMMM_ewald%rcut * QMMM_ewald%rcut
    sqAlPi    = QMMM_ewald%sqalpha / sqrt( pi )

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(Rho,fmm,stress_atom,meshcrd,meshidx) &
!$OMP& PRIVATE(imesh,iat,de,dist,dist2,ew_pre)
!$OMP DO SCHEDULE(DYNAMIC,2)
    do iat = 1, na_mm
      ! IMPORTANT: The current implementation of this loop is
      ! done to reduce the accumulation of precision errors.
      ! avoid modifying or simplifying it unless you intend
      ! to test this rigorously.
      f_at(:)     = 0.0_dp
      str_at(:,:) = 0.0_dp

      do imesh = 1, nmesh
        drij(1) = meshcrd(1,imesh) - rmm(1,iat)
        drij(2) = meshcrd(2,imesh) - rmm(2,iat)
        drij(3) = meshcrd(3,imesh) - rmm(3,iat)
        call pbc_displ_vector( ucell, kcell, drij )
        dist2 = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)

        ! Forces exerted on point charges by electrons.
        De = 0.0_dp
        if ( dist2 > rcut2 ) cycle

        ! If our density point and point charge are too close,
        ! this might diverge. Thus, we avoid it. Be aware that
        ! with a low value for rhocut, this might be numerically
        ! unstable.
        dist  = sqrt( dist2 )
        dist  = dist * dist * dist + rcut_rho3
        dist  = dist ** (1.0_dp / 3.0_dp)
        dist2 = dist * dist

        ew_pre = 2.0_dp * sqAlPi *exp( -QMMM_ewald%alpha * dist2 ) + &
                 erfc( QMMM_ewald%sqalpha * dist ) / dist
        ew_pre = ew_pre / dist2

        De      = Rho( meshidx(imesh) ) * ew_pre
        f_at(:) = f_at(:) + De * drij(:)

        do icrd = 1, 3
          str_at(icrd,:) = str_at(icrd,:) - drij(icrd) * De * drij(:)
        enddo
      enddo     ! grid

      fmm(:,iat)           = fmm(:,iat) + pc(iat) * f_at(:)
      stress_atom(:,:,iat) = stress_atom(:,:,iat) + str_at(:,:)
    enddo   ! mm atom
!$OMP END DO
!$OMP END PARALLEL

#ifdef MPI
    nullify( frc_buf, str_buf )

    call re_alloc( frc_buf, 1, 3, 1, na_mm, 'frc_buf', 'mmforce_ewald' )
    call re_alloc( str_buf, 1, 3, 1, 3, 1, na_mm, 'str_buf', 'mmforce_ewald' )

    frc_buf(:,:)   = fmm(:,:)
    str_buf(:,:,:) = stress_atom(:,:,:)

    fmm(:,:)           = 0.0_dp
    stress_atom(:,:,:) = 0.0_dp
    call MPI_AllReduce( frc_buf(1,1), fmm(1,1), 3*na_mm, MPI_double_precision,&
                        MPI_sum, MPI_Comm_World, MPIerr )
    call MPI_AllReduce( str_buf(1,1,1), stress_atom(1,1,1), 9*na_mm, &
                        MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerr )

    call de_alloc( frc_buf, 'frc_buf', 'mmforce_ewald' )
    call de_alloc( str_buf, 'str_buf', 'mmforce_ewald' )
#endif
    fmm(:,:) = fmm(:,:) * dvol * 2.0_dp

    !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(stress_atom,stress) &
    !$OMP& PRIVATE(de)
    !$OMP DO SCHEDULE(DYNAMIC,2)
    do ivec = 1, 3
    do icrd = 1, 3
      De = 0.0_dp
      do iat = 1, na_mm
        De = De + stress_atom(icrd,ivec,iat) * pc(iat)
      enddo
      stress(icrd,ivec) = stress(icrd,ivec) + De * dVol * volInv
    enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call de_alloc( stress_atom, 'stress_atom', 'mmforce_ewald' )

    ! Reciprocal space part of the ewald summation.
    ! Calculate structure factors for QM atoms.
    ! We first initialize and then loop over mesh points.
    S_qm_real(:) = 0.0_dp
    S_qm_imag(:) = 0.0_dp

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(S_qm_real,S_qm_imag,Rho) &
!$OMP& PRIVATE(imesh,kr,Sqm_r,Sqm_i,ikew)
!$OMP DO SCHEDULE(DYNAMIC,2)
    do ikew = 1, QMMM_ewald%n_kpoints
      ! IMPORTANT: The current implementation of this loop is
      ! done to reduce the accumulation of precision errors.
      ! avoid modifying or simplifying it unless you intend
      ! to test this rigorously.
      ! Sqm_r and Sqm_i might seem superfluous but they are
      ! actually key to avoid numerical issues.
      Sqm_r = 0.0_dp
      Sqm_i = 0.0_dp

      do imesh = 1, nmesh
        kr = QMMM_ewald%kx(ikew) * meshcrd(1,imesh) &
           + QMMM_ewald%ky(ikew) * meshcrd(2,imesh) &
           + QMMM_ewald%kz(ikew) * meshcrd(3,imesh)
        Sqm_r = Sqm_r + Rho( meshidx(imesh) ) * cos(kr)
        Sqm_i = Sqm_i + Rho( meshidx(imesh) ) * sin(kr)
      enddo

      S_qm_real(ikew) = S_qm_real(ikew) + Sqm_r
      S_qm_imag(ikew) = S_qm_imag(ikew) + Sqm_i
    enddo
!$OMP END DO
!$OMP END PARALLEL

    call de_alloc( meshcrd, 'meshcrd', 'mmforce_ewald' )
    call de_alloc( meshidx, 'meshidx', 'mmforce_ewald' )

#ifdef MPI
    nullify( S_buffer )
    buff_size = QMMM_ewald%n_kpoints

    call re_alloc( S_buffer, 1, buff_size, 'S_buffer', 'mmforce_ewald' )

    S_buffer  = S_qm_real
    S_qm_real = 0.0_dp
    call MPI_AllReduce( S_buffer(1), S_qm_real(1), buff_size, &
                        MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerr )

    S_buffer  = S_qm_imag
    S_qm_imag = 0.0_dp
    call MPI_AllReduce( S_buffer(1), S_qm_imag(1), buff_size, &
                        MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerr )

    call de_alloc( S_buffer, 'S_buffer', 'mmforce_ewald' )
#endif
    ! Calculate the reciprocal part of the force on classical atoms
    ! due to QM grid points
    inv4a  = 1.0_dp / ( 4.0_dp * QMMM_ewald%alpha )
    !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(fmm) &
    !$OMP& PRIVATE(iat, ikew, krecip, kr, de, expterm, kronij, kpref)
    !$OMP DO SCHEDULE(DYNAMIC,2)
    do iat = 1, na_mm
      ! IMPORTANT: The current implementation of this loop is
      ! done to reduce the accumulation of precision errors.
      ! avoid modifying or simplifying it unless you intend
      ! to test this rigorously.
      f_at(:)     = 0.0_dp

      do ikew = 1, QMMM_ewald%n_kpoints
        krecip(1) = QMMM_ewald%kx(ikew)
        krecip(2) = QMMM_ewald%ky(ikew)
        krecip(3) = QMMM_ewald%kz(ikew)

        expterm = exp( -QMMM_ewald%kmod2(ikew) * inv4a )
        expterm = expterm / QMMM_ewald%kmod2(ikew)

        kr = krecip(1) * rmm(1,iat) + krecip(2) * rmm(2,iat) &
           + krecip(3) * rmm(3,iat)
        De = expterm * ( cos(kr) * S_qm_imag(ikew) &
                       - sin(kr) * S_qm_real(ikew) )

        f_at(:) = f_at(:) + De * krecip(:)
      enddo

      fmm(:,iat) = fmm(:,iat) + f_at(:) * pc(iat) * pi4vol * dvol * 2.0_dp
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! Calculate structure factors for all MM atoms. This is
    ! used only for Stress calculation
    nullify( S_mm_real, S_mm_imag )
    call re_alloc( S_mm_real, 1, QMMM_ewald%n_kpoints, 'S_mm_real', &
                   'mmforce_ewald' )
    call re_alloc( S_mm_imag, 1, QMMM_ewald%n_kpoints, 'S_mm_imag', &
                   'mmforce_ewald' )
    S_mm_real(:) = 0.0_dp
    S_mm_imag(:) = 0.0_dp

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(S_mm_real,S_mm_imag) &
!$OMP& PRIVATE(kr)
!$OMP DO SCHEDULE(DYNAMIC,2)
    do ikew = 1, QMMM_ewald%n_kpoints
    do iat = 1, na_mm
      kr = QMMM_ewald%kx(ikew) * rmm(1,iat) &
         + QMMM_ewald%ky(ikew) * rmm(2,iat) &
         + QMMM_ewald%kz(ikew) * rmm(3,iat)
      S_mm_real(ikew) = S_mm_real(ikew) + pc(iat) * cos(kr)
      S_mm_imag(ikew) = S_mm_imag(ikew) + pc(iat) * sin(kr)
    enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    str_at = 0.0_dp
    ! Stress calculation.
    do ikew = 1, QMMM_ewald%n_kpoints
      ! IMPORTANT: The current implementation of this loop is
      ! done to reduce the accumulation of precision errors.
      ! avoid modifying or simplifying it unless you intend
      ! to test this rigorously.
      krecip(1) = QMMM_ewald%kx(ikew)
      krecip(2) = QMMM_ewald%ky(ikew)
      krecip(3) = QMMM_ewald%kz(ikew)

      expterm = exp( -QMMM_ewald%kmod2(ikew) * inv4a )
      expterm = expterm / QMMM_ewald%kmod2(ikew)

      De = expterm * ( S_mm_imag(ikew) * S_qm_imag(ikew) &
                     + S_mm_real(ikew) * S_qm_real(ikew) )

      kpref = 1.0_dp + QMMM_ewald%kmod2(ikew) * inv4a
      kpref = 2.0_dp * kpref / QMMM_ewald%kmod2(ikew)

      do iVec = 1, 3
      do iCrd = 1, 3
        kronij = real( int( ( iCrd + iVec - abs( iCrd - iVec ) ) &
                          / ( iCrd + iVec + abs( iCrd - iVec ) ) ), kind=dp )
        str_at(iCrd,iVec) = str_at(iCrd,iVec) + De * &
                          ( kronij - kpref * krecip(iCrd) * krecip(iVec) )
      enddo
      enddo
    enddo
    stress(:,:) = stress(:,:) + str_at(:,:) * pi4vol * dvol * volInv

    call de_alloc( S_qm_real, 'S_qm_real', 'mmforce_ewald' )
    call de_alloc( S_qm_imag, 'S_qm_imag', 'mmforce_ewald' )
    call de_alloc( S_mm_real, 'S_mm_real', 'mmforce_ewald' )
    call de_alloc( S_mm_imag, 'S_mm_imag', 'mmforce_ewald' )
  end subroutine mmforce_ewald

  subroutine QMMM_mmforce_vna( na_qm, r_qm, isa, ucell, fal, stressl, nodes )
    !! Adds the Vna contributions to forces over MM atoms.
    use atmfuncs      , only : rcut, phiatm
    use precision     , only : dp
    use QMMM_helper   , only : pbc_displ_vector
    use QMMM_structure, only : mm_charges

    implicit none
    integer, intent(in)   :: na_qm
      !! Number of QM atoms (i.e. SIESTA number of atoms).
    real(dp), intent(in)  :: r_qm(3,na_qm)
      !! Positions of the QM atoms.
    integer , intent(in)  :: isa(na_qm)
      !! Atomic species of QM atoms.
    real(dp), intent(in)  :: ucell(3,3)
      !! Unit cell vectors.
    real(dp), intent(inout) :: fal(3,na_qm)
      !! Forces over QM atoms in current node.
    real(dp), intent(inout) :: stressl(3,3)
      !! Stress contribution from current node.
    integer , intent(in)  :: nodes
      !! Number of parallel processes.

    real(dp) :: kcell(3,3), drij(3), qmcut2, dist2, va, grva(3), volinv, &
                f_at(3), str_at(3,3)
    integer  :: iat, jat, is, icrd

    external :: reclat
    real(dp), external :: volcel

    call reclat( ucell, kcell, 0 )
    if ( mm_charges%n < 1 ) return

    volinv = 1.0_dp / volcel( ucell )

    do iat = 1, mm_charges%n
      f_at(:)     = 0.0_dp
      str_at(:,:) = 0.0_dp

      do jat = 1, na_qm
        is     = isa(jat)
        qmcut2 = rcut(is,0)
        qmcut2 = qmcut2 * qmcut2

        drij(:) = mm_charges%r(:,iat) - r_qm(:,jat)
        call pbc_displ_vector( ucell, kcell, drij )

        dist2 = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)
        if ( dist2 > qmcut2 ) cycle

        call phiatm( is, 0, drij, va, grva )

        f_at(:)    = f_at(:) + grva(:)
        ! Here we divide over nodes, since fal should be reduced afterwards.
        ! The fal array is local to the node, but here we are calculating
        ! forces over ALL qm atoms.
        fal(:,jat) = fal(:,jat) - grva(:) * mm_charges%pc(iat) / nodes
        do icrd = 1, 3
          str_at(:,icrd) = str_at(:,icrd) + drij(:) * grva(icrd)
        enddo
      enddo

      mm_charges%f(:,iat) = mm_charges%f(:,iat) &
                          + f_at(:) * mm_charges%pc(iat)
      str_at(:,:)  = str_at(:,:) * mm_charges%pc(iat) * volInv

      ! Same reasoning as with fal applies to stressl.
      stressl(:,:) = stressl(:,:) + str_at(:,:) / nodes
      mm_charges%stress(:,:) = mm_charges%stress(:,:) + str_at(:,:)
    enddo
  end subroutine QMMM_mmforce_vna

end module QMMM_mmforces_m
