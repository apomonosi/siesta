module QMMM_pcpot_m
   !! This module contains subroutines to calculate the potential
   !! due to solvent molecules at each point of the mesh.
   !
   ! The main method to generate this potential is via the Ewald
   ! summation method. This essentially separates the potential
   ! into a real-space and reciprocal-space contributions, which
   ! follow these equations for a given point j in the grid:
   !
   ! Vreal(j) = sum_i q_i erfc( sqrt(a) * d_ij) / d_ij
   !
   ! Vrecip(j) = (1/V) * sum_k (4pi / k^2) exp(-k^2 / 4a) * S_kr
   !   ->> S_kr = cos(k * r_j) * S_real - sin(k * r_j) * S_imag
   !         ->> S_real = sum_i q_i * cos(k * r_ij)
   !         ->> S_imag = sum_i q_i * sin(k * r_ij)
   !
   ! Where a (alpha) is the inverse of the square of the gaussian
   ! widths for the fictitious charges.
   !
   ! These sums are rapidly convergent and include both a real
   ! space cut-off and a reciprocal space cut-off. How both
   ! of these are related to alpha is explained in the Ewald module.
   !
   ! After everything is calculated, the entire potential is
   ! multiplied by -2. Presumably, this is to avoid issues with
   ! Rydberg units and due to the electron density being a positive
   ! number despite being a negative charge ( which appear later
   ! in things like rho * Vpc ).
   !
   ! An alternate, direct-cutoff method is provided just for
   ! debugging purposes.
   implicit none
   public :: QMMM_pcpot
   public :: QMMM_set_ewald_grid

   private
contains

  subroutine QMMM_pcpot( qmcell, ntpl, ntml, ntm, rhoatm )
    !! Calculates the potential for a set of partial charges.
    !! This subroutine acts as an interface for the main methods of calculation.
    !
    ! Dev Note: An MPI barrier is called at the end of the routine to
    ! properly compute timings in cases where grid distribution is unbalanced.
    use alloc         , only : re_alloc, de_alloc
    use mesh          , only : meshLim, nsm
    use precision     , only : dp, grid_p
    use QMMM_core     , only : doing_QMMM, QMMM_coulomb_type, &
                               COULOMB_EWALD, COULOMB_REAL_CUTOFF
    use QMMM_structure, only : mm_charges
#ifdef MPI
    use mpi_siesta    , only : MPI_Barrier, MPI_Comm_World
#endif

    implicit none
    integer     , intent(in) :: ntpl
      !! Total number of grid points.
    integer     , intent(in) :: ntml(3)
      !! Number of mesh points along each axis for this node.
    integer     , intent(in) :: ntm(3)
      !! Number of mesh points along each axis.
    real(dp)    , intent(in) :: qmcell(3,3)
      !! Unit cell vectors.
    real(dp)    , intent(in) :: rhoatm(ntpl)
      !! Atomic density matrix, used to fill the potential grid.

    integer :: box_offset(3)
#ifdef MPI
    integer :: MPIerr
#endif

    if ( .not. doing_QMMM() ) return
    if ( mm_charges%n < 1 ) return

    if ( associated( mm_charges%Vpc ) ) &
      call de_alloc( mm_charges%Vpc, 'Vpc', 'QMMM_pcpot' )

    nullify( mm_charges%Vpc )
    call re_alloc( mm_charges%Vpc, 1, ntpl, 'Vpc', 'QMMM_pcpot' )

    box_offset(:) = ( meshLim(1,:) -1 ) * nsm

    select case ( QMMM_coulomb_type() )
    case ( COULOMB_EWALD )
      call pcpot_ewald( mm_charges%n, mm_charges%r, mm_charges%pc, qmcell, &
                        ntpl, ntml, ntm, box_offset, rhoatm, mm_charges%Vpc)
    case ( COULOMB_REAL_CUTOFF )
      call pcpot_direct( mm_charges%n, mm_charges%r, mm_charges%pc, qmcell, &
                         ntpl, ntml, ntm, box_offset, rhoatm, mm_charges%Vpc)
    case default
      call de_alloc( mm_charges%Vpc, 'Vpc', 'QMMM_pcpot' )
      call die( "ERROR: QMMM_mmforces - Wrong type of Coulomb interaction." )
    end select

#ifdef MPI
    call MPI_Barrier( MPI_COMM_WORLD, MPIerr)
#endif
  end subroutine QMMM_pcpot

  subroutine QMMM_set_ewald_grid( qmcell )
    use precision   , only : dp
    use QMMM_ewald_m, only : QMMM_ewald
    use QMMM_helper , only : get_lattice_type

    implicit none
    real(dp), intent(in)  :: qmcell(3,3)
      !! Unit cell vectors.
    real(dp)  :: kcell(3,3)
    character :: lattice_type

    ! We calculate the reciprocal lattice vectors
    lattice_type = get_lattice_type( qmcell )
    call reclat( qmcell, kcell, 0 )

    call QMMM_ewald%set_kgrid( kcell, lattice_type )
  end subroutine QMMM_set_ewald_grid

  subroutine pcpot_ewald( na_mm, rmm, pc, qmcell, ntpl, ntml, ntm, box0, &
                          Rho, Vpc )
    !! Calculates the potential using the Ewald summation method.
    use alloc         , only : re_alloc, de_alloc
    use precision     , only : dp, grid_p
    use QMMM_core     , only : QMMM_density_rcut
    use QMMM_ewald_m  , only : QMMM_ewald, smooth_potential
    use QMMM_helper   , only : get_lattice_type, pbc_displ_vector
    use units         , only : pi2
    use sys           , only : die

    implicit none
    integer     , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer     , intent(in)  :: ntpl
      !! Total number of grid points.
    integer     , intent(in)  :: ntml(3)
      !! Number of mesh points along each axis for this node.
    integer     , intent(in)  :: ntm(3)
      !! Number of mesh points along each axis.
    integer     , intent(in)  :: box0(3)
      !! The starting mesh cell index for this node.
    real(dp)    , intent(in)  :: rmm(3,na_mm)
      !! Atomic positions.
    real(dp)    , intent(in)  :: pc(na_mm)
      !! Atomic (classical) partial charges.
    real(dp)    , intent(in)  :: qmcell(3,3)
      !! Unit cell vectors.
    real(grid_p), intent(in)  :: Rho(ntpl)
      !! Electron density in the grid.
    real(grid_p), intent(out) :: Vpc(ntpl)
      !! Electrostatic potential for the QM region.

    character :: lattice_type
    integer   :: ix, iy, iz, iat, imesh, ix0, iy0, iz0, ikew, &
                 ntot(3), Len_to_Mesh_ele, Len_to_Mesh_smth

    real(dp)  :: ym, xm, zm, dist, rcut_rho3, drij(3), rcut2,       &
                 pi4vol, const_sin, kr, kcell(3,3), Mesh_to_Length, &
                 smooth_plus_electrode, pref, vpci, alpha4inv, scos, ssin
    real(dp), pointer   :: S_real(:), S_imag(:)
    real(dp), external  :: volcel

    rcut_rho3 = QMMM_density_rcut() ** 3
    rcut2     = QMMM_ewald%rcut * QMMM_ewald%rcut
    pi4vol    = 2.0_dp * pi2 / volcel( qmcell )
    alpha4inv = 1.0_dp / ( 4.0_dp * QMMM_ewald%alpha )

    ! We calculate the reciprocal lattice vectors
    lattice_type = get_lattice_type( qmcell )
    call reclat( qmcell, kcell, 0 )

    ntot(:) = QMMM_ewald%n_kbox
    Vpc(:) = 0.0_grid_p

    nullify( S_real, S_imag )
    call re_alloc( S_real, 1, QMMM_ewald%n_kpoints, 'S_real', &
                   'pcpot_ewald' )
    call re_alloc( S_imag, 1, QMMM_ewald%n_kpoints, 'S_imag', &
                   'pcpot_ewald' )

    ! Real part of the Ewald summation
    ! Loop over the mesh points
!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(Vpc,S_real,S_imag,Rho) &
!$OMP& PRIVATE(dist,imesh,ix,iy,iz,xm,ym,zm,vpci)
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC,2)
    do iz0 = 0, ntml(3) -1
    do iy0 = 0, ntml(2) -1
    do ix0 = 0, ntml(1) -1
      imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

      ! We use only grid points where there is density of charge.
      if ( .not. (abs(Rho(imesh)) > 0.0_dp) ) cycle
      ix = ix0 + box0(1)
      iy = iy0 + box0(2)
      iz = iz0 + box0(3)

      xm = qmcell(1,1) * ix / ntm(1) + qmcell(1,2) * iy / ntm(2) &
         + qmcell(1,3) * iz / ntm(3)
      ym = qmcell(2,1) * ix / ntm(1) + qmcell(2,2) * iy / ntm(2) &
         + qmcell(2,3) * iz / ntm(3)
      zm = qmcell(3,1) * ix / ntm(1) + qmcell(3,2) * iy / ntm(2) &
         + qmcell(3,3) * iz / ntm(3)

      ! Loop over MM atoms
      Vpci = 0.0_dp
      do iat = 1, na_mm
        drij(1) = xm - rmm(1,iat)
        drij(2) = ym - rmm(2,iat)
        drij(3) = zm - rmm(3,iat)

        call pbc_displ_vector( qmcell, kcell, drij )
        dist = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)

        if ( dist > rcut2 ) cycle
        ! Calculation of the external potential due to point charges plus
        ! Gaussian distributions.

        ! Real-space sum
        dist = sqrt( dist )
        dist = dist * dist * dist + rcut_rho3
        dist = dist ** (1.0_dp / 3.0_dp)

        Vpci = Vpci + pc(iat) * &
                      erfc( QMMM_ewald%sqalpha * dist ) / dist
      enddo ! at sv

      Vpc(imesh) = Vpc(imesh) + Vpci
    enddo     ! grid X
    enddo     ! grid Y
    enddo     ! grid Z
!$OMP END DO
!$OMP END PARALLEL

    ! Reciprocal part of the Ewald summation.
    ! Calculate structure factors for atoms.
    S_real(:) = 0.0_dp
    S_imag(:) = 0.0_dp

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(S_real,S_imag) &
!$OMP& PRIVATE(kr,iat,ikew,scos,ssin)
!$OMP DO SCHEDULE(DYNAMIC,2)
    do ikew = 1, QMMM_ewald%n_kpoints
      ! Sqm_r and Sqm_i might seem superfluous but they are
      ! actually key to avoid numerical issues.
      scos = 0.0_dp
      ssin = 0.0_dp

      do iat = 1, na_mm
        kr = QMMM_ewald%kx(ikew) * rmm(1,iat) &
           + QMMM_ewald%ky(ikew) * rmm(2,iat) &
           + QMMM_ewald%kz(ikew) * rmm(3,iat)

        scos = scos + pc(iat) * cos(kr)
        ssin = ssin + pc(iat) * sin(kr)
      enddo

      S_real(ikew) = S_real(ikew) + scos
      S_imag(ikew) = S_imag(ikew) + ssin
    enddo
!$OMP END DO
!$OMP END PARALLEL

    ! Calculate the reciprocal part of the potential on QM grid points.
    ! Loop over the mesh points.
!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(Vpc,S_real,S_imag,Rho) &
!$OMP& PRIVATE(imesh,kr,ix,iy,iz,xm,ym,zm,vpci,ssin,scos,pref)
!$OMP DO SCHEDULE(DYNAMIC,2)
    do iz0 = 0, ntml(3) -1
    do iy0 = 0, ntml(2) -1
    do ix0 = 0, ntml(1) -1
      imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

      ! We use only grid points where there is density of charge.
      if ( .not. (abs(Rho(imesh)) > 0.0_dp) ) cycle
      ix = ix0 + box0(1)
      iy = iy0 + box0(2)
      iz = iz0 + box0(3)

      xm = qmcell(1,1) * ix / ntm(1) + qmcell(1,2) * iy / ntm(2) &
         + qmcell(1,3) * iz / ntm(3)
      ym = qmcell(2,1) * ix / ntm(1) + qmcell(2,2) * iy / ntm(2) &
         + qmcell(2,3) * iz / ntm(3)
      zm = qmcell(3,1) * ix / ntm(1) + qmcell(3,2) * iy / ntm(2) &
         + qmcell(3,3) * iz / ntm(3)

      Vpci = 0.0_dp
      do ikew = 1, QMMM_ewald%n_kpoints
        kr = QMMM_ewald%kx(ikew) * xm + QMMM_ewald%ky(ikew) * ym + &
             QMMM_ewald%kz(ikew) * zm

        pref = exp( - QMMM_ewald%kmod2(ikew) * alpha4inv )
        pref = pref / QMMM_ewald%kmod2(ikew)

        scos = cos(kr) * S_real(ikew)
        ssin = sin(kr) * S_imag(ikew)

        Vpci = Vpci + pref * ( scos + ssin )
      enddo
      Vpc(imesh) = Vpc(imesh) + Vpci * pi4vol
    enddo     ! grid X
    enddo     ! grid Y
    enddo     ! grid Z
!$OMP END DO
!$OMP END PARALLEL
    call de_alloc( S_real, 'S_real', 'pcpot_ewald' )
    call de_alloc( S_imag, 'S_imag', 'pcpot_ewald' )

    !  Potential smoothing around the electrodes, if present.
    ! The smoothing makes some assumptions about the system:
    !
    !  * iz is the number of mesh points in z direction, and in this case, 360.
    !  * The electrode is in the region between points 1-60 and 300-360.
    !  * The smoothing is done by multiplying the final mesh by a quarter of
    !    a trigonometric function (going from 0 to 1).
    !  * The potential is set strictly to zero in the electrode region.
    !  * The size of the smoothing region is half of the electrode.
    if ( smooth_potential%on ) then
      write(6,*) "siesta-QMMM: WARNING: To use the smoothing function "//&
                 "the electrodes must be in the extremes of the box."
      write(6,*)

      Len_to_Mesh_ele  = &
        INT( smooth_potential%lelectrode * ntm(3) / qmcell(3,3) )
      Len_to_Mesh_smth = &
        INT( smooth_potential%lsmooth    * ntm(3) / qmcell(3,3) )
      smooth_plus_electrode = smooth_potential%lelectrode + &
                              smooth_potential%lsmooth

      const_sin = pi2 / ( 4.0_dp * Len_to_Mesh_smth )
!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(Vpc) &
!$OMP& PRIVATE(Mesh_to_Length,imesh,iz)
!$OMP DO SCHEDULE(DYNAMIC,2)
      do iz0 = 0, ntml(3) -1
      do iy0 = 0, ntml(2) -1
      do ix0 = 0, ntml(1) -1
        imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

        iz = iz0 + box0(3)

        Mesh_to_Length = iz * qmcell(3,3) / ntm(3)

        if ( ( Mesh_to_Length < smooth_plus_electrode ) .and. &
             ( Mesh_to_Length > smooth_potential%lelectrode ) ) &
          Vpc(imesh) = Vpc(imesh) * ( sin( const_sin * &
                       (iz - Len_to_Mesh_ele) ) ) ** 2

        if ( ( Mesh_to_Length < (qmcell(3,3) - smooth_potential%lelectrode ) )&
          .and. ( Mesh_to_Length > (qmcell(3,3) - smooth_plus_electrode) ) )  &
          Vpc(imesh) = Vpc(imesh) * ( sin( const_sin * &
                       (iz - ntm(3) -1 + Len_to_Mesh_ele) ) ) ** 2

        if ( ( Mesh_to_Length < smooth_potential%lelectrode ) .or. &
             ( Mesh_to_Length > (qmcell(3,3) - smooth_potential%lelectrode) )) &
          Vpc(imesh) = 0.0_grid_p
      enddo
      enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      write(6,*) "siesta-QMMM: The Vext was smoothed out at the region"//&
                 " of the electrodes."
      write(6,*)
    endif ! smooth_potential%on

    ! Fix units
    Vpc(:) = -2.0_grid_p * Vpc(:)
  end subroutine pcpot_ewald

    subroutine pcpot_direct( na_mm, rmm, pc, qmcell, ntpl, ntml, ntm, box0, &
                           Rho, Vpc)
    !! Calculates the potential using a cut-off scheme.
    use precision     , only : dp, grid_p
    use QMMM_core     , only : QMMM_density_rcut
    use QMMM_helper   , only : get_lattice_type, pbc_displ_vector
    use QMMM_neighbour, only : num_mmvecs, grid_veclist, grid_nr

    implicit none
    integer     , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer     , intent(in)  :: ntpl
      !! Total number of grid points.
    integer     , intent(in)  :: box0(3)
      !! The starting mesh cell index for this node.
    integer     , intent(in)  :: ntml(3)
      !! Number of mesh points along each axis for this node.
    integer     , intent(in)  :: ntm(3)
      !! Number of mesh points along each axis.
    real(dp)    , intent(in)  :: rmm(3,na_mm)
      !! Atomic positions.
    real(dp)    , intent(in)  :: pc(na_mm)
      !! Atomic (classical) partial charges.
    real(dp)    , intent(in)  :: qmcell(3,3)
      !! Unit cell vectors.
    real(grid_p), intent(in)  :: Rho(ntpl)
      !! Electron density in the grid.
    real(grid_p), intent(out) :: Vpc(ntpl)
      !! Electrostatic potential for the QM region.

    character :: lattice_type
    integer   :: ix, iy, iz, iat, js, icrd, ivec, imesh, ix0, iy0, iz0
    real(dp)  :: ym, xm, zm, dist, rcut_rho2, rcut_rho_A, drij(3), kcell(3,3)

    Vpc(:) = 0.0_grid_p
    if ( num_mmvecs == 0 ) return

    rcut_rho_A = QMMM_density_rcut()
    rcut_rho2  = rcut_rho_A * rcut_rho_A
    rcut_rho_A = 1.0_dp / rcut_rho_A
    lattice_type = get_lattice_type( qmcell )
    call reclat( qmcell, kcell, 0 )

    ! Loop over the mesh points
    do iz0 = 0, ntml(3) -1
    do iy0 = 0, ntml(2) -1
    do ix0 = 0, ntml(1) -1
      imesh = 1 + ix0 + ntml(1) * iy0 + ntml(1) * ntml(2) * iz0

      ! Only for grid points where there is density of charge
      if ( abs( Rho(imesh) ) > 0.0_dp ) then
        ix = ix0 + box0(1)
        iy = iy0 + box0(2)
        iz = iz0 + box0(3)

        xm = qmcell(1,1) * ix / ntm(1) + qmcell(1,2) * iy / ntm(2) &
           + qmcell(1,3) * iz / ntm(3)
        ym = qmcell(2,1) * ix / ntm(1) + qmcell(2,2) * iy / ntm(2) &
           + qmcell(2,3) * iz / ntm(3)
        zm = qmcell(3,1) * ix / ntm(1) + qmcell(3,2) * iy / ntm(2) &
           + qmcell(3,3) * iz / ntm(3)

        do iat = 1, num_mmvecs ! Loop over MM atoms
          js = grid_veclist(iat)

          if ( lattice_type == 'D' ) then
            drij(1) = xm - rmm(1,js) + grid_nr(1,iat) * qmcell(1,1)
            drij(2) = ym - rmm(2,js) + grid_nr(2,iat) * qmcell(2,2)
            drij(3) = zm - rmm(3,js) + grid_nr(3,iat) * qmcell(3,3)
          else
            drij(1) = xm - rmm(1,js)
            drij(2) = ym - rmm(2,js)
            drij(3) = zm - rmm(3,js)

            do icrd = 1, 3
            do ivec = 1, 3
              drij(icrd) = drij(icrd) + grid_nr(ivec,iat) * qmcell(icrd,ivec)
            enddo
            enddo
          endif
          call pbc_displ_vector( qmcell, kcell, drij )
          dist = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)

          if ( dist > rcut_rho2 ) then
            ! If our density point and point charge are too close,
            ! this might diverge. Thus, we avoid it.
            dist = 1.0_dp / sqrt( dist )
            Vpc(imesh) = Vpc(imesh) + pc(js) * dist
          else
            Vpc(imesh) = Vpc(imesh) + pc(js) * rcut_rho_A
          endif
        enddo ! MM atoms
      endif   ! abs(Rho(imesh)) > 0

    enddo     ! grid X
    enddo     ! grid Y
    enddo     ! grid Z

    ! Change units to Ry
    Vpc(:) = -2.0_grid_p * Vpc(:)
  end subroutine pcpot_direct
end module QMMM_pcpot_m
