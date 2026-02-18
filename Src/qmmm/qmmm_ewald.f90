module QMMM_ewald_m
  !! This module contains data structures and general options regarding
  !! Ewald summation in QM/MM calculations using SIESTA as the QM code.
  !
  ! See for more details:
  ! https://gitlab.com/siesta-project/analysis-tools/qmmm-driver/-/blob/master/Docs/Documentation-QMMM-CarlosSanz-2011.pdf?ref_type=heads
  !
  ! And the book "Understanding Molecular Dynamics", by Frenkel-Smit.
  !
  use precision, only: dp

  implicit none
  public :: QMMM_set_ewald

  private
  type QMMM_ewald_t
    !! Contains all the information regarding QM-MM Ewald
    !! summation for Coulomb interactions.
    real(dp) :: rcut  = 20.0_dp
      !! Ewald radius that separates direct and reciprocal
      !! summations. Used to set alpha and k-cutoff defaults.
    real(dp) :: alpha = 0.0_dp
      !! Ewald squared inverse of the gaussian width.
    real(dp) :: sqalpha = 0.0_dp
      !! Inverse of the gaussian widths.
    real(dp) :: kcut  = 0.0_dp
      !! Ewald cut-off radius in the reciprocal space.
    integer  :: k_nmax = 30
      !! Maximum number of points in the reciprocal space.
    integer  :: n_kbox(3)  = 0
      !! Total number of points in the reciprocal space.
    integer  :: n_kpoints  = 0
      !! Total number of non-zero points for the reciprocal summation.
    integer, pointer :: point_x(:)
      !! X coordinates (integer) for the reciprocal points.
    integer, pointer :: point_y(:)
      !! Y coordinates (integer) for the reciprocal points.
    integer, pointer :: point_z(:)
      !! Z coordinates (integer) for the reciprocal points.
    real(dp), pointer :: kx(:)
      !! X coordinates for the reciprocal points in reciprocal space.
    real(dp), pointer :: ky(:)
      !! Y coordinates for the reciprocal points in reciprocal space.
    real(dp), pointer :: kz(:)
      !! Z coordinates for the reciprocal points in reciprocal space.
    real(dp), pointer :: kmod2(:)
      !! Squared modulus of the current kpoint.
  contains
    procedure :: set_kgrid => set_ewald_grid
  end type QMMM_ewald_t

  type QMMM_smooth_t
    !! Smooth potential options for electrodes.
    logical  :: on = .false.
      !! Whether potential smoothing is activated.
    real(dp) :: lelectrode = 0.0_dp
      !! Electrode length
    real(dp) :: lsmooth    = 0.0_dp
      !! Smoothing length
  end type QMMM_smooth_t

  ! NOTE: The smoothing makes some assumptions about the system:
  !  * iz is the number of mesh points in z direction, and in this case, 360.
  !  * The electrode is in the region between points 1-60 and 300-360.
  !  * The smoothing is done by multiplying the final mesh by a quarter of
  !    a trigonometric function (going from 0 to 1).
  !  * The potential is set strictly to zero in the electrode region.
  !  * The size of the smoothing region is half of the electrode.

  type(QMMM_ewald_t) , public :: QMMM_ewald
  type(QMMM_smooth_t), public :: smooth_potential
contains

  subroutine QMMM_set_ewald( )
    !! Initialises Ewald summation variables and reads inputs.
    !
    ! The error of Ewald summation follows a shape of the type
    ! exp( -s^2 ) / s^2, for both reciprocal and real space summations.
    ! If we want to make both errors equal (for simplicity), we
    ! can arrive to the following relations (see bibliography):
    !
    ! * alpha = s^2 / rcut ^2
    ! * kcut = 2 * s * sqrt(alpha)
    !
    ! If we set s to 2.6, we guarantee an error of about 10^-4
    ! in both sums.
    use fdf      , only : fdf_get
    use precision, only : dp
    use units    , only : pi

    implicit none
    real(dp) :: sfactor

    QMMM_ewald%rcut = fdf_get( 'QMMM.Ewald.rcut', QMMM_ewald%rcut, 'Bohr' )

    sfactor    = 2.6_dp
    QMMM_ewald%alpha = sfactor * sfactor / ( QMMM_ewald%rcut * QMMM_ewald%rcut )
    QMMM_ewald%alpha = fdf_get( 'QMMM.Ewald.alpha', QMMM_ewald%alpha )

    QMMM_ewald%sqalpha = sqrt( QMMM_ewald%alpha )
    QMMM_ewald%kcut    = 2.0_dp * sfactor * QMMM_ewald%sqalpha
    QMMM_ewald%kcut    = fdf_get( 'QMMM.Ewald.kcut', QMMM_ewald%kcut )

    smooth_potential%on = fdf_get( 'QMMM.SmoothElectrode', smooth_potential%on )
    if ( smooth_potential%on ) then
      smooth_potential%lelectrode = &
        fdf_get( 'QMMM.SmoothElectrode.Electrode', smooth_potential%lelectrode,&
                 'Bohr' )
      smooth_potential%lsmooth = &
        fdf_get( 'QMMM.SmoothElectrode.Smooth', smooth_potential%lsmooth, &
                 'Bohr' )
    endif
  end subroutine QMMM_set_ewald

  subroutine set_ewald_grid( this, kcell, lattice_type )
    !! Calculates those points that belong to the reciprocal grid and stores
    !! their integer coordinates. The output is stored within this same
    !! module.
    use alloc      , only : de_alloc, re_alloc
    use parallel   , only : IOnode
    use precision  , only : dp
    use QMMM_helper, only : vec_norm
    use units      , only : pi

    implicit none
    class(QMMM_ewald_t), intent(inout) :: this

    real(dp) , intent(in) :: kcell(3,3)
    character, intent(in) :: lattice_type

    integer  :: n1, n2, n3, n1m, n2m, n3m, npoint, ni
    logical  :: skip_k
    real(dp) :: krecip(3), drij(3), twopi, kmod2, kcut2, k_neigh(6,3), &
                krecip_displaced(3), kmod2b

    kcut2 = this%kcut * this%kcut
    twopi = 2.0_dp * pi

    drij(:) = kcell(:,1)
    n1m = INT( QMMM_ewald%kcut / ( twopi * vec_norm(drij, 3) ) )
    if ( n1m > this%k_nmax ) call message( "WARNING", &
      'Ewald n1 > 30. Increasing the Ewald cut-off might be advisable.' )

    drij(:) = kcell(:,2)
    n2m = INT( QMMM_ewald%kcut / ( twopi * vec_norm(drij, 3) ) )
    if ( n2m > this%k_nmax ) call message( "WARNING", &
      'Ewald n2 > 30. Increasing the Ewald cut-off might be advisable.' )

    drij(:) = kcell(:,3)
    n3m = INT( QMMM_ewald%kcut / ( twopi * vec_norm(drij, 3) ) )
    if ( n3m > this%k_nmax ) call message( "WARNING", &
      'Ewald n3 > 30. Increasing the Ewald cut-off might be advisable.' )

    if ( IOnode ) &
      write(*,*) "Number of Ewald reciprocal points: ", n1m, n2m, n3m
    this%n_kbox(1) = n1m
    this%n_kbox(2) = n2m
    this%n_kbox(3) = n3m

    ! Fix for blindspots when lattice type is not D, i.e. not
    ! orthogonal. Essentially we store the 6 nearest neighbour cells,
    ! one in each direction. See later for how this is used.
    do ni = 1, 3
      k_neigh(ni, :) = 2 * twopi * (/ n1m * kcell(1, ni), &
                       n2m * kcell(2, ni), n3m * kcell(3, ni) /)
    enddo
    do ni = 4, 6
      k_neigh(ni, :) = -k_neigh(ni-3, :)
    enddo

    npoint = 0
    do n1 = -n1m, n1m
    do n2 = -n2m, n2m
    do n3 = -n3m, n3m
      if ( (n1 == 0) .and. (n2 == 0) .and. (n3 == 0) ) cycle

      if ( lattice_type == 'D' ) then
        krecip(1) = n1 * twopi * kcell(1,1)
        krecip(2) = n2 * twopi * kcell(2,2)
        krecip(3) = n3 * twopi * kcell(3,3)

        kmod2 = krecip(1) * krecip(1) + krecip(2) * krecip(2) &
              + krecip(3) * krecip(3)

        if ( kmod2 > kcut2 ) cycle
      else
        krecip(1) = twopi * ( n1 * kcell(1,1) + n2 * kcell(1,2) &
                            + n3 * kcell(1,3) )
        krecip(2) = twopi * ( n1 * kcell(2,1) + n2 * kcell(2,2) &
                            + n3 * kcell(2,3) )
        krecip(3) = twopi * ( n1 * kcell(3,1) + n2 * kcell(3,2) &
                            + n3 * kcell(3,3) )

        kmod2 = krecip(1) * krecip(1) + krecip(2) * krecip(2) &
              + krecip(3) * krecip(3)

        if ( kmod2 > kcut2 ) then
          ! Fix for "blindspots" in reciprocal space for hexagonal cells.
          ! There might be some k-points that end up outside the cut-off
          ! due to how the calculation above is done. Thus, we need to check
          ! whether the same k-point falls under the cut-off when displaced
          ! an entire (reciprocal) unit cell.
          skip_k = .true.
          do ni = 1, 6
            krecip_displaced = krecip - k_neigh(ni, :)
            kmod2 = krecip_displaced(1) * krecip_displaced(1) + &
                    krecip_displaced(2) * krecip_displaced(2) + &
                    krecip_displaced(3) * krecip_displaced(3)
            if ( kmod2 > kcut2 ) cycle
            krecip = krecip_displaced
            skip_k = .false.
            exit
          enddo

          if ( skip_k ) cycle
        endif
      endif
      npoint = npoint +1
    enddo
    enddo
    enddo

    this%n_kpoints = npoint

    call de_alloc(this%point_x, 'point_x', 'set_ewald_grid')
    call de_alloc(this%point_y, 'point_y', 'set_ewald_grid')
    call de_alloc(this%point_z, 'point_z', 'set_ewald_grid')
    call de_alloc(this%kx, 'kx', 'set_ewald_grid')
    call de_alloc(this%ky, 'ky', 'set_ewald_grid')
    call de_alloc(this%kz, 'kz', 'set_ewald_grid')
    call de_alloc(this%kmod2, 'kmod2', 'set_ewald_grid')

    nullify( this%point_x, this%point_y, this%point_z )
    nullify( this%kx, this%ky, this%kz, this%kmod2 )
    call re_alloc(this%point_x, 1, this%n_kpoints, 'point_x', 'set_ewald_grid')
    call re_alloc(this%point_y, 1, this%n_kpoints, 'point_y', 'set_ewald_grid')
    call re_alloc(this%point_z, 1, this%n_kpoints, 'point_z', 'set_ewald_grid')
    call re_alloc(this%kx, 1, this%n_kpoints, 'kx', 'set_ewald_grid')
    call re_alloc(this%ky, 1, this%n_kpoints, 'ky', 'set_ewald_grid')
    call re_alloc(this%kz, 1, this%n_kpoints, 'kz', 'set_ewald_grid')
    call re_alloc(this%kmod2, 1, this%n_kpoints, 'kmod2', 'set_ewald_grid')

    npoint = 0
    do n1 = -n1m, n1m
    do n2 = -n2m, n2m
    do n3 = -n3m, n3m
      if ( (n1 == 0) .and. (n2 == 0) .and. (n3 == 0) ) cycle

      if ( lattice_type == 'D' ) then
        krecip(1) = n1 * twopi * kcell(1,1)
        krecip(2) = n2 * twopi * kcell(2,2)
        krecip(3) = n3 * twopi * kcell(3,3)

        kmod2 = krecip(1) * krecip(1) + krecip(2) * krecip(2) &
              + krecip(3) * krecip(3)

        if ( kmod2 > kcut2 ) cycle
      else
        krecip(1) = twopi * ( n1 * kcell(1,1) + n2 * kcell(1,2) &
                            + n3 * kcell(1,3) )
        krecip(2) = twopi * ( n1 * kcell(2,1) + n2 * kcell(2,2) &
                            + n3 * kcell(2,3) )
        krecip(3) = twopi * ( n1 * kcell(3,1) + n2 * kcell(3,2) &
                            + n3 * kcell(3,3) )

        kmod2 = krecip(1) * krecip(1) + krecip(2) * krecip(2) &
              + krecip(3) * krecip(3)

        if ( kmod2 > kcut2 ) then
          ! Fix for "blindspots" in reciprocal space for hexagonal cells.
          skip_k = .true.
          do ni = 1, 6
            krecip_displaced = krecip - k_neigh(ni, :)
            kmod2b = krecip_displaced(1) * krecip_displaced(1) + &
                     krecip_displaced(2) * krecip_displaced(2) + &
                     krecip_displaced(3) * krecip_displaced(3)
            if ( kmod2b > kcut2 ) cycle
            krecip = krecip_displaced
            kmod2  = kmod2b
            skip_k = .false.
            exit
          enddo

          if ( skip_k ) cycle
        endif
      endif

      npoint = npoint +1
      this%point_x(npoint) = n1
      this%point_y(npoint) = n2
      this%point_z(npoint) = n3
      this%kx(npoint) = krecip(1)
      this%ky(npoint) = krecip(2)
      this%kz(npoint) = krecip(3)
      this%kmod2(npoint) = kmod2
    enddo
    enddo
    enddo

  end subroutine set_ewald_grid
end module QMMM_ewald_m
