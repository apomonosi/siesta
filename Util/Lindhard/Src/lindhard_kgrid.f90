module lindhard_kgrid_m
  !! Contains routines to deal with the kpoint-grid.
  use precision, only : dp

  implicit none
  private
  public :: setup_kgrid
  public :: evaluate
  public :: kpoint_crd

contains

  subroutine setup_kgrid( ucell, gridk, dkg, igmin, ng )
    use precision, only : dp

    implicit none
    real(dp), intent(in)  :: ucell(3,3)
      !! Unit cell vectors.
    real(dp), intent(out) :: gridk(3,3)
      !! Vectors defining the k-grid parallelepiped.
    real(dp), intent(out) :: dkg(3)
      !! K grid origin.
    integer , intent(out) :: igmin(3)
      !! Minimun k grid index in each direction.
    integer , intent(out) :: ng(3)
      !! Number of grid divisions of the i reciprocal lattice vector (RLV)

    logical  :: firm_displ
    integer  :: kscell(3,3)
    real(dp) :: kdispl(3)

    firm_displ = .false.
    call setup_kscell( ucell, kscell, kdispl, firm_displ )
    call find_kgrid( ucell, kscell, kdispl, firm_displ, gridk, dkg, igmin, ng )

  end subroutine setup_kgrid

  subroutine setup_kscell( cell, kscell, kdispl, firm_displ )
    !! The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
    !! If both are present, kgrid_Monkhorst_Pack has priority. If none is
    !! present, the cutoff default is zero, producing only the gamma point.
    !
    ! Examples of fdf data specifications:
    !  kgrid_cutoff  50. Bohr
    !  %block kgrid_Monkhorst_Pack  # Defines kscell and kdispl
    !     4  0  0   0.50               # (kscell(i,1),i=1,3), kdispl(1)
    !     0  4  0   0.50               # (kscell(i,2),i=1,3), kdispl(2)
    !     0  0  4   0.50               # (kscell(i,3),i=1,3), kdispl(3)
    !  %endblock kgrid_Monkhorst_Pack
    use fdf      , only : fdf_block, fdf_bline, fdf_bintegers, fdf_bnvalues, &
                          fdf_bvalues, fdf_get, block_fdf, parsed_line
    use m_minvec , only : minvec
    use precision, only : dp
    use sys      , only : die

    implicit none
    real(dp), intent(in)  :: cell(3,3)
      !! Unit cell vectors in real space cell(ixyz,ivec)
    integer , intent(out) :: kscell(3,3)
      !! Supercell in terms of k-grid unit cell.
    real(dp), intent(out) :: kdispl(3)
      !! Grid origin in k-grid-vector coordinates.
    logical , intent(out) :: firm_displ
      !! Whether we use a displacement provided via Monkhorst-Pack block.

    integer           i, j, factor(3,3), expansion_factor
    real(dp)          scmin(3,3),  vmod, cutoff
    real(dp)          ctransf(3,3)
    logical           mp_input
    type(block_fdf)   :: kblock
    type(parsed_line),pointer :: pline

    real(dp), parameter :: defcut = 0.0_dp
    integer, dimension(3,3), parameter :: unit_matrix =  &
                         reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

    mp_input = fdf_block('kgrid_Monkhorst_Pack',kblock)
    if ( mp_input ) then
      do i = 1,3
        if (.not. fdf_bline(kblock,pline)) &
          call die('kgridinit: ERROR in kgrid_Monkhorst_Pack block')

        kscell(1,i) = fdf_bintegers(pline,1)
        kscell(2,i) = fdf_bintegers(pline,2)
        kscell(3,i) = fdf_bintegers(pline,3)
        if ( fdf_bnvalues(pline) > 3 ) then
          kdispl(i) = mod(fdf_bvalues(pline,4), 1._dp)
        else
          kdispl(i) = 0._dp
        end if
      enddo
      firm_displ = .true.

    else
      cutoff = fdf_get('kgrid_cutoff',defcut,'Bohr')

      kdispl(1:3) = 0.0_dp  ! Might be changed later
      firm_displ = .false.

      ! Find equivalent rounded unit-cell
      call minvec( cell, scmin, ctransf )

      expansion_factor = 1
      do j = 1,3
        factor(j,1:3) = 0
        vmod = sqrt(dot_product(scmin(1:3,j),scmin(1:3,j)))
        factor(j,j) = int(2.0_dp*cutoff/vmod) + 1
        expansion_factor = expansion_factor * factor(j,j)
      enddo

      ! Generate actual supercell skeleton
      kscell = matmul(int(ctransf), factor)

      ! Avoid confusing permutations
      if (expansion_factor == 1) kscell = unit_matrix
    endif
  end subroutine setup_kscell

  subroutine find_kgrid( cell, kscell, displ, firm_displ, gridk, dkg, igmin, ng)
    !! Finds origin, size and indexes for k-point grid.
    use units,    only : pi
    use m_minvec, only : minvec

    implicit none
    integer, intent(in)     :: kscell(3,3)
      !! Supercell in terms of k-grid unit cell.
      !! scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
    real(dp), intent(in)    :: cell(3,3)
      !! Unit cell vectors in real space cell(ixyz,ivec)
    logical , intent(in)    :: firm_displ
      !! Whether we use a displacement provided via Monkhorst-Pack block.
    real(dp), intent(inout) :: displ(3)
      !! Grid origin in k-grid-vector coordinates:
      !! origin(ix) = sum_j gridk(ix,j)*displ(j)
    real(dp), intent(out)   :: gridk(3,3)
      !! Vectors defining the k-grid parallelepiped.
    real(dp), intent(out)   :: dkg(3)
      !! K grid origin.
    integer , intent(out)   :: igmin(3)
      !! Minimun k grid index in each direction.
    integer , intent(out)   :: ng(3)
      !! Number of grid divisions of the i reciprocal lattice vector (RLV)


    external :: idiag, reclat
    integer  :: i, j, kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3), nktot, proj(3,3)
    real(dp) :: dkx(3), dscell(3,3), gscell(3,3), scell(3,3)

    real(dp), parameter :: tiny = 1.e-12_dp

    gridk(:,:) = 0.0_dp
    dkg(:)     = 0.0_dp
    igmin(:)   = 0
    ng(:)      = 0

    !Find total number of points (determinant of kscell)
    nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) + &
                 kscell(2,1) * kscell(3,2) * kscell(1,3) + &
                 kscell(3,1) * kscell(1,2) * kscell(2,3) - &
                 kscell(1,1) * kscell(3,2) * kscell(2,3) - &
                 kscell(2,1) * kscell(1,2) * kscell(3,3) - &
                 kscell(3,1) * kscell(2,2) * kscell(1,3) )

    ! Find k-grid supercell
    do i = 1,3
      do j = 1,3
        scell(j,i) = cell(j,1) * kscell(1,i) + cell(j,2) * kscell(2,i) +&
                     cell(j,3) * kscell(3,i)
      enddo
    enddo

    ! Equivalent supercell DA with the property that there exists a primitive
    ! cell (pa') such that DA_i = N_i*pa'_i (See Moreno and Soler)
    call idiag( 3, kscell, kdsc, ml, mr, maux )

    proj(:,:) = 0  ! Possible sign changes
    do i = 1, 3
      proj(i,i) = 1
      if (kdsc(i,i) < 0) proj(i,i) = -1
    enddo
    kdsc = matmul( kdsc, proj )
    mr   = matmul( mr  , proj )

    ! Set the displacements if not firm (i.e., specified by the user). Even if
    ! firm, warn if a better choice is actually possible.
    do j = 1, 3
      if (mod(kdsc(j,j),2) == 0) then
        if ( firm_displ .and. (abs(displ(j) - 0.5_dp) > tiny) ) then
          write(6,"(a,i4,a,2f8.2)") "k-point displ. along", j, &
                                    " input, could be: ", displ(j), 0.5_dp
        else
          displ(j) = 0.5_dp
        endif
      else
        if ( firm_displ .and. (abs(displ(j) - 0.0_dp) > tiny) ) then
          write(6,"(a,i4,a,2f8.2)") "k-point displ. along", j, &
                                    " input, could be: ", displ(j), 0.0d0
        else
          displ(j) = 0.0_dp
        endif
      endif
    enddo

    dscell = matmul(scell,mr)
    ! Find k-grid unit vectors
    call reclat( dscell, gridk, 1 )

    ! Find grid origin in cartesian coordinates
    call reclat( scell, gscell, 1 )
    do j = 1,3
      dkx(j) = gscell(j,1) * displ(1) + gscell(j,2) * displ(2) + &
               gscell(j,3) * displ(3)
    enddo

    ! Find grid origin in gridk coordinates
    do i = 1,3
      dkg(i) = ( dkx(1) * dscell(1,i) + dkx(2) * dscell(2,i) + &
                 dkx(3) * dscell(3,i) ) / (2.0_dp*pi)
    enddo

    ! Find total range of grid indexes
    do j = 1,3
      ng(j)    = kdsc(j,j)
      igmin(j) = -( (ng(j)-1) / 2)
    enddo
  end subroutine find_kgrid

  subroutine evaluate( d, nkx, nky, nkz, xfrac, val )
    !! Interpolates energy on grid points.
    use precision, only : dp

    implicit none
    integer , intent(in)  :: nkx
      !! Number of kpoints in the x direction.
    integer , intent(in)  :: nky
      !! Number of kpoints in the y direction.
    integer , intent(in)  :: nkz
      !! Number of kpoints in the z direction.
    real(dp), intent(in)  :: d(0:nkx-1,0:nky-1,0:nkz-1)
      !! Energy grid.
    real(dp), intent(in)  :: xfrac(3)
      !! Reduced coordinates of point.
    real(dp), intent(out) :: val
      !! Interpolated value of the energy.

    integer  :: n(3), lo(3), hi(3), k
    real(dp) :: r(3), x(3), y(3), nk

    n(1) = nkx
    n(2) = nky
    n(3) = nkz
    ! Find the right 3D "grid cube" and the reduced coordinates of the point in
    ! it. The double mod assures that negative numbers are well treated (the
    ! idea is to bring the coordinates to the [0,n(k)) interval).
    do k = 1, 3
      nk    = real( n(k), kind=dp )
      r(k)  = modulo( n(k)*xfrac(k), nk )
      lo(k) = int( r(k) )
      hi(k) = mod ( lo(k)+1, n(k) )
      x(k)  = r(k) - lo(k)
      y(k)  = 1 - x(k)
    enddo

    ! Compute charge density by linear interpolation.
    val = d(lo(1),lo(2),lo(3)) * y(1) * y(2) * y(3) + &
          d(lo(1),lo(2),hi(3)) * y(1) * y(2) * x(3) + &
          d(lo(1),hi(2),lo(3)) * y(1) * x(2) * y(3) + &
          d(lo(1),hi(2),hi(3)) * y(1) * x(2) * x(3) + &
          d(hi(1),lo(2),lo(3)) * x(1) * y(2) * y(3) + &
          d(hi(1),lo(2),hi(3)) * x(1) * y(2) * x(3) + &
          d(hi(1),hi(2),lo(3)) * x(1) * x(2) * y(3) + &
          d(hi(1),hi(2),hi(3)) * x(1) * x(2) * x(3)
  end subroutine evaluate

  subroutine kpoint_crd( kv, kvmin, nkx, nky, nkz, ix, iy, iz, k )
    !! Computes k-point coordinates from indices in x,y,z directions.
    use precision, only : dp

    implicit none
    integer , intent(in)  :: nkx
      !! Number of kpoints in the x direction.
    integer , intent(in)  :: nky
      !! Number of kpoints in the y direction.
    integer , intent(in)  :: nkz
      !! Number of kpoints in the z direction.
    integer , intent(in)  :: ix
      !! Coordinate in the x direction.
    integer , intent(in)  :: iy
      !! Coordinate in the y direction.
    integer , intent(in)  :: iz
      !! Coordinate in the z direction.
    real(dp), intent(in)  :: kvmin(3)
      !! Minimum value of k-grid.
    real(dp), intent(in)  :: kv(3,3)
      !! Cartesian vectors.
    real(dp), intent(out) :: k(3)
      !! Output kpoint.

    integer  :: jcrd
    real(dp) :: v(3)

    v(1) = real( ix-1, kind=dp ) / real( nkx, kind=dp )
    v(2) = real( iy-1, kind=dp ) / real( nky, kind=dp )
    v(3) = real( iz-1, kind=dp ) / real( nkz, kind=dp )

    k(:) = kvmin(:)
    do jcrd = 1, 3
      k(:) = k(:) + kv(:,jcrd) * v(jcrd)
    enddo
  end subroutine kpoint_crd
end module lindhard_kgrid_m