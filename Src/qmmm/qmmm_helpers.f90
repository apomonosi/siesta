module QMMM_helper
  !! This module contains auxiliary subroutines and functions used
  !! in QM/MM calculations.
  public :: vec_norm
  public :: scalar_vec
  public :: get_lattice_type
  public :: get_pbc_vectors
  public :: pbc_displ_vector

contains

  function vec_norm( v1, vsize ) result( vnorm )
    !! Calculates the cartesian norm of a vector.
    use precision, only : dp
    implicit none
    integer , intent(in) :: vsize
      !! Vector sizes of v1 and v2.
    real(dp), intent(in) :: v1(vsize)
      !! The vector.

    real(dp) :: vnorm

    vnorm = sqrt( scalar_vec( v1, v1, vsize ) )
  end function vec_norm

  function scalar_vec( v1, v2, vsize ) result( scalar_t )
    !! Calculates the scalar product between two vectors of any dimension.
    use precision, only : dp
    implicit none
    integer , intent(in) :: vsize
      !! Vector sizes of v1 and v2.
    real(dp), intent(in) :: v1(vsize)
      !! The first vector.
    real(dp), intent(in) :: v2(vsize)
      !! The second vector.

    integer  :: iCrd
    real(dp) :: scalar_t

    scalar_t = 0.0_dp
    do iCrd = 1, vsize
      scalar_t = scalar_t + v1(iCrd) * v2(iCrd)
    enddo
  end function scalar_vec

  function get_lattice_type( cell ) result ( cell_t )
    !! Returns "D" if the lattice vectors lie along the X, Y and Z axis
    !! respectively, i.e. if the cell matrix is diagonal. Returns "G"
    !! otherwise
    use precision, only: dp

    implicit none
    real(dp), intent(in) :: cell(3,3)
      !! Cell vectors.
    character            :: cell_t
    integer              :: iCrd, jCrd

    cell_t = 'D'
    do iCrd = 1, 3
    do jCrd = 1, 3
      if ( iCrd /= jCrd ) then
        if ( abs(cell(iCrd,jCrd)) > 1e-14_dp ) cell_t = 'G'
      endif
    enddo
    enddo
  end function get_lattice_type

  subroutine get_pbc_vectors( cell, kcell, drij, nr )
    !! Gets the number of reciprocal cells scanned in
    !! each direction by a given drij vector.
    use precision, only: dp

    implicit none
    real(dp), intent(in)  :: cell(3,3)
      !! The direct cell vectors.
    real(dp), intent(in)  :: kcell(3,3)
      !! The reciprocal cell vectors.
    real(dp), intent(in)  :: drij(3)
      !! The input vector used for scanning.
    integer , intent(out) :: nr(3)
      !! The number of cells scanned by the input vector.

    character :: lattice_type
    real(dp)  :: nr_r(3)

    lattice_type = get_lattice_type( cell )
    if ( lattice_type == 'D' ) then
      nr_r(1) = ANINT( drij(1) / cell(1,1) )
      nr_r(2) = ANINT( drij(2) / cell(2,2) )
      nr_r(3) = ANINT( drij(3) / cell(3,3) )
    else
      ! Here it is assumed that kcell is matrix of the
      ! reciprocal vectors without 2*pi factor.
      nr_r(1) = ANINT( drij(1) * kcell(1,1) + drij(2) * kcell(2,1) + &
                       drij(3) * kcell(3,1) )
      nr_r(2) = ANINT( drij(1) * kcell(1,2) + drij(2) * kcell(2,2) + &
                       drij(3) * kcell(3,2) )
      nr_r(3) = ANINT( drij(1) * kcell(1,3) + drij(2) * kcell(2,3) + &
                       drij(3) * kcell(3,3) )
    endif
    nr(:) = - INT( nr_r(:) )
  end subroutine get_pbc_vectors

  subroutine pbc_displ_vector( cell, kcell, drij )
    !! Readequates a distance vector drij taking into account periodic
    !! boundary conditions.

    use precision, only: dp
    implicit none
    real(dp), intent(in)    :: cell(3,3)
      !! The direct cell vectors.
    real(dp), intent(in)    :: kcell(3,3)
      !! The reciprocal cell vectors.
    real(dp), intent(inout) :: drij(3)
      !! The input vector used for scanning.

    character :: lattice_type
    integer   :: iCrd, iVec, nk

    lattice_type = get_lattice_type( cell )
    if ( lattice_type == 'D' ) then
      do iCrd = 1, 3
        drij(iCrd) = drij(iCrd) - ANINT( drij(iCrd) / cell(iCrd,iCrd)) &
                                * cell(iCrd,iCrd)
      enddo
    else
      do iCrd = 1, 3
        ! Here it is assumed that kcell is matrix of the reciprocal vectors
        ! without 2*pi factor.
        nk = INT( ANINT( drij(1) * kcell(1,iCrd) + drij(2) * kcell(2,iCrd) &
                  + drij(3) * kcell(3,iCrd) ) )
        if ( nk == 0 ) cycle

        do iVec = 1,3
          drij(iCrd) = drij(iCrd) - cell(iVec,iCrd) * nk
        enddo
      enddo
    endif
  end subroutine pbc_displ_vector

end module QMMM_helper