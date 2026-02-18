module qmmm_pbc
  !! Subroutines used to deal with the PBC
  implicit none
  public :: pbc_displ_vector
  public :: get_pbc_vectors
  public :: get_lattice_type
  public :: reccel

contains

  subroutine pbc_displ_vector( lattice_type, cell, kcell, dr )
    use precision, only: dp

    implicit none
    character, intent(in)    :: lattice_type
    real(dp) , intent(in)    :: cell(3,3)
    real(dp) , intent(in)    :: kcell(3,3)
    real(dp) , intent(inout) :: dr(3)

    real(dp) :: icell(3,3)
    integer k, l, nk

    if ( lattice_type == 'D' ) then
      do k = 1, 3
        dr(k) = dr(k) - ANINT(dr(k) / cell(k,k)) * cell(k,k)
      enddo
    else
      ! We actually need the inverse of the unit cell vectors.
      ! We assume here that the reciprocal vectors matrix does not contain the
      ! 2*pi factor.
      do k = 1, 3
        nk = INT( ANINT(dr(1)*kcell(1,k) + dr(2)*kcell(2,k) + dr(3)*kcell(3,k)))
        if ( nk == 0 ) cycle

        do l = 1, 3
          dr(l) = dr(l) - cell(l,k) * nk
        enddo
      enddo
    endif
  end subroutine pbc_displ_vector

  subroutine get_pbc_vectors( lattice_type, cell, kcell, drij, nr )
    use precision, only: dp
    implicit none

    character, intent(in)  :: lattice_type
    real(dp) , intent(in)  :: cell(3,3)
    real(dp) , intent(in)  :: kcell(3,3)
    real(dp) , intent(in)  :: drij(3)
    integer  , intent(out) :: nr(3)

    integer :: k
    if ( lattice_type == 'D' ) then
      do k = 1, 3
        nr(k) = -INT( ANINT(drij(k) / cell(k,k)) )
      enddo
    else
      ! We actually need the inverse of the unit cell vectors.
      ! We assume here that the reciprocal vectors matrix does not contain the
      ! 2*pi factor.

      do k = 1, 3
        ! Here it is assumed that kcell is matrix of the reciprocal vectors
        ! without 2*pi factor
        nr(k) = -INT( ANINT(drij(1)*kcell(1,k) + drij(2)*kcell(2,k) + &
                            drij(3)*kcell(3,k)) )
      enddo
    endif
  end subroutine get_pbc_vectors

  character function get_lattice_type( cell )
    !! RETURN 'D' IF THE LATTICE VECTORS LIE ALONG THE X-, Y- AND Z-AXIS
    !! RESPECTIVELY. IN OTHER WORDS, THE CELL MATRIX IS DIAGONAL.
    !! RETURN 'G' OTHERWISE.
    use precision, only: dp

    implicit none
    real(dp), intent(in) :: cell(3,3)
    integer :: i, j

    get_lattice_type = 'D'
    do i = 1, 3
    do j = 1, 3
      if ( i /= j ) then
        if ( abs(cell(i,j)) > 1e-10_dp ) get_lattice_type = 'G'
      endif
    enddo
    enddo
  end function get_lattice_type

  subroutine reccel( N, A, B, IOPT )
    !! CALCULATES RECIPROCAL LATTICE VECTORS B.
    !! THEIR PRODUCT WITH DIRECT LATTICE VECTORS A IS 1 (IF IOPT=0) OR
    !! 2*PI (IF IOPT=1). N IS THE SPACE DIMENSION. WRITTEN BY J.M.SOLER.
    use precision, only : dp
    use sys      , only : die
    use units    , only : pi

    implicit none
    integer , intent(in)  :: n
    integer , intent(in)  :: iopt
    real(dp), intent(in)  :: A(N,N)
    real(dp), intent(out) :: B(N,N)
    real(dp) :: c, ci
    integer  :: i

    C = 1.0_dp
    if ( IOPT == 1 ) C = 2.0_dp * pi

    select case ( N )
    case ( 1 )
      B(1,1) = C / A(1,1)
    case ( 2 )
      C = C / (A(1,1) * A(2,2) - A(1,2) * A(2,1))
      B(1,1) =  A(2,2) * C
      B(1,2) = -A(2,1) * C
      B(2,1) = -A(1,2) * C
      B(2,2) =  A(1,1) * C
    case ( 3 )
      B(1,1) = A(2,2) * A(3,3) - A(3,2) * A(2,3)
      B(2,1) = A(3,2) * A(1,3) - A(1,2) * A(3,3)
      B(3,1) = A(1,2) * A(2,3) - A(2,2) * A(1,3)
      B(1,2) = A(2,3) * A(3,1) - A(3,3) * A(2,1)
      B(2,2) = A(3,3) * A(1,1) - A(1,3) * A(3,1)
      B(3,2) = A(1,3) * A(2,1) - A(2,3) * A(1,1)
      B(1,3) = A(2,1) * A(3,2) - A(3,1) * A(2,2)
      B(2,3) = A(3,1) * A(1,2) - A(1,1) * A(3,2)
      B(3,3) = A(1,1) * A(2,2) - A(2,1) * A(1,2)
      do i = 1, 3
        CI = C / ( A(1,I) * B(1,I) + A(2,I) * B(2,I) + A(3,I) * B(3,I) )
        B(:,I) = B(:,I) * CI
      enddo
    case default
      call die( 'RECCEL: NOT PREPARED FOR N>3.' )
    end select
  end subroutine reccel
end module qmmm_pbc