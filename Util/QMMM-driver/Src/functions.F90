module functions
  !! This module contains several common geometry-related functions
  !! used by the program, namely distances, angles and dihedrals.
  use precision, only : dp

  implicit none
  public :: dist2_v2
  public :: angle_v2
  public :: scalar_v2
  public :: norm_v2
  public :: scalar_vec
  public :: vec_norm
  public :: vec_norm2
  public :: dihedral_v
  public :: dihedral_rescaled
  public :: dihedral_gradient

  private
  real(dp), parameter :: TINY_VAL = 1.0e-14_dp
    !! Used when comparing a result to zero.

contains

  function dist2_v2( r1, r2 ) result( dv2 )
    !! Calculates the squared norm of a distance between two points.
    use precision, only : dp

    implicit none
    real(dp) :: r1(3)
      !! Point 1.
    real(dp) :: r2(3)
      !! Point 2.
    real(dp) :: dv2
    real(dp) :: dr(3)

    dr  = r1 - r2
    dv2 = scalar_v2( dr, dr )
  end function dist2_v2

  function angle_v2( dr12, dr32 ) result( angle_t )
    !! Calculates the angle between two vectors sharing a point.
    use precision, only : dp
    use mm_units , only : rads
    implicit none

    real(dp) :: dr12(3)
      !! Vector 1.
    real(dp) :: dr32(3)
      !! Vector 2.
    real(dp) :: angle_t, scalarp

    scalarp = dr12(1) * dr32(1) + dr12(2) * dr32(2) + dr12(3) * dr32(3)
    angle_t = norm_v2( dr12 ) * norm_v2( dr32 )

    angle_t = scalarp / angle_t
    if ( angle_t >  1.0_dp ) angle_t =  1.0_dp
    if ( angle_t < -1.0_dp ) angle_t = -1.0_dp

    angle_t = acos( angle_t ) * rads
  end function angle_v2

  function scalar_v2( r1, r2 ) result( scalar_t )
    !! Calculates the scalar product between two vectors.
    use precision, only : dp
    implicit none

    real(dp), intent(in) :: r1(3)
      !! First vector.
    real(dp), intent(in) :: r2(3)
      !! Second vector.

    real(dp) :: scalar_t
    scalar_t = r1(1) * r2(1) + r1(2) * r2(2) + r1(3) * r2(3)
  end function scalar_v2

  function norm_v2( r1 ) result( scalar_t )
    !! Calculates the norm of a vector.
    use precision, only : dp
    implicit none

    real(dp), intent(in) :: r1(3)
      !! Vector to calculate the norm.

    real(dp) :: scalar_t
    scalar_t = sqrt( scalar_v2(r1, r1) )
  end function norm_v2

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

  function vec_norm2( v1, vsize ) result( vnorm )
    !! Calculates the squared cartesian norm of a vector.
    use precision, only : dp
    implicit none
    integer , intent(in) :: vsize
      !! Vector sizes of v1 and v2.
    real(dp), intent(in) :: v1(vsize)
      !! The vector.

    real(dp) :: vnorm

    vnorm = scalar_vec( v1, v1, vsize )
  end function vec_norm2

  function dihedral_v( r12, r32, r43, mvec_o, nvec_o, mxn_o, mnorm_o, &
                       nnorm_o ) result( dihe_t )
    !! Calculates the dihedral angle taking the three vectors between
    !! four points. Optionally, it may also return the rotors, their
    !! norms, and their absolute values.
    use mm_units , only : rads
    use precision, only : dp
    use sys      , only : die

    implicit none
    real(dp), intent(in)  :: r12(3)
      !! Distance vector between atoms 1 and 2.
    real(dp), intent(in)  :: r32(3)
      !! Distance vector between atoms 3 and 2.
    real(dp), intent(in)  :: r43(3)
      !! Distance vector between atoms 4 and 3.
    real(dp), optional, intent(out) :: mvec_o(3)
      !! M, rotor r12 x r32.
    real(dp), optional, intent(out) :: nvec_o(3)
      !! N, rotor r43 x r32.
    real(dp), optional, intent(out) :: mxn_o
      !! Scalar product between M and N.
    real(dp), optional, intent(out) :: mnorm_o
      !! Norm of M.
    real(dp), optional, intent(out) :: nnorm_o
      !! Norm of N.
    real(dp) :: dihe_t
    real(dp) :: cos_dihe, mvec(3), nvec(3), mxn, mnorm, nnorm

    mvec(1) = r12(2) * r32(3) - r12(3) * r32(2)
    mvec(2) = r12(3) * r32(1) - r12(1) * r32(3)
    mvec(3) = r12(1) * r32(2) - r12(2) * r32(1)

    nvec(1) = r32(3) * r43(2) - r32(2) * r43(3)
    nvec(2) = r32(1) * r43(3) - r32(3) * r43(1)
    nvec(3) = r32(2) * r43(1) - r32(1) * r43(2)

    mxn   = scalar_v2( mvec, nvec )
    mnorm = norm_v2( mvec )
    nnorm = norm_v2( nvec )

    cos_dihe = mxn / ( mnorm * nnorm )
    if ( abs( cos_dihe ) > 1.00001_dp ) &
      call die( 'dihedral_energy: abs(cos_dihe) > 1.00001' )

    if ( cos_dihe >  1.0_dp ) cos_dihe =  1.0_dp
    if ( cos_dihe < -1.0_dp ) cos_dihe = -1.0_dp

    dihe_t = acos( cos_dihe ) * rads

    if ( present( mvec_o  ) ) mvec_o  = mvec
    if ( present( nvec_o  ) ) nvec_o  = nvec
    if ( present( mxn_o   ) ) mxn_o   = mxn
    if ( present( mnorm_o ) ) mnorm_o = mnorm
    if ( present( nnorm_o ) ) nnorm_o = nnorm
  end function dihedral_v

  function dihedral_rescaled( r12, r32, r43, r31, r41, mvec_o, nvec_o, &
                              mxn_o, mnorm_o, nnorm_o ) result( dihe_t )
    !! Calculates the dihedral angle given by four points, taking into
    !! account periodic boundary conditions and additional corrections.
    use mm_units , only : rads
    use precision, only : dp
    use sys      , only : die

    implicit none
    real(dp), intent(in)  :: r12(3)
      !! Distance vector between atoms 1 and 2.
    real(dp), intent(in)  :: r32(3)
      !! Distance vector between atoms 3 and 2.
    real(dp), intent(in)  :: r43(3)
      !! Distance vector between atoms 4 and 3.
    real(dp), intent(in)  :: r31(3)
      !! Distance vector between atoms 3 and 1.
    real(dp), intent(in)  :: r41(3)
      !! Distance vector between atoms 3 and 1.
    real(dp), optional, intent(out) :: mvec_o(3)
      !! M, rotor r12 x r32.
    real(dp), optional, intent(out) :: nvec_o(3)
      !! N, rotor r43 x r32.
    real(dp), optional, intent(out) :: mxn_o
      !! Scalar product between M and N.
    real(dp), optional, intent(out) :: mnorm_o
      !! Norm of M.
    real(dp), optional, intent(out) :: nnorm_o
      !! Norm of N.
    real(dp) :: dihe_t
    real(dp) :: cos_dihe, mvec(3), nvec(3), mxn, mnorm, nnorm, plane(3), l1, l2

    mvec(1) = r12(2) * r32(3) - r12(3) * r32(2)
    mvec(2) = r12(3) * r32(1) - r12(1) * r32(3)
    mvec(3) = r12(1) * r32(2) - r12(2) * r32(1)

    nvec(1) = r32(3) * r43(2) - r32(2) * r43(3)
    nvec(2) = r32(1) * r43(3) - r32(3) * r43(1)
    nvec(3) = r32(2) * r43(1) - r32(1) * r43(2)

    mxn   = scalar_v2( mvec, nvec )
    mnorm = norm_v2( mvec )
    nnorm = norm_v2( nvec )

    cos_dihe = mxn / ( mnorm * nnorm )
    if ( abs( cos_dihe ) > 1.00001_dp ) &
      call die( 'dihedral_energy: abs(cos_dihe) > 1.00001' )

    if ( cos_dihe >  1.0_dp ) cos_dihe =  1.0_dp
    if ( cos_dihe < -1.0_dp ) cos_dihe = -1.0_dp

    dihe_t = acos( cos_dihe ) * rads


    plane(1) = r12(3) * r31(2) - r12(2) * r31(3)
    plane(2) = r12(1) * r31(3) - r12(3) * r31(1)
    plane(3) = r12(2) * r31(1) - r12(1) * r31(2)

    ! Distance(l) from atom4 to the ABCD = 0 plane
    l1 = plane(1) * r41(1) + plane(2) * r41(2) + plane(3) * r41(3)
    l2 = sqrt( plane(1) * plane(1) + plane(2) * plane(2) &
             + plane(3) * plane(3) )

    if ( ( l1 / l2 ) < 0.0_dp ) dihe_t = 360.0_dp - dihe_t

    if ( present( mvec_o  ) ) mvec_o  = mvec
    if ( present( nvec_o  ) ) nvec_o  = nvec
    if ( present( mxn_o   ) ) mxn_o   = mxn
    if ( present( mnorm_o ) ) mnorm_o = mnorm
    if ( present( nnorm_o ) ) nnorm_o = nnorm
  end function dihedral_rescaled

  subroutine dihedral_gradient( dr12, dr32, dr43, dr31, dr42, &
                                atom, prefactor, fce )
    !! Calculates the force over a dihedral angle in a more general case.
    use precision, only : dp

    implicit none
    real(dp), intent(in)  :: dr12(3)
      !! Distance vector between atoms 1 and 2.
    real(dp), intent(in)  :: dr32(3)
      !! Distance vector between atoms 3 and 2.
    real(dp), intent(in)  :: dr43(3)
      !! Distance vector between atoms 4 and 3.
    real(dp), intent(in)  :: dr31(3)
      !! Distance vector between atoms 3 and 1.
    real(dp), intent(in)  :: dr42(3)
      !! Distance vector between atoms 4 and 2.
    integer , intent(in)  :: atom
      !! Atom of the dihedral we are interested in for the derivative.
      !! Goes from 1 to 4.
    real(dp), intent(in)  :: prefactor
      !! Externally calculated prefactor. Must contain information
      !! k * ( ang - ang_eq ).
    real(dp), intent(out) :: fce(3)
      !! Forces the chosen atom.

    integer  :: iCrd
    real(dp) :: mvec(3), nvec(3), mxn, nnorm, mnorm, &
                dmr(3), dnr(3), dmn, prue, dtot, dscalar, dihe_t

    dihe_t = dihedral_v( dr12, dr32, dr43, mvec, nvec, mxn, mnorm, nnorm )

    fce(:) = 0.0_dp
    prue = mxn / ( mnorm * nnorm )
    prue = ( 1.0_dp - prue * prue )
    if ( abs(prue) < TINY_VAL ) return

    dtot = - prefactor / sqrt(prue)

    do iCrd = 1, 3
      dmr = 0.0_dp
      dnr = 0.0_dp

      select case ( atom )
      case ( 1 )
        select case ( icrd )
        case ( 1 )
          dmr(2) = -dr32(3) ! dz23
          dmr(3) =  dr32(2) ! dy32
        case ( 2 )
          dmr(1) =  dr32(3) ! dz32
          dmr(3) = -dr32(1) ! dx23
        case ( 3 )
          dmr(1) = -dr32(2) ! dy23
          dmr(2) =  dr32(1) ! dx32
        end select

      case ( 2 )
        select case ( icrd )
        case ( 1 )
          dmr(2) =  dr31(3) ! dz31
          dmr(3) = -dr31(2) ! dy13
          dnr(2) = -dr43(3) ! dz34
          dnr(3) =  dr43(2) ! dy43
        case ( 2 )
          dmr(1) = -dr31(3) ! dz13
          dmr(3) =  dr31(1) ! dx31
          dnr(1) =  dr43(3) ! dz43
          dnr(3) = -dr43(1) ! dx34
        case ( 3 )
          dmr(1) =  dr31(2) ! dy31
          dmr(2) = -dr31(1) ! dx13
          dnr(1) = -dr43(2) ! dy34
          dnr(2) =  dr43(1) ! dx43
        end select

      case ( 3 )
        select case ( icrd )
        case ( 1 )
          dmr(2) =  dr12(3) ! dz12
          dmr(3) = -dr12(2) ! dy21
          dnr(2) =  dr42(3) ! dz42
          dnr(3) = -dr42(2) ! dy24
        case ( 2 )
          dmr(1) = -dr12(3) ! dz21
          dmr(3) =  dr12(1) ! dx12
          dnr(1) = -dr42(3) ! dz24
          dnr(3) =  dr42(1) ! dx42
        case ( 3 )
          dmr(1) =  dr12(2) ! dy12
          dmr(2) = -dr12(1) ! dx21
          dnr(1) =  dr42(2) ! dy42
          dnr(2) = -dr42(1) ! dx24
        end select

      case ( 4 )
        select case ( icrd )
        case ( 1 )
          dnr(2) = -dr32(3) ! dz23
          dnr(3) =  dr32(2) ! dy32
        case ( 2 )
          dnr(1) =  dr32(3) ! dz32
          dnr(3) = -dr32(1) ! dx23
        case ( 3 )
          dnr(1) = -dr32(2) ! dy23
          dnr(2) =  dr32(1) ! dx32
        end select

      end select

      dmn = ( mnorm / nnorm ) * scalar_v2( nvec, dnr ) &
          + ( nnorm / mnorm ) * scalar_v2( mvec, dmr )

      dscalar = scalar_v2( nvec, dmr ) + scalar_v2( mvec, dnr )

      fce(iCrd)  = dtot * ( dscalar * mnorm * nnorm - dmn * mxn ) &
                        / ( nnorm * nnorm * mnorm * mnorm )
    enddo
  end subroutine dihedral_gradient

end module functions
