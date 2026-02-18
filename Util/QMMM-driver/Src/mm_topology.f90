module mm_topology
  !! This module contains all of the data structures related to
  !! classical MM parameters in general.
  use precision, only: dp

  implicit none
  private

  type :: ffbond_t
    character(len=5) :: type
    real(dp)         :: k
    real(dp)         :: r_eq

  contains
    procedure, pass :: energy => bond_energy
    procedure, pass :: grad   => bond_gradient
  end type ffbond_t

  type :: ffangle_t
    character(len=8) :: type
    real(dp)         :: k
    real(dp)         :: r_eq

  contains
    procedure, pass :: energy => angle_energy
    procedure, pass :: grad_e => angle_gradient_extreme
    procedure, pass :: grad_m => angle_gradient_middle
  end type ffangle_t

  type :: ffdihedral_t
    character(len=11) :: type
    real(dp)          :: k
    real(dp)          :: r_eq
    real(dp)          :: per
    integer           :: multi

  contains
    procedure, pass :: energy => dihedral_energy
    procedure, pass :: grad   => dihedral_gradient
  end type ffdihedral_t

  type, public :: fftopology
    integer :: nbonds
    integer :: nangles
    integer :: ndihe
    integer :: nimp

    type(ffbond_t)    , allocatable :: bonds(:)
    type(ffangle_t)   , allocatable :: angles(:)
    type(ffdihedral_t), allocatable :: dihedrals(:)
    type(ffdihedral_t), allocatable :: impropers(:)
  end type fftopology

  type, public :: atom_connect_t
    integer :: nbonds  = 0
    integer :: nangl_e = 0
    integer :: nangl_m = 0
    integer :: ndihe_e = 0
    integer :: ndihe_m = 0
    integer :: nimp    = 0

    integer :: bond_at(6) = 0
    integer, pointer :: angm_at(:,:)
    integer, pointer :: ange_at(:,:)
    integer, pointer :: dihm_at(:,:)
    integer, pointer :: dihe_at(:,:)
    integer, pointer :: imp_at(:,:)
  end type atom_connect_t

  type, public :: atom_t
    !! Generic atom datatype for both QM and MM atoms.
    character(len=4) :: attype = ""
      !! Atom type.
    character(len=4) :: atname = ""
      !! Atom name.
    integer          :: z      = 0
      !! Atomic number.
    real(dp)         :: mass   = 0.0_dp
      !! Atomic mass (in amu).
    real(dp)         :: r(3)   = 0.0_dp
      !! Atom position.
    real(dp)         :: lj_Em  = 0.0_dp
      !! Lennard-Jones energy minimum value.
    real(dp)         :: lj_Rm  = 0.0_dp
      !! Lennard-Jones energy minimum position.
  end type atom_t

  type, public, extends(atom_t) :: qm_atom_t
    !! Datatype specific for QM atoms.
    integer          :: spec = 0
      !! QM atomic species in siesta.
  end type qm_atom_t

  type, public, extends(atom_t) :: mm_atom_t
    !! Datatype specific for MM atoms.
    character(len=4) :: aaname          = ""
      !! Name of the residue that contains this atom.
    integer          :: aanum           = 0
      !! Index of the residue that contains this atom.
    real(dp)         :: pc              = 0.0_dp
      !! Classical partial charge of the atom.
    integer          :: graph_layer     = 0
      !! Graphite layer number (used only in graphite-like structures)
    logical          :: is_blocked      = .false.
      !! Check whether the atom is blocked.
    logical          :: is_qm_neighbour = .false.
      !! Check whether this atom is included within the QMMM cut-off.
  end type mm_atom_t

contains

  function bond_energy( self, dr ) result( bond_e )
    !! Calculates bond energy for the current distance.
    use functions, only : norm_v2
    use precision, only : dp

    implicit none
    class(ffbond_t), intent(in) :: self
      !! Current bond.
    real(dp)       , intent(in) :: dr(3)
      !! Distance vector between atoms 1 and 2.
    real(dp) :: bond_e
    real(dp) :: drij

    drij   = norm_v2( dr ) - self%r_eq
    bond_e = self%k * drij * drij
  end function bond_energy

  subroutine bond_gradient( self, dr, bond_g )
    !! Calculates bond energy gradients for the current distance.
    use functions, only : norm_v2
    use precision, only : dp

    implicit none
    class(ffbond_t), intent(in)  :: self
      !! Current bond.
    real(dp)       , intent(in) :: dr(3)
      !! Distance vector between atoms 1 and 2.
    real(dp)       , intent(out) :: bond_g(3)
    real(dp) :: drij

    drij  = norm_v2( dr )
    bond_g(:) = dr(:) * 2.0_dp * self%k * ( drij - self%r_eq ) / drij
  end subroutine bond_gradient

  function angle_energy( self, dr12, dr32 ) result( angle_e )
    !! Calculates angle energy for the current angle.
    use functions, only : angle_v2
    use mm_units , only : rads
    use precision, only : dp

    implicit none
    class(ffangle_t), intent(in) :: self
      !! Current bond.
    real(dp)        , intent(in)  :: dr12(3)
      !! Distance vector between atoms 1 and 2.
      !! Atom 2 is the middle atom.
    real(dp)        , intent(in)  :: dr32(3)
      !! Distance vector between atoms 3 and 2.
      !! Atom 2 is the middle atom.
    real(dp) :: angle_e
    real(dp) :: ang_v

    ang_v = ( angle_v2( dr12, dr32 ) - self%r_eq ) / rads

    angle_e = self%k * ang_v * ang_v
  end function angle_energy

  subroutine angle_gradient_extreme( self, dr12, dr32, angle_g )
    !! Calculates angle energy gradients for an atom in one of the extremes.
    use functions, only : norm_v2, scalar_v2, angle_v2
    use mm_units , only : rads
    use precision, only : dp

    implicit none
    class(ffangle_t), intent(in)  :: self
      !! Current angle.
    real(dp)        , intent(in)  :: dr12(3)
      !! Distance vector between atoms 1 and 2.
      !! Atom 1 is the current atom, atom 2 is the middle atom.
    real(dp)        , intent(in)  :: dr32(3)
      !! Distance vector between atoms 3 and 2.
      !! Atom 2 is the middle atom.
    real(dp)        , intent(out) :: angle_g(3)

    real(dp) :: r12, r32, scal, ang_v

    ang_v = angle_v2( dr12, dr32 )
    r12  = norm_v2( dr12 )
    r32  = norm_v2( dr32 )
    scal = scalar_v2( dr12, dr32 ) / ( r12 * r12 * r32 * r32 )

    angle_g(:) = scal * r32 * dr12(:) / r12 - dr32(:) / (r12 * r32)
    angle_g(:) = angle_g(:) / ( sqrt( 1.0_dp - ( scal * r12 * r32 ) ** 2 ))
    angle_g(:) = angle_g(:) * 2.0_dp * self%k * ( ang_v - self%r_eq ) / rads
  end subroutine angle_gradient_extreme

  subroutine angle_gradient_middle( self, dr12, dr32, angle_g )
    !! Calculates angle energy gradients for an atom in the vertex.
    use functions, only : norm_v2, scalar_v2, angle_v2
    use mm_units , only : rads
    use precision, only : dp

    implicit none
    class(ffangle_t), intent(in)  :: self
      !! Current angle.
    real(dp)        , intent(in)  :: dr12(3)
      !! Distance vector between atoms 1 and 2.
      !! Atom 1 is the current atom, atom 2 is the middle atom.
    real(dp)        , intent(in)  :: dr32(3)
      !! Distance vector between atoms 3 and 2.
      !! Atom 2 is the middle atom.
    real(dp)        , intent(out) :: angle_g(3)

    real(dp) :: r12, r32, scal, ang_v

    ang_v = angle_v2( dr12, dr32 )

    r12  = norm_v2( dr12 )
    r32  = norm_v2( dr32 )
    scal = scalar_v2( dr12, dr32 ) / ( r12 * r12 * r32 * r32 )

    angle_g(:) =  ( dr12(:) + dr32(:) ) / (r12 * r32) - &
                    scal * ( dr12(:) * r32 / r12 + dr32(:) * r12 / r32 )
    angle_g(:) = angle_g(:) / ( sqrt( 1.0_dp - ( scal * r12 * r32 ) ** 2 ) )
    angle_g(:) = angle_g(:) * 2.0_dp * self%k * ( ang_v - self%r_eq ) / rads
  end subroutine angle_gradient_middle

  function dihedral_energy( self, dr12, dr32, dr43 ) result( dihe_e )
    use functions, only : dihedral_v
    use mm_units , only : rads
    use precision, only : dp

    implicit none
    class(ffdihedral_t), intent(in)  :: self
      !! Current dihedral angle.
    real(dp)           , intent(in)  :: dr12(3)
      !! Distance vector between atoms 1 and 2.
    real(dp)           , intent(in)  :: dr32(3)
      !! Distance vector between atoms 3 and 2.
    real(dp)           , intent(in)  :: dr43(3)
      !! Distance vector between atoms 4 and 3.
    real(dp) :: dihe_e
    real(dp) :: dihe_t

    dihe_t = dihedral_v( dr12, dr32, dr43 ) * abs( self%per ) - self%r_eq
    dihe_e = ( self%k / self%multi ) * ( 1.0_dp + cos( dihe_t / rads ) )
  end function dihedral_energy

  subroutine dihedral_gradient( self, fce, iatom, dr12, dr32, dr43, dr31, dr42 )
    !! Calculates the force over a dihedral angle using the AMBER force field.
    use functions, only : dihedral_v, scalar_v2
    use mm_units , only : rads
    use precision, only : dp
    use qmmm_pbc , only : pbc_displ_vector

    implicit none
    class(ffdihedral_t), intent(in)  :: self
      !! Current dihedral angle.
    real(dp)           , intent(in)  :: dr12(3)
      !! Distance vector between atoms 1 and 2.
    real(dp)           , intent(in)  :: dr32(3)
      !! Distance vector between atoms 3 and 2.
    real(dp)           , intent(in)  :: dr43(3)
      !! Distance vector between atoms 4 and 3.
    real(dp)           , intent(in)  :: dr31(3)
      !! Distance vector between atoms 3 and 1.
    real(dp)           , intent(in)  :: dr42(3)
      !! Distance vector between atoms 4 and 2.
    integer            , intent(in)  :: iatom
      !! Atom of the dihedral we are interested in for the derivative.
      !! Goes from 1 to 4.
    real(dp)           , intent(out) :: fce(3)
      !! Forces for the chosen atom.

    integer  :: iCrd
    real(dp) :: mvec(3), nvec(3), mxn, nnorm, mnorm, &
                dmr(3), dnr(3), dmn, prue, dtot, dscalar, dihe_t

    dihe_t = dihedral_v( dr12, dr32, dr43, mvec, nvec, mxn, mnorm, nnorm )

    fce(:) = 0.0_dp

    prue = mxn / ( mnorm * nnorm )
    prue = ( 1.0_dp - prue * prue )
    if ( abs(prue) < 1e-14_dp ) return

    dtot = sin( ( self%per * dihe_t - self%r_eq ) / rads )
    dtot = ( self%k / self%multi ) * dtot * self%per / sqrt(prue)

    do iCrd = 1, 3
      dmr = 0.0_dp
      dnr = 0.0_dp

      select case ( iatom )
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

end module mm_topology