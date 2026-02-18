!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module grimme_dispersion_m
  !! This module acts as an interface for Grimme's DFT-D3 library, which
  !! itself comprises two main sections: the D3 model itself (dftd3 module)
  !! and the MCTC library, which handles geometry and coordinates setup.
  !!
  !! We only need to call initialization once, and then energies and
  !! forces are calculated after SCF convergence.
  !!
  !! Even though the interfaced library holds different datatypes
  !! for parameters, we also keep our own dftd3_data_t so as to ease
  !! I/O from the fdf and avoid some problems in case the interface is
  !! changed in the future. This structure is initialized during module
  !! initialization in dftd3_initialize.
  !!
  !! Initialization can be called at any point before DFT-D3 forces
  !! contribution calculation, but it only needs to be called once. It
  !! also requires XC data to be already initialized in the code.
  !!
  !! dftd3_energy_forces, the subroutine which calculates forces, energy,
  !! and stress contributions, can also be called at any point after atom
  !! data is set up, since it only requires coordinates and atom identities.
  !! This is a post-SCF correction, so it does not enter the SCF loop.
  !!
  !! dftd3_get_periodic attemps to read the DFTD3.Periodic input variable to
  !! get the periodicity for D3. Otherwise, it relies on the periodicity
  !! guessed by SIESTA.
  !!
  !!
  !! THEORY BACKGROUND
  !! -----------------
  !!
  !! For ease of reading, we will present here three of the equations
  !! present in Grimme's work (DOI: 10.1063/1.3382344). Equation numbers
  !! correspond to those in the paper.
  !!
  !! (2) E(D3) = E(2-body) + E(3-body)
  !!
  !! (3) E(2-body) = sum_{A,B} ( s6 * C6AB / (rAB)^6 ) * f6(rAB)
  !!               + sum_{A,B} ( s8 * C8AB / (rAB)^8 ) * f8(rAB)
  !!
  !! Where A,B are atoms, sn and CnAB are parameters, rAB the distance
  !! between atoms, and fn a damping function:
  !!
  !! (4) fn(rAB) = 1 / ( 1 + 6 * ( rAB / (Srn * R0AB))^(-alp_n) )
  !!
  !! Where Srn, R0AB and alp_n are more parameters. Only CnAB and R0AB
  !! depend on the atoms, while Srn, Sn, and alp_n depend only on the DFT
  !! functional chosen.
  !!
  !! The 3-body term is a bit more complicated but at the same time
  !! requires less parametrization:
  !!
  !! (14) E(3-body) = Sum_{A,B,C} f3(rABC) * E(ABC)
  !!
  !! (11) E(ABC) = C9ABC * (3* cos(Ta) * cos(Tb) * cos(Tc) + 1) /
  !!                       (rAB * rBC * rCA)^3
  !!
  !! (13) C9ABC = - sqrt( C6AB * C6BC * C6AC )
  !!
  !! Where the value of C9ABC is, in truth, an approximation. Ta, Tb
  !! and Tc are the angles of the triangle formed by A, B and C, and
  !! the damping function f3 uses a value of 16 for alp and 4/3 for Sr.
#ifdef SIESTA__DFTD3
  use dftd3    , only : d3_param, rational_damping_param, zero_damping_param

  use precision, only : dp
  implicit none

  public :: dftd3_initialize
  public :: dftd3_energy_forces
  public :: dftd3_finalize

  private

  type dftd3_data_t
    !! Datatype to contain all of the SIESTA
    !! input data for D3. s8 and rs8 are kept with those names
    !! (instead of s8 and rs8) for consistency with the D3 library.
    real(dp) :: s6       = 1.0_dp
      !! S6 coefficient that pre-multiplies all
      !! of the two-body C6 terms. (eq 3)
    real(dp) :: rs6      = 1.0_dp
      !! S_{r6} prefactor in the damping factor quotient. (eq 4)
    real(dp) :: s8       = 1.0_dp
      !! S8 coefficient that pre-multiplies all
      !! of the two-body C8 terms. (eq 3)
    real(dp) :: rs8      = 1.0_dp
      !! S_{r8} prefactor in the damping factor quotient.
      !! Usually set to 1.0. (eq 4)
    real(dp) :: s9       = 1.0_dp
      !! Weight for 3-body interaction terms. When zero, no 3-body
      !! terms are calculated.
    real(dp) :: a1       = 0.4_dp
      !! First parameter in Becke-Johnson damping
    real(dp) :: a2       = 5.0_dp
      !! First parameter in Becke-Johnson damping
    real(dp) :: alp      = 14.0_dp
      !! Exponent for the C6 damping factor. The C8 factor is alp + 2.0.
    logical  :: BJ_damp  = .true.
      !! Whether we use Becke-Johnson damping or zero damping variants.
    real(dp) :: cutoff_cn = 10.0_dp
      !! Cut-off for coordination.
    real(dp) :: cutoff_2b = 60.0_dp
      !! Cut-off for 2-body interactions.
    real(dp) :: cutoff_3b = 40.0_dp
      !! Cut-off for 3-body interactions.
    logical  :: periodic(3) = .true.
      !! Whether the cell is periodic in the direction of each of the lattice
      !! vectors.
    logical :: override_periodic = .false.
      !! Whether we set the periodicity directions from the FDF options.
  end type dftd3_data_t

  type dftd3_geom_t
    !! This datatype is used to filter out all atoms that are not with
    !! an atomic number between 1 and 94; mainly ghost atoms and synthetics.
    integer :: nat    = 0
      !! Number of atoms.
    real(dp), pointer :: crd(:,:) => null()
      !! Atomic coordinates.
    real(dp), pointer :: frc(:,:) => null()
      !! Atomic forces.
    integer, pointer  :: z(:) => null()
      !! Atomic number.
    character(len=4), pointer :: sym(:) => null()
      !! Atomic symbols.
    integer, pointer  :: sidx(:) => null()
      !! Original atom index within SIESTA.
  contains
    procedure, pass :: setup      => dftd3_geom_setup
    procedure, pass :: update     => dftd3_geom_update
    procedure, pass :: get_forces => dftd3_geom_get_forces
    procedure, pass :: reset      => dftd3_geom_reset
  end type dftd3_geom_t

  type(d3_param)               :: grimme_d3_dat
    !! A datatype within the external library that contains all relevant D3
    !! parameters.
  type(rational_damping_param) :: grimme_d3_dparam
    !! Contains damping information in the case of rational BJ damping.
  type(zero_damping_param)     :: grimme_d3_zparam
    !! Contains damping information in the case of zero damping.
  type(dftd3_data_t)           :: siesta_d3_dat
    !! All of the SIESTA input data that is passed to the D3 library.
  type(dftd3_geom_t)           :: siesta_d3_atoms
    !! Contains pre-screened atoms that have actual D3 components.

  logical :: first_call = .true.
    !! Used mainly to print warning messages.

contains

  subroutine dftd3_initialize( )
    !! Reads input data and initializes the D3 library.
    use fdf      , only : fdf_get
    use sys      , only : die

    use dftd3    , only : get_rational_damping, get_zero_damping, &
                          new_rational_damping, new_zero_damping
    use mctc_env , only : error_type
    use gridXC   , only : getXC=>gridxc_getXC
    use sys      , only : message
    use m_cite   , only : add_citation

    implicit none
    integer           :: n_xc
    type(error_type) , allocatable :: d3_error
    character(len=20), allocatable :: fun_type(:), fun_author(:)

    siesta_d3_dat%BJ_damp = fdf_get( 'DFTD3.BJdamping', siesta_d3_dat%BJ_damp )

    call getXC( n = n_xc )
    if ( n_xc > 1 ) &
      call die( 'DFT-D3 functional initialization is not supported'&
              &' for mixed/cocktail functionals.'  )

    allocate( fun_type(n_xc), fun_author(n_xc) )
    call getXC( n = n_xc, func = fun_type, auth = fun_author )

    if ( fdf_get("DFTD3.UseXCDefaults", .true.) ) then
      if ( fun_type(1) == 'LDA' ) &
        call die( 'DFT-D3 functional initialization is not supported'&
                &' for LDA.' )
      if ( fun_type(1) == 'VDW' ) &
        call die( 'DFT-D3 functional initialization is not supported'&
                &' for VDW functionals.' )

      select case ( fun_author(1) )
      case ( 'PBE', 'pbe', 'PBESOL', 'pbesol', 'PBEsol', 'REVPBE', &
             'revpbe', 'revPBE', 'BLYP', 'LYP', 'blyp', 'lyp',    &
             'rpbe', 'RPBE', 'hse6', 'HSE6', 'pbe0', 'PBE0' )
        if ( (fun_author(1) == 'hse6') .or. (fun_author(1) == 'HSE6') ) &
          fun_author(1) = 'hse06'
          call message( 'INFO',  'DFT-D3: loading default parameters for' //&
                        ' functional ' // trim(fun_author(1)) // '.' )
      case default
        if ( allocated(fun_type)   ) deallocate( fun_type   )
        if ( allocated(fun_author) ) deallocate( fun_author )
        call die( 'This functional is not available with DFT-D3 '&
                &'functional initialization. Set DFTD3.UseXCDefaults '&
                &'to false and add custom parameters for s6, rs6 and s8.')
      end select

      siesta_d3_dat%s9 = fdf_get( 'DFTD3.s9', siesta_d3_dat%s9 )
      if ( siesta_d3_dat%BJ_damp ) then
        call add_citation("10.1002/jcc.21759")
        call get_rational_damping( grimme_d3_dat, fun_author(1), d3_error, &
                                   s9 = siesta_d3_dat%s9 )
      else
        call add_citation("10.1063/1.3382344")
        call get_zero_damping( grimme_d3_dat, fun_author(1), d3_error, &
                               s9 = siesta_d3_dat%s9 )
      endif
      if ( allocated( d3_error ) ) then
         call message('WARNING',"d3_error in dftd3_initialize")
         return
      endif
    endif


    if ( fun_type(1) == 'VDW' ) &
      call message( 'WARNING', 'D3 corrections should ideally be used in '//&
                    'combination with vdW functionals.' )

    if ( allocated(fun_type)   ) deallocate( fun_type   )
    if ( allocated(fun_author) ) deallocate( fun_author )

    ! We overwrite parameters with custom values if present in the fdf.
    siesta_d3_dat%s6  = fdf_get( 'DFTD3.s6'   , siesta_d3_dat%s6  )
    siesta_d3_dat%rs6 = fdf_get( 'DFTD3.rs6'  , siesta_d3_dat%rs6 )
    siesta_d3_dat%s8  = fdf_get( 'DFTD3.s8'   , siesta_d3_dat%s8  )
    siesta_d3_dat%rs8 = fdf_get( 'DFTD3.rs8'  , siesta_d3_dat%rs8 )
    siesta_d3_dat%a1  = fdf_get( 'DFTD3.a1'   , siesta_d3_dat%a1  )
    siesta_d3_dat%a2  = fdf_get( 'DFTD3.a2'   , siesta_d3_dat%a2  )
    siesta_d3_dat%alp = fdf_get( 'DFTD3.alpha', siesta_d3_dat%alp )
    siesta_d3_dat%s9  = fdf_get( 'DFTD3.s9'   , siesta_d3_dat%s9  )

    grimme_d3_dat%s6  = siesta_d3_dat%s6
    grimme_d3_dat%rs6 = siesta_d3_dat%rs6
    grimme_d3_dat%s8  = siesta_d3_dat%s8
    grimme_d3_dat%s9  = siesta_d3_dat%s9
    grimme_d3_dat%rs8 = siesta_d3_dat%rs8
    grimme_d3_dat%a1  = siesta_d3_dat%a1
    grimme_d3_dat%a2  = siesta_d3_dat%a2
    grimme_d3_dat%alp = siesta_d3_dat%alp
    if ( siesta_d3_dat%BJ_damp ) then
      call new_rational_damping( grimme_d3_dparam, grimme_d3_dat )
    else
      call new_zero_damping( grimme_d3_zparam, grimme_d3_dat )
    endif

    ! These should be in atomic units.
    siesta_d3_dat%cutoff_cn = fdf_get( 'DFTD3.CoordinationCutoff', &
                                       siesta_d3_dat%cutoff_cn , 'Bohr' )
    siesta_d3_dat%cutoff_2b = fdf_get( 'DFTD3.2BodyCutoff', &
                                       siesta_d3_dat%cutoff_2b , 'Bohr' )
    siesta_d3_dat%cutoff_3b = fdf_get( 'DFTD3.3BodyCutoff', &
                                       siesta_d3_dat%cutoff_3b , 'Bohr' )

  end subroutine dftd3_initialize

  subroutine dftd3_energy_forces( natoms, atm_crd, iza, lattice_vecs, &
                                  periodic_i, Edisp, Frc, stress )
    !! Calculates D3 corrections to energies, forces and stress.
    use alloc         , only : re_alloc, de_alloc
    use dftd3         , only : new_d3_model, get_dispersion, d3_model, &
                               realspace_cutoff
    use mctc_io       , only : structure_type, new_structure

    implicit none
    integer , intent(in)    :: natoms
      !! Number of atoms in unit cell.
    real(dp), intent(in)    :: atm_crd(3, natoms)
      !! Atom coordinates.
    integer , intent(in)    :: iza(natoms)
      !! Atomic numbers.
    real(dp), intent(in)    :: lattice_vecs(3,3)
      !! Cell vectors.
    logical , intent(in)    :: periodic_i(3)
      !! Cell vectors that are periodic, as detected by SIESTA.
    real(dp), intent(out)   :: Edisp
      !! Dispersion correction to energies.
    real(dp), intent(inout) :: Frc(3, natoms)
      !! Dispersion correction to forces.
    real(dp), intent(inout) :: stress(3,3)
      !! Dispersion correction to stress.

    external                      :: timer
    real(dp), external            :: volcel

    integer                       :: ii, jj
    real(dp)                      :: cell_vol
    real(dp), pointer, contiguous :: Gdisp(:,:), Gstress(:,:)

    type(realspace_cutoff)        :: cuts
      ! Stores cut-off data.
    type(structure_type)          :: molec
      ! Stores molecule information and PBC.
    type(d3_model)                :: grimme_d3_model
      ! Stores all coefficients and pre-calculations needed for the D3 model.

    call timer( 'dftd3', 1 )

    nullify( Gdisp, Gstress )
    call re_alloc( Gdisp  , 1, 3, 1, natoms, 'Gdisp' , &
                   'dftd3_energy_forces' )
    call re_alloc( Gstress, 1, 3, 1,      3, 'Gstress', &
                   'dftd3_energy_forces' )

    call siesta_d3_atoms%setup( natoms, atm_crd, iza )
    call siesta_d3_atoms%update( natoms, atm_crd )

    if ( .not. siesta_d3_dat%override_periodic ) &
      siesta_d3_dat%periodic = periodic_i
    call dftd3_get_periodic( siesta_d3_dat%periodic, &
                             siesta_d3_dat%override_periodic )
    call new_structure( molec, siesta_d3_atoms%z, siesta_d3_atoms%sym, &
                        siesta_d3_atoms%crd, lattice = lattice_vecs,   &
                        periodic = siesta_d3_dat%periodic )

    call new_d3_model( grimme_d3_model, molec )

    cuts%cn    = siesta_d3_dat%cutoff_cn
    cuts%disp2 = siesta_d3_dat%cutoff_2b
    cuts%disp3 = siesta_d3_dat%cutoff_3b

    if ( siesta_d3_dat%BJ_damp ) then
      call get_dispersion( molec, grimme_d3_model, grimme_d3_dparam,        &
                           cuts, Edisp, gradient = siesta_d3_atoms%frc,     &
                           sigma = Gstress )
    else
      call get_dispersion( molec, grimme_d3_model, grimme_d3_zparam,        &
                           cuts, Edisp, gradient = siesta_d3_atoms%frc,     &
                           sigma = Gstress )
    endif

    ! The *2 factor is due to the conversion from Hartree to Rydberg
    Edisp = Edisp * 2.0_dp

    cell_vol = volcel( lattice_vecs )
    do jj = 1, 3
    do ii = 1, 3
      stress(ii,jj) = stress(ii,jj) - Gstress(ii,jj) * 2.0_dp / cell_vol
    enddo
    enddo

    call siesta_d3_atoms%get_forces( natoms, Gdisp )
    do jj = 1, natoms
    do ii = 1, 3
      Frc(ii,jj) = Frc(ii,jj) - Gdisp(ii,jj) * 2.0_dp
    enddo
    enddo

    call de_alloc( Gdisp  , 'Gdisp'  , 'dftd3_energy_forces' )
    call de_alloc( Gstress, 'Gstress', 'dftd3_energy_forces' )

    call timer( 'dftd3', 2 )
  end subroutine dftd3_energy_forces

  subroutine dftd3_get_periodic( is_periodic, fdfoverride )
    !! Sets system periodicity from the one guessed by SIESTA. This can also be
    !! manually set up via the DFTD3.Periodic input variable.
    use fdf      , only : fdf_list, fdf_islist
    use m_io     , only : io_getout
    use sys      , only : message, die

    implicit none
    logical, intent(inout) :: is_periodic(3)
      !! Whether a given lattice vector is a periodic direction.
    logical, intent(inout) :: fdfoverride
      !! Whether we get the periodicity from the fdf.
    integer          :: nvec, ivec, per_list(3), stdout

    if (.not. first_call ) return
    first_call = .false.

    call io_getout( stdout )

    if ( fdf_islist('DFTD3.Periodic') ) then
      nvec = -1
      call fdf_list( 'DFTD3.Periodic', nvec, per_list )

      if ( nvec > 3 ) then
        call message( 'WARNING', 'DFTD3.Periodic has more than 3 elements.'//&
                        ' Only the first 3 elements will be read.')
        nvec = 3
      endif
      fdfoverride = .true.

      if ( nvec > 0 ) then
        is_periodic(1:3) = .false.
        call fdf_list( 'DFTD3.Periodic', nvec, per_list )

        do ivec = 1, nvec
          if ( (per_list(ivec) > 3) .or. (per_list(ivec) < 0) ) &
            call die( 'DFTD3.Periodic indices must be 1, 2 or 3.' )

          is_periodic( per_list(ivec) ) = .true.
        enddo

      else
        call die( 'DFTD3.Periodic declared but no elements found in list.' )

      endif

    else
      call message( 'INFO', 'Using SIESTA-detected periodicity for D3.'//&
                            ' Use DFTD3.Periodic to override.' )
      if (is_periodic(1)) &
        write(stdout,'(A)') 'Periodicity found for lattice vector 1.'
      if (is_periodic(2)) &
        write(stdout,'(A)') 'Periodicity found for lattice vector 2.'
      if (is_periodic(3)) &
        write(stdout,'(A)') 'Periodicity found for lattice vector 3.'
    endif
  end subroutine dftd3_get_periodic

  subroutine dftd3_finalize( )
    !! Deallocates internal datastructures.

    call siesta_d3_atoms%reset()
  end subroutine dftd3_finalize

  subroutine dftd3_geom_setup( self, natoms, atom_crd, atomZ )
    !! Generates a subset of atoms with only those that have an atomic number
    !! within 1 and 94, since the rest are not included in the D3 model.
    use alloc         , only : re_alloc
    use precision     , only : dp
    use periodic_table, only : symbol

    implicit none
    class(dftd3_geom_t), intent(inout) :: self
      !! Data structure that will contain the atoms filtered.
    integer           , intent(in)    :: natoms
      !! Total number of atoms.
    real(dp)          , intent(in)    :: atom_crd(3,natoms)
      !! Atomic coordinates.
    integer           , intent(in)    :: atomZ(natoms)
      !! Atomic numbers.

    integer :: iat, filtat

    ! If nat > 0, then the structure has already been set up.
    if ( self%nat > 0 ) return

    ! First loop to get atoms and then filter them.
    filtat = 0
    do iat = 1, natoms
      if ( (atomZ(iat) < 95) .and. (atomZ(iat) > 0)) filtat = filtat +1
    enddo
    self%nat = filtat

    nullify( self%crd )
    call re_alloc( self%crd, 1, 3, 1, self%nat, 'crd', 'dftd3' )

    nullify( self%z )
    call re_alloc( self%z, 1, self%nat, 'z', 'dftd3' )

    nullify( self%sidx )
    call re_alloc( self%sidx, 1, self%nat, 'sidx', 'dftd3' )

    nullify( self%sym )
    call re_alloc( self%sym, 1, self%nat, 'sym', 'dftd3' )

    nullify( self%frc )
    call re_alloc( self%frc, 1, 3, 1, self%nat, 'frc', 'dftd3' )

    filtat = 0
    self%frc(:,:) = 0.0_dp
    do iat = 1, natoms
      if ( (atomZ(iat) > 94) .or. (atomZ(iat) < 1)) cycle

      filtat = filtat +1
      self%sidx(filtat) = iat
      self%z(filtat)    = atomZ(iat)
      self%sym(filtat)  = adjustl( symbol( atomZ(iat) ) )

      self%crd(:,filtat) = atom_crd(:,iat)
    enddo
  end subroutine dftd3_geom_setup

  subroutine dftd3_geom_update( self, natoms, atom_crd )
    !! Updates the coordinates in the structure data.
    use precision     , only : dp

    implicit none
    class(dftd3_geom_t), intent(inout) :: self
      !! Data structure containing the filtered atoms
    integer           , intent(in)    :: natoms
      !! Total number of atoms.
    real(dp)          , intent(in)    :: atom_crd(3,natoms)
      !! Atomic coordinates.

    integer :: iat

    do iat = 1, self%nat
      self%crd(:,iat) = atom_crd(:,self%sidx(iat))
    enddo
  end subroutine dftd3_geom_update

  subroutine dftd3_geom_get_forces( self, natoms, atom_frc )
    !! Gets the forces for the filtered atoms and expands
    !! it into the full range of atoms.
    use precision     , only : dp

    implicit none
    class(dftd3_geom_t), intent(inout) :: self
      !! Data structure containing the filtered atoms.
    integer           , intent(in)    :: natoms
      !! Total number of atoms.
    real(dp)          , intent(out)   :: atom_frc(3,natoms)
      !! Atomic coordinates.

    integer :: iat

    atom_frc(:,:) = 0.0_dp
    do iat = 1, self%nat
      atom_frc(:,self%sidx(iat)) = self%frc(:,iat)
    enddo
  end subroutine dftd3_geom_get_forces

  subroutine dftd3_geom_reset( self )
    !! Resets the datastructure.
    use alloc, only : de_alloc

    implicit none
    class(dftd3_geom_t), intent(inout) :: self
      !! Data structure to be destroyed.

    self%nat    = 0
    call de_alloc( self%crd, 'crd', 'dftd3' )
    nullify( self%crd )
    call de_alloc( self%z, 'z', 'dftd3' )
    nullify( self%z )
    call de_alloc( self%sidx, 'sidx', 'dftd3' )
    nullify( self%sidx )
    call de_alloc( self%sym, 'sym', 'dftd3' )
    nullify( self%sym )
    call de_alloc( self%frc, 'frc', 'dftd3' )
    nullify( self%frc )

  end subroutine dftd3_geom_reset
#endif
end module grimme_dispersion_m
