module electric_field
  !! This module deals with the interaction between classical atomic
  !! partial charges and an external constant electric field.
  use precision  , only : dp

  implicit none

  public :: elecfield
  private

  logical  :: first_call  = .true.
    !! Whether it's the first time calling this routines.
  logical  :: field_found = .false.
    !! If we found an Electric field input.
  logical  :: TS_pot_found = .false.
    !! If we are doing a transiesta run, we can get an electric field.
  real(dp) :: efield(3)
    !! The magnitude of the electric field, if present.

contains
  subroutine elecfield( nac, rclas, ucell, mm_atoms, fdummy, Eelec )
    !! This subroutine sets up the external electric field and
    !! then calculates the interaction between the clasiscal
    !! atoms and said electric field, if required.
    use precision , only : dp
    use mm_topology, only : mm_atom_t

    implicit none
    integer , intent(in)    :: nac
      !! Total number of atoms.
    real(dp), intent(in)    :: rclas(3,nac)
      !! Atomic positions.
    real(dp), intent(in)    :: ucell(3,3)
      !! Size of the simulation box.
    type(mm_atom_t), intent(in) :: mm_atoms(nac)
      !! Atomic partial (classical) charges.
    real(dp), intent(inout) :: fdummy(3,nac)
      !! Atomic forces.
    real(dp), intent(inout) :: Eelec
      !! Energy contribution of the electric field.

    Eelec = 0.0_dp
    ! Read the electric field from the fdf file
    if ( first_call ) then
      efield(:) = 0.0_dp

      call read_transiesta_field( efield, TS_pot_found, ucell )
      call read_custom_efield( efield, field_found, TS_pot_found )

      first_call  = .false.
    endif

    if ( .not. (field_found .or. TS_pot_found) ) return

    ! Calculation of the force on the MM atoms due to the field
    call electricfield( nac, rclas, mm_atoms, efield, fdummy, Eelec )

  end subroutine elecfield

  subroutine read_custom_efield( elefield, custom_found, ts_found )
    !! Reads an arbitrary external field.
    use fdf      , only : block_fdf, parsed_line
    use fdf      , only : fdf_block, fdf_bline, fdf_bintegers, fdf_bmatch, &
                          fdf_bvalues, fdf_bnames, fdf_bclose, fdf_convfac, fdf_bnnames
    use mm_units , only : kcal_mol_eV
    use precision, only : dp
    use sys      , only : die

    implicit none
    real(dp), intent(inout) :: elefield(3)
      !! The electric field.
    logical , intent(out)   :: custom_found
      !! Whether we found the data in the fdf.
    logical , intent(in)    :: ts_found
      !! Whether we already found a transiesta potential.

    character(len=20) :: units
    real(dp)          :: cfactor
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    if ( fdf_block( 'ExternalElectricField', bfdf ) ) then
      if ( ts_found ) &
        write(*,*) 'WARNING: Using a custom electric field will override '//&
                   'the field coming from a transiesta potential.'
      elefield(:) = 0.0_dp

      do while ( fdf_bline(bfdf,pline) )
        if ( .not. fdf_bmatch(pline,"vvvn") ) &
          call die("Wrong format in ElectricField block")
        units   = trim(fdf_bnames(pline,1))
        cfactor = fdf_convfac(units,'V/Ang')

        elefield(1) = fdf_bvalues(pline,1) * cfactor
        elefield(2) = fdf_bvalues(pline,2) * cfactor
        elefield(3) = fdf_bvalues(pline,3) * cfactor
      end do

      call fdf_bclose( bfdf )

      write(*,*)
      write(*,'(a)') 'qmmm: Running with an External Electric Field of:'
      write(*,'(3F8.4,3x,A)') elefield(1:3),' V/Ang'
      write(*,'(a)') 'qmmm: WARNING: System must be neutral.'

      elefield     = elefield * kcal_mol_eV
      custom_found = .true.
    endif

  end subroutine read_custom_efield

  subroutine read_transiesta_field( elefield, found_ts, ucell )
    !! Reads an arbitrary external field.
    use fdf       , only : fdf_get, leqi
    use mm_units  , only : kcal_mol_eV
    use precision , only : dp
    use sys       , only : die
    use units     , only : eV, Ang

    implicit none
    real(dp), intent(in)    :: ucell(3,3)
      !! Size of the simulation box.
    real(dp), intent(inout) :: elefield(3)
      !! The electric field.
    logical , intent(out)   :: found_ts
      !! Whether we found the data in the fdf.

    real(dp)          :: potential
    integer           :: divide_vector, icrd
    character(len=80) :: solv_method

    potential = 0.0_dp
    found_ts  = .false.

    solv_method = fdf_get( 'SolutionMethod', 'diagon' )

    if ( .not. leqi(solv_method, 'transiesta') ) return

    potential     = fdf_get('TS.Voltage', potential, 'eV')
    divide_vector = fdf_get('QMMM.ElectrodeCellVector', 3)
    found_ts      = .true.

    do icrd = 1, 3
      if ( abs(ucell(icrd, divide_vector)) < 1.d-6 ) cycle

      elefield(icrd) = - Ang * potential / ucell(icrd, divide_vector)
    enddo

    write(*,'(A)') "TS electric field over MM charges (V/Ang):"
    write(*,'(5x,F14.6,F14.6,F14.6)') elefield(1:3)

    elefield(:) = elefield(:) * kcal_mol_eV
  end subroutine read_transiesta_field

  subroutine electricfield( nac, r, mm_atoms, elefield, force, Eelec )
    !! Calculates electric field contributions to energies and
    !! forces.
    use precision, only: dp
    use units    , only: Ang
    use mm_topology, only : mm_atom_t

    implicit none
    integer , intent(in)    :: nac
      !! Total number of atoms.
    real(dp), intent(in)    :: r(3,nac)
      !! Atomic positions.
    type(mm_atom_t), intent(in) :: mm_atoms(nac)
      !! Atomic partial (classical) charges.
    real(dp), intent(in)    :: elefield(3)
      !! Electric field magnitude.
    real(dp), intent(inout) :: force(3,nac)
      !! Atomic forces.
    real(dp), intent(inout) :: Eelec
      !! Energy contribution of the electric field.

    integer :: iatom

    ! Calculation of the energies and forces due to the field to be
    ! added on fdummy.
    do iatom = 1, nac
      Eelec = Eelec - mm_atoms(iatom)%pc * &
              ( r(1,iatom) * elefield(1) + r(2,iatom) * elefield(2) &
              + r(3,iatom) * elefield(3) ) / Ang

      force(1:3,iatom) = force(1:3,iatom) + &
                         mm_atoms(iatom)%pc * elefield(1:3)
    enddo
  end subroutine electricfield

end module electric_field
