! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module basis_enthalpy
  implicit none
  private
  public :: write_basis_enthalpy

CONTAINS

  subroutine write_basis_enthalpy( FE, FEharris, nspec, nat_u, species, &
                                   lastorb, norb_u, indorb )
    use files    , only : slabel
    use m_io     , only : io_getout
    use precision, only : dp
    use units    , only : eV

    real(dp), intent(in) :: FE
      !! Electronic (Free)Energy
    real(dp), intent(in) :: FEharris
      !! Electronic (Free)HarrisEnergy
    integer, intent(in) :: nspec
      !! Total number of atomic species.
    integer, intent(in) :: nat_u
      !! Total number of atoms in unit cell.
    integer, intent(in) :: species(nat_u)
      !! Speciex index for each atom.
    integer, intent(in) :: lastorb(0:nat_u)
      !! Global orbital index for the last orbital in an atom.
    integer, intent(in) :: norb_u
      !! Total number of orbitals in unit cell.
    integer, intent(in) :: indorb(norb_u)
      !! For an orbital, the index within a given atom.

    real(dp) :: basis_enthalpy, basis_harris_enthalpy, tot_vol, avg_pressure, &
                orb_enth
    integer  :: iu, fout, is, ia

    real(dp), allocatable :: species_pres(:), species_vol(:)

    allocate( species_pres(nspec) )

    call io_getout( fout )

    write(fout,'(A)') ""
    write(fout,'(A)') "Basis Enthalpy Calculation:"

    call read_species_pressure( species_pres, nspec, fout )
    write(fout,'(A)') ""

    allocate( species_vol(nspec) )
    call orb_volumes( tot_vol, species_vol, nspec, nat_u, species, &
                      lastorb, norb_u, indorb )

    orb_enth = 0.0_dp
    do is = 1, nspec
      orb_enth = orb_enth + species_pres(is) * species_vol(is)
    enddo
    basis_harris_enthalpy = FEHarris + orb_enth
    basis_enthalpy        = FE       + orb_enth

    write(fout,"(4x,a37,f18.6,a3)") &
      "Orbital volume contribution        = ", orb_enth / eV, ' eV'
    write(fout,"(4x,a37,f18.6,a3)") &
      "(Free)E + p_basis*V_orbitals       = ", basis_enthalpy / eV, ' eV'
    write(fout,"(4x,a37,f18.6,a3)") &
      "(Free)Eharris+ p_basis*V_orbitals  = ", basis_harris_enthalpy / eV, ' eV'

    avg_pressure = 0.0_dp
    do ia = 1, nat_u
      is = species(ia)
      avg_pressure = avg_pressure + species_pres(is)
    enddo
    avg_pressure = avg_pressure / nat_u

    write(fout,'(A)') "WARNING: BASIS_ENTHALPY and BASIS_HARRIS_ENTHALPY"//&
                      " files are deprecated. They will be removed in future"//&
                      " releases."
    write(fout,'(A)') "Please use system_label.BASIS_ENTHALPY in your scripts instead."
    call write_old_output_file( basis_enthalpy, FE, avg_pressure, tot_vol, &
                            "BASIS_ENTHALPY" )
    call write_old_output_file( basis_harris_enthalpy, FEHarris, avg_pressure, &
                            tot_vol, "BASIS_HARRIS_ENTHALPY" )
    call write_output_file( basis_enthalpy, FE, basis_harris_enthalpy, &
                            FEHarris, avg_pressure, tot_vol, &
                            trim(slabel)//".BASIS_ENTHALPY" )
  end subroutine write_basis_enthalpy

  subroutine read_species_pressure( species_pres, nspec, fout )
    !! Reads the basis set pressure for each atomic species. If only the
    !! global basis_pressure is present, it applies the same pressure to
    !! all species.
    use chemical , only : species_label, atomic_number
    use precision, only : dp
    use units    , only : GPa
    use fdf      , only : fdf_get, fdf_convfac, block_fdf, parsed_line
    use fdf      , only : fdf_block, fdf_bline, fdf_bmatch, fdf_bclose, &
                          fdf_bintegers, fdf_bvalues, fdf_bnames, leqi

    implicit none
    integer , intent(in)  :: fout
      !! Unit for standard output.
    integer , intent(in)  :: nspec
      !! Number of species.
    real(dp), intent(out) :: species_pres(nspec)
      !! Basis pressure for each atomic species.

    real(dp)                   :: cfactor, basis_pressure
    integer                    :: spec_i, spec_Z, spec_Z2
    character(len=20)          :: specname, upress
    logical                    :: spec_found
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    ! We read the input as GPa and then convert to SIESTA internal units.
    basis_pressure = fdf_get( "BasisPressure", 0.2_dp, "GPa" )
    species_pres(:) = basis_pressure

    if ( fdf_block("BasisPressure.Specs", bfdf) ) then
      do while ( fdf_bline(bfdf,pline) )

        if ( fdf_bmatch(pline,"nvn") ) then
          ! species name - pressure
          specname = trim(fdf_bnames(pline,1))

          upress  = trim(fdf_bnames(pline,2))
          cfactor = fdf_convfac( upress, 'GPa' )

          spec_found = .false.
          do spec_i = 1, nspec
            if ( leqi( specname, species_label(spec_i) ) ) then
              spec_found = .true.
              species_pres(spec_i) = fdf_bvalues(pline, 1) * cfactor
            endif
          enddo
          if ( .not. spec_found ) &
            write(fout,'(A)') "WARNING - BasisPressure.Specs block: Species "//&
                              trim(specname)//" not found."

        else if ( fdf_bmatch(pline,"nivn") ) then
          ! z - number - pressure
          if ( .not. leqi(fdf_bnames(pline,1), 'Z') ) &
            call die("Wrong format in BasisPressure.Specs block.")

          spec_Z  = fdf_bintegers(pline, 1)
          upress  = trim(fdf_bnames(pline,2))
          cfactor = fdf_convfac( upress, 'GPa' )

          spec_found = .false.
          do spec_i = 1, nspec
            if ( spec_Z == atomic_number(spec_i) ) then
              spec_found = .true.
              species_pres(spec_i) = fdf_bvalues(pline, 2) * cfactor
            endif
          enddo
          if ( .not. spec_found ) &
            write(fout,'(A, I0)') "WARNING - BasisPressure.Specs block: No"//&
                                  " species found with atomic number ", spec_Z            

        else if ( fdf_bmatch(pline,"nninivn") ) then  
          ! z from - number - to - number - pressure
          if ( .not. leqi(fdf_bnames(pline,1), 'Z') ) &
            call die("Wrong format in BasisPressure.Specs block.")

          spec_Z  = fdf_bintegers(pline, 1)
          spec_Z2 = fdf_bintegers(pline, 2)
          upress  = trim(fdf_bnames(pline,4))
          cfactor = fdf_convfac( upress, 'GPa' )
          do spec_i = 1, nspec
            if ( (spec_Z2 >= atomic_number(spec_i)) .and. &
                 (spec_Z  <= atomic_number(spec_i)) ) then
              spec_found = .true.
              species_pres(spec_i) = fdf_bvalues(pline, 3) * cfactor
            endif
          enddo

        else
          call die("Wrong format in BasisPressure.Specs block.")
        endif

      enddo
    endif

    do spec_i = 1, nspec
      specname = species_label(spec_i)
      write(fout,'(4x,A,I0,3A,F14.7,A)') &
        "Basis Pressure for species ", spec_i, "(",trim(specname),") : ", &
        species_pres(spec_i), " GPa"
    enddo

    ! Convert to SIESTA internal units (Ry/Bohr**3)
    basis_pressure  = GPa * basis_pressure
    species_pres(:) = GPa * species_pres(:)
  end subroutine read_species_pressure

  subroutine orb_volumes( orb_vol, species_vol, nspec, nat_u, species, &
                          lastorb, norb_u, indorb )
    !! Computes the total volume of the basis orbitals in the unit cell. It
    !! outputs the total orbital volume and the volume per species. Only a
    !! single representative of each "nlz" shell is used (e.g., just one and
    !! not 5 for a 'd' shell).
    use atmfuncs   , only: lofio, rcut
    use precision  , only: dp
    use units      , only: pi

    integer , intent(in)  :: nspec
      !! Total number of atomic species.
    integer , intent(in)  :: nat_u
      !! Total number of atoms in unit cell.
    integer , intent(in)  :: species(nat_u)
      !! Speciex index for each atom.
    integer , intent(in)  :: lastorb(0:nat_u)
      !! Global orbital index for the last orbital in an atom.
    integer , intent(in)  :: norb_u
      !! Total number of orbitals in unit cell.
    integer , intent(in)  :: indorb(norb_u)
      !! For an orbital, the index within a given atom.
    real(dp), intent(out) :: orb_vol
      !! Total orbital volume.
    real(dp), intent(out) :: species_vol(nspec)
      !! Orbital volume per species.

    integer  :: ia, ioa, is, l, io
    real(dp) :: r

    orb_vol     = 0.0_dp
    species_vol = 0.0_dp

    do ia = 1, nat_u
      is = species(ia)

      do io = lastorb(ia-1)+1, lastorb(ia)
        ioa = indorb(io)
        l   = lofio(is,ioa)
        r   = rcut(is,ioa)

        ! The 1 / (2l+1) factor is to only take into account only one
        ! orbital per (n,l) shell.
        species_vol(is) = species_vol(is) &
                        + ( (4.0_dp/3.0_dp)*pi*r**3 ) / ( 2 * l + 1 )
       enddo
    enddo

    do is = 1, nspec
      orb_vol = orb_vol + species_vol(is)
    enddo
  end subroutine orb_volumes

  subroutine write_old_output_file( enthalpy, freeE, avg_pressure, tot_vol, fname )
    !! Prints output file with basis enthalpy info.
    use m_io     , only : io_assign, io_close
    use precision, only : dp
    use units    , only : eV, GPa

    implicit none
    real(dp), intent(in) :: enthalpy
      !! Basis enthalpy (or Harris).
    real(dp), intent(in) :: avg_pressure
      !! Average basis pressure.
    real(dp), intent(in) :: tot_vol
      !! Total orbital volume.
    real(dp), intent(in) :: freeE
      !! Free Energy.
    character(len=*), intent(in) :: fname
      !! Name for outputfile.

    integer :: iu

    call io_assign( iu )
    open( iu , file=trim(fname), form = "formatted" , status = "replace" )
    write( iu, * ) enthalpy / eV
    ! TO-DO: Deprecate these four lines in future releases.
    write( iu, * ) "The above number is the electronic (free)energy:", freeE/eV
    write( iu, * ) "Plus the pressure : ", avg_pressure,  &
                                    " ( ", avg_pressure/GPa, " GPa)"
    write( iu, * ) "     times the orbital volume (in Bohr**3): ", tot_vol
    write( iu, * ) "WARNING: This file is deprecated. It will be removed in"//&
                   " future releases."
    call io_close( iu )

  end subroutine write_old_output_file

  subroutine write_output_file( enthalpy, freeE, enthalpy_h, freeE_h, &
                                avg_pressure, tot_vol, fname )
    !! Prints output file with basis enthalpy info.
    use m_io     , only : io_assign, io_close
    use precision, only : dp
    use units    , only : eV, GPa, Ang

    implicit none
    real(dp), intent(in) :: enthalpy
      !! Basis enthalpy.
    real(dp), intent(in) :: enthalpy_h
      !! Basis Harris enthalpy.
    real(dp), intent(in) :: freeE
      !! Free Energy.
    real(dp), intent(in) :: freeE_h
      !! Harris Energy.
    real(dp), intent(in) :: avg_pressure
      !! Average basis pressure.
    real(dp), intent(in) :: tot_vol
      !! Total orbital volume.
    character(len=*), intent(in) :: fname
      !! Name for outputfile.

    integer :: iu

    call io_assign( iu )
    open( iu , file=trim(fname), form = "formatted" , status = "replace" )

    write( iu, * ) "Basis enthalpy [eV] : ", enthalpy / eV
    write( iu, * ) "Basis Harris enthalpy [eV] : ", enthalpy_h / eV
    write( iu, * ) "Free energy [eV] : ", freeE / eV
    write( iu, * ) "Harris energy [eV] : ", freeE_h / eV
    write( iu, * ) "Orbital volume [Ang**3] : ", tot_vol / (Ang * Ang * Ang)
    write( iu, * ) "Average basis pressure [eV/Ang**3] : ", &
                   avg_pressure * Ang * Ang * Ang / eV
    write( iu, * ) "Average basis pressure [GPa] : ", avg_pressure / GPa

    call io_close( iu )

  end subroutine write_output_file
end module basis_enthalpy
