module assign_qmmm
  !! This module contains a single routine that assings masses and species
  !! to all atoms.

  implicit none
  public :: assign_atoms
contains

  subroutine assign_atoms( na_qm, na_mm, mm_atoms, qm_atoms )
    !! Subroutine that assigns masses and atomic symbols (sym).
    use fdf           , only : block_fdf, parsed_line
    use fdf           , only : fdf_block, fdf_bline, fdf_bnnames, fdf_bnvalues,&
                               fdf_bnames, fdf_bvalues, fdf_bclose
    use mm_topology   , only : mm_atom_t, qm_atom_t
    use periodic_table, only : ATMASS, SYMBOL
    use precision     , only : dp
    use sys           , only : die

    implicit none
    integer         , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    type(mm_atom_t) , intent(inout) :: mm_atoms(na_mm)
      !! MM atom names.
    type(qm_atom_t) , intent(inout) :: qm_atoms(na_qm)
      !! QM atomic numbers.

    integer          :: iatom, katom
    real(dp)         :: mass
    character(len=1) :: atn1
    character(len=2) :: atn2, name
    character(len=4) :: atn4

    character(len=2), allocatable :: sym(:)
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    allocate( sym( na_qm+na_mm ) )
    ! Atoms in the QM region.
    do iatom = 1, na_qm
      qm_atoms(iatom)%mass = ATMASS( qm_atoms(iatom)%z )
      sym(iatom)           = SYMBOL( qm_atoms(iatom)%z )
      if ( abs(qm_atoms(iatom)%mass) < 1.0d-10 ) &
        call die( "qmmm assign: There are QM atoms without mass." )
    enddo

    ! Atoms in the MM region.
    katom        = na_qm +1
    do iatom = 1, na_mm
      atn4 = mm_atoms(iatom)%atname
      atn1 = atn4(1:1)
      atn2 = atn4(1:2)

      ! We need to first check atn1 in cases there are overlaps with two-letter
      ! atom names.
      if ( atn1 == 'H'  ) mm_atoms(iatom)%z =  1
      if ( atn1 == 'B'  ) mm_atoms(iatom)%z =  5
      if ( atn1 == 'C'  ) mm_atoms(iatom)%z =  6
      if ( atn1 == 'N'  ) mm_atoms(iatom)%z =  7
      if ( atn1 == 'O'  ) mm_atoms(iatom)%z =  8
      if ( atn1 == 'F'  ) mm_atoms(iatom)%z =  9
      if ( atn1 == 'P'  ) mm_atoms(iatom)%z = 15
      if ( atn1 == 'S'  ) mm_atoms(iatom)%z = 16
      if ( atn1 == 'K'  ) mm_atoms(iatom)%z = 19
      if ( atn1 == 'V'  ) mm_atoms(iatom)%z = 23
      if ( atn1 == 'Y'  ) mm_atoms(iatom)%z = 39
      if ( atn1 == 'I'  ) mm_atoms(iatom)%z = 53
      if ( atn1 == 'W'  ) mm_atoms(iatom)%z = 74
      if ( atn1 == 'U'  ) mm_atoms(iatom)%z = 92

      if ( atn2 == 'He' ) mm_atoms(iatom)%z =  2
      if ( atn2 == 'Li' ) mm_atoms(iatom)%z =  3
      if ( atn2 == 'Be' ) mm_atoms(iatom)%z =  4
      if ( atn2 == 'Ne' ) mm_atoms(iatom)%z = 10
      if ( atn2 == 'Na' ) mm_atoms(iatom)%z = 11
      if ( atn2 == 'Mg' ) mm_atoms(iatom)%z = 12
      if ( atn2 == 'Al' ) mm_atoms(iatom)%z = 13
      if ( atn2 == 'Si' ) mm_atoms(iatom)%z = 14
      if ( atn2 == 'Cl' ) mm_atoms(iatom)%z = 17
      if ( atn2 == 'Ar' ) mm_atoms(iatom)%z = 18
      if ( atn2 == 'Ca' ) mm_atoms(iatom)%z = 20
      if ( atn2 == 'Sc' ) mm_atoms(iatom)%z = 21
      if ( atn2 == 'Ti' ) mm_atoms(iatom)%z = 22
      if ( atn2 == 'Cr' ) mm_atoms(iatom)%z = 24
      if ( atn2 == 'Mn' ) mm_atoms(iatom)%z = 25
      if ( atn2 == 'Fe' ) mm_atoms(iatom)%z = 26
      if ( atn2 == 'Co' ) mm_atoms(iatom)%z = 27
      if ( atn2 == 'Ni' ) mm_atoms(iatom)%z = 28
      if ( atn2 == 'Cu' ) mm_atoms(iatom)%z = 29
      if ( atn2 == 'Zn' ) mm_atoms(iatom)%z = 30
      if ( atn2 == 'Ga' ) mm_atoms(iatom)%z = 31
      if ( atn2 == 'Ge' ) mm_atoms(iatom)%z = 32
      if ( atn2 == 'As' ) mm_atoms(iatom)%z = 33
      if ( atn2 == 'Se' ) mm_atoms(iatom)%z = 34
      if ( atn2 == 'Br' ) mm_atoms(iatom)%z = 35
      if ( atn2 == 'Kr' ) mm_atoms(iatom)%z = 36
      if ( atn2 == 'Rb' ) mm_atoms(iatom)%z = 37
      if ( atn2 == 'Sr' ) mm_atoms(iatom)%z = 38
      if ( atn2 == 'Zr' ) mm_atoms(iatom)%z = 40
      if ( atn2 == 'Nb' ) mm_atoms(iatom)%z = 41
      if ( atn2 == 'Mo' ) mm_atoms(iatom)%z = 42
      if ( atn2 == 'Tc' ) mm_atoms(iatom)%z = 43
      if ( atn2 == 'Ru' ) mm_atoms(iatom)%z = 44
      if ( atn2 == 'Rh' ) mm_atoms(iatom)%z = 45
      if ( atn2 == 'Pd' ) mm_atoms(iatom)%z = 46
      if ( atn2 == 'Ag' ) mm_atoms(iatom)%z = 47
      if ( atn2 == 'Cd' ) mm_atoms(iatom)%z = 48
      if ( atn2 == 'In' ) mm_atoms(iatom)%z = 49
      if ( atn2 == 'Sn' ) mm_atoms(iatom)%z = 50
      if ( atn2 == 'Sb' ) mm_atoms(iatom)%z = 51
      if ( atn2 == 'Te' ) mm_atoms(iatom)%z = 52
      if ( atn2 == 'Xe' ) mm_atoms(iatom)%z = 54
      if ( atn2 == 'Cs' ) mm_atoms(iatom)%z = 55
      if ( atn2 == 'Ba' ) mm_atoms(iatom)%z = 56
      if ( atn2 == 'La' ) mm_atoms(iatom)%z = 57
      if ( atn2 == 'Ce' ) mm_atoms(iatom)%z = 58
      if ( atn2 == 'Pr' ) mm_atoms(iatom)%z = 59
      if ( atn2 == 'Nd' ) mm_atoms(iatom)%z = 60
      if ( atn2 == 'Pm' ) mm_atoms(iatom)%z = 61
      if ( atn2 == 'Sm' ) mm_atoms(iatom)%z = 62
      if ( atn2 == 'Eu' ) mm_atoms(iatom)%z = 63
      if ( atn2 == 'Gd' ) mm_atoms(iatom)%z = 64
      if ( atn2 == 'Tb' ) mm_atoms(iatom)%z = 65
      if ( atn2 == 'Dy' ) mm_atoms(iatom)%z = 66
      if ( atn2 == 'Ho' ) mm_atoms(iatom)%z = 67
      if ( atn2 == 'Er' ) mm_atoms(iatom)%z = 68
      if ( atn2 == 'Tm' ) mm_atoms(iatom)%z = 69
      if ( atn2 == 'Yb' ) mm_atoms(iatom)%z = 70
      if ( atn2 == 'Lu' ) mm_atoms(iatom)%z = 71
      if ( atn2 == 'Hf' ) mm_atoms(iatom)%z = 72
      if ( atn2 == 'Ta' ) mm_atoms(iatom)%z = 73
      if ( atn2 == 'Re' ) mm_atoms(iatom)%z = 75
      if ( atn2 == 'Os' ) mm_atoms(iatom)%z = 76
      if ( atn2 == 'Ir' ) mm_atoms(iatom)%z = 77
      if ( atn2 == 'Pt' ) mm_atoms(iatom)%z = 78
      if ( atn2 == 'Au' ) mm_atoms(iatom)%z = 79
      if ( atn2 == 'Hg' ) mm_atoms(iatom)%z = 80
      if ( atn2 == 'Tl' ) mm_atoms(iatom)%z = 81
      if ( atn2 == 'Pb' ) mm_atoms(iatom)%z = 82
      if ( atn2 == 'Bi' ) mm_atoms(iatom)%z = 83
      if ( atn2 == 'Po' ) mm_atoms(iatom)%z = 84
      if ( atn2 == 'At' ) mm_atoms(iatom)%z = 85
      if ( atn2 == 'Rn' ) mm_atoms(iatom)%z = 86
      if ( atn2 == 'Fr' ) mm_atoms(iatom)%z = 87
      if ( atn2 == 'Ra' ) mm_atoms(iatom)%z = 88
      if ( atn2 == 'Ac' ) mm_atoms(iatom)%z = 89
      if ( atn2 == 'Th' ) mm_atoms(iatom)%z = 90
      if ( atn2 == 'Pa' ) mm_atoms(iatom)%z = 91
      if ( atn2 == 'Np' ) mm_atoms(iatom)%z = 93
      if ( atn2 == 'Pu' ) mm_atoms(iatom)%z = 94

      if ( mm_atoms(iatom)%z == 0 ) &
        call die( "qmmm assign: There are MM atoms without atomic number." )

      mm_atoms(iatom)%mass = ATMASS( mm_atoms(iatom)%z )
      sym(katom)           = SYMBOL( mm_atoms(iatom)%z )

      if ( abs(mm_atoms(iatom)%mass) < 1.0d-10 ) &
        call die( "qmmm assign: There are MM atoms without mass." )
      katom = katom +1
    enddo

    ! Read atomic masses of different atoms from fdf block
    if ( fdf_block( 'NewMasses', bfdf ) ) then
      ! This is a very basic fdf block reading, it should be updated for
      ! newer versions of SIESTA block reading.
      do while ( fdf_bline(bfdf,pline) )
        if ( fdf_bnnames(pline) < 1 ) &
          call die("qmmm assign: Problem reading name from 'NewMasses' block.")
        name = trim(fdf_bnames(pline, 1))

        if ( fdf_bnvalues(pline) < 1 ) &
          call die("qmmm assign: Problem reading mass from 'NewMasses' block.")
        mass = fdf_bvalues(pline,1)

        write(6,"(/,a, A2, a, f6.2)") "qmmm assign: Read atomic mass of:  ", &
                                      name, "as", mass
        ! Assigns new masses
        do iatom = 1, na_qm
          if ( name == sym(iatom) ) qm_atoms(iatom)%mass = mass
        enddo
        do iatom = 1, na_mm
          if ( name == sym(na_qm+iatom) ) mm_atoms(iatom)%mass = mass
        enddo
      enddo
      call fdf_bclose( bfdf )
    endif

    deallocate( sym )
    return
  end subroutine assign_atoms
end module assign_qmmm