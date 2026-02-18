module mm_ene_frc
  !! This module contains subroutines that calculate the classical (MM-MM)
  !! interactions: bonds, angles, dihedrals, improper torsions, and Lennard-Jones.
  !
  ! To-do: Properly allocate all temporary arrays to avoid stack overflow issues.
  use mm_topology, only : fftopology, atom_connect_t, mm_atom_t
  use precision  , only : dp

  implicit none
  public :: mm_ene_fce
  public :: mm_dealloc

  ! Storage arrays.
  integer, allocatable, save :: ng1type(:,:)
    !! Type of bond for first neighbours.
  integer, allocatable, save :: angetype(:,:)
    !! Angle types when atom is in one of the extrema.
  integer, allocatable, save :: angmtype(:,:)
    !! Angle types when atom is in the middle.
  integer, allocatable, save :: dihety(:,:)
    !! Dihedral types  when atom is in one of the extrema.
  integer, allocatable, save :: dihmty(:,:)
    !! Dihedral types  when atom is in the middle.
  integer, allocatable, save :: evaldihe(:,:,:)
    !!
  integer, allocatable, save :: evaldihm(:,:,:)
    !!
  logical, allocatable, save :: evaldihelog(:,:)
    !!
  logical, allocatable, save :: evaldihmlog(:,:)
    !!
  integer, allocatable, save :: impty(:,:)
    !! Improper torsion types.
  integer, allocatable, save :: nglist(:)
    !! Neighbour list indexes.
  integer, allocatable, save :: nglistxat(:)
    !! Neighbours for a given atom.
  integer, allocatable, save :: nr(:,:)
    !! Neighbours across N number of cells.
  integer, allocatable, save :: nonbonded(:,:)
    !! Nonbonded pair interaction lists.
  integer, allocatable, save :: nonbondedxat(:)
    !! Nonbonded pair interactions for a given atom.
  integer, allocatable, save :: scaled(:,:)
    !! Scaled 1-4 pair interaction lists.
  integer, allocatable, save :: scalexat(:)
    !! Scaled 1-4 pair interactions for a given atom.

contains
  subroutine mm_dealloc( )
    !! Deallocates all module arrays.
    implicit none

    if ( allocated(ng1type     ) ) deallocate( ng1type      )
    if ( allocated(angetype    ) ) deallocate( angetype     )
    if ( allocated(angmtype    ) ) deallocate( angmtype     )
    if ( allocated(dihety      ) ) deallocate( dihety       )
    if ( allocated(dihmty      ) ) deallocate( dihmty       )
    if ( allocated(evaldihe    ) ) deallocate( evaldihe     )
    if ( allocated(evaldihm    ) ) deallocate( evaldihm     )
    if ( allocated(evaldihelog ) ) deallocate( evaldihelog  )
    if ( allocated(evaldihmlog ) ) deallocate( evaldihmlog  )
    if ( allocated(impty       ) ) deallocate( impty        )
    if ( allocated(nglist      ) ) deallocate( nglist       )
    if ( allocated(nglistxat   ) ) deallocate( nglistxat    )
    if ( allocated(nr          ) ) deallocate( nr           )
    if ( allocated(nonbonded   ) ) deallocate( nonbonded    )
    if ( allocated(nonbondedxat) ) deallocate( nonbondedxat )
    if ( allocated(scaled      ) ) deallocate( scaled       )
    if ( allocated(scalexat    ) ) deallocate( scalexat     )
  end subroutine mm_dealloc

  subroutine mm_ene_fce( natot, na_qm, na_mm, rclas, Etot_amber, fcetot_amber, &
                         stress_amber, mm_atoms, mm_connectivity, mm_top,      &
                         upd_ngb, rcut_mm, sfc, water, masst, cell, latt_typ,  &
                         coulombtype, ewald )
    !! This subroutine acts as a hub, joining all MM forces and energy
    !! contributions for output with the main driver.
    use coulomb_m, only : ewald_data_t
    use mm_units , only : kcal_mol, kcal_mol_eV
    use qmmm_pbc , only : reccel
    use units    , only : Ang

    implicit none
    integer          , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer          , intent(in)    :: na_mm
      !! Number of MM atoms.
    integer          , intent(in)    :: natot
      !! Total number of atoms (QM+MM).
    real(dp)         , intent(in)    :: rclas(3,natot)
      !! Positions of the MM atoms.
    type(mm_atom_t) , intent(in)     :: mm_atoms(na_mm)
      !! MM atom information.
    type(fftopology), intent(in)     :: mm_top
      !! Forcefield topology containing all classical bonded parameters.
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
    !! Data structure containing the connectivity data for each MM atom.

    real(dp)         , intent(in)    :: rcut_mm
      !! Distance cut-off for MM-MM interactions.
    real(dp)         , intent(in)    :: masst(natot)
      !! Atomic masses.
    real(dp)         , intent(in)    :: sfc
      !! Smoothing function cut-off. It is technically constant...
    logical          , intent(in)    :: water
      !! Whether we are using water restrain potentials.
    real(dp)         , intent(in)    :: cell(3,3)
      !! Periodic cell vectors.
    character(len=1) , intent(in)    :: latt_typ
      !! Type of PBC lattice.
    character(len=10), intent(in)    :: coulombtype
      !! When we are using Ewald or cut-off schemes.
    type(ewald_data_t), intent(in)   :: ewald
      !! Ewald summation parameters.

    real(dp)         , intent(out)   :: Etot_amber
      !! MM contributions to energy.
    real(dp)         , intent(out)   :: fcetot_amber(3,na_mm)
      !! MM contributions to forces.
    real(dp)         , intent(out)   :: stress_amber(3,3)
      !! MM contributions to cell stress.
    logical          , intent(inout) :: upd_ngb
      !! Whether we need to update neighbour lists.

    integer  :: iat
    real(dp) :: Ebond_amber, Eangle_amber, Edihe_amber, Eimp_amber, &
                Elj_amber, Eelec_amber, Elj_amber14, Eelec_amber14, &
                ewat, amber_cell(3,3), amber_kcell(3,3), cell_v

    real(dp), allocatable :: ramber(:,:), EmA(:), RmA(:)
    real(dp), allocatable :: fcebond_amber(:,:), fceangle_amber(:,:), &
                             fcedihe_amber(:,:), fceimp_amber(:,:),   &
                             fcelj_amber(:,:)  , fceelec_amber(:,:) , fwat(:,:)
    real(dp), external    :: volcel

    allocate( ramber(3,na_mm), EmA(na_mm), RmA(na_mm) )
    allocate( fcebond_amber(3,na_mm), fceangle_amber(3,na_mm), &
              fcedihe_amber(3,na_mm), fceimp_amber(3,na_mm),   &
              fcelj_amber(3,na_mm)  , fceelec_amber(3,na_mm) , fwat(3,na_mm) )


    ! Initializes forces and energies.
    Etot_amber    = 0.0_dp ; Ebond_amber   = 0.0_dp ; Eangle_amber   = 0.0_dp
    Edihe_amber   = 0.0_dp ; Eimp_amber    = 0.0_dp ; Elj_amber      = 0.0_dp
    Eelec_amber   = 0.0_dp ; Elj_amber14   = 0.0_dp ; Eelec_amber14  = 0.0_dp
    fcetot_amber  = 0.0_dp ; fcebond_amber = 0.0_dp ; fceangle_amber = 0.0_dp
    fcedihe_amber = 0.0_dp ; fceimp_amber  = 0.0_dp ; fcelj_amber    = 0.0_dp
    fceelec_amber = 0.0_dp ;
    stress_amber  = 0.0_dp ; ewat = 0.0_dp ; fwat = 0.0_dp ;

    ! Readapts arrays and converts to the necessary units.
    do iat = 1, na_mm
      ramber(1:3,iat) = rclas(1:3,na_qm+iat) / Ang
      EmA(iat) = mm_atoms(iat)%lj_Em * 2.0_dp * kcal_mol
      RmA(iat) = mm_atoms(iat)%lj_Rm / Ang * (2.0_dp ** (-5.0_dp / 6.0_dp))
    enddo
    amber_cell(1:3,1:3) = cell(1:3,1:3) / Ang
    call reccel( 3, amber_cell, amber_kcell, 0 )
    cell_v = volcel( amber_cell )

    ! Calculates bond energies and forces.
    call amber_bonds( na_mm, mm_connectivity, ramber, Ebond_amber, mm_atoms, &
                      mm_top, fcebond_amber, stress_amber, amber_cell,       &
                      amber_kcell, cell_v, latt_typ )

    ! Angle energies and forces.
    call amber_angles( na_mm, ramber, Eangle_amber, mm_atoms, mm_top, &
                       mm_connectivity, fceangle_amber, stress_amber, &
                       amber_cell, amber_kcell, cell_v, latt_typ )

    ! Dihedral energies and forces.
    call amber_dihes( na_mm, ramber, Edihe_amber, mm_atoms, mm_top,  &
                      mm_connectivity, fcedihe_amber, stress_amber,  &
                      amber_cell, amber_kcell, cell_v, latt_typ )

    ! Energies and forces from improper torsions.
    call amber_improper( na_mm, ramber, Eimp_amber, mm_atoms, mm_top, &
                         mm_connectivity, fceimp_amber, stress_amber, &
                         amber_cell, amber_kcell, cell_v, latt_typ )

    ! Energies and forces from nonbonded terms (i.e. LJ and coulomb)
    call amber_nonbonded( na_mm, mm_connectivity, ramber, Elj_amber,    &
                          Eelec_amber, Elj_amber14, Eelec_amber14, EmA, &
                          RmA, fceelec_amber, fcelj_amber, stress_amber,&
                          coulombtype, upd_ngb, rcut_mm, mm_atoms, sfc, &
                          amber_cell, ewald, amber_kcell, cell_v, latt_typ )

    ! Calculates the water restraint potential if needed.
    if ( water ) &
      call waters( na_qm, na_mm, natot, rclas, masst, mm_atoms, &
                   ewat, fwat, amber_cell, amber_kcell, latt_typ )

    ! Total energy and forces.
    Etot_amber = Ebond_amber + Eangle_amber + Edihe_amber + Eimp_amber + &
                 Elj_amber + Eelec_amber + Elj_amber14 + Eelec_amber14 + ewat

    write(*,'(A)') "siesta-qmmm: AMBER MM Region Energy Components: /eV"
    write(*,'(A20,F14.7)') "Bonds             = ", Ebond_amber / kcal_mol_eV
    write(*,'(A20,F14.7)') "Angles            = ", Eangle_amber/ kcal_mol_eV
    write(*,'(A20,F14.7)') "Dihedrals         = ", Edihe_amber / kcal_mol_eV
    write(*,'(A20,F14.7)') "Improper Torsions = ", Eimp_amber  / kcal_mol_eV
    write(*,'(A20,F14.7)') "Lennard-Jones     = ", Elj_amber   / kcal_mol_eV
    write(*,'(A20,F14.7)') "Lennard-Jones 1-4 = ", Elj_amber14 / kcal_mol_eV
    write(*,'(A20,F14.7)') "Non-bonded elec.  = ", Eelec_amber / kcal_mol_eV
    write(*,'(A20,F14.7)') "Non-bonded e. 1-4 = ", Eelec_amber14 / kcal_mol_eV
    if (water) write(*,'(A20,F14.7)') "Water terms       = ", eWat / kcal_mol_eV

    do iat = 1, na_mm
      fcetot_amber(1:3,iat) = fcebond_amber(1:3,iat) + fceangle_amber(1:3,iat) &
                            + fcedihe_amber(1:3,iat) + fceimp_amber(1:3,iat)   &
                            + fcelj_amber(1:3,iat)   + fceelec_amber(1:3,iat)  &
                            + fwat(1:3,iat)
      fcetot_amber(1:3,iat) = -1.0_dp * fcetot_amber(1:3,iat)
    enddo

    deallocate( fcebond_amber, fceangle_amber, fcedihe_amber, fceimp_amber, &
                fcelj_amber, fceelec_amber, fwat, ramber, EmA, RmA )
  end subroutine mm_ene_fce

  subroutine amber_bonds( na_mm, mm_connectivity, ramber, Ebond_amber,  &
                          mm_atoms, mm_top, fcebond_amber, stress_amber,&
                          cell, kcell, cell_v, latt_typ )
    !! Calculates energies and forces resulting from harmonic bonds.
    use qmmm_pbc , only : pbc_displ_vector

    implicit none
    integer         , intent(in)    :: na_mm
      !! Number of MM atoms.
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
    !! Data structure containing the connectivity data for each MM atom.

    real(dp)        , intent(in)    :: ramber(3,na_mm)
      !! Positions of the MM atoms.
    type(mm_atom_t), intent(in)     :: mm_atoms(na_mm)
      !! MM atom types.
    type(fftopology), intent(in)    :: mm_top
      !! Forcefield topology containing all classical bonded parameters.

    real(dp)        , intent(in)    :: cell(3,3)
      !! Periodic cell vectors.
    real(dp)        , intent(in)    :: kcell(3,3)
      !! Periodic cell vectors in the reciprocal space.
    real(dp)        , intent(in)    :: cell_v
      !! Unit cell volume.
    character(len=1), intent(in)    :: latt_typ
      !! Type of PBC lattice.
    real(dp)        , intent(inout) :: Ebond_amber
      !! Bond energies.
    real(dp)        , intent(inout) :: fcebond_amber(3,na_mm)
      !! Bond contributions to forces.
    real(dp)        , intent(inout) :: stress_amber(3,3)
      !! Bond contributions to cell stress.

    character(len=4) :: ty1, ty2, ty3
    character(len=5) :: tybond
    integer          :: iat, jat, ibnd, icrd
    real(dp)         :: stress_fact, dr(3), dg(3)

    logical, save :: first = .true.

    ! In the first run we assign bond types.
    if ( first ) then
      allocate( ng1type(na_mm,6) )
      do iat = 1, na_mm
      do jat = 1, mm_connectivity(iat)%nbonds
        do ibnd = 1, mm_top%nbonds
          tybond = mm_top%bonds(ibnd)%type
          ty1    = tybond(1:2)
          ty2    = tybond(4:5)
          ty3    = mm_atoms( mm_connectivity(iat)%bond_at(jat) )%attype

          if ( (mm_atoms(iat)%attype == ty1) .and. (ty3 == ty2) ) then
            ng1type(iat,jat) = ibnd
          elseif ( (mm_atoms(iat)%attype == ty2) .and. (ty3 == ty1) ) then
            ng1type(iat,jat) = ibnd
          endif
        enddo
      enddo
      enddo

      first = .false.
    endif

    ! Proper forces and energies calculation.
    stress_fact = 1.0_dp / cell_v
    do iat = 1, na_mm
    do jat = 1, mm_connectivity(iat)%nbonds
      dr(1:3) = ramber(1:3,iat) - ramber(1:3,mm_connectivity(iat)%bond_at(jat))
      call pbc_displ_vector( latt_typ, cell, kcell, dr )

      Ebond_amber = Ebond_amber + mm_top%bonds( ng1type(iat,jat) )%energy( dr )

      call mm_top%bonds( ng1type(iat,jat) )%grad( dr, dg )

      fcebond_amber(:,iat) = fcebond_amber(:,iat) + dg(:)
      do icrd = 1, 3
        stress_amber(1:3,icrd)  = stress_amber(1:3,icrd) + &
                                  stress_fact * ramber(1:3,iat) * dg(icrd)
      enddo
    enddo
    enddo

    Ebond_amber = 0.5_dp * Ebond_amber
  end subroutine amber_bonds

  subroutine amber_angles( na_mm, ramber, Eangle_amber, mm_atoms, mm_top, &
                           mm_connectivity, fceangle_amber, stress_amber, &
                           cell, kcell, cell_v, latt_typ )
    !! Calculates energies and forces resulting from harmonic angles.
    use qmmm_pbc , only : pbc_displ_vector

    implicit none
    integer         , intent(in)    :: na_mm
      !! Number of MM atoms.
    real(dp)        , intent(in)    :: ramber(3,na_mm)
      !! Positions of the MM atoms.
    type(mm_atom_t) , intent(in)    :: mm_atoms(na_mm)
      !! MM atom types.
    type(fftopology), intent(in)    :: mm_top
      !! Forcefield topology containing all classical bonded parameters.
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.
    real(dp)        , intent(in)    :: cell(3,3)
      !! Periodic cell vectors.
    real(dp)        , intent(in)    :: kcell(3,3)
      !! Periodic cell vectors in the reciprocal space.
    real(dp)        , intent(in)    :: cell_v
      !! Unit cell volume.
    character(len=1), intent(in)    :: latt_typ
      !! Type of PBC lattice.
    real(dp)        , intent(inout) :: Eangle_amber
      !! Angle energies.
    real(dp)        , intent(inout) :: fceangle_amber(3,na_mm)
      !! Angle contributions to forces.
    real(dp)        , intent(inout) :: stress_amber(3,3)
      !! Angle contributions to cell stress.

    real(dp)         :: dr(3), dr12(3), dr32(3), stress_fact
    integer          :: iat, iang, jang, icrd
    character(len=4) :: ty1, ty2, ty3
    character(len=8) :: tyangle

    real(dp), allocatable :: fext(:,:), fmid(:,:)

    logical, save :: first = .true.

    allocate( fext(3,na_mm), fmid(3,na_mm) )
    fext  = 0.0_dp
    fmid  = 0.0_dp

    ! First we assign angle types.
    if ( first ) then
      allocate( angetype(na_mm,25), angmtype(na_mm,25) )
      do iat  = 1, na_mm
      do jang = 1, mm_connectivity(iat)%nangl_e
      do iang = 1, mm_top%nangles
        tyangle = mm_top%angles(iang)%type
        ty1 = tyangle(1:2)
        ty2 = tyangle(4:5)
        ty3 = tyangle(7:8)

        if ( (mm_atoms(iat)%attype == ty1) .and. &
             (mm_atoms(mm_connectivity(iat)%ange_at(jang,1))%attype == ty2) &
             .and. &
             (mm_atoms(mm_connectivity(iat)%ange_at(jang,2))%attype == ty3) ) then
          angetype(iat,jang) = iang

        elseif ( (mm_atoms(iat)%attype == ty3) .and. &
                   (mm_atoms(mm_connectivity(iat)%ange_at(jang,1))%attype == ty2)&
             .and. (mm_atoms(mm_connectivity(iat)%ange_at(jang,2))%attype == ty1) ) then
          angetype(iat,jang) = iang
        endif
      enddo
      enddo
      enddo

      do iat  = 1, na_mm
      do jang = 1, mm_connectivity(iat)%nangl_m
      do iang = 1, mm_top%nangles
        tyangle = mm_top%angles(iang)%type
        ty1     = tyangle(1:2)
        ty2     = tyangle(4:5)
        ty3     = tyangle(7:8)

        if ( (mm_atoms(iat)%attype == ty2) .and. &
             (mm_atoms(mm_connectivity(iat)%angm_at(jang,1))%attype == ty1) &
             .and. (mm_atoms(mm_connectivity(iat)%angm_at(jang,2))%attype == ty3) ) then
          angmtype(iat,jang) = iang

        elseif ( (mm_atoms(iat)%attype == ty2) .and. &
                 (mm_atoms(mm_connectivity(iat)%angm_at(jang,1))%attype == ty3) &
              .and. (mm_atoms(mm_connectivity(iat)%angm_at(jang,2))%attype == ty1) ) then
          angmtype(iat,jang) = iang
        endif
      enddo
      enddo
      enddo

      first = .false.
    endif

    stress_fact = 1.0_dp / cell_v

    ! Atoms in the angle extreme.
    do iat  = 1, na_mm
    do jang = 1, mm_connectivity(iat)%nangl_e
      dr12(1:3) = ramber(1:3,iat) &
                - ramber(1:3,mm_connectivity(iat)%ange_at(jang,1))
      dr32(1:3) = ramber(1:3,mm_connectivity(iat)%ange_at(jang,2)) &
                - ramber(1:3,mm_connectivity(iat)%ange_at(jang,1))

      call pbc_displ_vector( latt_typ, cell, kcell, dr12 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr32 )

      Eangle_amber = Eangle_amber + &
                     mm_top%angles( angetype(iat,jang) )%energy( dr12, dr32 )

      call mm_top%angles( angetype(iat,jang) )%grad_e( dr12, dr32, dr )

      fext(:,iat) = fext(:,iat) + dr(:)
      do icrd = 1, 3
        stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + &
                                 stress_fact * ramber(1:3,iat) * dr(icrd)
      enddo
    enddo
    enddo

    ! Atoms in the angle center.
    do iat  = 1, na_mm
    do jang = 1, mm_connectivity(iat)%nangl_m
      dr12(:) = ramber(:,mm_connectivity(iat)%angm_at(jang,1)) &
              - ramber(:,iat)
      dr32(:) = ramber(:,mm_connectivity(iat)%angm_at(jang,2)) &
              - ramber(:,iat)

      call pbc_displ_vector( latt_typ, cell, kcell, dr12 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr32 )

      Eangle_amber = Eangle_amber + &
                     mm_top%angles( angmtype(iat,jang) )%energy( dr12, dr32 )

      call mm_top%angles( angmtype(iat,jang) )%grad_m( dr12, dr32, dr )

      fmid(:,iat) = fmid(:,iat) + dr(:)
      do icrd = 1, 3
        stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + &
                                 stress_fact * ramber(1:3,iat) * dr(icrd)
      enddo
    enddo
    enddo

    Eangle_amber = Eangle_amber / 3.0_dp

    do iat = 1, na_mm
      fceangle_amber(1:3,iat) = fext(1:3,iat) + fmid(1:3,iat)
    enddo

    deallocate( fext, fmid )
  end subroutine amber_angles

  subroutine amber_dihes( na_mm, ramber, Edihe_amber, mm_atoms, mm_top,     &
                          mm_connectivity, fcedihe_amber,&
                          stress_amber, cell, kcell, cell_v, latt_typ )
    !! Calculates energy and forces contributions from dihedral angles.
    use qmmm_pbc , only : pbc_displ_vector
    use precision, only : dp

    implicit none
    integer          , intent(in)    :: na_mm
      !! Number of MM atoms.
    real(dp)         , intent(in)    :: ramber(3,na_mm)
      !! Positions of the MM atoms.
    type(mm_atom_t) , intent(in)     :: mm_atoms(na_mm)
      !! MM atom types.
    type(fftopology), intent(in)     :: mm_top
      !! Forcefield topology containing all classical bonded parameters.
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.

    real(dp)         , intent(in)    :: cell(3,3)
      !! Periodic cell vectors.
    real(dp)         , intent(in)    :: kcell(3,3)
      !! Periodic cell vectors in the reciprocal space.
    real(dp)         , intent(in)    :: cell_v
      !! Unit cell volume.
    character(len=1) , intent(in)    :: latt_typ
      !! Type of PBC lattice.
    real(dp)         , intent(inout) :: Edihe_amber
      !! Dihedral energies.
    real(dp)         , intent(inout) :: fcedihe_amber(3,na_mm)
      !! Dihedral contributions to forces.
    real(dp)         , intent(inout) :: stress_amber(3,3)
      !! Dihedral contributions to cell stress.


    integer           :: iat, idihe, jdihe, ldihe, icrd, at_dih(3)
    character(len=4)  :: ty1, ty2, ty3, ty4
    character(len=11) :: tydihe
    real(dp)          :: dr12(3), dr32(3), dr43(3), dr31(3), dr42(3), &
                         fce(3), stress_fact

    real(dp), allocatable :: fext(:,:), fmid(:,:), perdihe(:)

    logical , save :: first = .true.

    allocate( fext(3,na_mm), fmid(3,na_mm), perdihe(mm_top%ndihe) )
    do idihe = 1, mm_top%ndihe
      perdihe(idihe) = mm_top%dihedrals(idihe)%per
    enddo
    fext  = 0.0_dp
    fmid  = 0.0_dp

    ! We assign dihedral types in the first run.
    if ( first ) then
      allocate( dihety(na_mm,100)     , dihmty(na_mm,100), &
                evaldihe(na_mm,100,5) , evaldihm(na_mm,100,5) )
      allocate( evaldihelog(na_mm,100), evaldihmlog(na_mm,100) )

      evaldihelog = .false.
      evaldihmlog = .false.
      evaldihm    = 0
      evaldihe    = 0

      ! Atom on one extreme of the dihedral.
      do iat   = 1, na_mm
      do ldihe = 1, mm_connectivity(iat)%ndihe_e
        dihety(iat,ldihe) = 1
        jdihe = 0

        do idihe = 1, mm_top%ndihe
          tydihe = mm_top%dihedrals(idihe)%type
          ty1 = tydihe(1:2)
          ty2 = tydihe(4:5)
          ty3 = tydihe(7:8)
          ty4 = tydihe(10:11)

          at_dih(1:3) = mm_connectivity(iat)%dihe_at(ldihe,1:3)
          if ( any(at_dih(1:3) == 0) ) cycle

          if ( ty1 == 'X ' ) then
            if ( (mm_atoms(at_dih(1))%attype == ty2) .and. &
                 (mm_atoms(at_dih(2))%attype == ty3) ) then
              dihety(iat,ldihe) = idihe

            elseif ( (mm_atoms(at_dih(1))%attype == ty3) .and. &
                     (mm_atoms(at_dih(2))%attype == ty2) ) then
              dihety(iat,ldihe) = idihe
            endif

          else
            if ( (mm_atoms(iat)%attype       == ty1) .and. &
                 (mm_atoms(at_dih(1))%attype == ty2) .and. &
                 (mm_atoms(at_dih(2))%attype == ty3) .and. &
                 (mm_atoms(at_dih(3))%attype == ty4) ) then

              dihety(iat,ldihe) = idihe
              if ( perdihe(idihe) < 0.0_dp ) then
                evaldihelog(iat,ldihe) = .true.
                jdihe = jdihe +1
                evaldihe(iat,ldihe,jdihe)   = idihe
                evaldihe(iat,ldihe,jdihe+1) = idihe +1
              endif

            elseif ( (mm_atoms(iat)%attype       == ty4) .and. &
                     (mm_atoms(at_dih(1))%attype == ty3) .and. &
                     (mm_atoms(at_dih(2))%attype == ty2) .and. &
                     (mm_atoms(at_dih(3))%attype == ty1) ) then

              dihety(iat,ldihe) = idihe
              if ( perdihe(idihe) < 0.0_dp ) then
                evaldihelog(iat,ldihe) = .true.
                jdihe = jdihe +1
                evaldihe(iat,ldihe,jdihe)   = idihe
                evaldihe(iat,ldihe,jdihe+1) = idihe +1
              endif
            endif
          endif
        enddo ! idihe
      enddo   ! ldihe
      enddo   ! iat

      ! Atom in the middle of the dihedral.
      do iat   = 1, na_mm
      do ldihe = 1, mm_connectivity(iat)%ndihe_m
        dihmty(iat,ldihe) = 1
        jdihe = 0

        do idihe = 1, mm_top%ndihe
          tydihe = mm_top%dihedrals(idihe)%type
          ty1 = tydihe(1:2)
          ty2 = tydihe(4:5)
          ty3 = tydihe(7:8)
          ty4 = tydihe(10:11)

          at_dih(1:3) = mm_connectivity(iat)%dihm_at(ldihe,1:3)
          if ( any(at_dih(1:3) == 0) ) cycle

          if ( ty1 == 'X ' ) then
            if ( (mm_atoms(iat)%attype       == ty2) .and. &
                 (mm_atoms(at_dih(2))%attype == ty3) ) then
              dihmty(iat,ldihe) = idihe

            elseif ( (mm_atoms(iat)%attype       == ty3) .and. &
                     (mm_atoms(at_dih(2))%attype == ty2) ) then
              dihmty(iat,ldihe) = idihe
            endif

          else
            if ( (mm_atoms(at_dih(1))%attype == ty1) .and. &
                 (mm_atoms(iat)%attype       == ty2) .and. &
                 (mm_atoms(at_dih(2))%attype == ty3) .and. &
                 (mm_atoms(at_dih(3))%attype == ty4) ) then

              dihmty(iat,ldihe) = idihe
              if ( perdihe(idihe) < 0.0_dp ) then
                evaldihmlog(iat,ldihe) = .true.
                jdihe = jdihe +1
                evaldihm(iat,ldihe,jdihe)   = idihe
                evaldihm(iat,ldihe,jdihe+1) = idihe +1
              endif

            elseif ( (mm_atoms(at_dih(1))%attype == ty4) .and. &
                     (mm_atoms(iat)%attype       == ty3) .and. &
                     (mm_atoms(at_dih(2))%attype == ty2) .and. &
                     (mm_atoms(at_dih(3))%attype == ty1) ) then

              dihmty(iat,ldihe) = idihe
              if ( perdihe(idihe) < 0.0_dp ) then
                evaldihmlog(iat,ldihe) = .true.
                jdihe = jdihe +1
                evaldihm(iat,ldihe,jdihe)   = idihe
                evaldihm(iat,ldihe,jdihe+1) = idihe +1
              endif
            endif
          endif
        enddo ! idihe
      enddo   ! ldihe
      enddo   ! iat

      first = .false.
    endif ! First run.

    ! Energy and forces contributions from dihedrals with the atom
    ! in a extreme.
    do iat = 1, mm_top%ndihe
      if ( perdihe(iat) < 0.0_dp ) perdihe(iat) = -perdihe(iat)
    enddo
    stress_fact = 1.0_dp / cell_v

    do iat   = 1, na_mm
    do ldihe = 1, mm_connectivity(iat)%ndihe_e
      at_dih(1:3) = mm_connectivity(iat)%dihe_at(ldihe,1:3)
      if ( any(at_dih(1:3) == 0) ) cycle

      dr12 = ramber(:,iat)       - ramber(:,at_dih(1))
      dr32 = ramber(:,at_dih(2)) - ramber(:,at_dih(1))
      dr31 = ramber(:,at_dih(2)) - ramber(:,iat)
      dr43 = ramber(:,at_dih(3)) - ramber(:,at_dih(2))
      dr42 = ramber(:,at_dih(3)) - ramber(:,at_dih(1))
      call pbc_displ_vector( latt_typ, cell, kcell, dr12 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr32 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr31 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr43 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr42 )

      if ( evaldihelog(iat,ldihe) ) then
        do idihe = 1, 5
          if ( evaldihe(iat,ldihe,idihe) /= 0 ) then
            jdihe = evaldihe(iat,ldihe,idihe)
            Edihe_amber = Edihe_amber + &
                          mm_top%dihedrals(jdihe)%energy( dr12, dr32, dr43 )
            call mm_top%dihedrals(jdihe)%grad( fce, 1, dr12, dr32, dr43, &
                                                             dr31, dr42 )

            do icrd = 1, 3
              fext(icrd,iat) = fext(icrd,iat) + fce(icrd)
              stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact &
                                     * ramber(1:3,iat) * fce(icrd)
            enddo
          endif
        enddo
      else
        idihe = dihety(iat,ldihe)
        Edihe_amber = Edihe_amber + &
                      mm_top%dihedrals(idihe)%energy( dr12, dr32, dr43 )
        call mm_top%dihedrals(idihe)%grad( fce, 1, dr12, dr32, dr43, &
                                                         dr31, dr42 )

        do icrd = 1, 3
          fext(icrd,iat) = fext(icrd,iat) + fce(icrd)
          stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact &
                                 * ramber(1:3,iat) * fce(icrd)
        enddo
      endif
    enddo
    enddo

    ! Energy and forces contributions from dihedrals with the atom
    ! in the middle.
    do iat   = 1, na_mm
    do ldihe = 1, mm_connectivity(iat)%ndihe_m
      at_dih(1:3) = mm_connectivity(iat)%dihm_at(ldihe,1:3)
      if ( any(at_dih(1:3) == 0) ) cycle

      dr12 = ramber(:,at_dih(1)) - ramber(:,iat)
      dr32 = ramber(:,at_dih(2)) - ramber(:,iat)
      dr31 = ramber(:,at_dih(2)) - ramber(:,at_dih(1))
      dr43 = ramber(:,at_dih(3)) - ramber(:,at_dih(2))
      dr42 = ramber(:,at_dih(3)) - ramber(:,iat)
      call pbc_displ_vector( latt_typ, cell, kcell, dr12 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr32 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr31 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr43 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr42 )

      if ( evaldihmlog(iat,ldihe) ) then
        do idihe = 1, 5
          if ( evaldihm(iat,ldihe,idihe) /= 0 ) then
            jdihe = evaldihm(iat,ldihe,idihe)
            Edihe_amber = Edihe_amber + &
                          mm_top%dihedrals(jdihe)%energy( dr12, dr32, dr43 )
            call mm_top%dihedrals(jdihe)%grad( fce, 2, dr12, dr32, dr43, &
                                               dr31, dr42 )

            do icrd = 1, 3
              fmid(icrd,iat) = fmid(icrd,iat) + fce(icrd)
              stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact &
                                     * ramber(1:3,iat) * fce(icrd)
            enddo
          endif
        enddo
      else
        idihe = dihmty(iat,ldihe)
        Edihe_amber = Edihe_amber + &
                      mm_top%dihedrals(idihe)%energy( dr12, dr32, dr43 )
        call mm_top%dihedrals(idihe)%grad( fce, 2, dr12, dr32, dr43, &
                                           dr31, dr42 )

        do icrd = 1, 3
          fmid(icrd,iat) = fmid(icrd,iat) + fce(icrd)
          stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact &
                                 * ramber(1:3,iat) * fce(icrd)
        enddo
      endif
    enddo
    enddo

    Edihe_amber = 0.25_dp * Edihe_amber

    do iat = 1, na_mm
      fcedihe_amber(1:3,iat) = fext(1:3,iat) + fmid(1:3,iat)
    enddo

    deallocate( fext, fmid, perdihe )
  end subroutine amber_dihes

  subroutine amber_improper( na_mm, ramber, Eimp_amber, mm_atoms, mm_top, &
                             mm_connectivity, fimp, stress_amber, cell, &
                             kcell, cell_v, latt_typ )
    use qmmm_pbc , only : pbc_displ_vector
    use precision, only : dp

    implicit none
    integer          , intent(in)    :: na_mm
      !! Number of MM atoms.
    real(dp)         , intent(in)    :: ramber(3,na_mm)
      !! Positions of the MM atoms.
    type(mm_atom_t) , intent(in)     :: mm_atoms(na_mm)
      !! MM atom types.
    type(fftopology), intent(in)     :: mm_top
      !! Forcefield topology containing all classical bonded parameters.
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.

    real(dp)         , intent(in)    :: cell(3,3)
      !! Periodic cell vectors.
    real(dp)         , intent(in)    :: kcell(3,3)
      !! Periodic cell vectors in the reciprocal space.
    real(dp)         , intent(in)    :: cell_v
      !! Unit cell volume.
    character(len=1) , intent(in)    :: latt_typ
      !! Type of PBC lattice.
    real(dp)         , intent(inout) :: Eimp_amber
      !! Improper torsions' energies.
    real(dp)         , intent(inout) :: fimp(3,na_mm)
      !! Improper torsion contributions to forces.
    real(dp)         , intent(inout) :: stress_amber(3,3)
      !! Improper torsion contributions to cell stress.

    character(len=4 ) :: ty1, ty2, ty3, ty4
    character(len=11) :: tyimp
    integer           :: iat, iimp, jimp, icrd, at_i(4), atin
    real(dp)          :: dr12(3), dr32(3), dr43(3), dr31(3), dr42(3), &
                         fce(3), stress_fact

    logical, save :: first = .true.

    ! Assign improper types in the first call.
    if (first) then
      allocate( impty(na_mm,25) )
      do iat  = 1, na_mm
      do jimp = 1, mm_connectivity(iat)%nimp
        impty(iat,jimp) = 1

        do iimp = 1, mm_top%nimp
          tyimp = mm_top%impropers(iimp)%type
          ty1   = tyimp(1:2)
          ty2   = tyimp(4:5)
          ty3   = tyimp(7:8)
          ty4   = tyimp(10:11)
          at_i(1:4) = mm_connectivity(iat)%imp_at(jimp,1:4)

          if ( (ty1 == 'X ') .and. (ty2 == 'X ') .and. (ty4 == 'X ') ) then
            if ( mm_atoms(at_i(3))%attype == ty3 ) then
              impty(iat,jimp) = iimp
            elseif ( mm_atoms(at_i(2))%attype == ty3 ) then
              impty(iat,jimp) = iimp
            endif

          elseif ( (ty1 == 'X ') .and. (ty2 == 'X ') ) then
            if ( (mm_atoms(at_i(3))%attype == ty3) .and. (mm_atoms(at_i(4))%attype == ty4) ) then
              impty(iat,jimp) = iimp
            elseif ( (mm_atoms(at_i(1))%attype == ty4) .and. &
                     (mm_atoms(at_i(2))%attype == ty3) ) then
              impty(iat,jimp) = iimp
            endif

          elseif ( ty1 == 'X ' ) then
            if ( (mm_atoms(at_i(2))%attype == ty2) .and. (mm_atoms(at_i(3))%attype == ty3) .and.&
                 (mm_atoms(at_i(4))%attype == ty4) ) then

              impty(iat,jimp) = iimp
            elseif ( (mm_atoms(at_i(3))%attype == ty2) .and. (mm_atoms(at_i(2))%attype == ty3) &
                     .and. (mm_atoms(at_i(1))%attype == ty4) ) then

              impty(iat,jimp) = iimp
            endif
          else
            if ( (mm_atoms(at_i(1))%attype == ty1) .and. (mm_atoms(at_i(2))%attype == ty2) .and.&
                 (mm_atoms(at_i(3))%attype == ty3) .and. (mm_atoms(at_i(4))%attype == ty4) ) then

              impty(iat,jimp) = iimp
            elseif ( (mm_atoms(at_i(4))%attype == ty1) .and. (mm_atoms(at_i(3))%attype == ty2) &
               .and. (mm_atoms(at_i(2))%attype == ty3) .and. (mm_atoms(at_i(1))%attype == ty4) &
               ) then

              impty(iat,jimp) = iimp
            endif
          endif
        enddo ! iimp
      enddo   ! jimp
      enddo   ! iat

      first = .false.
    endif ! First call.

    stress_fact = 1.0_dp / cell_v

    ! Calculates forces and energy contributions from improper torsions.
    do iat  = 1, na_mm
    do jimp = 1, mm_connectivity(iat)%nimp
      at_i(1:4) = mm_connectivity(iat)%imp_at(jimp,1:4)
      iimp      = impty(iat,jimp)

      dr12 = ramber(:,at_i(1)) - ramber(:,at_i(2))
      dr32 = ramber(:,at_i(3)) - ramber(:,at_i(2))
      dr31 = ramber(:,at_i(3)) - ramber(:,at_i(1))
      dr43 = ramber(:,at_i(4)) - ramber(:,at_i(3))
      dr42 = ramber(:,at_i(4)) - ramber(:,at_i(2))
      call pbc_displ_vector( latt_typ, cell, kcell, dr12 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr32 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr31 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr43 )
      call pbc_displ_vector( latt_typ, cell, kcell, dr42 )

      Eimp_amber = Eimp_amber + mm_top%impropers(iimp)%energy(dr12, dr32, dr43)

      ! Searches which atom of the improper is the current atom.
      if ( iat == at_i(1) ) atin = 1
      if ( iat == at_i(2) ) atin = 2
      if ( iat == at_i(3) ) atin = 3
      if ( iat == at_i(4) ) atin = 4

      call mm_top%impropers(iimp)%grad( fce, atin, dr12, dr32, dr43, &
                                        dr31, dr42 )

      do icrd = 1, 3
        fimp(icrd,iat) = fimp(icrd,iat) + fce( icrd )
        stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                 ramber(1:3,iat) * fce( icrd )
      enddo
    enddo
    enddo

    Eimp_amber = 0.25_dp * Eimp_amber

  end subroutine amber_improper

  subroutine amber_nonbonded( na_mm, mm_connectivity, ramber, Elj_amber, &
                              Eelec_amber, Elj_amber14, Eelec_amber14,   &
                              Em, Rm, felec, flj, stress_amber, coulombtype,&
                              upd_ngb, rcut_mm, mm_atoms, sfc, cell,     &
                              ewald, kcell, cell_v, latt_typ )
    !! Calculates non-bonded contributions to energies and forces, i.e.
    !! coulomb and Lennard-Jones terms.
    use coulomb_m     , only : ewald_data_t
    use functions     , only : scalar_v2
    use mm_units      , only : kcal_mol_eV, eps0
    use qmmm_mm_neighbour, only : mm_mneighb, mm_jan, mm_r2ij, mm_xij, maxnna
    use qmmm_pbc      , only : get_pbc_vectors, pbc_displ_vector
    use sys           , only : die
    use units         , only : pi, pi2, Ang, eV

    implicit none
    integer          , intent(in)    :: na_mm
      !! Number of MM atoms.
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.
    real(dp)         , intent(in)    :: ramber(3,na_mm)
      !! Positions of the MM atoms.
    real(dp)         , intent(in)    :: Em(na_mm)
      !! Lennard-Jones epsilon.
    real(dp)         , intent(in)    :: Rm(na_mm)
      !! Lennard-Jones Rmin.

    type(mm_atom_t)  , intent(in)    :: mm_atoms(na_mm)
      !! MM atoms info.
    real(dp)         , intent(in)    :: sfc
      !! Smoothing function cut-off. It is technically constant...

    real(dp)         , intent(in)    :: rcut_mm
      !! Distance cut-off for MM-MM interactions.

    real(dp)         , intent(inout) :: cell(3,3)
      !! Periodic cell vectors.
    real(dp)         , intent(in)    :: kcell(3,3)
      !! Periodic cell vectors in the reciprocal space.
    real(dp)         , intent(in)    :: cell_v
      !! Unit cell volume.
    character(len=1) , intent(in)    :: latt_typ
      !! Type of PBC lattice.
    character(len=10), intent(in)    :: coulombtype
      !! When we are using Ewald or cut-off schemes.
    type(ewald_data_t), intent(in)   :: ewald
      !! Ewald summation parameters.

    logical          , intent(inout) :: upd_ngb
      !! Whether we update neighbour lists.
    real(dp)         , intent(inout) :: Elj_amber
      !! Lennard-Jones energies.
    real(dp)         , intent(inout) :: Elj_amber14
      !! Lennard-Jones energies for 1-4 scaled interactions.
    real(dp)         , intent(inout) :: Eelec_amber
      !! Coulomb interaction energies.
    real(dp)         , intent(inout) :: Eelec_amber14
      !! Coulomb interaction energies for 1-4 scaled interactions.
    real(dp)         , intent(inout) :: felec(3,na_mm)
      !! Coulomb contributions to forces.
    real(dp)         , intent(inout) :: flj(3,na_mm)
      !! Lennard-Jones contributions to forces.
    real(dp)         , intent(inout) :: stress_amber(3,3)
      !! Coulomb and LJ contributions to cell stress.

    integer  :: dimvec, n_pointer, nr_indx(3), n1, n2, n3, n1max, n2max, n3max,&
                iat, jat, kat, lat, mat, icrd, ivec, inei, nnei, ni
    logical  :: skip_at, skip_k
    real(dp) :: Rij, Eij, drdist, unit_f, E2, dr2, Rij6, dr2_3, fac, ca, cb, &
                cd, rinn, rout, drij(3), fel, const2, const3, const4, const5,&
                const6, const7, kmod2, kronij, De, De2, S_real, S_imag, kr,  &
                k_neigh(6,3), krecip_displaced(3), SS, krecip(3), stress_fact,&
                lattice_vector_len, vector_len_max

    integer , allocatable      :: nglistemp(:), nrtemp(:,:)
    real(dp), allocatable      :: pc_cos_kr(:), pc_sin_kr(:)

    integer , save :: nna     = 200
    logical , save :: first   = .true.
    logical , save :: bigcell = .true.
    real(dp), save :: rcoor   = 0.0_dp

    allocate( pc_cos_kr(na_mm), pc_sin_kr(na_mm) )

    ! These units include a permitivity (epsilon) of 1.0.
    unit_f = kcal_mol_eV / ( 2.0_dp * pi2 * eps0 * Ang * eV )

    fel    = 0.0_dp
    flj    = 0.0_dp
    felec  = 0.0_dp

    dimvec = na_mm * 3000
    const2      = 2.0_dp * pi2 / cell_v
    stress_fact = 1.0_dp / cell_v

    drij(:) = kcell(:,1)
    n1max = int( ewald%kcut / ( pi2 * sqrt( scalar_v2( drij, drij ) ) ) )

    drij(:) = kcell(:,2)
    n2max = int( ewald%kcut / ( pi2 * sqrt( scalar_v2( drij, drij ) ) ) )

    drij(:) = kcell(:,3)
    n3max = int( ewald%kcut / ( pi2 * sqrt( scalar_v2( drij, drij ) ) ) )

    ! Assigns atoms to the nonbonded interactions list.
    if ( first ) then
      allocate( nonbonded(na_mm,100), scaled(na_mm,100), scalexat(na_mm), &
                nonbondedxat(na_mm) )

      do iat = 1, na_mm
        kat = 1

        do jat=1, mm_connectivity(iat)%nbonds
          if ( iat < mm_connectivity(iat)%bond_at(jat) ) then
            nonbonded(iat,kat) = mm_connectivity(iat)%bond_at(jat)
            kat = kat +1
          endif
        enddo
        do jat = 1, mm_connectivity(iat)%nangl_e
          if ( iat < mm_connectivity(iat)%ange_at(jat,2) ) then
            nonbonded(iat,kat) = mm_connectivity(iat)%ange_at(jat,2)
            kat = kat +1
          endif
        enddo
        nonbondedxat(iat) = kat -1

        ! Looks for 1-4 bonded (scaled)
        kat = 1
        do jat = 1, mm_connectivity(iat)%ndihe_e
          if ( iat < mm_connectivity(iat)%dihe_at(jat,3) ) then
            skip_at = .false.

            do lat = 1, kat-1 ! Skips repeated dihedrals.
              if ( mm_connectivity(iat)%dihe_at(jat,3) == scaled(iat,lat) ) &
                skip_at = .true.
            enddo
            do lat = 1, nonbondedxat(iat)
              ! Skips if already accounted for as 1-3 interaction.
              if ( mm_connectivity(iat)%dihe_at(jat,3) == nonbonded(iat,lat) ) &
                skip_at = .true.
            enddo

            if ( skip_at ) cycle
            scaled(iat,kat) = mm_connectivity(iat)%dihe_at(jat,3)
            kat = kat +1
          endif
        enddo
        scalexat(iat) = kat -1
      enddo

      vector_len_max = 0.0_dp
      do icrd = 1, 3
        lattice_vector_len = sqrt( cell(icrd,1)*cell(icrd,1)  + &
            cell(icrd,2)*cell(icrd,2)+cell(icrd,3)*cell(icrd,3) )

        if ( lattice_vector_len > vector_len_max ) &
          vector_len_max = lattice_vector_len
      enddo

      ! sfc: smooth-function cut-off skin in the neighbor list of 2Ang.
      rcoor = rcut_mm + sfc + 2.0_dp
      bigcell = .false.
      if ( (2.0_dp * rcoor) < vector_len_max ) bigcell = .true.

      ! Initialise neighbour search to create the cell arrays.
      call mm_mneighb( cell, rcoor, na_mm, ramber, 0, 1, nna )
      if ( associated(mm_jan ) ) deallocate( mm_jan )
      if ( associated(mm_r2ij) ) deallocate( mm_r2ij)
      if ( associated(mm_xij ) ) deallocate( mm_xij )
      allocate( mm_jan(maxnna), mm_r2ij(maxnna), mm_xij(3,maxnna) )

      first = .false.
    endif ! First call

    if ( upd_ngb ) then
      ! Initialise neighbour search to create the cell arrays.
      call mm_mneighb( cell, rcoor, na_mm, ramber, 0, 1, nna )
      if ( associated(mm_jan ) ) deallocate( mm_jan )
      if ( associated(mm_r2ij) ) deallocate( mm_r2ij)
      if ( associated(mm_xij ) ) deallocate( mm_xij )
      allocate( mm_jan(maxnna), mm_r2ij(maxnna), mm_xij(3,maxnna) )

      if ( allocated(nglist)    ) deallocate( nglist    )
      if ( allocated(nr)        ) deallocate( nr        )
      if ( allocated(nglistxat) ) deallocate( nglistxat )
      allocate( nglistemp(dimvec), nrtemp(3,dimvec), nglistxat(na_mm) )


      nglistemp = 0
      nrtemp    = 0
      nglistxat = 0
      rcoor = rcut_mm + sfc + 2.0_dp

      nnei = 1
      if ( bigcell ) then
        do iat = 1, na_mm ! Look for neighbours of atom ia
          call mm_mneighb( cell, rcoor, na_mm, ramber, iat, 0, nna )
          if ( NNA > MAXNNA ) call die( 'NNA bigger than MAXNNA' )

          do inei = 1, nna
            if ( .not. (abs(mm_r2ij(inei)) > 0.0_dp) ) cycle

            jat = mm_jan(inei)
            if ( iat > jat ) cycle
            if ( (mm_atoms(jat)%aaname == 'HOH') .and. (mm_atoms(jat)%atname /= 'O') ) cycle

            skip_at = .false.
            ! Exclude atoms connected by 1 or 2 covalent bonds
            do kat = 1, nonbondedxat(iat)
              if ( nonbonded(iat,kat) == jat ) skip_at = .true.
            enddo
            if ( skip_at ) cycle

            ! Exclude atoms 1 and 4 in covalently-bonded chain 1-2-3-4.
            do kat = 1, scalexat(iat)
              if ( scaled(iat,kat) == jat ) skip_at = .true.
            enddo
            if ( skip_at ) cycle

            drij(1:3) = mm_xij(1:3,inei) - ramber(1:3,jat) + ramber(1:3,iat)
            call get_pbc_vectors( latt_typ, cell, kcell, drij, nr_indx )

            if ( mm_atoms(jat)%aaname == 'HOH' ) then
              nglistemp(nnei)    = jat
              nglistemp(nnei+1)  = jat +1
              nglistemp(nnei+2)  = jat +2
              nrtemp(1:3,nnei)   = nr_indx(1:3)
              nrtemp(1:3,nnei+1) = nr_indx(1:3)
              nrtemp(1:3,nnei+2) = nr_indx(1:3)

              nnei = nnei +3
            else if ( mm_atoms(jat)%aaname /= 'HOH' ) then
              nglistemp(nnei)  = jat
              nrtemp(1:3,nnei) = nr_indx(1:3)
              nnei = nnei +1
            endif
          enddo ! inei, neighbours

          nglistxat(iat) = nnei -1

          if ( (nnei-1) > dimvec ) then
            write( 6, * ) 'Dimension Neighbour list (required, used)=', &
                          nnei-1, dimvec
            call die( 'MM Energy and Forces: Stopping Program...' )
          endif
        enddo
      else ! bigcell = f
        do iat = 1, na_mm ! Look for neighbours of atom ia
          call mm_mneighb( cell, rcoor, na_mm, ramber, iat, 0, nna )
          if ( NNA > MAXNNA ) call die( 'NNA bigger than MAXNNA.' )

          do inei = 1, nna
            if ( .not. (abs(mm_r2ij(inei)) > 0.0_dp) ) cycle

            jat = mm_jan(inei)
            if ( iat > jat ) cycle
            if ( (mm_atoms(jat)%aaname == 'HOH') .and. (mm_atoms(jat)%atname /= 'O') ) cycle


            drij(1:3) = ramber(1:3,iat) - ramber(1:3,jat)
            call pbc_displ_vector( latt_typ, cell, kcell, drij )
            dr2 = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)


            if ( abs(mm_r2ij(inei) - dr2) < 0.1_dp ) then
              skip_at = .false.
              ! Exclude atoms connected by 1 or 2 covalent bonds
              do kat = 1, nonbondedxat(iat)
                if ( nonbonded(iat,kat) == jat ) skip_at = .true.
              enddo
              if ( skip_at ) cycle

              ! Exclude atoms 1 and 4 in covalently-bonded chain 1-2-3-4.
              do kat = 1, scalexat(iat)
                if ( scaled(iat,kat) == jat ) skip_at = .true.
              enddo
              if ( skip_at ) cycle
            endif

            drij(1:3) = mm_xij(1:3,inei) - ramber(1:3,jat) + ramber(1:3,iat)
            call get_pbc_vectors( latt_typ, cell, kcell, drij, nr_indx )

            if ( mm_atoms(jat)%aaname == 'HOH' ) then
              nglistemp(nnei)    = jat
              nglistemp(nnei+1)  = jat +1
              nglistemp(nnei+2)  = jat +2
              nrtemp(1:3,nnei)   = nr_indx(1:3)
              nrtemp(1:3,nnei+1) = nr_indx(1:3)
              nrtemp(1:3,nnei+2) = nr_indx(1:3)

              nnei = nnei +3
            else if ( mm_atoms(jat)%aaname /= 'HOH' ) then
              nglistemp(nnei)  = jat
              nrtemp(1:3,nnei) = nr_indx(1:3)
              nnei = nnei +1
            endif
          enddo ! inei, neighbours

          nglistxat(iat) = nnei -1

          if ( (nnei-1) > dimvec ) then
            write( 6, * ) 'Dimension Neighbour list (required, used)=', &
                          nnei-1, dimvec
            call die( 'MM Energy and Forces: Stopping Program' )
          endif
        enddo ! iat
      endif

      ! Allocates the true neighbour list.
      allocate( nglist(nnei-1), nr(3,nnei-1) )

      nglist(1:nnei-1) = nglistemp(1:nnei-1)
      nr(1:3,1:nnei-1) = nrtemp(1:3,1:nnei-1)

      deallocate( nrtemp, nglistemp )

      upd_ngb = .false.
    endif

    ! Here we start the proper Energy and forces calculations.
    ! Starts with loop for 1,4-scaled nonbonded interactions
    do iat = 1, na_mm
    do kat = 1, scalexat(iat)
      jat = scaled(iat,kat)
      if ( ((mm_atoms(iat)%aaname == 'GRAP') .and. &
            (mm_atoms(jat)%aaname == 'GRAP')) ) cycle

      drij(1:3) = ramber(1:3,iat) - ramber(1:3,jat)

      call pbc_displ_vector( latt_typ, cell, kcell, drij )
      dr2    = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)
      drdist = sqrt( dr2 )

      if ( drdist > rcut_mm ) cycle
      Rij   = Rm(iat) + Rm(jat)
      Eij   = sqrt( Em(iat) * Em(jat) )
      Rij6  = Rij * Rij * Rij * Rij * Rij * Rij
      dr2_3 = dr2 * dr2 * dr2
      Elj_amber14 = Elj_amber14 + 0.5_dp * Eij * Rij6 / dr2_3 * &
                                          ( Rij6 / dr2_3 - 2.0_dp )
      fel = -6.0_dp * Eij * Rij6 / ( dr2 * dr2_3 ) * &
                            ( Rij6 / dr2_3 - 1.0_dp )

      do icrd = 1, 3
        flj(icrd,iat) = flj(icrd,iat) + drij(icrd) * fel
        flj(icrd,jat) = flj(icrd,jat) - drij(icrd) * fel
        stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + &
                                 stress_fact * drij(1:3) * drij(icrd) * fel
      enddo

      if ( coulombtype == 'ewald' ) then
        call coulomb_ewald_real( mm_atoms(iat)%pc, mm_atoms(jat)%pc, drdist,&
                                 E2, fel, ewald )
      else
        call coulomb_cutoff( mm_atoms(iat)%pc, mm_atoms(jat)%pc, drdist, &
                             dr2, E2, fel )
      endif

      Eelec_amber14 = Eelec_amber14 + unit_f * E2 / 1.2_dp
      fel = fel / 1.2_dp * unit_f

      do icrd = 1, 3
        felec(icrd,iat) = felec(icrd,iat) + drij(icrd) * fel
        felec(icrd,jat) = felec(icrd,jat) - drij(icrd) * fel
        stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + &
                                 stress_fact * drij(1:3) * drij(icrd) * fel
      enddo

    enddo ! kat
    enddo ! iat, scaled nonbonded

    ! Starts with non-scaled nonbonded interactions.

    ! rcut_mm is x0, and +sfc is x1.
    if ( coulombtype == 'ewald' ) then
      cb = - ( sfc / rcut_mm ) * ( erfc(ewald%sqr_alpha * rcut_mm) / rcut_mm &
                + 2.0_dp * ewald%sqr_alpha_pi * &
                exp(- ewald%alpha * rcut_mm *rcut_mm) )
      ca = - 0.5_dp * cb
      cd = ca - erfc(ewald%sqr_alpha * rcut_mm) / rcut_mm
    else
      cb = - sfc / ( rcut_mm * rcut_mm )
      ca = - 0.5_dp * cb
      cd = ca - 1.0_dp / rcut_mm
    endif
    rinn = rcut_mm * rcut_mm
    rout = (rcut_mm + sfc)**2


    n_pointer = 1 ! n_pointer: points the first neb atom of i in the neb list.
    do iat = 1, na_mm
      do kat = n_pointer, nglistxat(iat)
        jat = nglist(kat)

        if ( (mm_atoms(iat)%graph_layer /= mm_atoms(jat)%graph_layer) .or. &
             ((mm_atoms(iat)%graph_layer * mm_atoms(jat)%graph_layer) == 0) )&
           then

          drij(1:3) = ramber(1:3,iat) - ramber(1:3,jat)
          if ( latt_typ == 'D' ) then
            drij(1) = drij(1) + nr(1,kat) * cell(1,1)
            drij(2) = drij(2) + nr(2,kat) * cell(2,2)
            drij(3) = drij(3) + nr(3,kat) * cell(3,3)
          else
            do icrd = 1, 3
              drij(1:3) = drij(1:3) + nr(icrd,kat) * cell(1:3,icrd)
            enddo
          endif
          dr2 = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)

          ! Non-bonded interaction is smoothly decayed from sqrt(rinn)
          ! to sqrt(rout)
          if ( dr2 <= rinn ) then
            drdist = sqrt(dr2)
            Rij = Rm(iat) + Rm(jat)
            Eij = sqrt( Em(iat) * Em(jat) )

            Rij6  = Rij * Rij * Rij * Rij * Rij * Rij
            dr2_3 = dr2 * dr2 * dr2

            Elj_amber = Eij * (Rij6 / dr2_3) * ( (Rij6 / dr2_3) - 2.0_dp ) + &
                        Elj_amber
            fel       = -12.0_dp * Eij * Rij6 * ( Rij6 / dr2_3 - 1.0_dp ) &
                                 / ( dr2 * dr2_3 )

            do icrd = 1, 3
              flj(icrd,iat) = flj(icrd,iat) + drij(icrd) * fel
              flj(icrd,jat) = flj(icrd,jat) - drij(icrd) * fel

              stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                       drij(1:3) * drij(icrd) * fel
            enddo

            if ( coulombtype == 'ewald' ) then
              call coulomb_ewald_real( mm_atoms(iat)%pc, mm_atoms(jat)%pc, &
                                       drdist, E2, fel, ewald )
            else
              call coulomb_cutoff( mm_atoms(iat)%pc, mm_atoms(jat)%pc, drdist,&
                                   dr2, E2, fel )
            endif
            E2 = E2 * unit_f

            ! Add a constant to E2 to make the energy go to
            ! E = q_i*q_j/rcut*ca = q_i*q_j/rcut*(sfc/(2*rcut^2)) when
            ! drdist = rcut and where sfc: smooth-function cut-off.
            Eelec_amber = Eelec_amber + E2 + ( mm_atoms(iat)%pc * &
                          mm_atoms(jat)%pc * unit_f * cd )
            fel = fel * unit_f

            do icrd = 1, 3
              felec(icrd,iat) = felec(icrd,iat) + drij(icrd) * fel
              felec(icrd,jat) = felec(icrd,jat) - drij(icrd) * fel

              stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                       drij(1:3) * drij(icrd) * fel
            enddo

          elseif ( (dr2 > rinn) .and. (dr2 < rout) ) then
            drdist = sqrt(dr2)
            Rij = Rm(iat) + Rm(jat)
            Eij = sqrt( Em(iat) * Em(jat) )

            Rij6  = Rij * Rij * Rij * Rij * Rij * Rij
            dr2_3 = dr2 * dr2 * dr2

            Elj_amber = Eij * (Rij6 / dr2_3) * ( (Rij6 / dr2_3) - 2.0_dp ) &
                      + Elj_amber
            fel       = -12.0_dp * Eij * Rij6 * ( Rij6 / dr2_3 - 1.0_dp ) &
                                 / ( dr2 * dr2_3 )
            do icrd = 1, 3
              flj(icrd,iat) = flj(icrd,iat) + drij(icrd) * fel
              flj(icrd,jat) = flj(icrd,jat) - drij(icrd) * fel

              stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                        drij(1:3) * drij(icrd) * fel
            enddo

            ! Add coulomb and cut-off switching.
            E2  = mm_atoms(iat)%pc * mm_atoms(jat)%pc * unit_f
            fac = ( drdist - rcut_mm ) / sfc

            Eelec_amber = Eelec_amber + E2 * ( ca * fac * fac + cb * fac + ca )
            fel = ( E2 / drdist ) * ( 2.0_dp * ca * fac + cb ) / sfc

            do icrd = 1, 3
              felec(icrd,iat) = felec(icrd,iat) + drij(icrd) * fel
              felec(icrd,jat) = felec(icrd,jat) - drij(icrd) * fel

              stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                       drij(1:3) * drij(icrd) * fel
            enddo

          endif
        endif ! Graphite layers.
      enddo ! Loop over list, k

      n_pointer = nglistxat(iat) + 1
    enddo ! Loop over atoms iat

    if ( coulombtype == 'ewald' ) then
      ! Part of Ewald sum in the reciprocal space.
      ! Erecip_amber = 0.0_dp

      ! Calculate structure factors for all MM atoms.
      const4 = 0.5_dp / cell_v
      const6 = 1.0_dp / ( 4.0_dp * ewald%alpha )

      if (latt_typ == 'D') then
        ! Reciprocal-space sum
        do n1 = -n1max, n1max
        do n2 = -n2max, n2max
        do n3 = -n3max, n3max
          if ( (n1 == 0) .and. (n2 == 0) .and. (n3 == 0) ) cycle

          krecip(1) = n1 * pi2 * kcell(1,1)
          krecip(2) = n2 * pi2 * kcell(2,2)
          krecip(3) = n3 * pi2 * kcell(3,3)
          kmod2     = krecip(1) * krecip(1) + krecip(2) * krecip(2) + &
                      krecip(3) * krecip(3)

          if ( kmod2 > ewald%kcut_sq ) cycle

          S_real = 0.0_dp
          S_imag = 0.0_dp

          do iat = 1, na_mm
            kr = krecip(1) * ramber(1,iat) + krecip(2) * ramber(2,iat) + &
                 krecip(3) * ramber(3,iat)
            pc_cos_kr(iat) = mm_atoms(iat)%pc * cos(kr)
            pc_sin_kr(iat) = mm_atoms(iat)%pc * sin(kr)
            S_real = S_real + pc_cos_kr(iat)
            S_imag = S_imag + pc_sin_kr(iat)
          enddo

          const3 = ( const2 / kmod2 ) * unit_f * &
                    exp( - kmod2 / (4.0_dp * ewald%alpha) )
          SS = 0.5_dp * const3 * ( S_real * S_real + S_imag * S_imag )

          Eelec_amber  = Eelec_amber  + SS
          !Erecip_amber = Erecip_amber + SS
          do iat = 1, na_mm
            De = const3 * ( pc_cos_kr(iat) * S_imag - pc_sin_kr(iat) * S_real )
            felec(1:3,iat) = felec(1:3,iat) + De * krecip(1:3)
          enddo

          const5 = const3 * const4
          De2    = const5 * ( S_imag * S_imag + S_real * S_real )
          const7 = 2.0_dp * ( 1.0_dp + kmod2 * const6 ) / kmod2

          do ivec = 1, 3
          do icrd = 1, 3
            kronij = real( int( ( (icrd+ivec) - abs(icrd-ivec) ) / &
                                ( (icrd+ivec) + abs(icrd-ivec) ) ), kind = dp )
            stress_amber(icrd,ivec) = stress_amber(icrd,ivec) + De2 * ( kronij &
                                    - const7 * krecip(icrd) * krecip(ivec) )
          enddo
          enddo
        enddo
        enddo
        enddo
      else
        ! Reciprocal-space sum
        ! Fix for "blindspots" in reciprocal space for hexagonal cells.
        do ni = 1, 3
          k_neigh(ni, :) = 2 * pi2 * (/ n1max * kcell(1, ni), &
                           n2max * kcell(2, ni), n3max * kcell(3,ni) /)
        enddo
        do ni = 4, 6
          k_neigh(ni, :) = -k_neigh(ni-3, :)
        enddo

        do n1 = -n1max, n1max
        do n2 = -n2max, n2max
        do n3 = -n3max, n3max
          if ( (n1 == 0) .and. (n2 == 0) .and. (n3 == 0) ) cycle

          krecip(1) = pi2 * ( n1 * kcell(1,1) + n2 * kcell(1,2) &
                              + n3 * kcell(1,3) )
          krecip(2) = pi2 * ( n1 * kcell(2,1) + n2 * kcell(2,2) &
                              + n3 * kcell(2,3) )
          krecip(3) = pi2 * ( n1 * kcell(3,1) + n2 * kcell(3,2) &
                              + n3 * kcell(3,3) )
          kmod2     = krecip(1) * krecip(1) + krecip(2) * krecip(2) + &
                      krecip(3) * krecip(3)

          if ( kmod2 > ewald%kcut_sq ) then
            ! Fix for "blindspots" in reciprocal space for hexagonal cells.
            skip_k = .true.
            do ni = 1, 6
              krecip_displaced = krecip - k_neigh(ni, :)
              kmod2 = krecip_displaced(1) * krecip_displaced(1) + &
                      krecip_displaced(2) * krecip_displaced(2) + &
                      krecip_displaced(3) * krecip_displaced(3)
              if ( kmod2 > ewald%kcut_sq ) cycle
              krecip = krecip_displaced
              skip_k = .false.
              exit
            enddo

            if ( skip_k ) cycle
          endif

          S_real = 0.0_dp
          S_imag = 0.0_dp

          do iat = 1, na_mm
            kr = krecip(1) * ramber(1,iat) + krecip(2) * ramber(2,iat) + &
                 krecip(3) * ramber(3,iat)
            pc_cos_kr(iat) = mm_atoms(iat)%pc * cos(kr)
            pc_sin_kr(iat) = mm_atoms(iat)%pc * sin(kr)
            S_real = S_real + pc_cos_kr(iat)
            S_imag = S_imag + pc_sin_kr(iat)
          enddo

          const3 = ( const2 / kmod2 ) * unit_f * &
                    exp( - kmod2 / (4.0_dp * ewald%alpha) )
          SS = 0.5_dp * const3 * ( S_real * S_real + S_imag * S_imag )

          Eelec_amber  = Eelec_amber  + SS
          !Erecip_amber = Erecip_amber + SS
          do iat = 1, na_mm
            De = const3 * ( pc_cos_kr(iat) * S_imag - pc_sin_kr(iat) * S_real )
            felec(1:3,iat) = felec(1:3,iat) + De * krecip(1:3)
          enddo

          const5 = const3 * const4
          De2    = const5 * ( S_imag * S_imag + S_real * S_real )
          const7 = 2.0_dp * ( 1.0_dp + kmod2 * const6 ) / kmod2

          do ivec = 1, 3
          do icrd = 1, 3
            kronij = real( int( ( (icrd+ivec) - abs(icrd-ivec) ) / &
                                ( (icrd+ivec) + abs(icrd-ivec) ) ), kind = dp )
            stress_amber(icrd,ivec) = stress_amber(icrd,ivec) + De2 * ( kronij &
                                    - const7 * krecip(icrd) * krecip(ivec) )
          enddo
          enddo
        enddo
        enddo
        enddo
      endif

      ! From the reciprocal sum, substract contributions to the electrostatic
      ! interactions already contained in the many-body energy terms (1-,2-,
      ! 3-body and 1,4-scaled interactions).
      do iat = 1, na_mm
        do kat = 1, nonbondedxat(iat)
          ! Consider atoms connected by 1 or 2 covalent bonds
          jat = nonbonded(iat,kat)

          ! Check that the atom i and j are not connected by a dihedral
          skip_at = .false.
          do mat = 1, mm_connectivity(iat)%ndihe_e
            if ( mm_connectivity(iat)%dihe_at(mat,3) == jat ) then
              skip_at = .true.
              exit
            endif
          enddo
          if ( skip_at ) cycle

          drij(1:3) = ramber(1:3,iat) - ramber(1:3,jat)

          call pbc_displ_vector( latt_typ, cell, kcell, drij )
          dr2    = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)
          drdist = sqrt( dr2 )

          ! Substract the real part of the energy added in the reciprocal-space
          ! sum. The bond and angle terms already have the contribution from
          ! electrostatics.
          Eelec_amber  = Eelec_amber - ( mm_atoms(iat)%pc * mm_atoms(jat)%pc / drdist ) * &
                         erf( ewald%sqr_alpha * drdist ) * unit_f
          !Erecip_amber = Erecip_amber - ( mm_atoms(iat)%pc * mm_atoms(jat)%pc / drdist )* &
          !               erf( ewald%sqr_alpha * drdist ) * unit_f

          De = ( mm_atoms(iat)%pc * mm_atoms(jat)%pc / dr2 ) * unit_f * &
               ( erf( ewald%sqr_alpha * drdist ) / drdist - 2.0_dp * &
               ewald%sqr_alpha_pi * exp( -ewald%alpha * dr2 ) )

          do icrd = 1, 3
            felec(icrd,iat) = felec(icrd,iat) + drij(icrd) * De
            felec(icrd,jat) = felec(icrd,jat) - drij(icrd) * De

            stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                     drij(1:3) * drij(icrd) * De
          enddo
        enddo

        ! consider atoms 1 and 4 in covalently-bonded chain 1-2-3-4
        do kat = 1, scalexat(iat)
          jat       = scaled(iat,kat)
          drij(1:3) = ramber(1:3,iat) - ramber(1:3,jat)

          call pbc_displ_vector( latt_typ, cell, kcell, drij )
          dr2    = drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3)
          drdist = sqrt( dr2 )

          ! Substract the real part of the energy added in the reciprocal-space
          ! sum. The bond and angle terms have already the contribution from
          ! electrostatics.
          Eelec_amber14 = Eelec_amber14 - &
                      ( mm_atoms(iat)%pc * mm_atoms(jat)%pc / drdist ) * &
                      erf( ewald%sqr_alpha * drdist ) * 0.2_dp / 1.2_dp * unit_f
          fel = ( mm_atoms(iat)%pc * mm_atoms(jat)%pc / dr2 ) * &
                ( erf( ewald%sqr_alpha * drdist) / drdist - &
                2.0_dp * ewald%sqr_alpha_pi * exp( -ewald%alpha * dr2 ) ) &
                * 0.2_dp / 1.2_dp * unit_f

          do icrd = 1, 3
            felec(icrd,iat) = felec(icrd,iat) + drij(icrd) * fel
            felec(icrd,jat) = felec(icrd,jat) - drij(icrd) * fel

            stress_amber(1:3,icrd) = stress_amber(1:3,icrd) + stress_fact * &
                                     drij(1:3) * drij(icrd) * fel
          enddo
        enddo
      enddo

      ! Subtract the self-energy.
      do iat = 1, na_mm
        Eelec_amber  = Eelec_amber  - &
               ewald%sqr_alpha_pi * mm_atoms(iat)%pc * mm_atoms(iat)%pc * unit_f
      enddo
    endif ! Ewald

    deallocate( pc_cos_kr, pc_sin_kr )
  end subroutine amber_nonbonded

  subroutine coulomb_ewald_real( pc_i, pc_j, r_ij, E2, fel, ewald )
    !! Calculates point charge interactions using the Ewald method; this is
    !! only the real part.
    use coulomb_m, only : ewald_data_t

    implicit none
    real(dp), intent(in)  :: pc_i
      !! Partial charge of atom 1.
    real(dp), intent(in)  :: pc_j
      !! Partial charge of atom 2.
    real(dp), intent(in)  :: r_ij
      !! Distance between atoms 1 and 2.
    type(ewald_data_t), intent(in) :: ewald
      !! Ewald summation parameters.
    real(dp), intent(out) :: fel
      !! Electrostatic force component.
    real(dp), intent(out) :: E2
      !! Electrostatic contribution to energy.

    real(dp) :: r_ij2

    r_ij2 = r_ij * r_ij
    E2    =  ( pc_i * pc_j ) / r_ij  *   erfc( ewald%sqr_alpha * r_ij )
    fel   = -( pc_i * pc_j ) / r_ij2 * ( erfc( ewald%sqr_alpha * r_ij ) / r_ij &
            + 2.0_dp * ewald%sqr_alpha_pi * exp( - ewald%alpha * r_ij2 ) )
  end subroutine coulomb_ewald_real

  subroutine coulomb_cutoff( pc_i, pc_j, r_ij, r_ij2, E2, fel )
    !! Calculates point charge interactions using a real-space cut-off distance.

    implicit none
    real(dp), intent(in)  :: pc_i
      !! Partial charge of atom 1.
    real(dp), intent(in)  :: pc_j
      !! Partial charge of atom 2.
    real(dp), intent(in)  :: r_ij
      !! Distance between atoms 1 and 2.
    real(dp), intent(in)  :: r_ij2
      !! Squared distance between atoms 1 and 2.
    real(dp), intent(out) :: fel
      !! Electrostatic force component.
    real(dp), intent(out) :: E2
      !! Electrostatic contribution to energy.

    E2  = pc_i * pc_j / r_ij
    fel = - E2 / r_ij2
  end subroutine coulomb_cutoff

  subroutine waters( na_qm, na_mm, natot, rclas, masst, mm_atoms, &
                     ewat, fwat, cell, kcell, latt_typ )
    !! Calculates the restraint potential for water molecules. There are
    !! no stress contributions here.
    use functions, only : norm_v2
    use qmmm_pbc , only : pbc_displ_vector
    use units    , only : Ang

    implicit none
    integer         , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)  :: natot
      !! Total number of atoms (QM+MM).
    real(dp)        , intent(in)  :: rclas(3,natot)
      !! Positions of the MM atoms.
    type(mm_atom_t) , intent(in)  :: mm_atoms(na_mm)
      !! Information for MM atoms.
    real(dp)        , intent(in)  :: masst(natot)
      !! Atomic masses.
    real(dp)        , intent(in)  :: cell(3,3)
      !! Periodic cell vectors.
    real(dp)        , intent(in)  :: kcell(3,3)
      !! Reciprocal periodic cell vectors.
    character       , intent(in)  :: latt_typ
    real(dp)        , intent(out) :: ewat
      !! Energy contribution of the water restraint potential.
    real(dp)        , intent(out) :: fwat(3,na_mm)
      !! Forces contribution of the water restraint potential.

    integer  :: watlistnum, iat, iwat
    real(dp) :: rwat, m_center(3), rij, kte, dist2, mdist, drw(3), mtot

    integer , allocatable :: watlist(:)
    real(dp), allocatable :: ramber(:,:)

    allocate( watlist(2000) )
    allocate( ramber(3,natot) )

    kte       = 200.0_dp
    dist2     = 0.0_dp
    mdist     = 0.0_dp
    ewat      = 0.0_dp
    fwat(:,:) = 0.0_dp

    ! Calculates the center of mass of the system.
    ramber(1:3,1:natot) = rclas(1:3,1:natot) / Ang
    m_center = 0.0_dp
    mtot     = 0.0_dp
    do iat = 1, natot
      m_center(1:3) = m_center(1:3) + masst(iat) * ramber(1:3,iat)
      mtot          = mtot + masst(iat)
    enddo
    m_center = m_center / mtot

    do iat = 1, natot
      dist2 = ( ramber(1,iat) - m_center(1) )**2 + &
              ( ramber(2,iat) - m_center(2) )**2 + &
              ( ramber(3,iat) - m_center(3) )**2
      if ( dist2 > mdist ) mdist = dist2
    enddo

   rwat = sqrt(mdist) - 2.5_dp
   write( 6, '(a,2x,f8.4)' ) 'Water Cutoff Radius:', rwat

    ! Calculate the distance matrix, using only the MM atoms.
    ramber = 0.0_dp
    ramber(1:3,1:na_mm) = rclas(1:3,na_qm+1:natot) / Ang

    iwat = 1
    do iat=1, na_mm
      if ( .not. ((mm_atoms(iat)%aaname == 'HOH') .and. (mm_atoms(iat)%atname == 'O')) ) cycle

      drw(1:3) = ramber(1:3,iat) - m_center(1:3)
      call pbc_displ_vector( latt_typ, cell, kcell, drw )
      rij = norm_v2( drw )

      ! We check if the water is in the buffer zone.
      if ( rij < rwat ) cycle

      watlist(iwat) = iat
      iwat = iwat +1
    enddo

    watlistnum = iwat -1
    ewat = 0.0_dp
    fwat = 0.0_dp

    ! Calculation of energy and forces contribution.
    do iwat = 1, watlistnum
      iat = watlist(iwat)

      drw(1:3) = ramber(1:3,iat) - m_center(1:3)
      call pbc_displ_vector( latt_typ, cell, kcell, drw )
      rij = norm_v2( drw )
      if ( rij < rwat ) cycle

      ewat = ewat + kte * ( ( rij - rwat )**2 )
      drw(1:3) = 2.0_dp * kte * (rij - rwat) * drw(1:3) / rij

      fwat(1:3,iat) = fwat(1:3,iat) - drw(1:3)
    enddo

    deallocate( watlist, ramber )
  end subroutine waters

end module mm_ene_frc