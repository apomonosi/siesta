program SIESTA_QMMM

!! Program for QM-MM calculations. It uses the SIESTA code to treat at
!! DFT leVel the QM subsystem, and an in-house implementation of the
!! Amber99 forcefield parametrization to treat the MM subsystem.
!!
!! Original idea by D. Scherlis, D. Estrin and P. Ordejon. 2001.
!! Original QM-MM interfase with Siesta by D. Scherlis and P. Ordejon. 2001.
!! Original subroutines by D. Scherlis. 2001/2002.
!!
!! Solvent implementation using the Amber99 force field parametrization
!! by A. Crespo and M. Marti. subroutines by M. Marti. 2003.
!! Further modifications of QM-MM Siesta by A. Crespo and P. Ordejon. 2004.
!! DeVeloping of SIESTA-QMMM: programf or QM-MM calculations using the SIESTA
!! as SUBROUTINE by A. Crespo and P. Ordejon. 2004.
!! Crespo, May 2005.
!!
!! C. F. Sanz-Navarro (2008-2012)
!!
!! Code re-habilitation and modernization by Ernane de Freitas Martins and
!! Federico Pedron, 2021-2023.

  !! Modules that come from SIESTA.
  use alloc               , only : re_alloc, de_alloc
  use atomlist            , only : amass, iza
  use fdf                 , only : fdf_shutdown, fdf_get, fdf_deprecated
  use files               , only : slabel
  use qmmm_fsiesta        , only : siesta_units, siesta_quit, create_fifos, &
                                   siesta_send_coordinates, &
                                   siesta_receive_forces
  use m_cite              , only : init_citation, reset_citations
  use m_energies          , only : Ekinion, Etot, FreeE
  use m_fixed             , only : fixed
  use m_forces            , only : fa, cfa, ntcon
  use m_intramol_pressure , only : remove_intramol_pressure
  use m_kinetic           , only : tempion
  use m_steps             , only : inicoor, fincoor, istp, final
  use m_stress            , only : cstress, kin_stress, mstress, tstress, stress
  use m_timestamp         , only : timestamp
  use parallel            , only : IOnode
  use precision           , only : dp
  use siesta_geom         , only : na_u, ucell, vcell, xa, va, isa, &
                                   volume_of_some_cell
  use sys                 , only : die, message
  use siesta_cml          , only : cml_p
  use siesta_options      , only : bulkm, charnet, dt, sname, ftol, strtol, &
                                   usesavexv, tp, tempinit, writef, idyn,   &
                                   RemoveIntraMolecularPressure, varcel, dx
  use units               , only : amu, Ang, eV
  use zmatrix             , only : lUseZmatrix, iofaZmat, &
                                   CartesianForce_to_ZmatForce
#ifdef MPI
  use parallel            , only : Node, Nodes
  use mpi_siesta          , only : MPI_Comm_World, MPI_Init, MPI_Comm_Rank, &
                                   MPI_Comm_Size, MPI_Finalize
#endif

  ! Modules from this QMMM implementation.
  use assign_qmmm         , only : assign_atoms
  use center_all          , only : centermol
  use electric_field      , only : elecfield
  use coulomb_m           , only : coulombtype, ewald
  use ioxvrestr_m         , only : ioxvrestr
  use linkatoms           , only : linkatom, link_atoms, numlink, elink, &
                                   num_resonances, resonance
  use linkatoms           , only : link1, link2, link3, get_link_ff_parameters,&
                                   get_numlink, allocate_link_arrays, &
                                   deallocate_link_arrays
  use m_libsiesta_init    , only : libsiesta_init
  use mm_assign_m         , only : mm_assign, mm_atom_to_species
  use mm_ene_frc          , only : mm_ene_fce, mm_dealloc
  use mm_topology         , only : fftopology, atom_connect_t, mm_atom_t, &
                                   qm_atom_t
  use qm_read             , only : read_qm, readcrd
  use qmmm_constraints    , only : qmmm_fixed1, qmmm_fixed2, &
                                   add_siesta_qmmm_ntcon

  use qmmm_files_m        , only : qmmm_files
  use qmmm_ioxv_m         , only : qmmm_ioxv
  use qmmm_list_block     , only : qmmm_lst_blk
  use qmmm_lj             , only : ljef
  use qmmm_mm_neighbour   , only : update_qmmm_mneighb
  use qmmm_pbc            , only : get_lattice_type
  use qmmm_reinit_m       , only : qmmm_reinit
  use qmmm_move_m         , only : qmmm_move
  use qmmm_timers         , only : qmmm_timer
  use qmmm_write          , only : write_energy, write_siesta_input, &
                                   write_react_crd, write_pdbcrd, write_restr, &
                                   qmmm_write_forces, qmmm_write_stress_pressure
  use reinsert_atoms_m    , only : reinserting_atoms_in_box
  use restr_subs          , only : restr_read, restr_calc, restr_update, &
                                   restr_steps, hasRestraints, nrestr, restr_r
  use siesta_qmmm_options , only : mneigh_freq, center_qm_system, opt_idyn,    &
                                   siesta_qmmm_usesavexv, fc_idyn, therm_idyn, &
                                   read_siesta_qmmm_options, mmsteps, wricoord,&
                                   launch_siesta_flag
  implicit none

  ! General variables.
  integer  :: i, ia, istep, ix, j, jx, nspec, step
  logical  :: foundxv, relaxd, has_qm, has_mm
  real(dp) :: alat

  ! These should probably have their initialisations elsewhere.
  integer :: restart_istep   = 0

  real(dp)  :: volcel
  external  :: memory, volcel

  ! MM-related general variables
  type(fftopology) :: mm_top
  type(atom_connect_t), allocatable :: mm_connectivity(:)
  type(mm_atom_t)     , allocatable :: mm_atoms(:)
  type(qm_atom_t)     , allocatable :: qm_atoms(:)

  character(len=2) , allocatable :: atsym(:)

  integer , pointer :: blocklist(:), qmstep(:), atxres(:)
  real(dp), pointer :: fce_amber(:,:)

  integer  :: na_qm, na_mm, natot, nfree, imm, nroaa, ncount, int_dummy
  logical  :: actualiz, water, foundcrd, foundvat, foundva, writeipl, debug, &
              last_restr
  real(dp) :: Etots, Etot_amber, sfc, kcal, rcorteqm, rcortemm, &
              Elj, Eelec, stress_amber(3,3), fovermp, kin, tempe, &
              mscell(3,3), scell(3,3)

  ! Restrained optimization variables
  integer  :: isteprestr, istep_mv

  ! QM-MM cut-off variables
  real(dp) :: rcorteqmmm, radbloqmmm

  ! Solvent-related external variables
  external :: qmmm_prversion

  ! Periodic boundary conditions
  character         :: lattice_type
  character(len=13) :: string_path

  ! Format strings
  character(len=50) :: fmt1, fmtIERR, fmt2

#ifdef MPI
  integer  :: MPIerror  ! Return error code in MPI routines

  ! Initialise MPI and set processor number
  call MPI_Init( MPIerror )
  call MPI_Comm_Rank( MPI_Comm_World, Node , MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#endif

  fmt1    = '(/,"siesta-qmmm: ",71(1h*))'
  fmt2    = '("siesta-qmmm: ",71(1h*))'
  fmtIERR = '("siesta-qmmm:                  INPUT ERROR")'

  ! Initialise some variables
  relaxd   = .false.; has_qm = .true.; has_mm = .true.
  ntcon    = 0
  tempion  = 0.0_dp; tp    = 0.0_dp; vcell = 0.0_dp;
  strtol   = 0.0_dp; bulkm = 0.0_dp
  lattice_type = 'D'

  ! So far only serial version of siesta-qmmm
  IOnode = .true.

  ! Print version information
  call qmmm_prversion
  write(6,'(/,a,i4,a)') 'siesta-qmmm: Running in serial mode'

  ! Start time counters
  call timestamp('Start of run')
  call qmmm_timer( 'all', 0 )
  call qmmm_timer( 'Siesta-qmmm', 1 )
  call qmmm_timer( 'Setup', 1 )

  ! Factors of conversion to internal units
  kcal = 1.602177E-19_dp * 6.02214E23_dp / 4184.0_dp

  call qmmm_reinit( sname )
  call qmmm_files%set_filenames( slabel )
  call init_citation( slabel )

  ! Read the number of QM atoms. Some of this options are re-read in SIESTA.
  na_qm = fdf_get( 'NumberOfAtoms', 0 )
  if ( na_qm == 0 ) then
    write(6,'(/a)') 'siesta-qmmm: Running with no QM atoms.'
    has_qm = .false.
  endif

  ! Read the number of species
  nspec = fdf_get( 'NumberOfSpecies', 0 )
  allocate( atsym(1) )
  if ( has_qm ) then
    if ( nspec == 0) &
      call die("siesta-qmmm: You must specify the number of QM species.")

    ! Read QM coordinates and labels
    write(6,*)
    write(6,"('read: ',73(1h*))")
    deallocate( atsym )
    allocate( atsym(nspec) )
    allocate( qm_atoms(na_qm) )

    write(6,*)
    write(6,"('read: ',73(1h*))")

    call read_qm( na_qm, nspec, qm_atoms, atsym )
    call memory( 'A', 'D', na_qm, 'siesta-qmmm' )
  else
    allocate( qm_atoms(1) )
  endif


  ! Read the number of MM atoms
  call fdf_deprecated( 'NumberOfSolventAtoms', 'NumberofMMAtoms' )
  na_mm = fdf_get( 'NumberOfSolventAtoms', 0 )
  na_mm = fdf_get( 'NumberOfMMAtoms', na_mm )
  if ( na_mm == 0 ) then
    write(6,'(/a)') 'siesta-qmmm: Running with no MM atoms'
    has_mm = .false.
  endif

  call redcel( alat, ucell, scell, mscell )
  lattice_type = get_lattice_type( ucell )

  ! Read and assign Solvent variables
  allocate( mm_connectivity(na_mm) )
  allocate( mm_atoms(na_mm) )

  nullify( fce_amber )
  call re_alloc( fce_amber  , 1, 3, 1, na_mm , 'fce_amber'  , 'siesta_qmmm' )

  nullify( atxres )
  if ( has_mm ) then
    call re_alloc( atxres, 1, 200000, 'atxres', 'siesta_qmmm')
    call mm_assign( na_qm, na_mm, nroaa, mm_atoms, mm_connectivity, &
                    qm_atoms, rcorteqmmm, rcorteqm, rcortemm, sfc,  &
                    radbloqmmm, atxres, ucell, lattice_type, mm_top )

    ! Substract the MM charge from the total charge in order to calculate the
    ! charge on the QM region.
    do i=1,na_mm
      charnet = charnet - mm_atoms(i)%pc
    enddo

    if ( abs(charnet) < 0.001 ) charnet = 0.0_dp
  endif ! mm

  debug = fdf_get( 'SolventDebug', .false. )
  linkatom = .false.

  call allocate_link_arrays() ! Allocates arrays as size-0.
  if ( has_qm .and. (.not. has_mm) ) then
    allocate( resonance(num_resonances) )
    resonance(1)%path   = './'
    resonance(1)%weight = 1.0_dp

  else if ( has_qm .and. has_mm ) then
    call get_numlink( na_qm, na_mm, mm_atoms, qm_atoms, &
                      mm_top, debug, ucell, lattice_type )
    allocate( resonance(num_resonances) )
    resonance(1)%path   = './'
    resonance(1)%weight = 1.0_dp

    if ( numlink > 0 ) then
      linkatom = .true.

      do i = 1, num_resonances
        allocate( resonance(i)%atype(numlink) )
        allocate( resonance(i)%isa(numlink) )
        allocate( resonance(i)%iza(numlink) )
      enddo
      if ( num_resonances > 1 ) then
        do i = 2, num_resonances
          write( string_path, '(a,i2.2,a)' ) 'RESONANCE-', i, '/'
          resonance(i)%path=trim(string_path)

          ! IS THIS A GOOD IDEA??
          call system( 'mkdir -p '//trim(resonance(i)%path) )
          call system( 'cp *.psf '//trim(resonance(i)%path) )
        enddo
      endif

      do i = 1, num_resonances
      do j = 1, numlink
        resonance(i)%iza(j) = 0
      enddo
      enddo

      ! Allocate LinkAtom variables
      call deallocate_link_arrays()
      call allocate_link_arrays()

      ! Sets parameters for all LinkAtoms
      call link1( mm_connectivity, na_qm, na_mm, mm_atoms, qm_atoms, &
                  mm_top, nspec, debug, ucell, lattice_type )
      call get_link_ff_parameters( na_mm, mm_atoms, mm_top )
      ! Sets to zero Em for HL and CQM
      do i = 1, numlink
      do j = 1, link_atoms(i)%num_qm
        qm_atoms(link_atoms(i)%qm(j))%lj_Em = 0.0_dp
      enddo
      enddo
    endif ! numlink > 0
  endif   ! has_qm .and. has_mm

  ! Allocation of solvent+solute variables
  natot = na_mm + na_qm

  ! More initializations.
  Elj   = 0.0_dp; Elink = 0.0_dp; Etots = 0.0_dp;
  Etot  = 0.0_dp; Eelec = 0.0_dp;
  nfree = natot ; step = 0; foundvat = .false.

  writeipl      = fdf_get( 'WriteIniParLas', .false. )

  ! Read simulation data
  call read_siesta_qmmm_options( )

  ! Assignation of masses and species
  call assign_atoms( na_qm, na_mm, mm_atoms, qm_atoms )

  ! Center system  to make the QM system to coincide with SIESTA QM grid.
  if ( has_qm .and. center_qm_system ) &
    call centermol( na_qm, na_mm, qm_atoms, mm_atoms, ucell )

  ! Read fixed atom constraints
  nullify( blocklist )
  call re_alloc( blocklist, 1, natot, 'blocklist', 'siesta_qmmm')
  blocklist = 0
  call qmmm_fixed1( na_qm, na_mm, natot, nroaa, blocklist, mm_atoms, water )

  call reinserting_atoms_in_box( lattice_type, natot, na_qm, na_mm, &
                                 qm_atoms, mm_atoms, blocklist, ucell )

  ! Sets LinkAtoms positions
  if ( linkatom ) then
    call link3( na_qm, na_mm, qm_atoms, mm_atoms, mm_top, ucell, lattice_type )
    ! Sets to zero pc and Em for MMLink
    do i = 1, numlink
      do j = 1, link_atoms(i)%num_mm
        link_atoms(i)%pc(j) = mm_atoms(link_atoms(i)%mm(j))%pc
        link_atoms(i)%Em(j) = mm_atoms(link_atoms(i)%mm(j))%lj_Em
        mm_atoms(link_atoms(i)%mm(j))%lj_Em = 0.0_dp
      enddo

      mm_atoms(link_atoms(i)%mm(1))%pc = 0.0_dp
      do j = 2, link_atoms(i)%num_mm
        mm_atoms(link_atoms(i)%mm(j))%pc = mm_atoms(link_atoms(i)%mm(j))%pc &
                                 + link_atoms(i)%pc(1) &
                                 / real(link_atoms(i)%num_mm-1,kind=dp)
      enddo
    enddo
  endif

  ! Read the Coulomb type cut-off (non-periodic) or ewald (periodic) for the MM
  ! potential and QM/MM interaction calculation.
  coulombtype = fdf_get( 'Coulombtype', 'ewald' )
  call ewald%set_coef( rcortemm )

  ! Initialize missing stuff needed for Siesta-embedded routines.
  call libsiesta_init( na_qm, na_mm, numlink, qm_atoms, mm_atoms, link_atoms )
  call mm_atom_to_species( na_qm, na_mm, numlink, nspec, mm_atoms, isa )

  if ( natot /= na_u-numlink ) then
    write(*,*)'natot= '  , natot
    write(*,*)'numlink= ', numlink
    write(*,*)'na_u= '   , na_u
    call die( 'natot /= na_u - numlink' )
  endif

  natot = na_u
  na_u  = na_qm + na_mm
  usesavexv = siesta_qmmm_usesavexv
  va = 0.0_dp

  if (has_qm) then
    do i = 1, na_qm
      amass(i) = qm_atoms(i)%mass
      iza(i)   = qm_atoms(i)%z
    enddo
  endif
  if (has_mm) then
    do i = 1, na_mm
      amass(i+na_qm) = mm_atoms(i)%mass
      iza(i+na_qm)   = mm_atoms(i)%z
    enddo
  endif

  foundxv  = .false.
  foundvat = .false.
  ! Read cell shape and atomic positions from a former run.
  if ( siesta_qmmm_usesavexv ) then
    call qmmm_ioxv( qmmm_files%slabel, 'r', int_dummy, restart_istep, na_u, &
                    ucell, vcell, xa, va, foundxv, foundvat )
    lattice_type = get_lattice_type(ucell)
  endif

  ! Read a crd file from a former run
  call readcrd( na_qm, na_mm, amass, xa, va, foundcrd, foundvat )

  if ( has_qm .and. has_mm ) then ! Calculate Rcut & block list QM-MM
    call qmmm_lst_blk( na_qm, na_mm, nroaa, atxres, xa, rcorteqmmm, &
                       radbloqmmm, mm_atoms, qmmm_files%slabel )
  endif

  if ( (.not. opt_idyn ) .and. ( .not. fc_idyn ) .and. (.not. foundvat) ) then
    ! Build initial velocities according to Maxwell-Bolzmann distribution
    call vmb( na_u, tempinit, amass, xa, isa, va )

    ! Impose constraints to atomic movements by chAnging forces
    call qmmm_fixed2( na_qm, na_mm, na_u, nfree, blocklist, mm_atoms, fa, &
                      cfa, va )
    call add_siesta_qmmm_ntcon( ucell, na_qm, na_mm, isa, amass, xa, blocklist,&
                                mm_atoms, ntcon )

    ncount = 0
    if ( has_qm .and. has_mm ) then
      ! Convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
      fovermp = 0.009579038_dp * Ang * Ang / eV

      ! Calculate kinetic energy and temperature
      kin    = 0.0_dp
      do ia = 1, na_u
        if ( sqrt( va(1,ia)**2 + va(2,ia)**2 + va(3,ia)**2 ) == 0.0_dp ) cycle
        ncount = ncount+1
        do i = 1, 3
          kin = kin + 0.5_dp * amass(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

      tempe = 2.0_dp * kin / (3.0_dp * ncount - 3.0_dp) / 8.617d-5 / eV
      do i = 1, na_u
        if ( sqrt( va(1,i)**2 + va(2,i)**2 + va(3,i)**2 ) < 1.0d-14 ) cycle
        do ix = 1, 3
          va(ix,i) = va(ix,i) * sqrt(tempinit/tempe)
        enddo
      enddo

      if ( (3*ncount) /= (3*natot-ntcon) ) then
        call message('WARNING', '3*ncount/=3*natot-ntcon.')
        call message('WARNING', &
          'Unless you are running MD at 0K, review your constraints.' )
      endif
    endif
  endif

  ! Dump initial coordinates after qmmm_ioxv for debuging in PDB format
  if ( debug ) then
    write(6,'(/a)') 'Dump of initial coordinates for debuging: '
    if ( has_qm ) write(6,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)') &
         ( 'ATOM', i, atsym(qm_atoms(i)%spec), 'STO ', 'A', i, xa(1,i) / Ang, &
           xa(2,i) / Ang, xa(3,i) / Ang, i = 1, na_qm )
    if ( has_mm ) write(6,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)') &
         ( 'ATOM', na_qm+i, mm_atoms(i)%atname, mm_atoms(i)%aaname, 'A', &
           mm_atoms(i)%aanum, xa(1,na_qm+i) / Ang, xa(2,na_qm+i) / Ang,  &
           xa(3,na_qm+i) / Ang, i = 1, na_mm )
    write(6,'(A3)') 'END'
  endif

  ! Sets LinkAtoms positions for new positions
  if ( linkatom .and. ( foundxv .or. foundcrd ) ) &
    call link3( na_qm, na_mm, qm_atoms, mm_atoms, mm_top, ucell, lattice_type )

  ! Dump initial coordinates to output
  if ( writeipl ) then
    write(6,*)
    write(6,"('siesta-qmmm: ',71(1h=))")
    if ( has_qm ) then
      write(6,'(a)') 'siesta-qmmm: Atomic coordinates (Ang) and species.'
      write(6,"(i6,2x,3f10.5,2x,i3)") &
         ( ia, xa(1,ia) / Ang, xa(2,ia) / Ang, xa(3,ia) / Ang, &
           qm_atoms(ia)%spec, ia = 1, na_qm )
    endif

    ! Dump initial MM coordinates to output
    if ( has_mm ) then
      write(6,*)
      write(6,"('siesta-qmmm: Number of MM atoms: ',i6)") na_mm
      write(6,"('siesta-qmmm: Number of residues: ',i4)") nroaa
      write(6,'(/a)')'siesta-qmmm: MM atom coordinates (Ang)'
      write(6,"(i6,2x,A4,2x,A4,2x,i5,2x,3f12.6)") &
       ( i, mm_atoms(i)%atname, mm_atoms(i)%aaname, mm_atoms(i)%aanum, &
         xa(1,i+na_qm) / Ang, xa(2,i+na_qm) / Ang, xa(3,i+na_qm) / Ang,&
         i = 1, na_mm )
      write(6,'(/a)') 'siesta-qmmm: Atoms parameters'
    endif

    if ( has_qm .and. has_mm ) &
      write(6,"(i6,2x,A4,2x,f12.7,2x,f12.7)") ( i, qm_atoms(i)%attype, &
           qm_atoms(i)%lj_Em, qm_atoms(i)%lj_Rm, i = 1, na_qm )
    if ( has_mm ) &
      write(6,"(i6,2x,A4,2x,f9.6,2x,f9.6,2x,f9.6)") ( i, mm_atoms(i)%attype, &
          mm_atoms(i)%lj_Em, mm_atoms(i)%lj_Rm, mm_atoms(i)%pc, i = 1, na_mm )
      write(6,"('siesta-qmmm: ',71(1h=))")
  endif
  call pxfflush(6)

  if ( has_qm ) then
    ! Launch siesta process

    if ( num_resonances > 1 )  then
      call system('rm -f directories')

      do i = 1, num_resonances
        call system( 'echo '//trim(resonance(i)%path)//' >> directories' )
      enddo
    endif

    ! Write the 'siesta_qmmmslabel.siesta.fdf' file to be read by siesta
    ! as subroutine
    if ( launch_siesta_flag ) &
      call write_siesta_input( qmmm_files%qmlabel, na_qm, qm_atoms, &
                               coulombtype, rcorteqmmm, rcorteqm, rcortemm )

    ! Set physical units for communication with siesta
    call siesta_units( 'Bohr', 'Ry' )
    call create_fifos( qmmm_files%qmlabel )

    write( 6,'(a)' ) &
      'siesta-qmmm: unit cell vectors (Ang) from siesta run:    '
    write( 6,'(3f12.6)' ) ucell(1,1) / Ang, ucell(2,1) / Ang, ucell(3,1) / Ang
    write( 6,'(3f12.6)' ) ucell(1,2) / Ang, ucell(2,2) / Ang, ucell(3,2) / Ang
    write( 6,'(3f12.6)' ) ucell(1,3) / Ang, ucell(2,3) / Ang, ucell(3,3) / Ang
  endif

  call qmmm_timer( 'Setup', 2 )

  ! Start loop over restrained optimization steps
  call restr_read( ) ! Sets value for hasRestraints and restr_steps.
  nullify( qmstep )
  call re_alloc( qmstep, 1, restr_steps+1, 'qmstep', 'siesta_qmmm' )
  qmstep = 0

  Etot   = 0.0_dp ; Etots  = 0.0_dp
  fa     = 0.0_dp ; stress = 0.0_dp
  lattice_type = get_lattice_type(ucell)

  last_restr = .false.
  do isteprestr = 1, restr_steps +1
    if ( hasRestraints ) then
      write(6,*)
      write(6,'(A)')    '*******************************'
      write(6,'(A,i5)') '  Restrained Step : ', isteprestr
      write(6,'(A)')    '*******************************'
    endif

    if ( isteprestr == restr_steps +1 ) last_restr = .true.

    ! Start loop over coordinate changes
    istp = 0
    if ( restart_istep > 0 ) inicoor = restart_istep +1
    if ( inicoor > fincoor ) relaxd = .true.

    do istep = inicoor, fincoor
      ! This check should not be necessary.
      if ( inicoor > fincoor ) exit

      cml_p = .false.
      call timer( 'IterMD', 1 )

      call qmmm_timer( 'IterMD', 1 )
      istp = istp + 1

      write(6,'(2/,7x,A)')    '********************************'
      write(6,'(9x,a,2x,i6)')'siesta-qmmm:  Begin Step= ',istep
      write(6,'(7x,A,2/)')    '********************************'

      ! Update the Neb list eVery mneigh_freq steps
      actualiz = .false.
      if ( (mod(istep,mneigh_freq) == 0) .or. (istep == inicoor) ) &
        actualiz = .true.

      ! Reinsert atoms in the box to have better images unless the Nose
      ! thermostat is used since the self-consistency in the Nose algorithm
      ! is not compatible with the reinsertion of atoms
      if ( actualiz .and. (.not. therm_idyn) ) &
        call reinserting_atoms_in_box( lattice_type, na_u, na_qm, na_mm, &
                                       qm_atoms, mm_atoms, blocklist, ucell )


      if ( has_qm .and. has_mm .and. actualiz ) &
        call update_qmmm_mneighb( na_qm, na_mm, na_u, mm_atoms, xa, &
                                  ucell, rcorteqmmm, actualiz,      &
                                  lattice_type )

      ! Calculate Energy and Forces using Siesta as Subroutine
      if ( has_qm ) then
        call qmmm_timer( 'Siesta', 1 )
        call siesta_send_coordinates( qmmm_files%qmlabel, na_qm, xa, ucell, &
                                      na_mm, mm_atoms, (istep == fincoor) .and.&
                                      last_restr )
        call siesta_receive_forces( qmmm_files%qmlabel, na_qm, na_mm, Etots, &
                                    fa, stress )
        call qmmm_timer( 'Siesta', 2 )
        qmstep(isteprestr) = qmstep(isteprestr) + 1
      endif ! qm

      ! Start MMxQM loop
      do imm = 1, mmsteps
        step = step +1

        if ( mmsteps > 1 ) then
          write(6,*)
          write(6,'(A)')    '*******************************'
          write(6,'(A,i5)') '   MM x QM Step : ', imm
          write(6,'(A)')    '*******************************'
        endif

        if ( has_qm .and. has_mm ) call qmmm_timer( 'QM-MMcoupl', 1 )

        ! Calculation of last QM-MM interaction: LJ Energy and Forces only
        if ( has_qm .and. has_mm ) then
          call ljef( na_qm, na_u, xa, qm_atoms, mm_atoms, fa, stress, Elj, &
                     rcorteqmmm, ucell, lattice_type )
          call qmmm_timer( 'QM-MMcoupl', 2 )
        endif !qm & mm

        ! LinkAtom: set again link mm atoms parameters
        if ( has_qm .and. has_mm .and. linkatom )  then
          do i = 1, numlink
          do j = 1, link_atoms(i)%num_mm
            mm_atoms( link_atoms(i)%mm(j) )%pc    = link_atoms(i)%pc(j)
            mm_atoms( link_atoms(i)%mm(j) )%lj_Em = link_atoms(i)%Em(j)
          enddo
          enddo
        endif

        ! Calculate pure Solvent energy and forces
        if ( has_mm ) then
          call qmmm_timer( 'MMenergy', 1 )
          call mm_ene_fce( na_u, na_qm, na_mm, xa, Etot_amber, fce_amber,  &
                           stress_amber, mm_atoms, mm_connectivity, mm_top,&
                           actualiz, rcortemm, sfc, water, amass, ucell,   &
                           lattice_type, coulombtype, ewald )
          call qmmm_timer( 'MMenergy', 2 )
        endif

        ! converts fa to Kcal/mol/Ang
        fa(1:3,1:na_u) = fa(1:3,1:na_u) * Ang / eV * kcal
        if ( numlink > 0 ) then
          do ia = 1, numlink
          do j  = 1, num_resonances
            link_atoms(ia)%reson(j)%fa(1:3) = &
            link_atoms(ia)%reson(j)%fa(1:3) * Ang/eV * kcal
          enddo
          enddo
        endif
        stress(1:3,1:3) = stress(1:3,1:3) * Ang * Ang * Ang / eV * kcal

        ! add famber to fa
        if ( has_mm ) then
          fa(1:3,na_qm+1:na_u) = fa(1:3,na_qm+1:na_u) &
                                + fce_amber(1:3,1:na_mm)
          stress(1:3,1:3) = stress(1:3,1:3) + stress_amber(1:3,1:3)
        endif

        ! Calculation of LinkAtom Energy and Forces
        if ( has_qm .and. has_mm .and. linkatom ) then
          call link2( xa, na_u, na_mm, fa, stress, mm_top, &
                      Elink, lattice_type, ucell )
          ! Set again link atmos parameters to zero for next step
          do i = 1, numlink
          do j = 1, link_atoms(i)%num_mm
            link_atoms(i)%pc(j)           = mm_atoms(link_atoms(i)%mm(j))%pc
            mm_atoms( link_atoms(i)%mm(j) )%pc    = 0.0_dp

            if ( j > 1 ) &
              mm_atoms( link_atoms(i)%mm(j) )%pc = &
                                         mm_atoms( link_atoms(i)%mm(j) )%pc &
                                         + link_atoms(i)%pc(1) &
                                         / real(link_atoms(i)%num_mm-1, kind=dp)

            mm_atoms( link_atoms(i)%mm(j) )%lj_Em = 0.0_dp
          enddo
          enddo
        endif

        ! Calculation of restrained Optimization Energy and Forces
        if ( (imm == 1) .and. hasRestraints ) &
         call restr_calc( na_u, xa, fa, istp, isteprestr, ucell, lattice_type )

        ! Calculate electric field.
        if ( has_mm ) call elecfield( na_mm, xa(1:3,na_qm+1:na_u), ucell, &
                                      mm_atoms, fa(1:3,na_qm+1:na_u), Eelec )

        ! Converts fa to Ry/Bohr since the SIESTA subroutines for
        ! qmmm_move has got targeted parameters in SIESTA internal
        ! units (Ry, Bohr, etc..)
        fa(1:3,1:na_u)  = fa(1:3,1:na_u)  / Ang * eV / kcal
        stress(1:3,1:3) = stress(1:3,1:3) / (Ang * Ang * Ang) * eV / kcal

        ! Writes final energy decomposition
        Etot = Etots + Elj + ( ( Etot_amber + Elink + Eelec ) / kcal*eV )

        write(6,*)
        write(6,'(/,a)') 'siesta-qmmm: Energy Decomposition (eV):'
        if ( has_qm ) write(6,'(a,2x,F16.6)')          'Esiesta:', Etots / eV
        if ( has_qm .and. has_mm ) write(6,'(a,2x,F16.6)') 'Elj:    ', Elj / eV
        if ( has_mm )          write(6,'(a,2x,F16.6)') 'Esolv:  ', &
                                                    Etot_amber / kcal
        if ( abs(Elink) > 0.0_dp ) write(6,'(a,2x,F16.6)') 'Elink:  ', &
                                                            Elink / kcal
        if ( abs(Eelec) > 0.0_dp ) write(6,'(a,2x,F16.6)') 'Eefield:', &
                                                            Eelec / kcal
        if ( has_qm .or. has_mm ) write(6,'(a,2x,F16.6)') 'Etot:  ', Etot / eV
        call pxfflush(6)

        ! Sets fa to zero inside mmxqm step
        if ( has_qm .and. has_mm .and. (imm /= 1) ) then
          fa(1:3,1:na_qm) = 0.0_dp
          va(1:3,1:na_qm) = 0.0_dp
          if ( linkatom ) then
            do i = 1, numlink
              fa(1:3,link_atoms(i)%mm(1)+na_qm) = 0.0_dp
              va(1:3,link_atoms(i)%mm(1)+na_qm) = 0.0_dp
            enddo
          endif
        endif

        ! Compute stress without internal molecular pressure
        call remove_intramol_pressure( ucell, stress, na_u, xa, fa, mstress)

        ! Impose constraints to atomic movements by chAnging forces
        if ( RemoveIntraMolecularPressure ) then
          ! Consider intramolecular pressure-removal as another
          ! kind of constraint
          call fixed( ucell, mstress, na_qm, isa, amass, xa, fa, cstress, &
                      cfa, ntcon )
        else
          call fixed( ucell, stress, na_qm, isa, amass, xa, fa, cstress, &
                      cfa, ntcon )
        endif

        ! Impose constraints to atomic movements by chAnging forces
        call qmmm_fixed2( na_qm, na_mm, na_u, nfree, blocklist, mm_atoms,&
                          fa, cfa, va )
        call add_siesta_qmmm_ntcon( ucell, na_qm, na_mm, isa, amass, xa, &
                                    blocklist, mm_atoms, ntcon )

        ! Calculate and output Zmatrix forces
        if ( lUseZmatrix .and. opt_idyn ) then
          call CartesianForce_to_ZmatForce( na_u, xa, fa )
          if ( IOnode ) call iofaZmat()
        endif

        ! Compute kinetic contribution to stress
        kin_stress(1:3,1:3) = 0.0_dp
        volume_of_some_cell = volcel(ucell)
        do ia = 1, na_u
        do jx = 1, 3
        do ix = 1, 3
          kin_stress(ix,jx) = kin_stress(ix,jx) - &
              amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume_of_some_cell
        enddo
        enddo
        enddo

        ! Add kinetic term to stress tensor
        tstress = stress + kin_stress

        ! Force output
        if ( IOnode ) then
          call qmmm_write_forces( na_u, fa, cfa, ftol, inicoor == istep,&
                                  final, writef, idyn == 6, dx )
          call qmmm_write_stress_pressure( idyn, final, varcel, tp, &
                      RemoveIntraMolecularPressure, FreeE, stress, kin_stress, &
                      mstress, cstress, tstress, volume_of_some_cell )
        endif

        ! Accumulate coordinates in PDB/CRD file for animation
        call write_pdbcrd( na_qm, xa, na_u, step, wricoord, na_mm, &
                           nspec, atsym, qm_atoms, mm_atoms, ucell )

        istep_mv = istep ! Used only since qmmm_move uses istep as "inout".
         ! Move atoms, i.e. update coordinates.
        call qmmm_move( istep_mv, relaxd, na_qm, na_mm,qm_atoms, mm_atoms )
        if ( (.not. opt_idyn) .and. (.not. fc_idyn) ) then
          write(6,'(/,a,f12.6)') 'siesta-qmmm: Kinetic Energy (eV):', &
                                  Ekinion / eV
          write(6,'(/,a,f12.6)') &
            'siesta-qmmm: Total Energy + Kinetic (eV):', (Etot+Ekinion) / eV
          write(6,'(/,a,f12.3,a)') 'siesta-qmmm: System Temperature:', &
                                    tempion, ' K'
        endif

        ! Write Energy in file
        call write_energy( istep, dt, Etot, Ekinion, tempion, na_u, cfa )

        ! Save last atomic positions and velocities
        call qmmm_ioxv( qmmm_files%slabel, 'w', istep, int_dummy, na_u, ucell,&
                        vcell, xa, va, foundxv, foundva )

        ! Write atomic restraints each step
        if ( nrestr > 0 ) call write_restr( restr_r, nrestr )

        ! Sets variables for next siesta cycle
        fa           = 0.0_dp
        stress       = 0.0_dp
        lattice_type = get_lattice_type(ucell)

        ! Calculation Hlink's New positions
        if ( has_qm .and. has_mm .and. linkatom ) &
          call link3( na_qm, na_mm, qm_atoms, mm_atoms, mm_top, ucell, &
                      lattice_type )

      enddo ! imm

      if ( ( mmsteps /= 1 ) .and. ( imm /= 1 ) ) relaxd = .false.

      call qmmm_timer( 'IterMD', 2 )
    enddo ! istep, end of coordinate-relaxation loop

    if ( hasRestraints ) then
      call restr_update( )
      call write_react_crd( Etot, restr_r, isteprestr, na_qm, na_mm, na_u, &
                            xa, qm_atoms, mm_atoms, nspec, atsym )
      call ioxvrestr( na_u, ucell, xa, va, isteprestr )
    endif
  enddo ! isteprestr, end of restrain optimization loop

  ! Stop siesta process
  if ( has_qm ) call siesta_quit( 'all' )

  ! Dump last coordinates to output
  if ( writeipl ) then
    if( has_qm ) then
      write(6,'(/a)')'siesta-qmmm: Last atomic coordinates (Ang) and'//&
                      ' species'
      write(6,"(i6,2x,3f10.5,2x,i3)") ( ia, xa(1,ia) / Ang, xa(2,ia) / Ang, &
                                        xa(3,ia) / Ang, qm_atoms(ia)%spec,  &
                                        ia = 1, na_qm )
    endif
    if ( has_mm ) then
      write(6,'(/a)') 'siesta-qmmm: Last solvent coordinates (Ang)'
      write(6,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)') &
        ( 'ATOM', i, mm_atoms(i)%atname, mm_atoms(i)%aaname, 'A', &
          mm_atoms(i)%aanum, xa(1,na_qm+i) / Ang, xa(2,na_qm+i) / Ang, &
          xa(3,na_qm+i) / Ang, i = 1, na_mm )
    endif
  endif
  call pxfflush(6)

  call mm_dealloc( )
  call deallocate_link_arrays()

  call de_alloc( atxres     , 'atxres'     , 'siesta_qmmm' )
  call de_alloc( blocklist  , 'blocklist'  , 'siesta_qmmm' )
  call de_alloc( qmstep     , 'qmstep'     , 'siesta_qmmm' )

  do ia = 1, na_mm
    call de_alloc( mm_connectivity(ia)%imp_at , 'imp_at' , 'siesta_qmmm' )
    call de_alloc( mm_connectivity(ia)%ange_at, 'ange_at', 'siesta_qmmm' )
    call de_alloc( mm_connectivity(ia)%angm_at, 'angm_at', 'siesta_qmmm' )
    call de_alloc( mm_connectivity(ia)%dihe_at, 'dihe_at', 'siesta_qmmm' )
    call de_alloc( mm_connectivity(ia)%dihm_at, 'dihm_at', 'siesta_qmmm' )
  enddo
  deallocate( atsym, mm_atoms, qm_atoms, mm_connectivity )
  call reset_citations()

  ! Stop time counter and print final date/time
  call qmmm_timer( 'Siesta-qmmm', 2 )
  call qmmm_timer( 'all', 3 )
  call timestamp('End of run')

#if MPI
  call MPI_Finalize( MPIerror )
#endif

end program SIESTA_QMMM
