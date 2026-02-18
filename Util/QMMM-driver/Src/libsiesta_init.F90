module m_libsiesta_init
  !! Initializes the 'siesta' data structures used by the qmmm-driver
  !! program.
  private
  public :: libsiesta_init

  ! Most of the functionality here is not needed, but a careful
  ! analysis is needed before pruning, due to interactions
  ! between sections. The work can be re-used for general refactoring.
contains

  subroutine libsiesta_init( na_qm, na_mm, numlink, qm_atoms, mm_atoms, &
                             link_atoms )
    !! Calls SIESTA data structure initializations.
    use alloc             , only : re_alloc
    use atomlist          , only : initatomlists, amass, qtot, iza
    use atm_types         , only : species, nspecies
    use chemical          , only : read_chemical_types, number_of_species
    use files             , only : slabel
    use fdf               , only : fdf_get
    use linkatoms         , only : link_atom_t
    use m_energies        , only : DUext, Eharrs, Eions, Ekinion,   &
                                   Emad, Emeta, Emm, Entropy, FreeE
    use m_forces          , only : fa, cfa
    use m_iostruct        , only : write_struct
    use m_ioxv            , only : xv_file_read
    use m_steps           , only : inicoor, fincoor, istp
    use memory_log        , only : memory_report
    use metaforce         , only : initmeta
    use mm_topology       , only : qm_atom_t, mm_atom_t
    use molecularmechanics, only : inittwobody
    use parallel          , only : Node, IOnode
    use precision         , only : dp
    use qmmm_struct_init_m, only : qmmm_struct_init
    use siesta_cmlsubs    , only : siesta_cml_init
    use siesta_geom       , only : na_u, ucell, shape, xa, va, isa
    use siesta_options    , only : idyn, tempinit, nmove, charnet, ia1, ia2, &
                                   ifinal, istart
    use siesta_qmmm_options,only : read_siesta_options
    use sys               , only : die, bye
    use zmatrix           , only : lUseZmatrix, &
                                   write_canonical_ucell_and_Zmatrix
#ifdef MPI
    use parallel          , only : Nodes
#endif

    implicit none
    integer, intent(in) :: na_qm
      !! Number of classical (MM) atoms.
    integer, intent(in) :: na_mm
      !! Number of classical (MM) atoms.
    integer, intent(in) :: numlink
      !! Number of link atoms.
    type(qm_atom_t), intent(in) :: qm_atoms(na_qm)
      !! QM atoms.
    type(mm_atom_t), intent(in) :: mm_atoms(na_mm)
      !! MM atoms.
    type(link_atom_t), intent(in) :: link_atoms(numlink)
      !! Positions for linkatoms.

    logical :: atmonly, struct_only

    IOnode = ( Node == 0 )

    ! Print version information
    if ( IOnode ) then
      write(6,"(/,a,/)") " --- Initializing libsiesta data structures --- "
#ifdef MPI
      if ( Nodes > 1 ) then
        write(6,'(/,a,i4,a)') "* Running on ", Nodes, " nodes in parallel."
      else
        write(6,'(/,a,i4,a)') "* Running in serial mode with MPI."
      endif
#else
      write(6,'(/,a,i4,a)') "* Running in serial mode."
#endif
    endif

    ! Start time counter
    call timer( 'siesta', 0 )
    call timer( 'siesta', 1 )
    call timer( 'Setup', 1 )

    ! Variable initialization.
    DUext = 0.0_dp ; Eharrs  = 0.0_dp ; Entropy = 0.0_dp
    Eions = 0.0_dp ; Ekinion = 0.0_dp ; FreeE   = 0.0_dp
    Emad  = 0.0_dp ; Emeta   = 0.0_dp ; Emm     = 0.0_dp

    call siesta_cml_init( )
    call memory_report( level     = fdf_get('alloc_report_level',0), &
                        file      = trim(slabel)//'.alloc', &
                        threshold = fdf_get('alloc_report_threshold', 0._dp), &
                        printNow  = .false. )

    ! Initialise force field component
    call inittwobody( )

    ! Initialize pseudopotentials and atomic orbitals
    nspecies = 1
    if ( na_qm > 0 ) then
      call read_chemical_types( )
      nspecies = number_of_species() +1
    endif
    allocate( species(nspecies) )

    call update_chemical_types( )

    atmonly = fdf_get( 'Atom-Setup-Only', .false. )
    if ( atmonly ) call bye("End of atom setup")

    ! Read geometry and initialize atom lists.
    ! Sets na_u, isa, ucell
    call qmmm_struct_init( na_qm, na_mm, numlink, qm_atoms, mm_atoms, &
                           link_atoms )
    call initatomlists()    ! Sets iza

    ! Early exit if only checking the structure.
    struct_only = fdf_get( 'Output-Structure-Only', .false. )

    if ( IOnode ) then
      call write_struct( ucell, na_u, isa, iza, xa )
      if ( lUseZmatrix ) &
        call write_canonical_ucell_and_Zmatrix( filename = &
                                                "OUT.UCELL.ZMATRIX.INITIAL" )
    endif
    if ( struct_only ) call bye( "End of structure processing." )
    ! End of Initial Structure Processing.

    if ( IOnode ) then
      write(6,'(/,a,20("*"),a,28("*"))') "siesta: ", " Simulation parameters "
      write(6,'(a)') 'siesta:'
      write(6,'(a)') 'siesta: The following are some of the parameters'//&
                     ' of the simulation.'
      write(6,'(a)') 'siesta: A complete list of the parameters '// &
                     'used, including default values,'
      write(6,'(a,a)') 'siesta: can be found in file out.fdf'
      write(6,'(a)') 'siesta:'
    endif

    ! Allocate other arrays based on read sizes
    nullify( fa, cfa )
    call re_alloc( fa , 1, 3, 1, na_u, "fa" , "siesta_init" )
    call re_alloc( cfa, 1, 3, 1, na_u, "cfa", "siesta_init" )
    call read_siesta_options( na_u, nspecies )

    qtot = qtot - charnet ! qtot set in initatomlists, charnet set in redata

    ! Warn the user: if not doing a direct optimization, the Zmatrix
    ! coordinates are no longer updated. Only coordinates are treated.
    if ( lUseZmatrix ) then
      if ( idyn /= 0 ) then
        write(6,"(a)") " WARNING: Zmatrix form will be used only for input !"
        write(0,"(a)") " WARNING: Zmatrix form will be used only for input !"
      endif
    endif

    ! Madelung correction for charged systems.
    if ( abs(charnet) > 1.0e-10_dp) &
      call madelung( ucell, shape, charnet, Emad )

    ! Initialise metadynamic forces if required
    call initmeta( )

    if ( idyn == 0) then
      inicoor = 0
      fincoor = nmove
    else if ( (idyn >= 1) .and. (idyn <= 5) ) then
      inicoor = istart
      fincoor = ifinal
    else if ( idyn == 6 ) then
      inicoor = 0
      fincoor = (ia2-ia1+1)*3*2
    else if ( idyn == 7 ) then
      call die( "'PHONON' support is deprecated" )

    else if ( idyn == 8 ) then
      inicoor = 0
      fincoor = huge(1)
    else
      call die('siesta: Wrong value for idyn.')
    endif

    ! Build initial velocities according to Maxwell-Bolzmann distribution.
    if ( (idyn /= 0) .and. (idyn /= 6) .and. (.not. xv_file_read) ) &
      call vmb( na_u, tempinit, amass, xa, isa, va )

    istp = 0
    call timer( 'Setup', 2 )
    if ( IOnode ) call pxfflush( 6 )

  end subroutine libsiesta_init

  subroutine update_chemical_types( )
    !! Adds the "MM" type to atomic species information.
    use chemical, only : chemical_type, chemical_list

    implicit none
    integer :: ispec, nspec_old
    type(chemical_type), allocatable :: backup_list(:)

    nspec_old = 0
    if ( allocated(chemical_list) ) then
      nspec_old = size(chemical_list,1)
      allocate( backup_list(nspec_old) )

      do ispec = 1, nspec_old
        backup_list(ispec)%ps_file_spec = chemical_list(ispec)%ps_file_spec
        backup_list(ispec)%spec_label   = chemical_list(ispec)%spec_label
        backup_list(ispec)%z            = chemical_list(ispec)%z
      enddo

      deallocate( chemical_list )
      allocate( chemical_list(nspec_old+1) )

      do ispec = 1, nspec_old
        chemical_list(ispec)%ps_file_spec = backup_list(ispec)%ps_file_spec
        chemical_list(ispec)%spec_label   = backup_list(ispec)%spec_label
        chemical_list(ispec)%z            = backup_list(ispec)%z
      enddo
      deallocate(backup_list)
    else
      allocate( chemical_list(1) )
    endif

    ispec = nspec_old +1
    chemical_list(ispec)%ps_file_spec = "MM.psf"
    chemical_list(ispec)%spec_label   = "MM"
    chemical_list(ispec)%z            = 0
  end subroutine update_chemical_types
end module m_libsiesta_init
