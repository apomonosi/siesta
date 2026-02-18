! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module qmmm_struct_init_m
  private
  public :: qmmm_struct_init

contains
  subroutine qmmm_struct_init( na_qm, na_mm, numlink, qm_atoms, &
                               mm_atoms, link_atoms )
    use alloc         , only : re_alloc, de_alloc
    use atomlist      , only : elem, iza
    use atm_types     , only : nspecies
    use fdf           , only : fdf_get, leqi
    use files         , only : slabel
    use linkatoms     , only : link_atom_t
    use m_iostruct    , only : read_struct
    use m_iotdxv      , only : iotdxv, tdxv_file_read
    use m_ioxv        , only : ioxv, xv_file_read
    use m_mpi_utils   , only : broadcast
    use m_steps       , only : inicoor, fincoor
    use mm_topology   , only : qm_atom_t, mm_atom_t
    use parallel      , only : IOnode
    use precision     , only : dp
    use siesta_cml    , only : cml_p, cmlStartModule, cmlAddMolecule, mainXML,&
                               cmlAddLattice, cmlAddProperty, cmlEndModule
    use siesta_geom   , only : na_u, ucell, va, xa, vcell, isa, cisa
    use siesta_options, only : dt, usesaveddata
    use sys           , only : die
    use units         , only : Ang
    use zmatrix       , only : luseZmatrix

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

    real(dp), external :: volcel
    external           :: iozm

    character(len=22) :: dyntyp
    character(len=60) :: restart_file
    integer           :: ia, ix, istart_tmp, ifinal_tmp
    logical           :: foundxv, foundzm, found_restart, step_back, &
                         use_struct_file, usesavexv, usesavezm, writic

    ! Read number of atoms and coordinates, and unit cell
    ! Do it twice for backwards compatibility.
    use_struct_file = fdf_get( 'useStructFile'   , .false.         )
    use_struct_file = fdf_get( 'MD.useStructFile', use_struct_file )

    ! Sets na_u, xa, and isa
    if ( use_struct_file ) then
      call read_struct( na_u, ucell )

      write( 6, '(/,a)' ) 'siestaqmmm: WARNING: If present, an XV file '//&
                          'prevails over previous structure input.'
    else
      na_u = na_mm + na_qm + numlink
      nullify(xa)
      call re_alloc( xa, 1, 3, 1, na_u, 'xa', 'struct_init' )
      call re_alloc( isa, 1, na_u, 'isa', 'struct_init' )

      if ( na_qm > 0 ) then
        do ia = 1, na_qm
          xa(1:3,ia) = qm_atoms(ia)%r(:)
        enddo
      endif
      if ( na_mm > 0 ) then
        do ia = 1, na_mm
          xa(1:3,ia+na_qm) = mm_atoms(ia)%r(:)
        enddo
      endif
      if ( numlink > 0 ) then
        do ia = 1, numlink
          xa(1:3,na_u-numlink+ia) = &
            link_atoms(ia)%reson(1)%r(1:3)
        enddo
      endif

      do ia = 1, na_qm
        isa(ia) = qm_atoms(ia)%spec
      enddo
      if ( (na_mm > 0) .or. (numlink > 0) ) isa(na_qm+1:na_u) = nspecies
    endif

    ! Prepare iza here: it might be needed by ioxv
    nullify( iza, va )
    call re_alloc( iza, 1, na_u, 'iza', 'struct_init' )
    call re_alloc( va , 1, 3, 1, na_u, 'va', 'struct_init' )
    va(1:3,1:na_u) = 0.0_dp
    vcell(1:3,1:3) = 0.0_dp

    ! Options read here instead of in siesta_options
    usesavexv = fdf_get( 'MD.useSaveXV', usesaveddata )
    usesavezm = fdf_get( 'MD.useSaveZM', usesaveddata )
    writic    = fdf_get( 'WriteCoorInitial', .true. )

    ! Read Z-matrix coordinates and forces from file
    if ( luseZmatrix ) then
      foundzm = .false.

      if ( usesavezm ) then
        call iozm( 'read', ucell, vcell, xa, foundzm )
        if ( IOnode .and. (.not. foundzm) ) &
          write( 6, '(/,a)' ) 'siestaqmmm: WARNING: ZM file not found.'
      endif
    endif

    ! For TDDFT start or restart, structure must be read from TDXV.
    foundxv = .false.
    dyntyp  = trim( fdf_get( 'MD.TypeOfRun', 'CG' ) )
    if ( leqi(dyntyp,'TDED') ) then
      call iotdxv( 'read', ucell, vcell, na_u, isa, iza, xa, va, foundxv )
      if ( (.not. foundxv) .and. IOnode ) &
        write(6,'(/,a)') 'siestaqmmm: TDXV file not found.'

      xv_file_read = tdxv_file_read

    elseif ( usesavexv ) then
      call ioxv('read', ucell, vcell, na_u, isa, iza, xa, va, foundxv)
      istart_tmp = fdf_get('MD.InitialTimeStep', 1)
      ifinal_tmp = fdf_get('MD.FinalTimeStep', 1)

      if ( (.not. foundxv) .and. IOnode ) &
        write(6,'(/,a)') 'siestaqmmm: XV file not found.'

      ! For a Verlet/Nose/PR/NPR run with more than one time step, if
      ! the RESTART file is not found, backward-propagate the atomic
      ! positions by one time step using the Euler method
      if ( foundxv .and. ( (ifinal_tmp - istart_tmp) > 0 ) ) then
        ! In isa we have first the QM and then the MM. Force the mm isa to
        ! nspecies.
        if ( (na_mm > 0) .or. (numlink > 0) ) isa(na_qm+1:na_u) = nspecies

        step_back = .true.
        if ( leqi(dyntyp, 'verlet') ) then
          restart_file = trim(slabel) // '.VERLET_RESTART'
        elseif ( leqi(dyntyp, 'nose') ) then
          restart_file = trim(slabel) // '.NOSE_RESTART'
        elseif ( leqi(dyntyp, 'parrinellorahman') ) then
          restart_file = trim(slabel) // '.PR_RESTART'
        elseif ( leqi(dyntyp, 'noseparrinellorahman') ) then
          restart_file = trim(slabel) // '.NPR_RESTART'
        else
          step_back = .false.
        endif

        if ( step_back ) then
          if ( IOnode ) inquire( file = restart_file, exist = found_restart )
          call broadcast( found_restart )

          if ( .not. found_restart ) then
            xa = xa - dt * va

            if ( IOnode ) &
              write(6,'(a)') 'WARNING: '//trim(restart_file)//&
                             ' not found--reading only XV file'//&
                             ' and moving back 1 time step using Euler.'
          endif
        endif
      endif
    endif

    ! Dump initial coordinates to output
    if ( writic .and. IOnode ) then
      write(6,'(/a)') 'siestaqmmm: Atomic coordinates (Bohr) and species'
      write(6,"('siestaqmmm: ',2x,3f10.5,i3,3x,i6)") &
        ( (xa(ix,ia), ix=1,3), isa(ia), ia, ia=1, na_u )
    endif

    ! Automatic cell generation
    if ( volcel(ucell) < 1.0d-8 ) &
      call die( "You must explicity set a simulation cell for QMMM simulations.")

    ! Output of initial system details:
    if ( cml_p ) then
      call cmlStartModule( xf = mainXML, title = 'Initial System' )
      call cmlAddMolecule( xf = mainXML, natoms = na_u, coords= xa / Ang,&
                           elements = elem, atomRefs = cisa )

      call cmlAddLattice( xf = mainXML, cell = ucell, &
        units = 'siestaUnits:angstrom', dictref = 'siesta:ucell' )
      call cmlEndModule( xf = mainXML )
    endif

  end subroutine qmmm_struct_init
end module qmmm_struct_init_m
