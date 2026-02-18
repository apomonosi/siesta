module qmmm_write
  !! This module contains subroutines that handle the output for the
  !! QMMM driver.
  use precision   , only : dp
  use qmmm_files_m, only : qmmm_files
  use units       , only : Ang, eV

  implicit none
  public :: write_energy
  public :: write_react_crd
  public :: write_siesta_input
  public :: write_restr
  public :: write_pdbcrd
  public :: qmmm_write_forces
  public :: qmmm_write_positions
  public :: qmmm_write_stress_pressure

  private
contains
  subroutine write_energy( istp, dt, Epot, Ekin, temper, natoms, cfa )
    !! Writes energy contributions and temperature to a *.ene file.

    use precision          , only : dp
    use siesta_qmmm_options, only : opt_idyn, fc_idyn

    implicit none
    integer , intent(in) :: istp
      !! Current MD step.
    real(dp), intent(in) :: dt
      !! Time step in fs.
    real(dp), intent(in) :: Epot
      !! Potential energy of the system.
    real(dp), intent(in) :: Ekin
      !! Kinetic energy of the system.
    real(dp), intent(in) :: temper
      !! Average temperature of the system.
    integer , intent(in) :: natoms
      !! Number of atoms in unit cell.
    real(dp), intent(in) :: cfa(3,natoms)
      !! Atomic (constrained) forces.

    integer  :: funit
    real(dp) :: cfmax
    external :: io_assign, io_close

    cfmax = maxval( abs(cfa) )

    call io_assign( funit )
    open( funit, file = qmmm_files%ene, form = 'formatted', position = 'append', &
          status = 'unknown' )

    if ( opt_idyn .or. fc_idyn ) then
      if ( istp == 0 ) &
        write( funit, * ) '#step   Energy(eV)   Max. Force (eV/Ang)   '
      write( funit, '(i5,2x,F14.7,2x,F14.7)' ) istp, Epot / eV, cfmax * Ang / eV
    else
      if ( istp == 1 ) &
        write( funit, * ) '# time(fs)      K.E.(eV)     P.E.(ev)      '//&
                          'T.E.(eV)    Temp(k)'
        write( funit, '(F9.2,2x,F12.4,2x,F12.4,2x,F12.4,2x,F7.2)' ) &
          istp * dt, Ekin / eV, Epot / eV, (Epot + Ekin) / eV, temper
    endif

    call io_close( funit )
  end subroutine write_energy

  subroutine write_pdbcrd( naqm, rclas, natot, istp, wricoord, namm, &
                           nesp, atsym, qm_atoms, mm_atoms, mmcell )

    !! Writes atomic coordinates to an output *.pdb file.
    use linkatoms  , only : numlink, link_atoms, resonance
    use mm_topology, only : mm_atom_t, qm_atom_t
    use precision  , only : dp
    use sys        , only : die

    implicit none
    integer         , intent(in) :: naqm
      !! Number of QM atoms.
    integer         , intent(in) :: namm
      !! Number of MM atoms.
    integer         , intent(in) :: natot
      !! Total number of QM+MM atoms.
    integer         , intent(in) :: istp
      !! Current MD step.
    integer         , intent(in) :: wricoord
      !! Every wricoord steps we write the MDCRD file too.
    real(dp)        , intent(in) :: rclas(3,natot)
      !! Atomic coordinates.
    integer         , intent(in) :: nesp
      !! Total number of different atomic species.
    character(len=2), intent(in) :: atsym(nesp)
      !! Atom symbols for each species.
    type(qm_atom_t) , intent(in) :: qm_atoms(naqm)
      !! Structure containing information for QM atoms.
    type(mm_atom_t) , intent(in) :: mm_atoms(namm)
      !! Structure containing information for MM atoms.
    real(dp)        , intent(in) :: mmcell(3,3)
      !! Simulation cell.

    character(len=30) :: title
    character(len=4)  :: ch
    character(len=1)  :: chain(natot)
    integer           :: unitp, unitc, unitl, iat, ios
    logical           :: foundc
    real(dp)          :: degtorad

    logical , save    :: frstme = .true.
    real(dp), save    :: cell_a, cell_b, cell_c, &
                         cell_alp, cell_bet, cell_gam

    chain(1:natot) = 'A'
    do iat = 1, namm
      if ( .not. mm_atoms(iat)%is_qm_neighbour ) chain(iat+naqm) = 'L'
      if ( .not. mm_atoms(iat)%is_blocked ) chain(iat+naqm) = 'B'
      if ( (.not. mm_atoms(iat)%is_qm_neighbour) .and. &
           (.not. mm_atoms(iat)%is_blocked) ) chain(iat+naqm) = 'C'
    enddo

    if ( frstme ) then
      degtorad = 4.0_dp * atan( 1.0_dp ) / 180.0_dp
      call uncell( mmcell, cell_a, cell_b, cell_c, cell_alp, &
                   cell_bet, cell_gam, degtorad )

      ! Writes initial PDB coords in .init.pdb file
      call io_assign( unitp )
      open( unitp, file = qmmm_files%pdbi, form = 'formatted', status = 'unknown' )
      write( unitp, '(A,3(F9.3),3(F7.2),A)' ) 'CRYST1', &
        cell_a / Ang, cell_b / Ang, cell_c / Ang, &
        cell_alp, cell_bet, cell_gam, ' P 1           1'

      do iat = 1, naqm
        write( unitp, '(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)' ) &
          'ATOM', iat, atsym(qm_atoms(iat)%spec), 'STO ', chain(iat), iat, &
          rclas(1,iat) / Ang, rclas(2,iat) / Ang, rclas(3,iat) / Ang
      enddo

      if ( namm > 9999 ) then
        do iat = 1, namm
          write( unitp, '(A4,I7,2x,A4,A4,I5,4x,3f8.3)' ) &
            'ATOM', naqm+iat, mm_atoms(iat)%atname, mm_atoms(iat)%aaname, &
            mm_atoms(iat)%aanum, rclas(1,iat+naqm) / Ang, &
            rclas(2,iat+naqm) / Ang, rclas(3,iat+naqm) / Ang
        enddo
      else
        do iat = 1, namm
          write( unitp, '(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)' ) &
            'ATOM', naqm+iat, mm_atoms(iat)%atname, mm_atoms(iat)%aaname, &
            chain(naqm+iat), mm_atoms(iat)%aanum, rclas(1,iat+naqm) / Ang,&
            rclas(2,iat+naqm) / Ang, rclas(3,iat+naqm) / Ang
        enddo
      endif

      do iat = 1, numlink
      write( unitp, '(A4,I7,2x,A2,2X,A4,A,I4,4x,3f8.3)' ) &
        'ATOM', naqm+namm+iat, atsym(resonance(1)%isa(iat)), 'LA  ', 'A', iat,&
        link_atoms(iat)%reson(1)%r(1) / Ang, &
        link_atoms(iat)%reson(1)%r(2) / Ang, &
        link_atoms(iat)%reson(1)%r(3) / Ang
      enddo
      write( unitp, '(A3)' ) 'END'
      call io_close( unitp )

      ! Writes in .crd file at begining
      inquire( file = qmmm_files%mcrd, exist = foundc )
      if ( foundc ) then
        call io_assign( unitc )
        open( unitc, file = qmmm_files%mcrd, status = 'old' )
        read( unitc, '(a30)', iostat = ios ) title
        if ( ios /= 0 ) &
          call die ('write_pdbcrd: Problem while reading from CRD file.' )

        call io_close( unitc )

        ch = title(1:4)
        if ( ch /= 'File' ) then
          call io_assign( unitc )
          open( unitc, file = qmmm_files%mcrd, form = 'formatted', &
                status = 'unknown' )
          write( unitc, '(20a)' ) 'File CRD'
          call io_close( unitc )
        endif
      else
        call io_assign( unitc )
        open( unitc, file = qmmm_files%mcrd, form = 'formatted', &
              status = 'unknown' )
        write( unitc, '(20a)' ) 'File CRD'
        call io_close( unitc )
      endif

      frstme = .false.
    endif ! frstme
    call pxfflush( 6 )

    ! Writes in the .crd file if needed.
    if ( mod(istp, wricoord) == 0 ) then
      call io_assign( unitc )
      open( unitc, file = qmmm_files%mcrd, form = 'formatted', position = 'append',&
            status = 'unknown' )
      write( unitc, '(10F8.3)' ) ( rclas(1:3,iat)   / Ang, iat = 1, natot ), &
                                 ( link_atoms(iat)%reson(1)%r(1:3), &
                                   iat = 1, numlink )
      call io_close( unitc )
    endif

    ! Writes current coordinates in PDB format on the .last.pdb file.
    call io_assign( unitl )
    open( unitl, file = qmmm_files%pdbl, form = 'formatted', status = 'unknown' )
    write( unitl, '(A,3(F9.3),3(F7.2),A)' ) 'CRYST1', &
      cell_a / Ang, cell_b / Ang, cell_c / Ang, &
      cell_alp, cell_bet, cell_gam, ' P 1           1'

    do iat = 1, naqm
      write( unitl, '(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)' ) &
        'ATOM', iat, atsym(qm_atoms(iat)%spec), 'STO ', chain(iat), iat, &
        rclas(1,iat) / Ang, rclas(2,iat) / Ang, rclas(3,iat) / Ang
    enddo

    if ( namm > 9999 ) then
      do iat = 1, namm
        write( unitl, '(A4,I7,2x,A4,A4,I5,4x,3f8.3)' ) &
          'ATOM', naqm+iat, mm_atoms(iat)%atname, mm_atoms(iat)%aaname, &
          mm_atoms(iat)%aanum, rclas(1,naqm+iat) / Ang, &
          rclas(2,naqm+iat) / Ang, rclas(3,naqm+iat) / Ang
      enddo
    else
      do iat = 1, namm
        write( unitl, '(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)' ) &
          'ATOM', naqm+iat, mm_atoms(iat)%atname, mm_atoms(iat)%aaname, &
          chain(naqm+iat), mm_atoms(iat)%aanum, rclas(1,naqm+iat) / Ang,&
          rclas(2,naqm+iat) / Ang, rclas(3,naqm+iat) / Ang
      enddo
    endif

    do iat = 1, numlink
      write( unitl, '(A4,I7,2x,A2,2X,A4,A,I4,4x,3f8.3)' ) &
        'ATOM', naqm+namm+iat, atsym(resonance(1)%isa(iat)), 'LA  ', 'A', &
        iat, link_atoms(iat)%reson(1)%r(1) / Ang, &
        link_atoms(iat)%reson(1)%r(2) / Ang, &
        link_atoms(iat)%reson(1)%r(3) / Ang
    enddo
    write( unitl, '(A3)' ) 'END'
    call io_close( unitl )
  end subroutine write_pdbcrd

  subroutine write_react_crd( Etots, rtot, irestrstep, naqm, namm, natot, &
                              rclas, qm_atoms, mm_atoms, nesp, atsym )
    !! Writes the reaction coordinate and its associated PDB file.
    use mm_topology, only : mm_atom_t, qm_atom_t
    use precision  , only : dp

    implicit none
    real(dp)        , intent(in) :: Etots
      !! Total energy.
    real(dp)        , intent(in) :: rtot(:)
      !! Total coordinate value.
    integer         , intent(in) :: naqm
      !! Number of QM atoms.
    integer         , intent(in) :: namm
      !! Number of MM atoms.
    integer         , intent(in) :: natot
      !! Total number of QM+MM atoms.
    integer         , intent(in) :: irestrstep
      !! Current restrained optimization step.
    real(dp)        , intent(in) :: rclas(3,natot)
      !! Atomic coordinates.
    type(qm_atom_t) , intent(in) :: qm_atoms(naqm)
      !! Residue indices for each atom in the QM region.
    type(mm_atom_t) , intent(in) :: mm_atoms(namm)
      !! Residue indices for each atom in the MM region.
    integer         , intent(in) :: nesp
      !! Total number of different atomic species.
    character(len=2), intent(in) :: atsym(nesp)
      !! Atom symbols for each species.

    integer :: iat, funit
    character(len=10) :: nrest
    character(len=100):: fmt

    ! Writes .rce file
    call io_assign( funit )
    open( funit, file = qmmm_files%cene, form = 'formatted', position = 'append', &
          status = 'unknown' )

    write( nrest, '(I10)' ) size(rtot,1)
    write( fmt  , '(A)'   ) '('//trim(nrest)//'F10.4,2x,F14.7)'
    write( funit, trim(fmt) ) rtot(:), Etots / eV
    call io_close( funit )

    ! Writes .rcg file
    call io_assign( funit )
    open( funit, file = qmmm_files%cpdb, form = 'formatted', position = 'append', &
          status = 'unknown' )
    write( funit, '(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)' )

    do iat = 1, naqm
      write( funit, '(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)' ) &
            'ATOM', iat, atsym(qm_atoms(iat)%spec), 'STO ', 'A', iat, &
            rclas(1,iat) / Ang, rclas(2,iat) / Ang, rclas(3,iat) / Ang
    enddo
    do iat = 1, namm
      write( funit, '(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)' ) &
             'ATOM', naqm+iat, mm_atoms(iat)%atname, mm_atoms(iat)%aaname, &
             'A', mm_atoms(iat)%aanum, rclas(1,naqm+iat) / Ang, &
             rclas(2,naqm+iat) / Ang, rclas(3,naqm+iat) / Ang
    enddo
    write( funit, '(A3,i3)' ) 'END', irestrstep
    call io_close( funit )
  end subroutine write_react_crd

  subroutine write_siesta_input( siestaslabel, na_qm, qm_atoms, coulombtype, &
                                 rcut_qmmm, rcut_density, rcut_ew )
    !! Writes an additional fdf for SIESTA to read as an interface.
    use fdf                , only : fdf_deprecated, fdf_get, leqi, fdf_defined
    use fdf_parse          , only : parsed_line, digest, nnames, names, destroy, labeleq
    use linkatoms          , only : numlink, link_atoms, &
                                    num_resonances, resonance
    use mm_topology        , only : qm_atom_t
    use precision          , only : dp
    use sys                , only : die

    implicit none
    character(len=*) , intent(in) :: siestaslabel
      !! Special filename for files read by siesta.
    integer          , intent(in) :: na_qm
      !! Number of QM atoms.
    type(qm_atom_t)  , intent(in) :: qm_atoms(na_qm)
      !! Atomic coordinates for the QM region.
    character(len=10), intent(in) :: coulombtype
      !! Either ewald or cut-off coulomb interaction types.
    real(dp)         , intent(in) :: rcut_qmmm
      !! Cut-off radius for QM-MM interactions.
    real(dp)         , intent(in) :: rcut_density
      !! Cut-off radius for partial charge-density point proximity.
    real(dp)         , intent(in) :: rcut_ew
      !! Cut-off radius for Ewald direct summations.

    character(len=4)   :: isEnd
    character(len=18)  :: selfname
    character(len=80)  :: acf
    character(len=131) :: lineaux
    character(len=255) :: fileout
    integer            :: iline, ivec, iat, nlines, unit, iunit, ireson, ios, &
                          tmplen
    logical            :: do_smooth
    real(dp)           :: fscale, electr_l, smooth_l

    type(parsed_line), pointer :: pline


    character(len=131), allocatable :: line(:)
    integer           , allocatable :: length(:)


    call io_assign( iunit )

    open( iunit, file = qmmm_files%input )

    iline = 1
    ios   = 0
    do while ( .true. )
      read( iunit, fmt = '(a)', iostat = ios ) lineaux
      if ( ios /= 0 ) exit

      tmplen = len_trim( lineaux )
      if ( tmplen /= 0 ) iline = iline +1
    enddo

    nlines = iline -1
    allocate( line(iline), length(iline) )

    iline = 1
    ios   = 0
    rewind( iunit )
    do while ( .true. )
      read( iunit, fmt = '(a)', iostat = ios ) line(iline)
      if ( ios /= 0 ) exit

      length(iline) = len_trim( line(iline) )
      if ( length(iline) /= 0 ) iline = iline +1
    enddo
    call io_close( iunit )

    ! Format of atomic coordinates
    acf = fdf_get( 'AtomicCoordinatesFormat', 'Ang' )

    fscale = 1.0_dp
    selfname = "write_siesta_input"
    if ( leqi(acf, 'NotScaledCartesianBohr') .or. leqi(acf, 'Bohr') ) then
      write( 6, '(A)' ) selfname//&
        ': Atomic-coordinates input format  =  Cartesian coordinates (in Bohr)'
    else if (leqi(acf, 'NotScaledCartesianAng') .or. leqi(acf, 'Ang' ) ) then
      fscale = 1.0_dp / Ang
      write( 6, '(a)' ) selfname//&
        ': Atomic-coordinates input format  =  Cartesian coordinates (in Ang)'
    else
      write( 6, "(/, 'write_siesta_input: ',73(1h*))" )
      write( 6, '(a)' ) selfname//':    INPUT ERROR'
      write( 6, '(a)' ) selfname//&
        ': You must use one of the following coordinate options:'
      write( 6, '(a)' ) selfname//':     - NotScaledCartesianBohr (or Bohr)'
      write( 6, '(a)' ) selfname//':     - NotScaledCartesianAng (or Ang) '
      write( 6, "('read: ',73(1h*))")
      call die( 'Wrong coordinates unit.' )
    endif

    ! Variables related to the smooth of the MM potential. This are kept here
    ! only for backwards compatibility.
    call fdf_deprecated('SmoothPotential','QMMM.SmoothElectrode')
    call fdf_deprecated('ElectrodeLength','QMMM.SmoothElectrode.Electrode')
    call fdf_deprecated('SmoothLength'   ,'QMMM.SmoothElectrode.Smooth')

    do_smooth = fdf_get( 'SmoothPotential', .false. )
    do_smooth = fdf_get( 'QMMM.SmoothElectrode', do_smooth )
    electr_l  = fdf_get( 'ElectrodeLength', 5.0_dp, 'Bohr' )
    electr_l  = fdf_get( 'QMMM.SmoothElectrode.Electrode', electr_l, 'Bohr' )
    smooth_l  = fdf_get( 'SmoothLength'   , 5.0_dp, 'Bohr' )
    smooth_l  = fdf_get( 'QMMM.SmoothElectrode.Smooth', smooth_l, 'Bohr' )

    do ireson = 1, num_resonances
      fileout = trim( trim(resonance(ireson)%path) // siestaslabel ) // '.fdf'

      call io_assign( unit )
      open( unit, file = fileout)

      iline = 0
      do while ( .true. )
        iline = iline +1
        if ( iline > nlines ) exit

        lineaux = line(iline)
        pline => digest(trim(lineaux))

        if ( pline%ntokens < 1 ) cycle
        if ( nnames(pline) < 1 ) cycle

        if ( labeleq(names(pline,1), 'SystemName') ) then
          write( unit, '(a,2x,a)' ) 'SystemName', siestaslabel
        elseif ( labeleq(names(pline,1), 'SystemLabel') ) then
          write( unit, '(a,2x,a)' ) 'SystemLabel', siestaslabel
        elseif ( labeleq(names(pline,1), 'MD.TypeOfRun') ) then
          write( unit, '(a)' ) 'MD.TypeOfRun qmmm'
        elseif ( labeleq(names(pline,1), 'WriteMDXmol') ) then
          write( unit, '(a)' ) 'WriteMDXmol .false.'
        elseif ( labeleq(names(pline,1), 'NumberOfAtoms') ) then
          write( unit, '(a,3x,I6)' ) 'NumberOfAtoms ', na_qm+numlink
        elseif ( leqi(names(pline,1), "%block")) then
          if ( labeleq(names(pline,2), "LatticeVectors") ) then
            write( unit, '(a)' ) '%block LatticeVectors'
            ivec = 0
            do while (.true.)
              iline   = iline +1
              lineaux = line(iline)

              call destroy(pline)
              pline => digest(trim(lineaux))
              if ( pline%ntokens < 1 ) cycle

              write( unit, '(a)' ) trim( lineaux )
              ivec = ivec+1
              if (ivec == 3) exit
            enddo

            do while (.true.)
              iline   = iline +1
              lineaux = line(iline)
              isEnd = lineaux(1:4)

              if (leqi(isEnd, "%end")) exit
            enddo
            write( unit, '(a,3x,I6)' ) '%endblock LatticeVectors'
          elseif (labeleq(names(pline,2), "AtomicCoordinatesAndAtomicSpecies")) &
            then

            write( unit, '(a)' ) '%block AtomicCoordinatesAndAtomicSpecies'
            do iat = 1, na_qm
              write( unit, '(1X,3F8.3,I4)' ) &
                fscale * qm_atoms(iat)%r(1:3), qm_atoms(iat)%spec
            enddo

            if ( numlink > 0 ) then
              do iat = 1, numlink
                write( unit, '(1X,3F8.3,I4)' ) fscale * &
                  link_atoms(iat)%reson(ireson)%r(1:3), &
                  resonance(ireson)%isa(iat)
              enddo
            endif

            do while (.true.)
              iline   = iline +1
              lineaux = line(iline)
              isEnd = lineaux(1:4)

              if (leqi(isEnd, "%end")) exit
            enddo
            write( unit, '(a)' ) '%endblock AtomicCoordinatesAndAtomicSpecies'

          elseif ( labeleq(names(pline,2), "SoluteAtomTypes") .or. &
                   labeleq(names(pline,2), "QM.AtomTypes")    .or. &
                   labeleq(names(pline,2), "SolventInput")    .or. &
                   labeleq(names(pline,2), "MM.Atoms")        .or. &
                   labeleq(names(pline,2), "CutOffRadius") ) then
            do while (.true.)
              iline   = iline +1
              lineaux = line(iline)
              isEnd = lineaux(1:4)

              if (leqi(isEnd, "%end")) exit
            enddo
          else
            write( unit, '(a)' ) trim(lineaux)

          endif

          ! This labels are also skipped.
        elseif ( labeleq(names(pline,1), 'NumberOfSolventAtoms') .or. &
                 labeleq(names(pline,1), 'NumberOfMMAtoms')      .or. &
                 labeleq(names(pline,1), 'NumberOfSolventAtoms') .or. &
                 labeleq(names(pline,1), 'MM.Cutoff')            .or. &
                 labeleq(names(pline,1), 'Block.Cutoff')         .or. &
                 labeleq(names(pline,1), 'SolventDebug')         .or. &
                 labeleq(names(pline,1), 'WriteIniParLas')       .or. &
                 labeleq(names(pline,1), 'MD.MMxQMsteps')        .or. &
                 labeleq(names(pline,1), 'LaunchSiesta')         .or. &
                 labeleq(names(pline,1), 'MPICommand')           .or. &
                 labeleq(names(pline,1), 'NebListFreq')          .or. &
                 labeleq(names(pline,1), 'CenterQMSystem')       .or. &
                 labeleq(names(pline,1), 'SmoothPotential')      .or. &
                 labeleq(names(pline,1), 'ElectrodeLength')      .or. &
                 labeleq(names(pline,1), 'DefaultLinkSymbol')    .or. &
                 labeleq(names(pline,1), 'MM.ParmFile')          .or. &
                 labeleq(names(pline,1), 'SmoothLength')         .or. &
                 labeleq(names(pline,1), 'QMMM.SmoothElectrode') .or. &
                 labeleq(names(pline,1), 'QMMM.SmoothElectrode.Electrode') .or.&
                 labeleq(names(pline,1), 'QMMM.SmoothElectrode.Smooth') ) then
          ! Do nothing
        else
          write( unit, '(a)' ) trim(lineaux)

        endif

        call destroy(pline)
      enddo

      if (.not. fdf_defined('MD.TypeOfRun') ) &
        write( unit, '(a)') 'MD.TypeOfRun qmmm'
      if (.not. fdf_defined('QMMM.Enabled') ) &
        write( unit, '(a)') 'QMMM.Enabled T'
      if (.not. fdf_defined('QMMM.Driver') ) &
        write( unit, '(a)') 'QMMM.Driver SIESTA_DRIVER'
      if (.not. fdf_defined('QMMM.Driver.Type') ) &
        write( unit, '(a)') 'QMMM.Driver.Type pipe'

      if (.not. fdf_defined('QMMM.Driver.QMRegionFile') ) &
        write( unit, '(a,a)') 'QMMM.Driver.QMRegionFile ', &
                              trim(siestaslabel)//'.coords'
      if (.not. fdf_defined('QMMM.Driver.MMChargeFile') ) &
        write( unit, '(a,a)') 'QMMM.Driver.MMChargeFile ', &
                              trim(siestaslabel)//'.pc'
      if (.not. fdf_defined('QMMM.Driver.ForceOutFile') ) &
        write( unit, '(a,a)') 'QMMM.Driver.ForceOutFile ', &
                              trim(siestaslabel)//'.forces'

      if (.not. fdf_defined('QMMM.Cutoff') ) &
        write( unit, '(a,F14.7,a)') 'QMMM.Cutoff      ', rcut_qmmm   , ' Ang'
      if (.not. fdf_defined('QMMM.DensityCut') ) &
        write( unit, '(a,F14.7,a)') 'QMMM.DensityCut  ', rcut_density, ' Ang'

      if (.not. fdf_defined('QMMM.CoulombType') ) then
        if ( trim(coulombtype) == 'ewald' ) then
          write( unit, '(a)') 'QMMM.CoulombType 1'
          write( unit, '(a,F14.7,a)') &
            'QMMM.Ewald.rcut', rcut_ew, ' Ang'
        else
          write( unit, '(a)') 'QMMM.CoulombType 2'
        endif
      endif

      if ( do_smooth ) then
        write( unit, '(a)') 'QMMM.SmoothElectrode T'
        write( unit, '(a,F14.7,a)') &
          'QMMM.SmoothElectrode.Electrode', electr_l , ' Bohr'
        write( unit, '(a,F14.7,a)') &
          'QMMM.SmoothElectrode.Smooth'   , smooth_l , ' Bohr'
      endif

      call io_close( unit )
    enddo

    deallocate( line, length )
  end subroutine write_siesta_input

  subroutine write_restr( rtot, nrestr )
    !! Writes the restrained reaction coordinate. See harmonic_restraints for
    !! more details.
    use precision , only : dp

    implicit none
    real(dp) , intent(in) :: rtot(nrestr)
      !! Restraints coordinate values.
    integer  , intent(in) :: nrestr
      !! Total number of restraints.

    integer :: unit
    character(len=10) :: n_rest
    character(len=100):: fmt

    if ( nrestr < 1 ) return
    call io_assign( unit )
    open( unit, file = qmmm_files%wrt, form = 'formatted', position = 'append',&
          status = 'unknown' )

    write( n_rest, '(I10)' ) nrestr
    write( fmt   , '(A)'   ) '('//trim(n_rest)//'F10.4,2x,F14.7)'
    write( unit  , trim(fmt) ) rtot(:)

    call io_close( unit )
  end subroutine write_restr

  subroutine qmmm_write_forces( na_u, fa, cfa, ftol, first_step, final, &
                                writef, doingfc, dx )
    use iofa_m    , only : iofa
    use precision , only : dp
    use siesta_cml, only : cml_p, mainXML, cmlStartPropertyList, cmlAddComment,&
                           cmlAddProperty, cmlEndPropertyList
    use units     , only : Ang, eV
    ! Also uses ofc as an external subroutine.

    implicit none
    logical, intent(in) :: first_step
      !! Initial MD step.
    logical, intent(in) :: final
      !! Whether we are in the last MD step.
    logical, intent(in) :: writef
      !! Whether we write the forces to the main output file.
    logical, intent(in) :: doingfc
      !! Whether we are doing a Force Constant matrix calculation.
    integer, intent(in) :: na_u
      !! Number of atoms in unit cell.
    real(dp), intent(in) :: fa(3,na_u)
      !! Forces over atoms.
    real(dp), intent(in) :: cfa(3,na_u)
      !! Forces over atoms, with constraints applied.
    real(dp), intent(in) :: ftol
      !! Maximum force tolerance.
    real(dp), intent(in) :: dx
      !! Atomic displacement for force constant matrix.

    integer  :: ia, ix
    real(dp) :: cfmax, fmax, fres
    real(dp) :: ftot(3), cftot(3)

    call analyze_force( na_u, cfa, cfmax, cftot, fres )
    call analyze_force( na_u,  fa,  fmax,  ftot, fres )

    ! Almost the same forces output whether during simulation
    ! or at the end. Unfortunately not quite, therefore slightly
    ! tortuous logic below. If we are content to change format
    ! of output file slightly, this can be simplified.
    if ( .not. final ) then
      ! print forces to xml every step.
      ! output forces to stdout depending on writef
      if (cml_p) then
        call cmlStartPropertyList(mainXML, title = 'Forces')
        call cmlAddComment( mainXML,"Output: matrix fa(1:3,1:na_u)")
        call cmlAddProperty( xf = mainXML, value = fa*Ang/eV, &
             dictref = 'qmmm:forces', units = 'siestaUnits:evpa')
        call cmlAddProperty( xf = mainXML, value = ftot*Ang/eV, &
             dictref = 'qmmm:ftot', units = 'siestaUnits:evpa')
        call cmlAddProperty( xf = mainXML, value = fmax*Ang/eV, &
             dictref = 'qmmm:fmax', units = 'siestaUnits:evpa')
        call cmlAddProperty( xf = mainXML, value = fres*Ang/eV, &
             dictref = 'qmmm:fres', units = 'siestaUnits:evpa')
        call cmlAddProperty( xf = mainXML, value = cfmax*Ang/eV, &
             dictref = 'qmmm:cfmax', units = 'siestaUnits:evpa')
        call cmlEndPropertyList(mainXML)
      endif

      write(6,'(/,a)') 'qmmm: Atomic forces (eV/Ang):'
      if (writef) then
        write(6,'(i6,3f12.6)')( ia, ( fa(ix,ia) * Ang/eV, ix=1,3 ), ia=1, na_u )
      endif

      call iofa( na_u, fa , 'FA')
      if ( any( fa /= cfa ) ) then
        ! if any one constraint is not the same
        ! as the real forces, the FAC file
        ! will also be created
        call iofa( na_u, cfa , 'FAC' )
      endif

      write( 6, '(40("-"),/,a6,3f12.6)' ) 'Tot',( ftot(ix) * Ang / eV, ix =1,3 )
      write( 6, '(40("-"),/,a6, f12.6)' ) 'Max', fmax * Ang / eV
      write( 6, '(a6,f12.6,a26)' )        'Res', fres * Ang / eV, &
                                       '    sqrt( Sum f_i^2 / 3N )'
      write( 6, '(40("-"),/,a6, f12.6,a)' ) 'Max', cfmax * Ang / eV, &
                                            '    constrained'

      ! Write Force Constant matrix if FC calculation ...
      if ( doingfc ) then
         ! If the istep is the first step, then it
         ! must be the first
         call ofc( fa, dx, na_u, .false., first_step )
         if ( any( fa /= cfa ) ) then
            call ofc( cfa, dx, na_u, .true., first_step )
         endif
      endif

    else ! not final
      ! When finalizing, only print forces if they are sufficiently large.
      if ( fmax > ftol ) then
        write(6,'(/,a)') 'qmmm: Atomic forces (eV/Ang):'
        write(6,'(a,i6,3f12.6)') &
             ('qmmm: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
        write(6,'(a,40("-"),/,a,a6,3f12.6)') &
             'qmmm: ','qmmm: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
        if ( cml_p ) then
          call cmlStartPropertyList( mainXML, title = 'Force Summary' )
          call cmlAddComment( mainXML, "Output: matrix fa(1:3,1:na_u)" )
          call cmlAddProperty( xf = mainXML, value = fa * Ang/eV, &
               dictref = 'qmmm:forces', units = 'siestaUnits:evpa' )
          call cmlAddProperty( xf = mainXML, value =  ftot * Ang/eV, &
               dictref = 'qmmm:ftot', units = 'siestaUnits:evpa')
          call cmlEndPropertyList( mainXML )
        endif !cml_p
      endif
      if ( any(cfa /= fa) ) then
        if ( cfmax > ftol ) then
          write(6,'(/,a)') 'qmmm: Constrained forces (eV/Ang):'
          write(6,'(a,i6,3f12.6)') &
               ('qmmm: ',ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
          write(6,'(a,40("-"),/,a,a4,3f12.6)') &
               'qmmm: ','qmmm: ','Tot',(cftot(ix)*Ang/eV,ix=1,3)
          if (cml_p) then
            call cmlStartPropertyList( mainXML, &
                 title = 'Constrained Force Summary')
            call cmlAddComment( mainXML, &
                 "Output: matrix cfa(1:3,1:na_u)")
            call cmlAddProperty( xf = mainXML, value = cfa * Ang / eV, &
                 dictref = 'qmmm:cforces', units = 'siestaUnits:evpa' )
            call cmlAddProperty(xf = mainXML, value = cftot * Ang / eV, &
                 dictref = 'qmmm:cftot', units = 'siestaUnits:evpa')
            call cmlEndPropertyList( mainXML )
          endif !cml_p
        endif
      endif
    endif ! final

  contains

    subroutine analyze_force( natu, frc, f_max, f_tot, f_res )
        !! Gets the maximum force value, the total forces over the system
        !! and the forces mean squared residual.
      use precision, only: dp

      implicit none
      integer , intent(in) :: natu
        !! Number of Atoms in unit cell.
      real(dp), intent(in) :: frc(3,natu)
        !! Forces over atoms.
      real(dp), intent(out) :: f_max
        !! Maximum force component for any direction
      real(dp), intent(out) :: f_tot(3)
        !! Sum of forces for each cartesian direction.
      real(dp), intent(out) :: f_res
        !! Residual of the forces.

      integer :: iat, icrd

      f_max = 0._dp
      f_tot = 0._dp
      f_res = 0._dp
      do iat = 1, natu
      do icrd = 1, 3
        f_max       = max(f_max, abs(frc(icrd,iat)))
        f_tot(icrd) = f_tot(icrd) + frc(icrd,iat)
        f_res       = f_res + frc(icrd,iat) ** 2
      enddo
      enddo
      f_res = sqrt(f_res / (natu * 3))
    end subroutine analyze_force
  end subroutine qmmm_write_forces

  subroutine qmmm_write_stress_pressure( idyn, final, varcel, target_P,   &
                                         remove_int_press, FreeE, stress, &
                                         kin_stress, mstress, cstress,    &
                                         tstress, cell_vol )
    !! Outputs pressure and stress information.
    use precision, only : dp
    use siesta_cml, only : cml_p, mainXML, cmlStartPropertyList, cmlAddComment,&
    cmlAddProperty, cmlEndPropertyList
    use units, only : kbar, eV, Ang

    implicit none
    logical, intent(in) :: final
      !! Whether we are in the last MD step.
    logical, intent(in) :: remove_int_press
      !! Whether we remove intramolecular pressure.
    logical, intent(in) :: varcel
      !! Whether we are doing a variable cell run.
    real(dp), intent(in) :: target_P
      !! Target pressure (for optimization.).
    integer, intent(in) :: idyn
      !! Type of MD.
    real(dp), intent(in) :: FreeE
      !! Total (static) free energy.
    real(dp), intent(in) :: cell_vol
      !! Cell volume.
    real(dp), intent(in) :: stress(3,3)
    real(dp), intent(in) :: kin_stress(3,3)
    real(dp), intent(in) :: mstress(3,3)
    real(dp), intent(in) :: cstress(3,3)
    real(dp), intent(in) :: tstress(3,3)


    integer :: jx, ix
    real(dp):: Pmol     ! Molecular pressure (discounting Virial term)
    real(dp):: Psol     ! Pressure of "solid"
    real(dp):: Press    ! Pressure
    real(dp):: ps(3,3)  ! Auxiliary array

    if ( .not. final ) then
      ! Write Voigt components of total stress tensor
      ps = stress + kin_stress
      write(6,'(/,a,6f12.2)') 'Stress tensor Voigt[x,y,z,yz,xz,xy] (kbar):', &
           (ps(jx,jx)/kbar,jx=1,3), ps(2,3)/kbar, ps(1,3)/kbar, ps(1,2)/kbar

      Press = - ((ps(1,1) + ps(2,2) + ps(3,3)) / 3.0_dp)
      write(6,"(a,f14.4)") "(Free)E + p*V (eV/cell)", &
           (FreeE + Press*cell_vol) / eV

      if ( remove_int_press ) then
        ps = mstress + kin_stress
        write(6,'(/,a,6f12.2)') &
          'Inter molecular stress Voigt[x,y,z,yz,xz,xy] (kbar):', &
          (ps(jx,jx)/kbar,jx=1,3), ps(2,3)/kbar, ps(1,3)/kbar, ps(1,2)/kbar

        Press = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
        write(6,"(a,f14.4)") "(Free)E + p_inter_molec * V  (eV/cell)", &
             (FreeE + Press*cell_vol) / eV
      endif

      ! This use of the volume is OK, as it is called from state_analysis,
      ! before possibly changing the cell.
      ! Write "target enthalpy" (E + pV, where p is the *target* pressure)
      write(6,"(a,f14.4)") "Target enthalpy (eV/cell)", &
           (FreeE + target_P*cell_vol) / eV

      ! Write stress to CML file always
      if (cml_p) &
        call cmlAddProperty( xf = mainXML, value = stress*Ang**3, &
          dictref = 'siesta:stress', title = 'Stress',            &
          units = 'siestaUnits:evpa3')

      ! Output depends on dynamics option
      if ( (.not. ((idyn == 0) .and. (.not.varcel))) .and. &
           ( (idyn < 6) .or. (idyn > 7) ) ) then

        write(6,'(/,a,3(/,a,3f12.6))') &
              'siesta: Stress tensor (static) (eV/Ang**3):', &
              ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

        Psol = - ((stress(1,1) + stress(2,2) + stress(3,3))/3.0_dp)
        write(6,'(/,a,f20.8,a)')  'siesta: Pressure (static):', &
                                  Psol / kBar, '  kBar'

        if ( cml_p ) &
            call cmlAddProperty( xf = mainXML, value = Psol, &
                 dictref = 'siesta:psol', title = 'Pressure (Static)', &
                 units = 'siestaUnits:kBar')

        write(6,'(/,a,3(/,a,3f12.6))') &
              'siesta: Stress tensor (total) (eV/Ang**3):', &
              ('     ',(tstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

        Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
        write(6,'(/,a,f20.8,a)')  'siesta: Pressure (total):', &
                                  Psol/kBar, '  kBar'
        if ( cml_p ) then
          call cmlAddProperty( xf = mainXML, value = tstress*Ang**3, &
                 dictref = 'siesta:tstress', title = 'Total Stress', &
                 units = 'siestaUnits:evpa3')
          call cmlAddProperty( xf = mainXML, value = Psol, &
                 dictref = 'siesta:tpsol', title = 'Pressure (Total)', &
                 units = 'siestaUnits:kBar')
        endif

        if ( .not. remove_int_press ) return ! Early return.

        ps = mstress
        write(6,'(/,a,3(/,a,3f12.6))') &
              'siesta: Stress tensor (nonmol) (eV/Ang**3):', &
              ('     ',(ps(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

        Psol = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
        write(6,'(/,a,f20.8,a)') 'siesta: Pressure (nonmol):',&
                                   Psol/kBar, '  kBar'
        if ( cml_p ) then
          call cmlAddProperty( xf = mainXML, value = ps*Ang**3, &
            dictref = 'siesta:mstress', title = 'Stress tensor (normal)',&
            units = 'siestaUnits:evpa3')
          call cmlAddProperty( xf = mainXML, value = Psol, &
            dictref = 'siesta:pmol', title = 'Pressure (Nonmol)', &
            units = 'siestaUnits:kBar')
          endif

        ps = mstress + kin_stress
        write(6,'(/,a,3(/,a,3f12.6))') &
             'siesta: Stress tensor (nonmol+kin) (eV/Ang**3):', &
              ('     ',(ps(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

        Psol = - ((ps(1,1)+ps(2,2) +ps(3,3))/3.0_dp)
        write(6,'(/,a,f20.8,a)') &
              'siesta: Pressure (nonmol+kin):', Psol/kBar, '  kBar'

        if (cml_p) then
          call cmlAddProperty( xf = mainXML, value = ps*Ang**3, &
            dictref = 'siesta:tmstress', title = 'Stress tensor (nonmol+kin)',&
            units = 'siestaUnits:evpa3')
          call cmlAddProperty( xf = mainXML, value = Psol, &
            dictref = 'siesta:tpmol', title = 'Pressure (Nonmol+Kin)', &
            units = 'siestaUnits:kBar')
        endif
      endif ! idyn
    else ! final

      ! Print stress tensor unconditionally
      write(6,'(/,a,3(/,a,3f12.6))') &
           'siesta: Stress tensor (static) (eV/Ang**3):', &
          ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

      if ( cml_p ) &
        call cmlAddProperty( xf = mainXML, value = stress*Ang**3/eV, &
             dictref = 'siesta:stress', units = 'siestaUnits:eV_Ang__3')

      ! Print constrained stress tensor if different from unconstrained
      if ( any(cstress /= stress ) ) then
        write(6,'(/,a,3(/,a,3f12.6))') &
           'siesta: Constrained stress tensor (static) (eV/Ang**3):', &
           ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

        if (cml_p) &
          call cmlAddProperty( xf = mainXML, value = cstress*Ang**3/eV, &
               dictref = 'siesta:cstress', &
               units = 'siestaUnits:eV_Ang__3')
      endif

      Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) )/3.0_dp)
      Pmol = - (( mstress(1,1) + mstress(2,2) + mstress(3,3) )/3.0_dp)

      write(6,'(/,a,f18.6,a)') 'siesta: Cell volume =', &
            cell_vol/Ang**3, ' Ang**3'
      write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))') &
           'siesta: Pressure (static):', &
           'siesta: ','Solid',        'Molecule',      '  Units', &
           'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',&
           'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',&
           'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'

      if ( .not. cml_p ) return
      call cmlStartPropertyList( mainXML, title = 'Final Pressure' )
      call cmlAddProperty( xf = mainXML, value = cell_vol/Ang**3,&
           title = 'cell volume', dictref = 'siesta:cellvol', &
           units = 'siestaUnits:Ang__3' )
      call cmlAddProperty( xf = mainXML, value = Psol/kBar, &
           title = 'Pressure of Solid', dictref = 'siesta:pressSol', &
           units = 'siestaUnits:kbar' )
      call cmlAddProperty( xf = mainXML, value = Pmol/kBar, &
           title = 'Pressure of Molecule', dictref = 'siesta:pressMol', &
           units = 'siestaUnits:kbar' )
      call cmlEndPropertyList( mainXML )

    endif !final for stress & pressure
  end subroutine qmmm_write_stress_pressure

  subroutine qmmm_write_positions( has_moved )
    !! Writes coordinates and cell during molecular dynamics.
    use atomlist      , only : elem, iza
    use m_iostruct    , only : write_struct
    use siesta_cml    , only : cml_p, mainXML, cmlAddMolecule, cmlAddLattice
    use siesta_geom   , only : ucell, na_u, isa, xa, cisa
    use siesta_options, only : idyn
    use units         , only : Ang
    use zmatrix       , only : lUseZmatrix, write_canonical_ucell_and_Zmatrix

    implicit none
    logical, intent(in) :: has_moved
      !! Whether the structure has already changed due to forces/stresses.
    character(len=30)   :: zmatname

    call write_struct( ucell, na_u, isa, iza, xa, has_moved )
    if ( lUseZmatrix .and. (idyn == 0) ) then
      zmatname = "OUT.UCELL.ZMATRIX"
      if ( has_moved ) zmatname = "NEXT_ITER.UCELL.ZMATRIX"
      call write_canonical_ucell_and_Zmatrix( filename = trim(zmatname) )
    endif

    ! Note that this information should be written only at the
    ! time of processing the current geometry (in state_init)

    if ( (.not. has_moved) .and. cml_p ) then
      call cmlAddMolecule( xf = mainXML, natoms = na_u, elements = elem, &
                           atomRefs = cisa, coords = xa / Ang )
      call cmlAddLattice( xf = mainXML, cell = ucell / Ang, &
                          units = 'siestaUnits:Ang', dictref = 'siesta:ucell' )
    endif
  end subroutine qmmm_write_positions
end module qmmm_write
