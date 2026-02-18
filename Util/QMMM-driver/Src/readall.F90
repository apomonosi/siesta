module qm_read
  !! Handles the input of QM system information and coordinates.
  implicit none
  public :: read_qm
  public :: readcrd

  private
contains

  subroutine read_qm( na_qm, nesp, qm_atoms, atsym )
    !! Reads the information for the QM subsystem.
    use alloc     , only : re_alloc, de_alloc
    use fdf       , only : fdf_get, leqi
    use fdf       , only : block_fdf, parsed_line
    use fdf       , only : fdf_block, fdf_bline, fdf_bnvalues, fdf_bnintegers,&
                           fdf_bintegers, fdf_bvalues, fdf_bclose, fdf_bnames,&
                           fdf_bnnames
    use mm_topology,only : qm_atom_t
    use precision , only : dp
    use sys       , only : die
    use units     , only : Ang

    implicit none
    integer         , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)  :: nesp
      !! Total number of species.
    type(qm_atom_t) , intent(out) :: qm_atoms(na_qm)
      !! Species index for each QM atom.
    character(len=2), intent(out) :: atsym(nesp)
      !! Atomic symbol for each QM atom species.

    character(len=80) :: acf
    integer           :: iat, isp, iscale
    integer, pointer  :: atnum(:)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    ! Unit format of atomic coordinates.
    iscale = 0
    acf    = 'Ang'
    acf    = fdf_get( 'AtomicCoordinatesFormat', acf )

    if ( leqi( acf, 'NotScaledCartesianBohr' ) .or. leqi(acf,'Bohr') ) then
      iscale = 0
      write( 6, '(a,a)' ) 'read_qm: Atomic-coordinates input format  = ', &
                          '    Cartesian coordinates (in Bohr)'

    elseif ( leqi( acf, 'NotScaledCartesianAng' ) .or. leqi( acf, 'Ang' ) ) then
      iscale = 1
      write( 6, '(a,a)' ) 'qmmm_read: Atomic-coordinates input format  = ', &
                          '    Cartesian coordinates (in Ang)'

    else
      write( 6, "(/,'qmmm_read: ',73(1h*))" )
      write( 6, "('qmmm_read:                  INPUT ERROR')" )
      write( 6, '(2a)' ) 'qmmm_read: You must use one of the following', &
                         'coordinate options:'
      write( 6, '(a)' ) 'qmmm_read:     - NotScaledCartesianBohr (or Bohr)'
      write( 6, '(a)' ) 'qmmm_read:     - NotScaledCartesianAng (or Ang)'
      write( 6, "('qmmm_read: ',73(1h*))" )
      call die( 'qmmm_read: Wrong input units.' )
    endif

    ! Read atomic coordinates and species
    if ( fdf_block( 'AtomicCoordinatesAndAtomicSpecies', bfdf ) ) then

      do iat = 1, na_qm
        if (.not. fdf_bline(bfdf,pline)) &
         call die('qmmm_read: ERROR in AtomicCoordinatesAndAtomicSpecies block')
        if ( fdf_bnvalues(pline) < 3 ) &
          call die( "qmmm_read: problem reading QM region coordinates." )

        qm_atoms(iat)%r(1) = fdf_bvalues(pline,1)
        qm_atoms(iat)%r(2) = fdf_bvalues(pline,2)
        qm_atoms(iat)%r(3) = fdf_bvalues(pline,3)

        if ( fdf_bnintegers(pline) < 1 ) &
          call die( "qmmm_read: problem reading QM region species." )

        qm_atoms(iat)%spec = fdf_bintegers(pline,1)
      enddo

      call fdf_bclose( bfdf )
    else
      call die( "qmmm_read: You must specify the atomic coordinates." )
    endif

    ! Change coordinates to Siesta's units.
    ! iscale = 0 => Do nothing
    ! iscale = 1 => Multiply by 1./0.529177 (Ang --> Bohr)
    if ( iscale == 1 ) then
      do iat = 1, na_qm
        qm_atoms(iat)%r(:) = qm_atoms(iat)%r(:) * Ang
      enddo
    endif

    ! Read the atomic labels
    nullify( atnum )
    call re_alloc( atnum, 1, nesp, 'atnum', 'read_qm')
    if ( fdf_block( 'ChemicalSpeciesLabel', bfdf ) ) then
      do iat = 1, nesp
        if (.not. fdf_bline(bfdf,pline)) &
         call die('qmmm_read: ERROR in ChemicalSpeciesLabel block')

        if ( fdf_bnintegers(pline) < 2 ) &
          call die( "qmmm_read: problem reading QM species info." )
        if ( fdf_bnnames(pline) < 1 ) &
          call die( "qmmm_read: problem reading QM species name." )

        isp        = fdf_bintegers(pline,1)
        atnum(iat) = fdf_bintegers(pline,2)
        atsym(iat) = trim(fdf_bnames(pline,1))

        ! change sign of atnum if there are ghost atoms
        if ( atnum(iat) < 0) atnum(iat) = - atnum(iat)
      enddo
      call fdf_bclose( bfdf )
    else
      call de_alloc( atnum, 'atnum', 'read_qm')
      call die("qmmm_read: You must specify the atomic labels.")
    endif

    ! Assigns atomic number (qm_atoms()%z)
    do iat = 1, na_qm
      qm_atoms(iat)%z = 0
      do isp = 1, nesp
        if ( qm_atoms(iat)%spec == isp ) qm_atoms(iat)%z = atnum(isp)
      enddo

      if ( qm_atoms(iat)%z == 0 ) then
        call de_alloc( atnum, 'atnum', 'read_qm')
        call die( "qmmm_read: There are atomos without atomic number.")
      endif
    enddo
    call de_alloc( atnum, 'atnum', 'read_qm')
  end subroutine read_qm

  subroutine readcrd( na_qm, na_mm, masst, rat, vat, foundcrd, foundvat )
    !! Reads coordinates from a previous run.
    use fdf      , only : fdf_string
    use precision, only : dp
    use sys      , only : die
    use units    , only : Ang

    implicit none
    integer , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer , intent(in)    :: na_mm
      !! Number of MM atoms.
    real(dp), intent(in)    :: masst(na_qm+na_mm)
      !! Atomic masses.
    real(dp), intent(inout) :: rat(3,na_qm+na_mm)
      !! Cartesian positions for atoms.
    real(dp), intent(inout) :: vat(3,na_qm+na_mm)
      !! Cartesian velocities for atoms.
    logical , intent(out)   :: foundcrd
      !! Whether we found a coordinate file.
    logical , intent(out)   :: foundvat
      !! Whether we found velocities within the file.

    integer            :: iat, ltop, natoms, iu, ios
    real(dp)           :: time
    character(len=164) :: sname
    character(len=255) :: fname
    character(len=30)  :: title
    character(len=4)   :: ch

    ! Externals
    external :: io_assign, io_close

    foundcrd = .false.
    sname = fdf_string( 'SystemLabel', 'siesta' )
    fname = trim(sname) // '.mdcrd'

    inquire( file = fname, exist = foundcrd )

    if ( foundcrd ) then
      ios = 0
      call io_assign( iu )
      open( iu, file = fname, status = 'old' )
      read( iu, iostat = ios, fmt = '(a30)' ) title
      call readcrd_check_ios( ios )

      ch = title(1:4)
      if ( ch == 'File' ) then
        foundcrd = .false.
      else
        read( iu, *, iostat = ios ) natoms, time
        call readcrd_check_ios( ios )

        if ( natoms /= (na_qm + na_mm) ) &
          call die( 'readcrd: Wrong number of atoms!!!' )


        ! AMBER rst files have the following structure, assumming an odd
        ! number of atoms:
        !
        ! r11 r12 r13 r21 r22 r23
        ! ...
        ! rn1 rn2 rn3
        ! v11 v12 v13 v21 v22 v23
        ! ...
        ! vn1 vn2 vn3
        !
        ! We should have special care whether the system has an odd or even
        ! number of atoms to read things properly.

        ltop = natoms
        if ( mod(natoms, 2) == 0 ) ltop = ltop -1
        do iat = 1, natoms, 2
          read( iu, iostat = ios, fmt = '(6f12.7)') &
            rat(1,iat)  , rat(2,iat)  , rat(3,iat)  ,&
            rat(1,iat+1), rat(2,iat+1), rat(3,iat+1)
          call readcrd_check_ios( ios )
        enddo

        if ( mod(natoms, 2) == 0 ) then
          read( iu, iostat = ios, fmt = '(3f12.7)') &
            rat(1,natoms), rat(2,natoms), rat(3,natoms)
          call readcrd_check_ios( ios )
        endif

        do iat = 1, natoms, 2
          read( iu, iostat = ios, fmt = '(6f12.7)') &
            vat(1,iat)  , vat(2,iat)  , vat(3,iat)  ,&
            vat(1,iat+1), vat(2,iat+1), vat(3,iat+1)
          call readcrd_check_ios( ios )
        enddo

        if ( mod(natoms, 2) == 0 ) then
          read( iu, iostat = ios, fmt = '(3f12.7)') &
            vat(1,natoms), vat(2,natoms), vat(3,natoms)
          call readcrd_check_ios( ios )
        endif

        ! Unit corrections.
        rat( 1:3, 1:natoms ) = rat( 1:3, 1:natoms ) * Ang
        vat( 1:3, 1:natoms ) = vat( 1:3, 1:natoms ) * Ang * 0.020455_dp

        write( 6, '(/a)' ) &
          'readcrd: Reading Coordinates and Velocities from CRD file.'

        ! Change velocities if using deuterium.
        do iat = 1, na_qm + na_mm
          if ( abs(masst(iat) - 2.0_dp) > 1e-14_dp ) &
            vat(1:3,iat) = vat(1:3,iat) / sqrt(2.0_dp)
        enddo
      endif
      call io_close( iu )
    endif !foundcrd

    ! Check if velocities are non-zero
    do iat = 1,na_qm+na_mm
    do iu  = 1, 3
      if ( abs(vat(iu,iat)) > 0.0_dp ) foundvat = .true.
    enddo
    enddo
  end subroutine readcrd

  subroutine readcrd_check_ios( ios )
    !! Dies if there are errors within input read.
    use sys, only : die
    implicit none
    integer, intent(in) :: ios
      !! Result of iostat when reading.
    if ( ios /= 0 ) call die('readcrd: problem reading coordinate file.' )
  end subroutine readcrd_check_ios

end module qm_read