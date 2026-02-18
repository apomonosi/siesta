module mm_assign_m
  !! This module handles the assignation of classical parameters
  !! for the MM regions.
  ! This is basically a reformat of the old "solv_assign" subroutine.
  use mm_topology, only : fftopology, atom_connect_t, mm_atom_t, qm_atom_t
  use precision  , only : dp

  implicit none
  public :: mm_assign
  public :: mm_atom_to_species

  private
  integer, parameter :: MAX_NUM_RES = 20000
    !! Maximum number of residues.
  integer, parameter :: MAX_TYPE_RES = 200
    !! Maximum number of different residue types.
  integer, parameter :: MAX_CONNEC = 1000
    !! Maximum number of connectivities between MM atoms.

contains
  subroutine mm_assign( na_qm, na_mm, nroaa, mm_atoms, mm_connectivity,&
                        qm_atoms, rcut_qmmm, rcut_qm, rcut_mm, sfc,    &
                        radbloqmmm, atsinres, mmcell, lattice_type, mm_top )
    !! Assigns MM parameters.
    use graphite_m, only : find_graphite_layers
    use fdf       , only : fdf_get, fdf_deprecated, fdf_defined
    use fdf       , only : block_fdf, parsed_line
    use fdf       , only : fdf_block, fdf_bline, fdf_bintegers,&
                           fdf_bnames, fdf_bvalues, fdf_bclose
    use qmmm_pbc  , only : get_lattice_type
    use sys       , only : die
    use units     , only : Ang

    implicit none
    integer          , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer          , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer          , intent(out)   :: nroaa
      !! Number of classical residues.

    real(dp)         , intent(out)   :: rcut_qmmm
      !! Distance cut-off for QM-MM interactions.
    real(dp)         , intent(out)   :: rcut_mm
      !! Distance cut-off for MM-MM interactions.
    real(dp)         , intent(out)   :: rcut_qm
      !! Proximity cut-off between points of the density grid and MM positions.
    real(dp)         , intent(out)   :: sfc
      !! Smoothing function cut-off. It is technically constant...
    real(dp)         , intent(out)   :: radbloqmmm
      !! Freeze all MM atoms beyond this distance from the QM region.

    integer          , intent(out)   :: atsinres(MAX_NUM_RES)
      !! Number of atoms in a given residue.
    real(dp)         , intent(inout) :: mmcell(3,3)
      !! Periodic cell vectors.
    character(len=1) , intent(out)   :: lattice_type
      !! Type of PBC lattice.

    type(qm_atom_t)     , intent(inout) :: qm_atoms(na_qm)
    type(mm_atom_t)     , intent(inout) :: mm_atoms(na_mm)
    type(fftopology)    , intent(inout) :: mm_top
    type(atom_connect_t), intent(inout) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.

    character(len=1) :: ch1
    character(len=4) :: ch4
    character(len=20):: blockname
    integer          :: iat, jat, ires, jres, nres, iconn, ncon, ivec, famber
    logical          :: debug
    real(dp)         :: rcut, xmin, xmax
    character(len=4), allocatable :: resname(:), atnamea(:,:), attypea(:,:), &
                                     aanamea(:)
    integer         , allocatable :: atnu(:,:), nataa(:,:), con(:,:)   , &
                                     con2(:,:), atxres(:) , atomsxaa(:), &
                                     resnum(:)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    ! Externals
    real(dp) , external :: volcel
    external :: io_close

    ! Initialize variables.
    nroaa = 0
    sfc   = 2.0_dp ; ncon   = 0

    rcut_mm = 100.0_dp  ; rcut_qmmm  = 100.0_dp
    rcut_qm = 1.E-06_dp ; radbloqmmm = 100.0_dp

    debug = fdf_get( 'SolventDebug' , .false. )

    if ( debug ) then
      write( 6, '(/,a,70(1h=))' ) 'mm_assign: '
      write( 6, "(a)" ) 'mm_assign: Debug mode enabled.'
      write( 6, "(a)" ) 'mm_assign: Remember to check input and amber.parm.'
    endif

    ! Read MM coordinates by atom
    allocate( resnum(na_mm) )
    blockname = "MM.Atoms"
    if ( fdf_defined('SolventInput') ) then
      call fdf_deprecated('SolventInput', 'MM.Atoms')
      blockname = 'SolventInput'
    endif

    if ( fdf_block(trim(blockname),bfdf) ) then
      do iat = 1, na_mm
        if (.not. fdf_bline(bfdf,pline)) &
          call die( 'mm_assign: Problem while reading MM coordinates.')

        ! We do need a regular read. If not, atomnames that have a # will
        ! be parsed as comments. And in any case, since we are following the PDB
        ! format, it should be character-centric and not token-centric.
        read( pline%line, fmt = '(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)' ) &
              ch4, jat, mm_atoms(iat)%atname, mm_atoms(iat)%aaname, ch1, &
              resnum(iat), mm_atoms(iat)%r(1:3)
      enddo
      call fdf_bclose( bfdf )
    else
      write( 6, "(a)" ) 'mm_assign: SolventInput block not found.'
      call die( "mm_assign: You must specify the MM atom coordinates." )
    endif

    ! Change coordinates to SIESTA units
    do iat = 1, na_mm
      mm_atoms(iat)%r(1:3) = mm_atoms(iat)%r(1:3) * Ang
    enddo

    ! Read cutoff radii for different sections of the code.

    call fdf_deprecated('CutoffRadius', 'X.Cutoff')
    if ( fdf_block( 'CutOffRadius', bfdf ) ) then
      if ( fdf_bline(bfdf,pline) ) rcut_qm    = fdf_bvalues(pline, 1)
      if ( fdf_bline(bfdf,pline) ) rcut_qmmm  = fdf_bvalues(pline, 1)
      if ( fdf_bline(bfdf,pline) ) rcut_mm    = fdf_bvalues(pline, 1)
      if ( fdf_bline(bfdf,pline) ) radbloqmmm = fdf_bvalues(pline, 1)
      call fdf_bclose( bfdf )
    endif
    rcut_qm    = fdf_get( 'QMMM.DensityCut', rcut_qm   , 'Bohr' )
    rcut_mm    = fdf_get( 'MM.Cutoff'   , rcut_mm   , 'Ang' )
    rcut_qmmm  = fdf_get( 'QMMM.Cutoff' , rcut_qmmm , 'Ang' )
    radbloqmmm = fdf_get( 'Block.Cutoff', radbloqmmm, 'Ang' )

    rcut = max( rcut_mm, rcut_qmmm )

    ! Read lattice vectors for periodic boundary conditions.
    ! We create a cubic cell if there is no previous cell.
    if ( abs(volcel( mmcell )) < 1.0e-8_dp ) then
      do ivec = 1, 3
        xmin = huge( 1.0_dp )
        xmax = -xmin

        do iat = 1, na_mm
          xmin = min( xmin, mm_atoms(iat)%r(ivec) - rcut )
          xmax = max( xmax, mm_atoms(iat)%r(ivec) + rcut )
        enddo

        mmcell(ivec,ivec) = xmax - xmin
      enddo
    endif

    lattice_type = get_lattice_type( mmcell )

    ! Assigns the number and indices of residues.
    ires = 1
    mm_atoms(1)%aanum = ires
    do iat = 2, na_mm
      if ( resnum(iat) == resnum(iat-1) ) then
        mm_atoms(iat)%aanum = mm_atoms(iat-1)%aanum
      else
        ires = ires +1
        mm_atoms(iat)%aanum = ires
      endif
    enddo
    deallocate( resnum )
    nroaa = mm_atoms(na_mm)%aanum

    if ( nroaa == 0 ) &
      call die( "mm_assign: Number of residues cannot be zero." )

    if ( na_mm == 0 ) nroaa = 0

    ! Read QM atom types
    if ( na_qm /= 0 ) then
      blockname = "QM.AtomTypes"
      if ( fdf_defined('SoluteAtomTypes') ) then
        call fdf_deprecated('SoluteAtomTypes', 'QM.AtomTypes')
        blockname = 'SoluteAtomTypes'
      endif
      if ( fdf_block( trim(blockname), bfdf ) ) then
        do iat = 1, na_qm
          if (.not. fdf_bline(bfdf,pline)) &
            call die( 'mm_assign: Not enough QM atom types.')
          read( pline%line, * ) ch4

          qm_atoms(iat)%attype = ch4
        enddo
        call fdf_bclose( bfdf )
      else
        write( 6, "(a)" ) 'mm_assign: QM.AtomTypes block not found.'
        call die('mm_assign: You must specify QM atom types.')
      endif
    endif

    ! Check the values for cut-off radii
    if ( rcut_qm < 1.e-8_dp ) &
      call die( 'mm_assign: QM grid cut-off radius too close to zero.' )
    if ( rcut_qmmm < 1.e-8_dp ) &
      call die( 'mm_assign: QM-MM cut-off radius too close to zero' )

    ! External solvent connectivity block
    allocate( con2(2, MAX_CONNEC) )
    blockname = "MM.Connectivity"
    if ( fdf_defined('SolventConnectivity') ) then
      call fdf_deprecated('SolventConnectivity', 'MM.Connectivity')
      blockname = 'SolventConnectivity'
    endif
    if ( fdf_block( trim(blockname), bfdf ) ) then
      iconn = 1

      do while ( fdf_bline(bfdf,pline) )
        if ( iconn > MAX_CONNEC ) &
          call die('mm_assign: Amount of MM connectivities > MAX_CONNEC.')

        con2(1,iconn) = fdf_bintegers(pline, 1)
        con2(2,iconn) = fdf_bintegers(pline, 2)

        iconn = iconn +1
      enddo
      call fdf_bclose( bfdf )

      ncon = iconn -1

      allocate( con( 2, ncon ) )

      con( 1:2, 1:ncon ) = con2( 1:2, 1:ncon )

      if ( ncon /= 0 ) &
        write( 6, '(a)' ) 'mm_assign: Read new connectivities block.'
    else
      write( 6, '(a)' ) 'mm_assign: SolventConnectivity block not found.'
      ncon = 0
      allocate( con(2, ncon+1) )
    endif
    deallocate( con2 )
    if ( debug ) write( 6, '(a)' ) 'mm_assign: Done reading from input.'

    ! Now we need the parameters. First we check if the amber.parm file exists.
    call open_amber_parmfile( famber )

    ! If found, we call the subroutine that reads the bond, angle, dihedral
    ! and improper angle parameters.
    call amber_union_parms( mm_top, famber )
    if ( debug ) write( 6, '(a)' ) 'mm_assign: Read conectivity parameters.'

    ! Now we start the real assignation, according to each atom.
    allocate( atnu(nroaa,100), atnamea(nroaa,100), resname(nroaa), &
              atxres(nroaa) )

    ! Gets the amount of atoms per residue.
    allocate( atomsxaa(MAX_TYPE_RES), aanamea(MAX_TYPE_RES) )
    call atxaa( nres, aanamea, atomsxaa, famber )

    jat = 1 ; atxres = 0
    do ires = 1, nroaa
      resname(ires) = mm_atoms(jat)%aaname

      do jres = 1, nres
        if ( resname(ires) == aanamea(jres) ) atxres(ires) = atomsxaa(jres)
      enddo
      do iat = 1, atxres(ires)
        atnu(ires,iat)    = jat
        atnamea(ires,iat) = mm_atoms(jat)%atname

        if ( mm_atoms(atnu(ires,iat))%aanum /= &
             mm_atoms(atnu(ires,1))%aanum ) then
          write( 6, * ) 'mm_assign: Atom missing in residue: ', &
                        resname(ires), ires
          call die( 'mm_assign: MM atoms missing.' )
        endif

        jat = jat +1
      enddo

      if ( atxres(ires) == 0 ) then
        write( 6, * ) 'mm_assign: Wrong residue name  :', resname(ires), ires
        call die('mm_assign:  Wrong residue name(s).' )
      endif
    enddo
    deallocate( atomsxaa, aanamea )
    if ( debug ) write( 6, '(a)' ) 'mm_assign: Done with variable exchanges.'

    ! The following call extracts AMBER classical monoatomic parameters.
    allocate( attypea(nroaa,100), nataa(nroaa,100) )
    call paramats( nroaa, resname, atnamea, attypea, &
                   nataa, na_mm, atxres, mm_atoms, famber )
    if ( debug ) write( 6, '(a)' ) 'mm_assign: Called paramats.'

    ! And here we obtain Lennard-Jones parameters.
    call read_lj( nroaa, attypea, na_qm, na_mm, qm_atoms, mm_atoms, atxres, &
                  famber )
    if ( debug ) write( 6, '(a)' ) 'mm_assign: Called lj().'

    ! We then assign first neighbours according to connectivity.
    call amber_first_ngb( na_mm, nataa, nroaa, atxres, atnu, resname,    &
                          mm_connectivity, atnamea, ncon, con, mm_atoms, &
                          mmcell, lattice_type, famber )
    if ( debug ) then
      write( 6, '(a)' ) 'mm_assign: Done with first neighbours.'
      write( 6, '(a)' ) 'mm_assign: Connectivity matrix:'

      do iat = 1, na_mm
        write( 6, "('mm_assign: ',7I6)" ) iat, mm_connectivity(iat)%bond_at(1:6)
      enddo
    endif

    ! Here we calculate all connectivity magnitudes: bonds, angles, dihedrals
    ! and improper angles.
    call bon_ang_dih_imp( na_mm, mm_connectivity, atnamea, nroaa, &
                          atxres, atnu, resname, famber )
    if ( debug ) &
      write( 6, '(a)' ) 'mm_assign: Calculated bonds, angles, dihedrals.'

    ! We do not need the parm file anymore.
    call io_close( famber )

    do ires = 1, nroaa
      if ( resname(ires) == 'GRAP' ) then
        call find_graphite_layers( na_mm, mm_atoms, mm_connectivity )
        exit
      endif
    enddo

    ! We convert the name WAT to HOH, for internal conventions.
    ! Then, we verify that water atoms are in the correct order (oxygen first).
    do iat = 1, na_mm
      if ( mm_atoms(iat)%aaname == 'WAT' ) mm_atoms(iat)%aaname = 'HOH'
    enddo
    do ires = 1, nroaa
      if ( resname(ires) == 'WAT' ) resname(ires) = 'HOH'
    enddo
    do ires = 1, nroaa
      if ( .not. (resname(ires) == 'HOH') ) cycle

      do iat = 1, atxres(ires)
        if (.not. (atnamea(ires,iat) == 'O') ) cycle
        if ( iat == 1) cycle

        write( 6, * ) 'mm_assign: Wrong atom order in water residue ', ires
        call die( 'Wrong atom ordering in water molecules' )
      enddo
    enddo

    ! Sets up atsinres array, the number of atoms in each residue.
    atsinres(:) = 0
    if ( nroaa > MAX_NUM_RES ) &
      call die( 'mm_assign: increase atsinres vector dimension (MAX_NUM_RES).' )
    atsinres(1:nroaa) = atxres(1:nroaa)

    deallocate( atnamea, atnu, resname, atxres )
    deallocate( attypea, nataa )

    ! Small verification of parameters.
    do iat = 1, na_qm
      if ( (qm_atoms(iat)%attype == 'HO') .or. &
           (qm_atoms(iat)%attype == 'HW') ) cycle
      if ( (abs(qm_atoms(iat)%lj_Rm) > 0.0_dp)  .and. &
           (abs(qm_atoms(iat)%lj_Em) > 0.0_dp) ) cycle

      write( 6, '(a,i6)' ) 'mm_assign: Null LJ parameter for QM atom ', iat
      write( 6, '(a,f14.7,a,f14.7)' ) 'mm_assign: Rm = ', qm_atoms(iat)%lj_Rm, &
                                      ', Em = ', qm_atoms(iat)%lj_Em
    enddo

    do iat = 1, na_mm
      if ( (mm_atoms(iat)%attype == 'HO') .or. &
           (mm_atoms(iat)%attype == 'HW') ) cycle
      if ( (abs(mm_atoms(iat)%lj_Rm) > 0.0_dp)  .and. &
           (abs(mm_atoms(iat)%lj_Em) > 0.0_dp) ) cycle

      write( 6, '(a,i6)' ) 'mm_assign: Null LJ parameter for MM atom ', iat
      write( 6, '(a,f14.7,a,f14.7)' ) 'mm_assign: Rm = ', mm_atoms(iat)%lj_Rm, &
                                      ', Em = ', mm_atoms(iat)%lj_Em
    enddo

    if( debug ) then
      write( 6, '(a)' ) 'mm_assign: Exiting routine successfully.'
      write( 6, '(a,70(1h=),/)' ) 'mm_assign: '
    endif
  end subroutine mm_assign

  subroutine open_amber_parmfile( funit )
    !! Searches for amber parm file and assigns it to a unit.
    use fdf, only : fdf_get
    use sys, only : die

    implicit none
    integer, intent(out) :: funit
      !! File unit for amber parm file.

    integer :: len, status
    logical :: foundamber

    character(len=255) :: fname
    character(len=:), allocatable :: env_path

    external :: io_assign

    fname = "amber.parm"

    ! We first check the environment.
    call get_environment_variable( "SIESTA_MM_PARM_FILE", length = len, &
                                   status = status, trim_name = .true. )
    if ( (status == 0) .and. (len > 0) ) then
      allocate( character(len) :: env_path )
      call get_environment_variable( "SIESTA_MM_PARM_FILE", value = env_path &
                                     , trim_name = .true. )

      fname = trim(env_path)
      deallocate( env_path )
    endif

    fname = fdf_get( "MM.ParmFile", fname )

    foundamber = .false.
    inquire( file = trim(fname), exist = foundamber )
    if ( .not. foundamber ) then
      call die( "mm_assign: 'amber.parm' file not found." )
    else
      write(6,'(A)') "mm_assign: Using MM parameter file at "//trim(fname)//"."

      call io_assign( funit )
      open( file = trim(fname), unit = funit )
    endif
  end subroutine open_amber_parmfile

  subroutine check_ios_mm( ios, bname, iopt, iidx_i )
    use sys, only : die

    implicit none
    integer, intent(in) :: ios
      !! File I/O status when reading.
    integer, intent(in) :: bname
      !! Blockname to add to message.
    integer, intent(in), optional :: iopt
      !! Additional options for printing.
    integer, intent(in), optional :: iidx_i
      !! Index to print in message if necessary.

    character(len=255) :: msg
    character(len=5)   :: iidx, ios_s

    if ( ios == 0 ) return
    msg = 'mm_assign: ERROR:'

    ! Convert index to string in order to concat later.
    write( ios_s, '(I5)' ) ios
    if ( present( iidx_i ) ) write( iidx, '(I5)' ) iidx_i

    select case( bname )
    case (1) !! Classical atom information.
      msg = trim(msg)//' Problems reading residues block in amber.parm.'

      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of residues not found.'
        if ( iopt == 2 ) msg = trim(msg)//' iaas = '//iidx//'.'
        if ( iopt == 3 ) msg = trim(msg)//' iat = ' //iidx//'.'
      endif

    case (2) !! Lennard-Jones
      msg = trim(msg)//' Problems reading Lennard-Jones block in amber.parm.'

      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of LJ types not found.'
        if ( iopt == 2 ) msg = trim(msg)//' LJtype = '//iidx//'.'
      endif

    case (3) !! Connectivity and neighbours
      msg = trim(msg)//' Problems reading connectivity block in amber.parm.'

      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of residues not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read resname'// &
                                 ' and bonds. ResID = '//iidx//'.'
        if ( iopt == 3 ) msg = trim(msg)//' Attempted to read ID'// &
                                 ' of neighbours. ResID = '//iidx//'.'
      endif

    case (4) !! Improper torsions
      msg = trim(msg)//' Problems reading impropers block in amber.parm.'

      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of residues not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read resname'// &
                                 ' and imps per res. ResID = '//iidx//'.'
        if ( iopt == 3 ) msg = trim(msg)//' Attempted to read ID'// &
                                 ' of atoms in torsion. ResID = '//iidx//'.'
      endif

    case (5) !! Residue names and atoms.
      msg = trim(msg)//' Problems reading residues block in amber.parm.'
      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of residues not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read resname'// &
                                 ' and number of atoms. ResID = '//iidx//'.'
        if ( iopt == 3 ) msg = trim(msg)//' Attempted to read atom'// &
                                 ' data. ResID = '//iidx//'.'
      endif

    case (6) !! Bond data.
      msg = trim(msg)//' Problems reading bonds block in amber.parm.'
      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of bonds not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read bond data.'// &
                                 ' Bond ID = '//iidx//'.'
      endif

    case (7) !! Angle data.
      msg = trim(msg)//' Problems reading angles block in amber.parm.'
      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of angles not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read angle data.'// &
                                 ' Angle ID = '//iidx//'.'
      endif

    case (8) !! Dihedral data.
      msg = trim(msg)//' Problems reading dihedrals block in amber.parm.'
      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of dihedrals not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read dihedral '// &
                                 'angle data. Dihedral ID = '//iidx//'.'
      endif

    case (9) !! Improper torsions data.
      msg = trim(msg)//' Problems reading improper torsions block in amber.parm.'
      if ( present( iopt ) ) then
        if ( iopt == 1 ) msg = trim(msg)//' Number of torsions not found.'
        if ( iopt == 2 ) msg = trim(msg)//' Attempted to read improper '// &
                                 'torsions data. Improper ID = '//iidx//'.'
      endif

    case default
      msg = trim(msg)//' Problems with parameter file reading.'
    end select

    msg = trim(msg)//' IO status = '//ios_s//'.'

    call die( msg )
  end subroutine check_ios_mm

  subroutine paramats( nroaa, resname, atnamea, attypea, &
                       nataa, na_mm, atxres, mm_atoms, ui )
    !! Returns the atomic types and charges.
    use alloc    , only : re_alloc, de_alloc
    use sys      , only : die

    implicit none
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)  :: nroaa
      !! Number of classical residues.
    character(len=4), intent(in)  :: resname(nroaa)
      !! Classical residue names.
    character(len=4), intent(in)  :: atnamea(nroaa,100)
      !! Atom names within each residue.
    character(len=4), intent(out) :: attypea(nroaa,100)
      !! Atom types within each residue.
    integer         , intent(out) :: nataa(nroaa,100)
      !! Atoms index in each residue.
    integer         , intent(out) :: atxres(nroaa)
      !! Number of atoms in each residue.
    type(mm_atom_t) , intent(inout) :: mm_atoms(na_mm)
      !! MM atom types.
    integer         , intent(in)    :: ui
      !! File unit for amber parm file.

    character(len=12) :: option
    integer           :: n1, n2, n3, pnaas, ios, iaas, iat, ires, jat
    logical           :: search
    real(dp), pointer :: qaa(:,:)

    ! These data structures are used to simplify I/O operations.
    type temp_atom
      character(len=4) :: name   = ""
      character(len=4) :: type   = ""
      integer          :: idx    = 0
      integer          :: z      = 0
      real(dp)         :: chrg   = 0.0_dp
    end type temp_atom

    type temp_residue
      integer          :: natoms = 0
      character(len=4) :: name   = ""
      type(temp_atom), allocatable :: atoms(:)
    endtype temp_residue
    type(temp_residue), allocatable :: ff_residues(:)

    rewind( ui )

    ! Reads the number of atoms and the number of residues.
    ! The first line has both magnitudes; then reads each aminoacid and their
    ! variables.

    search = .true.
    ios    = 0
    do while ( search )
      read( ui, *, iostat = ios ) option
      if ( ios /= 0 ) call check_ios_mm( ios, 1 )

      if ( option == 'residues' ) then
        read( ui, *, iostat = ios ) pnaas
        if ( ios /= 0 ) call check_ios_mm( ios, 1, 1 )

        allocate( ff_residues(pnaas) )
        do iaas = 1, pnaas
          read( ui, *, iostat = ios ) ff_residues(iaas)%name, &
                                      ff_residues(iaas)%natoms
          if ( ios /= 0 ) call check_ios_mm( ios, 1, 2, iaas )

          allocate( ff_residues(iaas)%atoms(ff_residues(iaas)%natoms) )

          do iat = 1, ff_residues(iaas)%natoms
            read( ui, *, iostat = ios ) ff_residues(iaas)%atoms(iat)%name, &
                                        ff_residues(iaas)%atoms(iat)%type, &
                                        n1, n2, n3, &
                                        ff_residues(iaas)%atoms(iat)%idx , &
                                        ff_residues(iaas)%atoms(iat)%z   , &
                                        ff_residues(iaas)%atoms(iat)%chrg
            if ( ios /= 0 ) call check_ios_mm( ios, 1, 3, iat )
          enddo
        enddo

        search = .false.
      endif
    enddo

    ! Assigns charges and atom types according to the AMBER forcefield.
    nullify( qaa )
    call re_alloc( qaa, 1, nroaa, 1 , 100, 'qaa', 'paramats' )
    qaa = 100.0_dp ! Used as an internal check.

    do iaas = 1, pnaas
      do ires = 1, nroaa
        if ( .not. ( trim(resname(ires)) == &
                     trim(ff_residues(iaas)%name)) ) cycle

        do iat =1, ff_residues(iaas)%natoms
        do jat =1, ff_residues(iaas)%natoms
          if ( .not. ( trim(atnamea(ires,iat)) == &
                       trim(ff_residues(iaas)%atoms(jat)%name)) )&
            cycle
          qaa(ires,iat)     = ff_residues(iaas)%atoms(jat)%chrg
          attypea(ires,iat) = ff_residues(iaas)%atoms(jat)%type
          nataa(ires,iat)   = ff_residues(iaas)%atoms(jat)%idx
        enddo
        enddo
      enddo
      deallocate( ff_residues(iaas)%atoms )
    enddo

    jat = 1
    do ires = 1, nroaa
    do iat  = 1, atxres(ires)
      mm_atoms(jat)%attype = attypea(ires,iat)
      jat = jat +1
    enddo
    enddo

    jat = 1
    do ires = 1, nroaa
    do iat  = 1, atxres(ires)
      mm_atoms(jat)%pc = qaa(ires,iat)

      if ( abs( mm_atoms(jat)%pc - 100.0_dp ) < 1e-14_dp ) then
        write( 6, '(a,i5)' ) 'paramats: Wrong atom name:  ', jat
        call die( 'Paramats: Error while assigning charges.' )
      endif

      jat = jat +1
    enddo
    enddo

    call de_alloc( qaa, 'qaa', 'paramats' )
  end subroutine paramats

  subroutine read_lj( nroaa, attypea, na_qm, na_mm, qm_atoms, mm_atoms, atxres,&
                      ui )
    !! Assigns Lennard-Jones parameters according to the atom types in attypea.
    use alloc    , only : re_alloc, de_alloc
    use mm_units , only : kcal_mol
    use sys      , only : die
    use units    , only : Ang

    implicit none
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)  :: nroaa
      !! Number of classical residues.
    character(len=4), intent(in)  :: attypea(nroaa,100)
      !! Atom types within each residue.
    integer         , intent(in)  :: atxres(nroaa)
      !! Number of atoms in each residue.
    type(qm_atom_t) , intent(inout)  :: qm_atoms(na_qm)
      !! QM atoms.
    type(mm_atom_t) , intent(inout)  :: mm_atoms(na_mm)
      !! MM atoms.
    integer         , intent(in)  :: ui
      !! File unit for amber parm file.

    integer           :: ilj, nlj, ios, ires, iat
    logical           :: search
    character(len=12) :: option

    real(dp), pointer :: Ema(:,:), Rma(:,:)
    type temp_lj
      real(dp)         :: Em
      real(dp)         :: Rm
      character(len=4) :: type
    end type temp_lj
    type(temp_lj), allocatable :: ff_ljparm(:)

    allocate( ff_ljparm(MAX_TYPE_RES) )

    ! Opens and reads the parameter files.
    rewind( ui )

    ios    = 0
    search = .true.
    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 2 )

      if ( option == 'ljs' ) then
        read( ui, *, iostat = ios ) nlj
        call check_ios_mm( ios, 2, 1 )

        if ( nlj > MAX_TYPE_RES ) &
          call die( 'lj: LJ parameters must not exeed MAX_TYPE_RES.' )
        do ilj = 1, nlj
          read ( ui, *, iostat = ios ) ff_ljparm(ilj)%type, &
                    ff_ljparm(ilj)%Rm, ff_ljparm(ilj)%Em
          call check_ios_mm( ios, 2, 2, ilj )
        enddo

        search = .false.
      endif
    enddo

    ! Converts units.
    do ilj = 1, nlj
      ff_ljparm(ilj)%Rm = ff_ljparm(ilj)%Rm * 2.0_dp * Ang / &
                          ( (2.0_dp)**(1.0_dp/6.0_dp) )
      ff_ljparm(ilj)%Em = 0.5_dp * ff_ljparm(ilj)%Em / kcal_mol
    enddo

    ! Assigns the LJ parameters corresponding to MM atom types.
    nullify( Ema, Rma )
    call re_alloc( Ema, 1, nroaa, 1, 100, 'Ema', 'read_lj'  )
    call re_alloc( Rma, 1, nroaa, 1, 100, 'Rma', 'read_lj'  )

    do ires = 1, nroaa
    do iat  = 1, atxres(ires)
      if ( (attypea(ires,iat) == 'C' ) .or.  &
           (attypea(ires,iat) == 'CA') .or. (attypea(ires,iat) == 'CM') .or. &
           (attypea(ires,iat) == 'CC') .or. (attypea(ires,iat) == 'CV') .or. &
           (attypea(ires,iat) == 'CW') .or. (attypea(ires,iat) == 'CR') .or. &
           (attypea(ires,iat) == 'CB') .or. (attypea(ires,iat) == 'C*') .or. &
           (attypea(ires,iat) == 'CN') .or. (attypea(ires,iat) == 'CK') .or. &
           (attypea(ires,iat) == 'CQ') .or. (attypea(ires,iat) == 'CX') .or. &
           (attypea(ires,iat) == 'CY') .or. (attypea(ires,iat) == 'CD') ) then

        do ilj = 1, nlj
          if ( .not. (ff_ljparm(ilj)%type == 'C') ) cycle
          Rma(ires,iat) = ff_ljparm(ilj)%Rm
          Ema(ires,iat) = ff_ljparm(ilj)%Em
        enddo

      elseif ((attypea(ires,iat) == 'N' ) .or. (attypea(ires,iat) == 'NA') .or.&
              (attypea(ires,iat) == 'NB') .or. (attypea(ires,iat) == 'NC') .or.&
              (attypea(ires,iat) == 'N*') .or. (attypea(ires,iat) == 'N2') .or.&
              (attypea(ires,iat) == 'NO') .or. (attypea(ires,iat) == 'NP')) then

        do ilj = 1, nlj
          if ( .not. (ff_ljparm(ilj)%type == 'N') ) cycle
          Rma(ires,iat) = ff_ljparm(ilj)%Rm
          Ema(ires,iat) = ff_ljparm(ilj)%Em
        enddo

      else
        do ilj = 1, nlj
        if ( attypea(ires,iat) == ff_ljparm(ilj)%type ) then
            Rma(ires,iat) = ff_ljparm(ilj)%Rm
            Ema(ires,iat) = ff_ljparm(ilj)%Em
        endif
        enddo

      endif
    enddo
    enddo

    ! Collapses MM parameters to each MM atom.
    ilj = 1
    do ires = 1, nroaa
    do iat  = 1, atxres(ires)
      mm_atoms(ilj)%lj_Em = Ema(ires,iat)
      mm_atoms(ilj)%lj_Rm = Rma(ires,iat)
      ilj = ilj +1
    enddo
    enddo
    call de_alloc( Ema, 'Ema', 'read_lj'  )
    call de_alloc( Rma, 'Rma', 'read_lj'  )

    ! Assigns the LJ parameters corresponding to QM atom types.
    do iat  = 1, na_qm
      if ( (qm_atoms(iat)%attype == 'C' ) .or.  &
           (qm_atoms(iat)%attype == 'CA') .or. (qm_atoms(iat)%attype == 'CM') .or. &
           (qm_atoms(iat)%attype == 'CC') .or. (qm_atoms(iat)%attype == 'CV') .or. &
           (qm_atoms(iat)%attype == 'CW') .or. (qm_atoms(iat)%attype == 'CR') .or. &
           (qm_atoms(iat)%attype == 'CB') .or. (qm_atoms(iat)%attype == 'C*') .or. &
           (qm_atoms(iat)%attype == 'CN') .or. (qm_atoms(iat)%attype == 'CK') .or. &
           (qm_atoms(iat)%attype == 'CQ') .or. (qm_atoms(iat)%attype == 'CX') .or. &
           (qm_atoms(iat)%attype == 'CY') .or. (qm_atoms(iat)%attype == 'CD') ) then

        do ilj = 1, nlj
          if ( .not. (ff_ljparm(ilj)%type == 'C') ) cycle
          qm_atoms(iat)%lj_Rm = ff_ljparm(ilj)%Rm
          qm_atoms(iat)%lj_Em = ff_ljparm(ilj)%Em
        enddo

      elseif ( (qm_atoms(iat)%attype == 'N' ) .or. (qm_atoms(iat)%attype == 'NA') &
          .or. (qm_atoms(iat)%attype == 'NB') .or. (qm_atoms(iat)%attype == 'NC') &
          .or. (qm_atoms(iat)%attype == 'N*') .or. (qm_atoms(iat)%attype == 'N2') &
          .or. (qm_atoms(iat)%attype == 'NO') .or. (qm_atoms(iat)%attype == 'NP') &
          ) then

        do ilj = 1, nlj
          if ( .not. (ff_ljparm(ilj)%type == 'N') ) cycle
          qm_atoms(iat)%lj_Rm = ff_ljparm(ilj)%Rm
          qm_atoms(iat)%lj_Em = ff_ljparm(ilj)%Em
        enddo

      else
        do ilj = 1, nlj
        if ( qm_atoms(iat)%attype == ff_ljparm(ilj)%type ) then
          qm_atoms(iat)%lj_Rm = ff_ljparm(ilj)%Rm
          qm_atoms(iat)%lj_Em = ff_ljparm(ilj)%Em
        endif
        enddo
      endif
    enddo
    deallocate( ff_ljparm )
  end subroutine read_lj

  subroutine amber_first_ngb( na_mm, nataa, nroaa, atxres, atnu, resname, &
                              mm_connectivity, atnamea, ncon, con, mm_atoms,&
                              mmcell, lattice_type, ui )
    !! Assigns first neighbours accodring to connectivity.
    use functions, only : norm_v2
    use qmmm_pbc , only : pbc_displ_vector, reccel
    use units    , only : Ang

    implicit none
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)  :: nroaa
      !! Number of classical residues.
    integer         , intent(in)  :: nataa(nroaa,100)
      !! Atom index within each residue.
    integer         , intent(in)  :: atxres(nroaa)
      !! Number of atoms in each residue.
    integer         , intent(in)  :: atnu(nroaa,100)
      !! Global index of an atom within a given residue.
    character(len=4), intent(in)  :: resname(nroaa)
      !! Atom types within each residue.
    character(len=4), intent(in)  :: atnamea(nroaa,100)
      !! Atom types within each residue.
    integer         , intent(in)  :: ncon
      !! Total number of connectivity entries.
    integer         , intent(in)  :: con(2,ncon)
      !! Indices for each pair of atoms in a given connection.
    type(mm_atom_t) , intent(in)  :: mm_atoms(na_mm)
      !! Cartesian atomic positions.
    real(dp)        , intent(in)  :: mmcell(3,3)
      !! Unit cell vectors.
    character(len=1), intent(in)  :: lattice_type
      !! Type of lattice used.
    type(atom_connect_t), intent(inout) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.
    integer         , intent(in)  :: ui
      !! File unit for amber parm file.

    character(len=12) :: option
    character(len=1)  :: c1
    character(len=2)  :: c2
    character(len=4)  :: c4
    integer           :: ires, jres, iat, jat, ios, ineig, icon, nresid
    logical           :: search
    real(dp)          :: kcell(3,3), dx, dr(3), rij

    type temp_bond
      integer :: atom1 = 0
      integer :: atom2 = 0
    end type temp_bond

    type temp_residue
      integer          :: nbonds = 0
      character(len=4) :: name   = ""
      type(temp_bond), allocatable :: bonds(:)
    endtype temp_residue
    type(temp_residue), allocatable :: ff_residues(:)

    rewind( ui )

    ios    = 0
    search = .true.
    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 3 )

      if ( option == 'connectivity' ) then
        read( ui, *, iostat = ios ) nresid
        call check_ios_mm( ios, 1 )

        allocate( ff_residues(nresid) )
        do ires = 1, nresid
          read( ui, *, iostat = ios ) ff_residues(ires)%name, &
                                      ff_residues(ires)%nbonds
          call check_ios_mm( ios, 2, ires )

          allocate( ff_residues(ires)%bonds(ff_residues(ires)%nbonds) )
          do iat = 1, ff_residues(ires)%nbonds

            read( ui, *, iostat = ios ) ff_residues(ires)%bonds(iat)%atom1,&
                                        ff_residues(ires)%bonds(iat)%atom2
            call check_ios_mm( ios, 3, ires )
          enddo
        enddo

        search = .false.
      endif
    enddo

    ! Assigns first neighbours for each atom.
    do jres = 1, nresid
    do ires = 1, nroaa
      if ( .not. (resname(ires) == ff_residues(jres)%name) ) cycle

      do iat = 1, atxres(ires)
        ineig = 1

        do icon = 1, ff_residues(jres)%nbonds
          if ( nataa(ires,iat) == ff_residues(jres)%bonds(icon)%atom1 ) then
            do jat = 1, atxres(ires)
              if ( .not. (nataa(ires,jat) == &
                   ff_residues(jres)%bonds(icon)%atom2 ) ) cycle
              mm_connectivity(atnu(ires,iat))%bond_at(ineig) = atnu(ires,jat)
              ineig = ineig +1
            enddo

          elseif ( nataa(ires,iat) == ff_residues(jres)%bonds(icon)%atom2 ) then
            do jat = 1, atxres(ires)
              if ( .not. (nataa(ires,jat) == &
                   ff_residues(jres)%bonds(icon)%atom1 ) ) cycle
              mm_connectivity(atnu(ires,iat))%bond_at(ineig) = atnu(ires,jat)
              ineig = ineig +1
            enddo
          endif
        enddo ! Connections for atom iat.
      enddo   ! Atom iat in a residue.
    enddo
      deallocate( ff_residues(jres)%bonds )
    enddo
    deallocate( ff_residues )

    ! Neighbour assignments for graphene slabs.
    call reccel( 3, mmcell, kcell, 0 )

    do ires = 1, nroaa
      if ( .not. ((resname(ires) == 'GRAP') .or. (resname(ires) == 'GRAH')) )&
        cycle
      ineig = 1
      do jres = 1, nroaa
        if ( ires == jres ) cycle
        if ( .not. ((resname(jres) == 'GRAP') .or. (resname(jres) == 'GRAH')) )&
          cycle

        dr(:) = mm_atoms(atnu(ires,1))%r(:) - mm_atoms(atnu(jres,1))%r(:)
        call pbc_displ_vector( lattice_type, mmcell, kcell, dr )
        rij = norm_v2( dr )

        ! dx is now a dummy array.
        if ( resname(ires) == resname(jres) ) then
          dx = 1.7_dp * Ang
        else
          dx = 1.4_dp * Ang
        endif

        if ( rij < dx ) then
          mm_connectivity(atnu(ires,1))%bond_at(ineig) = atnu(jres,1)
          ineig = ineig +1
        endif
      enddo
    enddo

    ! Finds neighbours of 2 contiguous residues.
    do ires = 2,nroaa
      c4 = resname(ires)
      c1 = c4(1:1)

      if ( c1 == 'N' ) then ! Might be a terminal N.
        if ( .not. ((c4 == 'NME') .or. (c4 == 'NALB')) ) cycle
      endif

      do iat = 1, atxres(ires)
        if ( atnamea(ires,iat) /= 'N' ) cycle

        do jat = 1, atxres(ires-1)
          if ( atnamea(ires-1,jat) /= 'C') cycle
          mm_connectivity(atnu(ires  ,iat))%bond_at(3) = atnu(ires-1,jat)
          mm_connectivity(atnu(ires-1,jat))%bond_at(3) = atnu(ires,iat)
        enddo
      enddo
    enddo

    ! Find neighbours between two consecutive nucleotides.
    do ires = 2, nroaa
      c4 = resname(ires)
      c1 = c4(3:3)
      if ( c1 == '5' ) cycle

      do jres = 1, atxres(ires)
        if ( atnamea(ires,jres) /= 'P' ) cycle

        do iat = 1, atxres(ires-1)
          c4 = atnamea(ires-1,iat)
          c2 = c4(1:2)
          if ( c2 /= 'O3' ) cycle

          mm_connectivity(atnu(ires  ,jres))%bond_at(4) = atnu(ires-1,iat)
          mm_connectivity(atnu(ires-1,iat ))%bond_at(2) = atnu(ires,jres)
        enddo
      enddo
    enddo

    ! Assigns any additional connections forced in the input.
    do icon = 1, ncon
    do iat  = 1, 6
      if ( mm_connectivity(con(1,icon))%bond_at(iat) /= 0 ) cycle
      search = .false.

      mm_connectivity(con(1,icon))%bond_at(iat) = con(2,icon)

      do ineig = 1, 6

        ! Here we use search as a dummy variable to control flow.
        if ( mm_connectivity(con(2,icon))%bond_at(ineig) /= 0 ) cycle
        mm_connectivity(con(2,icon))%bond_at(ineig) = con(1,icon)
        search = .true.
        exit
      enddo
      if ( search ) exit ! Go to the next icon.
    enddo
    enddo
  end subroutine amber_first_ngb

  subroutine bon_ang_dih_imp( na_mm, mm_connectivity, atnamea, nroaa, atxres, &
                              atnu, resname, ui )
    ! Assigns all connectivity magnitudes: bonds, angles, dihedrals and
    ! improper torsions. Also, reads impropers data.
    use alloc    , only : re_alloc, de_alloc
    use sys      , only : die

    implicit none
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)  :: nroaa
      !! Number of classical residues.
    integer         , intent(in)  :: atxres(nroaa)
      !! Number of atoms in each residue.
    integer         , intent(in)  :: atnu(nroaa,100)
      !! Global index of an atom within a given residue.
    character(len=4), intent(in)  :: resname(nroaa)
      !! Atom types within each residue.
    character(len=4), intent(in)  :: atnamea(nroaa,100)
      !! Atom types within each residue.
    type(atom_connect_t), intent(inout) :: mm_connectivity(na_mm)
    !! Data structure containing the connectivity data for each MM atom.
    integer         , intent(in)  :: ui
      !! File unit for amber parm file.

    integer           :: iat, jat, kat, lat, mat, ires, jres, iimp, jimp, &
                         ios, nresid, ntot, impmax
    logical           :: search
    character(len=10) :: option
    integer, pointer  :: impnum(:,:)

    type temp_residue
      integer          :: nimpr = 0
      character(len=4) :: name  = ""
      character(len=4), allocatable :: impr(:,:)
    endtype temp_residue
    type(temp_residue), allocatable :: ff_residues(:)

    ! Reads and calculates improper torsions.
    rewind( ui )

    ios    = 0
    search = .true.
    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 4 )

      if ( option == 'impropers' ) then
        read( ui, *, iostat = ios ) nresid
        call check_ios_mm( ios, 4, 1 )

        allocate( ff_residues(nresid) )
        do ires = 1, nresid
          read( ui, *, iostat = ios ) ff_residues(ires)%name, &
                                      ff_residues(ires)%nimpr
          call check_ios_mm( ios, 4, 2, ires )

          allocate( ff_residues(ires)%impr(ff_residues(ires)%nimpr,4) )
          do iat = 1, ff_residues(ires)%nimpr
            read( ui, *, iostat = ios ) ff_residues(ires)%impr(iat,1), &
                                        ff_residues(ires)%impr(iat,1), &
                                        ff_residues(ires)%impr(iat,3), &
                                        ff_residues(ires)%impr(iat,4)
            call check_ios_mm( ios, 4, 3, ires )
          enddo
        enddo

        search = .false.
      endif
    enddo

    ! Assigns impropers according to atoms, starting from the
    ! second residue.
    impmax = nroaa * 25

    nullify( impnum )
    call re_alloc( impnum, 1, impmax, 1, 4, 'impnum', 'bon_ang_dih_imp' )
    impnum = 0
    ntot   = 1

    do ires = 1, nresid
    do jres = 1, nroaa
      if ( ff_residues(ires)%name /= resname(jres) ) cycle
      do iimp = 1, ff_residues(ires)%nimpr
        do iat = 1, 4

        if ( (ff_residues(ires)%impr(iimp,iat) == '+M') .and. (jres /= nroaa) ) then
          do jat = 1, atxres(jres+1)
            if ( atnamea(jres+1,jat) /= 'N') cycle
            impnum(ntot,iat) = atnu(jres+1,jat)
          enddo

        elseif ( (ff_residues(ires)%impr(iimp,iat) == '-M') .and. (jres /= 1) ) then
          do jat = 1,atxres(jres-1)
            if ( atnamea(jres-1,jat) /= 'C') cycle
            impnum(ntot,iat) = atnu(jres-1,jat)
          enddo

        else
          do jat = 1, atxres(jres)
            if ( atnamea(jres,jat) /= ff_residues(ires)%impr(iimp,iat) ) cycle
            impnum(ntot,iat) = atnu(jres,jat)
          enddo
        endif
        enddo

        ntot = ntot +1
        if ( ntot >= impmax ) &
          call die( 'mm_assign: increase the size of the improper matrix.' )
      enddo
    enddo
      deallocate( ff_residues(ires)%impr )
    enddo

    do iimp = 1, ntot
    do iat  = 1, 4
      if (impnum(iimp,iat) /= 0) cycle
      impnum(iimp,:) = 0
    enddo
    enddo

    do iat = 1, na_mm
      jimp = 0
      do iimp = 1, ntot
        if ( (impnum(iimp,1) == iat) .or. (impnum(iimp,2) == iat) .or. &
             (impnum(iimp,3) == iat) .or. (impnum(iimp,4) == iat) ) then
          jimp = jimp +1
        endif
      enddo
      mm_connectivity(iat)%nimp = jimp
      nullify( mm_connectivity(iat)%imp_at )
      call re_alloc( mm_connectivity(iat)%imp_at, 1, &
                     mm_connectivity(iat)%nimp, 1, 4, &
                     'imp_at' , 'siesta_qmmm' )
    enddo

    do iat = 1, na_mm
      jimp = 0
      do iimp = 1, ntot
        if ( (impnum(iimp,1) == iat) .or. (impnum(iimp,2) == iat) .or. &
             (impnum(iimp,3) == iat) .or. (impnum(iimp,4) == iat) ) then
          jimp = jimp +1
          mm_connectivity(iat)%imp_at(jimp,:) = impnum(iimp,:)
        endif
      enddo
    enddo
    call de_alloc( impnum, 'impnum', 'bon_ang_dih_imp' )
    deallocate( ff_residues )
    ! Finished with impropers.

    ! Now we use connectivity data for bonds, angles and dihedrals.
    ! Here we assign bonds.
    do iat = 1, na_mm
      mm_connectivity(iat)%nbonds = 0
      do jat = 1, 6
        if ( mm_connectivity(iat)%bond_at(jat) /= 0 ) &
          mm_connectivity(iat)%nbonds = mm_connectivity(iat)%nbonds +1
      enddo
    enddo

    ! Searches for angles with atom iat in the extremes.
    do iat = 1, na_mm
      ntot = 1
      do jat = 1, mm_connectivity(iat)%nbonds
        kat = mm_connectivity(iat)%bond_at(jat)

        do lat = 1, mm_connectivity(kat)%nbonds
          if ( mm_connectivity(kat)%bond_at(lat) == iat ) cycle
          ntot = ntot +1
        enddo
      enddo
      mm_connectivity(iat)%nangl_e = ntot -1

      nullify( mm_connectivity(iat)%ange_at )
      call re_alloc( mm_connectivity(iat)%ange_at, 1, &
                     mm_connectivity(iat)%nangl_e, 1, 2, &
                     'ange_at' , 'siesta_qmmm' )
    enddo

    do iat = 1, na_mm
      ntot = 1
      do jat = 1, mm_connectivity(iat)%nbonds
        kat = mm_connectivity(iat)%bond_at(jat)

        do lat = 1, mm_connectivity(kat)%nbonds
          if ( mm_connectivity(kat)%bond_at(lat) == iat ) cycle

          mm_connectivity(iat)%ange_at(ntot,1) = mm_connectivity(iat)%bond_at(jat)
          mm_connectivity(iat)%ange_at(ntot,2) = mm_connectivity(kat)%bond_at(lat)
          ntot = ntot +1
        enddo
      enddo
    enddo

    ! Searches for angles with atom iat in the middle.
    do iat = 1, na_mm
      ntot = 1

      do jat = 1, mm_connectivity(iat)%nbonds
      do kat = 1, mm_connectivity(iat)%nbonds
        if ( mm_connectivity(iat)%bond_at(kat) <= &
             mm_connectivity(iat)%bond_at(jat) ) cycle
        ntot = ntot +1
      enddo
      enddo

      mm_connectivity(iat)%nangl_m = ntot -1
      nullify( mm_connectivity(iat)%angm_at )
      call re_alloc( mm_connectivity(iat)%angm_at, 1, &
                     mm_connectivity(iat)%nangl_m, 1, 2, &
                     'angm_at' , 'siesta_qmmm' )
    enddo

    do iat = 1, na_mm
      ntot = 1

      do jat = 1, mm_connectivity(iat)%nbonds
      do kat = 1, mm_connectivity(iat)%nbonds
        if ( mm_connectivity(iat)%bond_at(kat) <= &
             mm_connectivity(iat)%bond_at(jat) ) cycle

        mm_connectivity(iat)%angm_at(ntot,1) = mm_connectivity(iat)%bond_at(jat)
        mm_connectivity(iat)%angm_at(ntot,2) = mm_connectivity(iat)%bond_at(kat)
        ntot = ntot +1
      enddo
      enddo
    enddo

    ! Searches for dihedrals with atom iat in the extremes.
    do iat = 1, na_mm
      ntot = 1

      do jat = 1, mm_connectivity(iat)%nangl_e
        kat = mm_connectivity(iat)%ange_at(jat,2)

        do lat = 1, mm_connectivity(kat)%nbonds
          mat = mm_connectivity(kat)%bond_at(lat)

          if ( mat /= mm_connectivity(iat)%ange_at(jat,1) ) then
            ntot = ntot +1
          endif
        enddo
      enddo
      mm_connectivity(iat)%ndihe_e = ntot -1
      nullify( mm_connectivity(iat)%dihe_at )
      call re_alloc( mm_connectivity(iat)%dihe_at, 1, &
                     mm_connectivity(iat)%ndihe_e, 1, 3, &
                     'dihe_at' , 'siesta_qmmm' )
    enddo

    do iat = 1, na_mm
      ntot = 1

      do jat = 1, mm_connectivity(iat)%nangl_e
        kat = mm_connectivity(iat)%ange_at(jat,2)

        do lat = 1, mm_connectivity(kat)%nbonds
          mat = mm_connectivity(kat)%bond_at(lat)

          if ( mat /= mm_connectivity(iat)%ange_at(jat,1) ) then
            mm_connectivity(iat)%dihe_at(ntot,1) = mm_connectivity(iat)%ange_at(jat,1)
            mm_connectivity(iat)%dihe_at(ntot,2) = kat
            mm_connectivity(iat)%dihe_at(ntot,3) = mat
            ntot = ntot +1
          endif
        enddo
      enddo
    enddo

    ! Searches for dihedrals with atom iat in the middle, starting
    ! from angles that have that same atom in the extreme.
    do iat = 1, na_mm
      ntot = 1

      do jat = 1, mm_connectivity(iat)%nangl_e
        kat = mm_connectivity(iat)%ange_at(jat,1)

        do lat = 1, mm_connectivity(iat)%nbonds
          mat = mm_connectivity(iat)%bond_at(lat)
          if ( mat /= kat ) then
            ntot = ntot +1
          endif
        enddo
      enddo

      mm_connectivity(iat)%ndihe_m = ntot -1
      nullify( mm_connectivity(iat)%dihm_at )
      call re_alloc( mm_connectivity(iat)%dihm_at, 1, &
                     mm_connectivity(iat)%ndihe_m, 1, 3, &
                     'dihm_at' , 'siesta_qmmm' )
    enddo

    do iat = 1, na_mm
      ntot = 1
      do jat = 1, mm_connectivity(iat)%nangl_e
        kat = mm_connectivity(iat)%ange_at(jat,1)

        do lat = 1, mm_connectivity(iat)%nbonds
          mat = mm_connectivity(iat)%bond_at(lat)
          if ( mat /= kat ) then
            mm_connectivity(iat)%dihm_at(ntot,1) = mat
            mm_connectivity(iat)%dihm_at(ntot,2) = kat
            mm_connectivity(iat)%dihm_at(ntot,3) = mm_connectivity(iat)%ange_at(jat,2)
            ntot = ntot +1
          endif
        enddo
      enddo
    enddo
  end subroutine bon_ang_dih_imp

  subroutine atxaa( ntres, aanamea, atomsxaa, ui )
    !! Returns how many atoms are in each residue, and eaech residue's name.
    use sys      , only : die

    implicit none
    integer         , intent(out) :: ntres
      !! Number of residue types.
    character(len=4), intent(out) :: aanamea(MAX_TYPE_RES)
      !! Name of each residue.
    integer         , intent(out) :: atomsxaa(MAX_TYPE_RES)
      !! Number of atoms per residue.
    integer         , intent(in)  :: ui
      !! File unit for amber parm file.

    character(len=4)  :: dummychar
    character(len=10) :: option
    integer           :: ires, iat, ios
    logical           :: search

    rewind( ui )

    ios    = 0
    search = .true.
    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 5 )

      if ( option == 'residues' ) then
        read( ui, *, iostat = ios ) ntres
        call check_ios_mm( ios, 5 , 1 )

        if ( ntres > MAX_TYPE_RES ) &
          call die( 'mm_assign: The number of different residues must'//&
                    ' not exceed MAX_TYPE_RES.' )

        do ires = 1, ntres
          read( ui, *, iostat = ios ) aanamea(ires),atomsxaa(ires)
          call check_ios_mm( ios, 5 , 2, ires )

          if ( atomsxaa(ires) > 100 ) &
            call die( 'mm_assign: Number of atoms in a given residue must'//&
                      ' not exceed 100.' )

          ! Basically skips over the next nAtoms lines.
          do iat = 1, atomsxaa(ires)
            read( ui, *, iostat = ios ) dummychar
            call check_ios_mm( ios, 5 , 3, ires )
          enddo
        enddo

        search = .false.
        exit
      endif
    enddo
  end subroutine atxaa

  subroutine amber_union_parms( mm_top, ui )
    !! Reads the number, constants and other parameters for bonds,
    !! angles, dihedrals and improper torsions.
    use sys      , only : die

    implicit none
    type(fftopology), intent(inout) :: mm_top
      !! Data structure with topology.
    integer         , intent(in)    :: ui
      !! File unit for amber parm file.

    character(len=12) :: option
    integer           :: ios, ipar
    logical           :: search

    rewind( ui )
    ! Bond data.

    ios    = 0
    search = .true.
    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 6 )

      if ( option == 'bonds' ) then
        read( ui, *, iostat = ios ) mm_top%nbonds
        call check_ios_mm( ios, 6, 1 )

        allocate( mm_top%bonds(mm_top%nbonds) )
        do ipar = 1, mm_top%nbonds
          read( ui, fmt = '(A5,2x,F5.1,4x,F6.4)', iostat = ios ) &
            mm_top%bonds(ipar)%type, mm_top%bonds(ipar)%k,       &
            mm_top%bonds(ipar)%r_eq
          call check_ios_mm( ios, 6, 2, ipar )
        enddo

        search = .false.
        exit
      endif
    enddo

    ! Angle data
    rewind ( ui )
    search = .true.
    ios    = 0

    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 7 )

      if ( option == 'angles' ) then
        read( ui, *, iostat = ios ) mm_top%nangles
        call check_ios_mm( ios, 7, 1 )

        allocate( mm_top%angles(mm_top%nangles) )
        do ipar = 1, mm_top%nangles
          read( ui, fmt = '(A8,3x,F5.1,6x,F6.2)', iostat = ios ) &
            mm_top%angles(ipar)%type, mm_top%angles(ipar)%k,     &
            mm_top%angles(ipar)%r_eq

          call check_ios_mm( ios, 7, 2, ipar )
        enddo

        search = .false.
        exit
      endif
    enddo

    ! Dihedrals
    rewind( ui )
    search = .true.
    ios    = 0

    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 8 )

      if ( option == 'dihes' ) then
        read( ui, *, iostat = ios ) mm_top%ndihe
        call check_ios_mm( ios, 8, 1 )

        allocate( mm_top%dihedrals(mm_top%ndihe) )
        do ipar = 1, mm_top%ndihe
          read( ui, fmt='(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)', iostat = ios )&
            mm_top%dihedrals(ipar)%type, mm_top%dihedrals(ipar)%multi, &
            mm_top%dihedrals(ipar)%k   , mm_top%dihedrals(ipar)%r_eq , &
            mm_top%dihedrals(ipar)%per

          call check_ios_mm( ios, 8, 2, ipar )
        enddo
        search = .false.
        exit
      endif
    enddo

    ! Improper torsions.
    rewind( ui )
    search = .true.
    ios    = 0

    do while ( search )
      read ( ui, *, iostat = ios ) option
      call check_ios_mm( ios, 9 )

      if ( option == 'imps' ) then
        read( ui, *, iostat = ios ) mm_top%nimp
        call check_ios_mm( ios, 9, 1 )

        allocate( mm_top%impropers(mm_top%nimp) )
        do ipar = 1, mm_top%nimp
          read( ui, fmt='(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)', iostat = ios )&
            mm_top%impropers(ipar)%type, mm_top%impropers(ipar)%multi, &
            mm_top%impropers(ipar)%k   , mm_top%impropers(ipar)%r_eq , &
            mm_top%impropers(ipar)%per
          call check_ios_mm( ios, 9, 1, ipar )

          ! Must be always positive.
          mm_top%impropers(ipar)%per = abs( mm_top%impropers(ipar)%per )
        enddo

        search = .false.
        exit
      endif
    enddo
  end subroutine amber_union_parms


  subroutine mm_atom_to_species( na_qm, na_mm, na_link, nspec, mm_atoms, isa )
    !! Stores the different MM atom types into the ISA array; useful for the
    !! XV file printing. It also outputs the different MM atom types found.

    implicit none
    integer        , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer        , intent(in)    :: na_mm
      !! Number of MM atoms.
    integer        , intent(in)    :: na_link
      !! Number of Link atoms.
    integer        , intent(in)    :: nspec
      !! Number of QM atom species.
    type(mm_atom_t), intent(in)    :: mm_atoms(na_mm)
      !! Data structure with the MM atoms information.
    integer        , intent(inout) :: isa(na_qm+na_mm+na_link)
      !! Array containing the species index.

    integer :: ia, itype, ntype
    logical :: found_type
    character(len=4), allocatable :: at_type(:)

    if ( na_mm < 1 ) return
    allocate( at_type(na_mm) )
    at_type = 'XXXX'

    ntype = 1
    do ia = 1, na_mm

      found_type = .false.
      do itype = 1, ntype
        if ( trim(mm_atoms(ia)%attype) == trim(at_type(itype)) ) then
          found_type = .true.
          exit
        endif
      enddo
      if ( found_type ) cycle

      ! If we found a new type we haven't previously stored.
      at_type(ntype) = mm_atoms(ia)%attype
      ntype = ntype +1
    enddo
    ntype = ntype -1

    ! Now we transpose this data into the ISA array.
    do ia = 1, na_mm
      do itype = 1, ntype
        if ( trim(mm_atoms(ia)%attype) == trim(at_type(itype)) ) exit
      enddo
      isa(na_qm+ia) = nspec + itype
    enddo

    ! Finally, we print the info to output.
    write(6,'(a)') " "
    write(6,'(a41,I3)') "Number of different MM atom types found: ", ntype
    do itype = 1, ntype
      write(6,'(a5,i3,a11,i3,a4,a4)') &
        "Type ", itype, " (species: ", nspec+itype, ") - ", &
        at_type(itype)
    enddo
    write(6,'(a)') " "

    deallocate( at_type )

    ! In the case of link atoms, we just assign the last possible type.
    if ( na_link < 1 ) return
    do ia = 1, na_link
      isa(na_qm+na_mm+ia) = nspec + ntype +1
    enddo

  end subroutine mm_atom_to_species

end module mm_assign_m