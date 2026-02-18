module linkatoms

  use precision, only: dp

  implicit none

  character(len=5) :: h_link_scaling_type = 'splam'
  logical :: got_numlink = .false.
  character(len=2) :: default_link_symbol = 'H?'
  logical :: linkatom = .false.

  integer :: default_link_type = 0
  integer :: numlink = 0
  integer :: NSpecialLinkAtoms = 0
  real(dp) :: Elink

  type, public :: resonace_link_t
    real(dp) :: r(3)    = 0.0_dp
    real(dp) :: fa(3)   = 0.0_dp
    real(dp) :: bond_r  = 0.0_dp
    real(dp) :: bond_k  = 0.0_dp
    real(dp) :: angle_k(3)   = 0.0_dp
    integer  :: params(25,4) = 0
  end type resonace_link_t

  type, public :: link_atom_t
    integer          :: num_qm   = 0
    integer          :: num_mm   = 0
    integer          :: qm(4)    = 0
    integer          :: mm(4)    = 0
    integer          :: mm2(4,3) = 0
    real(dp)         :: pc(4)    = 0.0_dp
    real(dp)         :: Em(4)    = 0.0_dp
    character(len=4) :: qm_type(4) = ''

    type(resonace_link_t), allocatable :: reson(:)
  end type link_atom_t

  type(link_atom_t), allocatable :: link_atoms(:)

  ! C. F. Sanz-Navarro (Nov. 2009)
  ! this module deals with the case of having resonance structures in the QM
  ! part, which cannot be accurately described by 1 type of link atoms and
  ! needs the linear combination of different ways (single and double bonds
  ! with the link atoms) to be able to model the system

  integer :: num_resonances = 1
  type resonance_type
     real(dp) :: weight
     real(dp) :: sp_weight
       !! Weight for "special" link atoms.
     character(len=4), pointer :: atype(:)
     integer, pointer :: isa(:)
     integer, pointer :: iza(:)
     character(len=13) :: path
  end type resonance_type
  type(resonance_type), allocatable :: resonance(:)

  type, public :: sp_resonace_link_t
    integer          :: isa
    character(len=4) :: attype
  end type sp_resonace_link_t

  type, public :: sp_link_atom_t
    integer :: qm
    integer :: mm
    type(sp_resonace_link_t), allocatable :: reson(:)
  end type sp_link_atom_t
  type(sp_link_atom_t), allocatable :: sp_link_atoms(:)

contains
  subroutine allocate_link_arrays()
    use alloc, only: re_alloc
    implicit none
    integer :: ia

    allocate( link_atoms(numlink) )
    if ( numlink > 0 ) then
      do ia = 1, numlink
        allocate( link_atoms(ia)%reson(num_resonances) )
      enddo
    endif
  end subroutine allocate_link_arrays

  subroutine deallocate_link_arrays()
    use alloc, only: de_alloc
    implicit none
    integer :: ia

    if ( NSpecialLinkAtoms > 0 ) then
      do ia = 1, numlink
        deallocate( link_atoms(ia)%reson )
      enddo
    endif
    if ( allocated(link_atoms) ) deallocate( link_atoms )

    if ( NSpecialLinkAtoms > 0 ) then
      do ia = 1, numlink
        deallocate( sp_link_atoms(ia)%reson )
      enddo
    endif
    if ( allocated(sp_link_atoms) ) deallocate( sp_link_atoms )
  end subroutine deallocate_link_arrays

  subroutine read_special_link_atoms( bfdf )
    use fdf       , only : fdf_get, leqi
    use fdf       , only : block_fdf, parsed_line
    use fdf       , only : fdf_block, fdf_bline, fdf_bclose, &
                           fdf_bintegers, fdf_bvalues, fdf_bnames
    use sys       , only : die

    implicit none
    type(block_fdf), intent(in) :: bfdf

    integer  :: i, j
    real(dp) :: total_weight

    character(len=80) :: acf, acf_default
    integer           :: iscale
    type(parsed_line), pointer :: pline

    ! Format of atomic coordinates
    acf_default = 'Ang'
    acf = fdf_get('AtomicCoordinatesFormat',acf_default)

    if (leqi(acf,'NotScaledCartesianBohr') .or. leqi(acf,'Bohr') ) then
      iscale = 0
      write(6,'(a,a)') 'read: Atomic-coordinates input format  = ',&
                       '    Cartesian coordinates (in Bohr)'
    elseif (leqi(acf,'NotScaledCartesianAng') .or. leqi(acf,'Ang') ) then
      iscale = 1
      write(6,'(a,a)') 'read: Atomic-coordinates input format  = ',&
                       '    Cartesian coordinates (in Ang)'
    else
      write(6,"(/,'read: ',73(1h*))")
      write(6,"('read:                  INPUT ERROR')")
      write(6,'(a)') 'read: You must use one of the following'//&
                     'coordinate options:'
      write(6,'(a)') 'read:     - NotScaledCartesianBohr (or Bohr)'
      write(6,'(a)') 'read:     - NotScaledCartesianAng (or Ang) '
      write(6,"('read: ',73(1h*))")
      call die( 'Error reading link atoms - Wrong AtomicCoordinatesFormat.' )
    endif


    total_weight = 0.0_dp
    if ( fdf_bline( bfdf, pline ) ) then
      do i = 1, num_resonances
        resonance(i)%sp_weight = fdf_bvalues(pline,i)
        total_weight = total_weight + resonance(i)%sp_weight
      enddo
    else
      call die( 'special_linkatoms: No resonance weights found.' )
    endif

    if ( abs(1.0_dp-total_weight) > 0.0001_dp ) &
      write(*,*) 'WARNING: The sum of resonance weights is not zero.'//&
                 'Total sum = ', total_weight

    do j = 1, NSpecialLinkAtoms
      if ( fdf_bline( bfdf, pline ) ) then
        do i = 1, num_resonances
          sp_link_atoms(j)%reson(i)%isa    = fdf_bintegers(pline, i)
          sp_link_atoms(j)%reson(i)%attype = trim( fdf_bnames(pline, i) )
        enddo

        sp_link_atoms(j)%qm = fdf_bintegers(pline, num_resonances+1)
        sp_link_atoms(j)%mm = fdf_bintegers(pline, num_resonances+2)
      else
        call die( 'special_linkatoms: Missing special link atom information.' )
      endif
    enddo
  end subroutine read_special_link_atoms

  subroutine get_numlink( na_qm, na_mm, mm_atoms, qm_atoms, &
                          mm_top, deb, cell, lattice_type )
    use alloc      , only : de_alloc, re_alloc
    use functions  , only : norm_v2
    use fdf        , only : block_fdf, parsed_line
    use fdf        , only : fdf_block, fdf_bline, fdf_bclose, fdf_bnintegers,&
                            fdf_bintegers
    use mm_topology, only : fftopology, mm_atom_t, qm_atom_t
    use qmmm_pbc   , only : pbc_displ_vector, reccel
    use sys        , only : die
    use units      , only : Ang

    implicit none
    integer         , intent(in) :: na_mm
    integer         , intent(in) :: na_qm
    type(mm_atom_t) , intent(in) :: mm_atoms(na_mm)
    type(qm_atom_t) , intent(in) :: qm_atoms(na_qm)
    type(fftopology), intent(in) :: mm_top
    logical         , intent(in) :: deb
    real(dp)        , intent(in) :: cell(3,3)
    character       , intent(in) :: lattice_type

    character         :: ch1
    character(len=4)  :: ch4, chi, chj
    character(len=5)  :: tybond
    integer           :: i, j, l
    logical           :: parameter_set, outofrange
    real(dp)          :: r1, dr(3), kcell(3,3), qmmm_range
    integer , pointer :: qm2link(:), min_at(:)
    real(dp), pointer :: rmin(:)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    got_numlink = .true.

    if ( na_qm == 0 ) return
    qmmm_range = 0.0_dp
    call reccel( 3, cell, kcell, 0 )

    nullify( qm2link, min_at, rmin )
    call re_alloc( qm2link, 1, na_qm, 'qm2link', 'get_numlink' )
    call re_alloc( min_at , 1, na_qm, 'min_at' , 'get_numlink' )
    call re_alloc( rmin   , 1, na_qm, 'rmin'   , 'get_numlink' )
    qm2link = 0
    min_at  = 0
    rmin    = 0.0_dp

    if ( fdf_block( 'SpecialLinkAtoms', bfdf ) ) then
      call message( 'WARNING', 'The resonance-type of scheme is DEPRECATED. ' )
      write(6,*) 'WARNING: The resonance-type of scheme is DEPRECATED. '//&
                'Unless you know what you are doing, try a different way to'//&
                ' partition the QM and MM regions that does not involve'//&
                ' breaking resonant bonds/structures.'

      if (.not. fdf_bline(bfdf,pline)) &
        call die('link1: ERROR in SpecialLinkAtoms block')

      if ( fdf_bnintegers(pline) < 2 ) &
          call die( "link1: problem reading number of resonances or" //&
                    " number of atoms in SpecialLinkAtoms." )

      NSpecialLinkAtoms = fdf_bintegers(pline, 1)
      num_resonances    = fdf_bintegers(pline, 2)

      if ( NSpecialLinkAtoms > 0 ) then
        if ( allocated( sp_link_atoms ) ) then
          do i = 1, NSpecialLinkAtoms
            if ( allocated(sp_link_atoms(i)%reson) ) &
              deallocate( sp_link_atoms(i)%reson )
          enddo
          deallocate( sp_link_atoms )
        endif

        allocate( sp_link_atoms(NSpecialLinkAtoms) )
        do i = 1, NSpecialLinkAtoms
          allocate( sp_link_atoms(i)%reson(num_resonances) )
        enddo

        call read_special_link_atoms( bfdf )
        call fdf_bclose( bfdf )
        do i = 1, NSpecialLinkAtoms
          qm2link(sp_link_atoms(i)%qm) = 1
        enddo
      endif
    endif

    numlink = NSpecialLinkAtoms

    ! Count total number of link atoms
    do i = 1, na_qm
      if ( qm2link(i) == 1 ) cycle
      dr = qm_atoms(i)%r(:) - mm_atoms(1)%r(:)

      call pbc_displ_vector( lattice_type, cell, kcell, dr )
      rmin(i)   = norm_v2( dr )
      min_at(i) = 1

      do j = 2, na_mm
        dr = qm_atoms(i)%r(:) - mm_atoms(j)%r(:)
        call pbc_displ_vector( lattice_type, cell, kcell, dr )
        r1 = norm_v2( dr )
        if ( r1 <= rmin(i) ) then
           rmin(i)   = r1
           min_at(i) = j
        endif
      enddo

      ! Only consider pairs of atoms close enough
      if ( rmin(i) > (2.8_dp * Ang) ) cycle

      chi = qm_atoms(i)%attype
      if ( chi(1:1) == 'H' ) cycle

      j   = min_at(i)
      chj = mm_atoms(j)%attype
      if ( chj(1:1) == 'H' ) cycle

      parameter_set = .false.
      qmmm_range    = 0.0_dp
      do l = 1, mm_top%nbonds
        tybond = mm_top%bonds(l)%type

        if ( (chi(1:2) == tybond(1:2)) .and. (chj(1:2) == tybond(4:5)) ) then
          qmmm_range    = 1.2_dp * Ang * mm_top%bonds(l)%r_eq
          parameter_set = .true.
          exit
        elseif ((chi(1:2) == tybond(4:5)) .and. (chj(1:2) == tybond(1:2))) then
          qmmm_range    = 1.2_dp * Ang * mm_top%bonds(l)%r_eq
          parameter_set = .true.
          exit
        endif
      enddo

      if ( .not. parameter_set ) then
        outofrange = &
          check_out_of_range( i, j, qm_atoms(i)%attype, mm_atoms(j)%attype, rmin(i)/Ang)
        if ( outofrange ) cycle
      endif

      ch4 = qm_atoms(i)%attype
      ch1 = ch4(1:1)

      if ( rmin(i) < qmmm_range ) then
        if ( deb ) &
          write(6,*) "get_numlink: Checking link atoms for QM atom no. ", i

        if ( ch1 == 'O' ) then
          if ( deb ) write(6,*) i, min_at(i), rmin(i)
          numlink = numlink +1

        elseif ( (ch1 == 'C') .or. (ch1 == 'N') ) then
          if ( deb ) write(6,*) i, min_at(i), rmin(i)
          numlink = numlink +1

        endif
      endif
    enddo

    call de_alloc( qm2link, 'qm2link', 'get_numlink' )
    call de_alloc( min_at , 'min_at' , 'get_numlink' )
    call de_alloc( rmin   , 'rmin'   , 'get_numlink' )
  end subroutine get_numlink

  ! Reads Link Atoms parameters
  subroutine link1( mm_connectivity, na_qm, na_mm, mm_atoms, qm_atoms, mm_top, &
                    nspec, deb, cell, lattice_type )
    use alloc       , only : de_alloc, re_alloc
    use fdf         , only : fdf_get
    use fdf         , only : block_fdf, parsed_line
    use fdf         , only : fdf_block, fdf_bline, fdf_bclose, fdf_bnintegers,&
                             fdf_bintegers, fdf_bnames, fdf_bnnames
    use functions   , only : norm_v2
    use mm_topology , only : fftopology, atom_connect_t, mm_atom_t, qm_atom_t
    use qmmm_pbc    , only : pbc_displ_vector, reccel
    use sys         , only : die
    use units       , only : Ang

    type(fftopology)    , intent(in) :: mm_top
    type(atom_connect_t), intent(in) :: mm_connectivity(na_mm)
      !! Data structure containing the connectivity data for each MM atom.
    integer         , intent(in) :: na_mm
    integer         , intent(in) :: na_qm
    integer         , intent(in) :: nspec
    type(mm_atom_t) , intent(in) :: mm_atoms(na_mm)
    type(qm_atom_t) , intent(in) :: qm_atoms(na_qm)
    logical         , intent(in) :: deb
    real(dp)        , intent(in) :: cell(3,3)
    character       , intent(in) :: lattice_type

    character        :: ch1
    character(len=2) :: atype, atsym(nspec), ch2
    character(len=4) :: chi, chj, ch4
    character(len=5) :: tybond
    integer          :: i, j, k, l, m, atnum(nspec), &
                        default_link_isa, default_link_iza
    logical          :: parameter_set
    real(dp)         :: r1, kcell(3,3), dr(3), qmmm_range

    integer , pointer :: qm2link(:), min_at(:)
    real(dp), pointer :: rmin(:)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    if ( na_qm   == 0 ) return
    if ( numlink == 0 ) return

    qmmm_range = 0.0_dp
    call reccel( 3, cell, kcell, 0 )

    nullify( qm2link, min_at, rmin )
    call re_alloc( qm2link, 1, na_qm, 'qm2link', 'link1' )
    call re_alloc( min_at , 1, na_qm, 'min_at' , 'link1' )
    call re_alloc( rmin   , 1, na_qm, 'rmin'   , 'link1' )
    qm2link = 0
    min_at  = 0
    rmin    = 0.0_dp

    if ( .not. got_numlink ) &
      call die( "link1: Implementation error. Need to call get_numlink first." )

    default_link_symbol = trim( fdf_get('DefaultLinkSymbol','H?') )
    default_link_isa    = 0
    default_link_iza    = 0

    if ( fdf_block( 'ChemicalSpeciesLabel', bfdf ) ) then
      do i = 1, nspec
        if (.not. fdf_bline(bfdf,pline)) &
         call die('link1: ERROR in ChemicalSpeciesLabel block')

        if ( fdf_bnintegers(pline) < 2 ) &
          call die( "link1: problem reading QM species info." )
        if ( fdf_bnnames(pline) < 1 ) &
          call die( "link1: problem reading QM species name." )

        j        = fdf_bintegers(pline,1)
        atnum(i) = fdf_bintegers(pline,2)
        atsym(i) = trim(fdf_bnames(pline,1))
      enddo
      call fdf_bclose( bfdf )
    else
      call die("link1: You must specify the atomic labels.")
    endif

    do k = 1, nspec
      atype = atsym(k)
      if ( default_link_symbol(1:1) == atype(1:1) ) then
        default_link_isa = k
        default_link_iza = atnum(k)
      endif
    enddo

    if ( fdf_block( 'SpecialLinkAtoms', bfdf ) ) then
      if (.not. fdf_bline(bfdf,pline)) &
        call die('link1: ERROR in SpecialLinkAtoms block')

      if ( fdf_bnintegers(pline) < 2 ) &
          call die( "link1: problem reading number of resonances or" //&
                    " number of atoms in SpecialLinkAtoms." )

      NSpecialLinkAtoms = fdf_bintegers(pline, 1)
      num_resonances    = fdf_bintegers(pline, 2)

      if ( NSpecialLinkAtoms > 0 ) then
        if ( allocated( sp_link_atoms ) ) then
          do i = 1, NSpecialLinkAtoms
            if ( allocated(sp_link_atoms(i)%reson) ) &
              deallocate( sp_link_atoms(i)%reson )
          enddo
          deallocate( sp_link_atoms )
        endif

        allocate( sp_link_atoms(NSpecialLinkAtoms) )
        do i = 1, NSpecialLinkAtoms
          allocate( sp_link_atoms(i)%reson(num_resonances) )
        enddo

        call read_special_link_atoms( bfdf )
        call fdf_bclose( bfdf )

        do i = 1, num_resonances
          resonance(i)%weight = resonance(i)%sp_weight
        enddo

        ! Assignation of CQM and CMM for special link atoms
        do i = 1, NSpecialLinkAtoms
          link_atoms(i)%qm(1)          = sp_link_atoms(i)%qm
          qm2link(sp_link_atoms(i)%qm) = 1
          link_atoms(i)%mm(1)          = sp_link_atoms(i)%mm!-na_qm
        enddo

        do i = 1, num_resonances
          ! Assignates atomic number (iza)
          do j = 1, NSpecialLinkAtoms
            resonance(i)%atype(j) = sp_link_atoms(j)%reson(i)%attype
            resonance(i)%isa(j)   = sp_link_atoms(j)%reson(i)%isa
            do k = 1, nspec
              if ( resonance(i)%isa(j) == k ) resonance(i)%iza(j) = atnum(k)
            enddo
            if ( resonance(i)%iza(j) == 0 ) &
              call die( "link1: Atoms with atomic number 0 or no number." )
          enddo
        enddo
      endif
    endif

    if ( numlink > NSpecialLinkAtoms ) then
      do j = 1, num_resonances
      do i = NSpecialLinkAtoms+1, numlink
        resonance(j)%atype(i) = default_link_symbol
        resonance(j)%isa(i)   = default_link_isa
        resonance(j)%iza(i)   = default_link_iza
      enddo
      enddo
    endif

    k = NSpecialLinkAtoms
    do i = 1, na_qm ! assignation of CQM and CMM
      if ( qm2link(i) == 1 ) cycle

      dr = qm_atoms(i)%r(:) - mm_atoms(1)%r(:)
      call pbc_displ_vector( lattice_type, cell, kcell, dr )
      rmin(i)   = norm_v2( dr )
      min_at(i) = 1

      do j = 2, na_mm
        dr = qm_atoms(i)%r(:) - mm_atoms(j)%r(:)
        call pbc_displ_vector( lattice_type, cell, kcell, dr )
        r1 = norm_v2( dr )
        if ( r1 <= rmin(i) ) then
          rmin(i)   = r1
          min_at(i) = j
        endif
      enddo

      ! Only consider pairs of atoms close enough
      if ( rmin(i) > 2.8_dp * Ang ) cycle
      chi = qm_atoms(i)%attype
      if ( chi(1:1) == 'H' ) cycle

      j   = min_at(i)
      chj = mm_atoms(j)%attype
      if ( chj(1:1) == 'H' ) cycle

      parameter_set = .false.
      qmmm_range    = 0.0d0
      do l = 1, mm_top%nbonds
        tybond = mm_top%bonds(l)%type
        if ( (chi(1:2) == tybond(1:2)) .and. (chj(1:2) == tybond(4:5)) ) then
          qmmm_range    = 1.2_dp * Ang * mm_top%bonds(l)%r_eq
          parameter_set = .true.
          exit
        elseif ((chi(1:2) == tybond(4:5)) .and. (chj(1:2) == tybond(1:2))) then
          qmmm_range    = 1.2_dp * Ang * mm_top%bonds(l)%r_eq
          parameter_set = .true.
          exit
        endif
      enddo
      if ( .not. parameter_set ) then
        if ( check_out_of_range( i, j, chi, chj, rmin(i) / Ang ) ) cycle
      endif

      ch4 = qm_atoms(i)%attype
      ch1 = ch4(1:1)
      if ( rmin(i) < qmmm_range ) then
        if ( deb ) write(6,*) "link1: Checking link atoms for QM atom no. ", i
        if ( ch1 == 'O') then
          if ( deb ) write(6,*) i, min_at(i), rmin(i)
          rmin(i) = 20.0_dp

          k = k+1
          if ( k > numlink ) call die( 'link1: ERROR - k > numlink' )

          link_atoms(k)%qm(1) = i
          link_atoms(k)%mm(1) = min_at(i)

        elseif ( (ch1 == 'C') .or. (ch1 == 'N') ) then
          if ( deb ) write(6,*) i, min_at(i), rmin(i)
          rmin(i) = 20.0_dp

          k = k +1
          if ( k > numlink ) call die( 'link1: ERROR - k > numlink' )

          link_atoms(k)%qm(1) = i
          link_atoms(k)%mm(1) = min_at(i)
        endif
      endif
    enddo ! QM atoms

    if ( k /= numlink ) then
      write(6,'(a)') "Total amount of link atoms scanned different from numlink."
      write(6,'(a)') "Total k = ", k, ', numlink = ',numlink
      call die( 'link1: ERROR - k /= numlink' )
    endif

    ! assignation of QM 1st neighbors
    rmin   = 0.0_dp
    min_at = 0
    do i = 1, numlink
      m   = 0
      ch4 = qm_atoms(link_atoms(i)%qm(1))%attype
      ch1 = ch4(1:1)

      ! sp oxygen or sp2 carbon/nitrogen
      if ( ch1 == 'O' ) then
        m = 2
      elseif ( (ch1 == 'C') .or. (ch1 == 'N') ) then
        m = 3
      endif

      ! sp3 carbon
      ch4 = qm_atoms(link_atoms(i)%qm(1))%attype
      ch2 = ch4(1:2)
      if ( ch2 == 'CT' ) m = 4
      if ( m == 0 ) then
        write(6,*) 'link1: Wrong link atom QM atom type.... Check parameters.'
        write(6,*) 'link1: QM atom no. ', link_atoms(i)%qm(1)
        call die( 'link1: ERROR - Wrong link atom QM atom type.' )
      endif

      ! loop over CQM neighbors
      do k = 2, m
        rmin(i)   = 10.0_dp
        min_at(i) = 0

        do j = 1, na_qm
          if ( (j == link_atoms(i)%qm(1)) .or. (j == link_atoms(i)%qm(2)) .or. &
               (j == link_atoms(i)%qm(3)) ) cycle

          dr = qm_atoms(link_atoms(i)%qm(1))%r(:) - qm_atoms(j)%r(:)
          call pbc_displ_vector( lattice_type, cell, kcell, dr )
          r1 = norm_v2( dr )

          if ( r1 <= rmin(i) ) then
            rmin(i)   = r1
            min_at(i) = j
          endif
        enddo

        if ( rmin(i) > 4.0_dp ) then
          write(6,*) 'link1: Wrong Link atom QM neighbour... Check geometry.'
          write(6,*) 'link1: rmin(i)= ', j, rmin(i)
          call die( 'link1: ERROR - Wrong link atom QM neighbour.' )
        endif

        link_atoms(i)%qm(k) = min_at(i)
      enddo ! k neighbours
    enddo ! link atoms

    do i = 1, numlink ! checking CQM neighbors
      ! sp2 carbon
      ch4 = qm_atoms(link_atoms(i)%qm(1))%attype
      ch2 = ch4(1:2)
      ch1 = ch4(1:1)
      if ( ch1 == 'O' ) then
        if ( (link_atoms(i)%qm(2) == 0) .and. (link_atoms(i)%qm(3) == 0) .and. &
             (link_atoms(i)%qm(4) == 0) ) &
          call die( 'link1: Wrong QM neighbor number for a sp oxygen.' )
      elseif ( (ch1 == 'C')) then
        if ( ((link_atoms(i)%qm(2) == 0) .and. (link_atoms(i)%qm(3) == 0)) .or. &
             ((link_atoms(i)%qm(2) == 0) .and. (link_atoms(i)%qm(4) == 0)) .or. &
             ((link_atoms(i)%qm(3) == 0) .and. (link_atoms(i)%qm(4) == 0)) ) &
          call die( 'link1: Wrong QM neighbor number for a sp2 carbon.' )
      endif

      ! sp3 carbon
      if ( ch2 == 'CT' ) then
        if ( (link_atoms(i)%qm(2) == 0) .or. (link_atoms(i)%qm(3) == 0) .or. &
             (link_atoms(i)%qm(4) == 0) ) &
          call die( 'link1: Wrong QM neighbor number for a sp3 carbon.' )
      endif
    enddo


    do i = 1, numlink

      ! asignation of QM link atoms types
      do j = 1, 4

        if ( link_atoms(i)%qm(j) == 0 ) then
          link_atoms(i)%qm_type(j) = 'XX'
        else
          link_atoms(i)%qm_type(j) = qm_atoms( link_atoms(i)%qm(j) )%attype
        endif
      enddo

      ! assignation of MM 1st neighbors
      do k = 2, 4
        link_atoms(i)%mm(k) = mm_connectivity(link_atoms(i)%mm(1))%bond_at(k-1)
      enddo

      ! assignation of MM 2nd neighbors
      do j = 2, 4
        if ( link_atoms(i)%mm(j) == 0 ) cycle
        m = 1
        do k = 1, 4
          if ( mm_connectivity(link_atoms(i)%mm(j))%bond_at(k) /= link_atoms(i)%mm(1) ) then
            link_atoms(i)%mm2(j,m) = mm_connectivity(link_atoms(i)%mm(j))%bond_at(k)
            m = m +1
          endif
        enddo
      enddo
    enddo

    link_atoms(:)%num_mm = 0
    link_atoms(:)%num_qm = 0
    do i = 1, numlink
    do j = 1, 4
      if ( link_atoms(i)%qm(j) /= 0 ) link_atoms(i)%num_qm = link_atoms(i)%num_qm +1
      if ( link_atoms(i)%mm(j) /= 0 ) link_atoms(i)%num_mm = link_atoms(i)%num_mm +1
    enddo
    enddo

    ! write LA params in file
    write(6,'(/,a)') 'siesta-qmmm: Link atom parameters:'
    if ( h_link_scaling_type == 'splam' ) then
      write(6,*)'  The link atom will be scaled by the modified version '//&
                'of the SPLAM method.'
    else
      write(6,*)'  The link atom will be located at a fix distance from '//&
                'the QM atom.'
    endif

    write(6,*)
    do i = 1, numlink
       write(6,'(a,2x,1I6)') 'linkatom:', na_qm+na_mm+i
       write(6,'(a,2x,4I6)') &
         'qmatoms: ', ( link_atoms(i)%qm(j)        , j=1, link_atoms(i)%num_qm )
       write(6,'(a,2x,4A4)') &
         'qmtypes: ', ( link_atoms(i)%qm_type(j)    , j=1, 4 )
       write(6,'(a,2x,4I6)') &
         'mmatoms: ', ( link_atoms(i)%mm(j)        , j=1, link_atoms(i)%num_mm )
       write(6,'(a,2x,4A4)') &
         'mmtypes: ', ( mm_atoms(link_atoms(i)%mm(j))%attype, j=1, link_atoms(i)%num_mm )
    enddo
    write(6,*)

    call de_alloc( qm2link, 'qm2link', 'link1' )
    call de_alloc( min_at , 'min_at' , 'link1' )
    call de_alloc( rmin   , 'rmin'   , 'link1' )
  end subroutine link1

  subroutine get_link_ff_parameters( n_mm, mm_atoms, mm_top )
    use alloc      , only : re_alloc, de_alloc
    use mm_topology, only : fftopology, mm_atom_t
    use precision  , only : dp
    use sys        , only : die

    integer         , intent(in) :: n_mm
    type(mm_atom_t) , intent(in) :: mm_atoms(n_mm)
    type(fftopology), intent(in) :: mm_top


    character(len=4)  :: ty1_h, ty2_h, ty3_h, tyl_h, ty1, ty2, ty3, ty4
    character(len=5)  :: tybond_h, tybond
    character(len=8)  :: tyangle_h, tyangle
    character(len=11) :: tydihe
    integer           :: resonance_id, i, j, jm, jq, k, l, m, parameter_bond_set
    logical           :: parameter_set

    real(dp), pointer :: perdihe2(:)

    nullify( perdihe2 )
    call re_alloc( perdihe2, 1, mm_top%ndihe, 'perdihe2', 'get_link_ff' )
    do i = 1, mm_top%ndihe
      perdihe2(i) = mm_top%dihedrals(i)%per
    enddo

    do resonance_id=1, num_resonances
      do i = 1, numlink ! asignation for E and F
        parameter_set = .false.
        do k = 1, mm_top%nbonds     !	bond Cqm -- Cmm parameter 1
          tybond = mm_top%bonds(k)%type
          ty1    = tybond(1:2)
          ty2    = tybond(4:5)

          if ((link_atoms(i)%qm_type(1) == ty1) .and. &
              (mm_atoms(link_atoms(i)%mm(1))%attype == ty2)) then
            link_atoms(i)%reson(resonance_id)%params(1,1) = k
            parameter_set = .true.
          elseif ( (link_atoms(i)%qm_type(1)     == ty2) .and. &
                   (mm_atoms(link_atoms(i)%mm(1))%attype == ty1) ) then
            link_atoms(i)%reson(resonance_id)%params(1,1) = k
            parameter_set = .true.
          endif
        enddo

        if ( .not. parameter_set ) then
          write(6,*) 'get_link_ff_parameters: Missing parameters for bond: ', &
                      link_atoms(i)%qm_type(1), '-',mm_atoms(link_atoms(i)%mm(1))%attype
          call die( 'No suitable parameter for Cqm-Cmm bond' )
        endif

        parameter_set = .false.  ! 	angles Cqm -- Cqm -- Cmm parameters 2 to 4
        do j = 2, 4
          if ( link_atoms(i)%qm(j) == 0 ) cycle

          do k = 1, mm_top%nangles
            tyangle = mm_top%angles(k)%type
            ty1     = tyangle(1:2)
            ty2     = tyangle(4:5)
            ty3     = tyangle(7:8)

            if ( (link_atoms(i)%qm_type(j) == ty1) .and. (link_atoms(i)%qm_type(1) == ty2) .and. &
                 (mm_atoms(link_atoms(i)%mm(1))%attype == ty3) ) then
              link_atoms(i)%reson(resonance_id)%params(j,1) = k
              parameter_set = .true.
            elseif ( (link_atoms(i)%qm_type(j) == ty3) .and. &
                     (link_atoms(i)%qm_type(1) == ty2) .and. &
                     (mm_atoms(link_atoms(i)%mm(1))%attype == ty1) ) then
              link_atoms(i)%reson(resonance_id)%params(j,1) = k
              parameter_set = .true.
            endif
          enddo
        enddo

        if ( .not. parameter_set ) then
          write(6,*) 'get_link_ff_parameters: Missing parameters for angle: ',&
               link_atoms(i)%qm_type(j), '-', link_atoms(i)%qm_type(1), '-', &
               mm_atoms(link_atoms(i)%mm(1))%attype
          call die('No suitable parameter for Cqm-Cqm-Cm angle' )
        endif
      enddo ! numlink

      do i = 1, numlink
        tyl_h = resonance(resonance_id)%atype(i)

        if ( '?' == tyl_h(2:2) ) then
          parameter_bond_set = 0
          do k = 1, mm_top%nbonds         !	bond Cqm -- Cmm parameter 1
            parameter_set = .false.
            tybond_h      = mm_top%bonds(k)%type
            ty1_h         = tybond_h(1:2)
            ty2_h         = tybond_h(4:5)

            if ( (link_atoms(i)%qm_type(1) == ty1_h) .and. &
                 (tyl_h(1:1) == ty2_h(1:1)) ) then
              link_atoms(i)%reson(resonance_id)%bond_r = mm_top%bonds(k)%r_eq
              link_atoms(i)%reson(resonance_id)%bond_k  = &
                 sqrt( mm_top%bonds( link_atoms(i)%reson(resonance_id)%params(1,1) )%k &
                                                   / mm_top%bonds(k)%k )
              parameter_set = .true.
              parameter_bond_set = k
            elseif ( (link_atoms(i)%qm_type(1) == ty2_h) .and. &
                     (tyl_h(1:1) == ty1_h(1:1)) ) then
              link_atoms(i)%reson(resonance_id)%bond_r = mm_top%bonds(k)%r_eq
              link_atoms(i)%reson(resonance_id)%bond_k  = &
                 sqrt( mm_top%bonds( link_atoms(i)%reson(resonance_id)%params(1,1) )%k &
                                                   / mm_top%bonds(k)%k )
              parameter_set = .true.
              parameter_bond_set = k
            endif
            if ( .not. parameter_set ) cycle

            do j = 2, 4
              if ( link_atoms(i)%qm(j) == 0 ) cycle
              parameter_set = .false.

              do l = 1, mm_top%nangles
                tyangle_h = mm_top%angles(l)%type
                ty1_h = tyangle_h(1:2)
                ty2_h = tyangle_h(4:5)
                ty3_h = tyangle_h(7:8)
                tyl_h = resonance(resonance_id)%atype(i)

                if ( (link_atoms(i)%qm_type(j) == ty1_h) .and. &
                     (link_atoms(i)%qm_type(1) == ty2_h) .and. &
                     (tyl_h(1:1) == ty3_h(1:1)) ) then
                  link_atoms(i)%reson(resonance_id)%angle_k(j-1) =&
                    mm_top%angles(l)%k
                  parameter_set = .true.
                  exit
                elseif ( (link_atoms(i)%qm_type(j) == ty3_h) .and. &
                         (link_atoms(i)%qm_type(1) == ty2_h) .and. &
                         (tyl_h(1:1) == ty1_h(1:1)) ) then
                  link_atoms(i)%reson(resonance_id)%angle_k(j-1) =&
                    mm_top%angles(l)%k
                  parameter_set = .true.
                  exit
                endif
              enddo
              if ( .not. parameter_set ) exit
            enddo
            if (parameter_set) exit
          enddo

          if ( .not. parameter_set ) then
            write(6,*)
            write(6,*) 'WARNING - get_link_ff_parameters'
            write(6,*) 'In Cqm-Cmm bond: ', link_atoms(i)%qm_type(1), &
                       '-', mm_atoms(link_atoms(i)%mm(1))%attype

            if ( parameter_bond_set == 0 ) then
              write(6,*) 'get_link_ff_parameters:'//&
                 'Default parameters set for Cqm-LA bond: ', &
                 link_atoms(i)%qm_type(1), '-', resonance(resonance_id)%atype(i)
              call set_default_H_link_bond_parameters( &
                 link_atoms(i)%qm_type(1), link_atoms(i)%reson(resonance_id)%bond_r,  &
                 mm_top%bonds(link_atoms(i)%reson(resonance_id)%params(1,1))%k, &
                 link_atoms(i)%reson(resonance_id)%bond_k )
            endif

            write(6,*) 'get_link_ff_parameters: '//&
                       'Default parameters are set for angles:'
            do j = 2, 4
              if ( link_atoms(i)%qm(j) == 0 ) cycle
              link_atoms(i)%reson(resonance_id)%angle_k(j-1) = &
                 mm_top%angles(link_atoms(i)%reson(resonance_id)%params(j,1))%k
              write(6,*) 'get_link_ff_parameters: Angle: ',  &
                 link_atoms(i)%qm_type(j), '-', link_atoms(i)%qm_type(1), '-', &
                 resonance(resonance_id)%atype(i)
            enddo
          endif

        else ! if ('?'/=tyl_h(2))
          parameter_set = .false.
          do l = 1, mm_top%nbonds
            tybond_h = mm_top%bonds(l)%type
            ty1_h    = tybond_h(1:2)
            ty2_h    = tybond_h(4:5)
            tyl_h    = resonance(resonance_id)%atype(i)

            if ( (link_atoms(i)%qm_type(1) == ty1_h) .and. &
                 (tyl_h(1:2) == ty2_h(1:2)) ) then
              link_atoms(i)%reson(resonance_id)%bond_r = mm_top%bonds(l)%r_eq
              link_atoms(i)%reson(resonance_id)%bond_k  = &
                  sqrt( mm_top%bonds(k)%k / mm_top%bonds(l)%k)
              parameter_set = .true.
              exit
            elseif ( (link_atoms(i)%qm_type(1) == ty2_h) .and. &
                     (tyl_h(1:2) == ty1_h(1:2)) ) then
              link_atoms(i)%reson(resonance_id)%bond_r = mm_top%bonds(l)%r_eq
              link_atoms(i)%reson(resonance_id)%bond_k  = &
                  sqrt(mm_top%bonds(k)%k/mm_top%bonds(l)%k)
              parameter_set = .true.
              exit
            endif
          enddo

          if ( .not. parameter_set ) then
            write(6,*) 'get_link_ff_parameters: Missing parameters for bond: ',&
                       link_atoms(i)%qm_type(1), '-', resonance(resonance_id)%atype(i)
            call die( 'ERROR: No suitable parameter for Cqm-LA bond.' )
          endif

          do j = 2, 4
            if ( link_atoms(i)%qm(j) == 0 ) cycle

            parameter_set = .false.
            do l = 1, mm_top%nangles  !angles Cqm -- Cqm -- LA parameters 2 to 4
              tyangle_h = mm_top%angles(l)%type
              ty1_h = tyangle_h(1:2)
              ty2_h = tyangle_h(4:5)
              ty3_h = tyangle_h(7:8)
              tyl_h = resonance(resonance_id)%atype(i)

              if ( (link_atoms(i)%qm_type(j) == ty1_h) .and. &
                   (link_atoms(i)%qm_type(1) == ty2_h) .and. &
                   (tyl_h(1:2) == ty3_h(1:2)) ) then
                link_atoms(i)%reson(resonance_id)%angle_k(j-1) = mm_top%angles(l)%k
                parameter_set = .true.
                exit
              elseif ( (link_atoms(i)%qm_type(j) == ty3_h) .and. &
                       (link_atoms(i)%qm_type(1) == ty2_h) .and. &
                       (tyl_h(1:2) == ty1_h(1:2)) ) then
                link_atoms(i)%reson(resonance_id)%angle_k(j-1) = mm_top%angles(l)%k
                parameter_set = .true.
                exit
              endif
            enddo

            if ( .not. parameter_set ) then
              write(6,*) 'get_link_ff_parameters: Missing parameters'//&
                ' for angle: ', link_atoms(i)%qm_type(j), '-', link_atoms(i)%qm_type(1), &
                '-', resonance(resonance_id)%atype(i)
              call die( 'No suitable parameter for Cqm-Cqm-LA angle.' )
            endif
          enddo  ! j=2,4
        endif  !  if ('?'==tyl_h(2))
      enddo ! numlink

      do i = 1, numlink ! 	angles Cqm -- Cmm -- X parameters 5 to 7
      do j = 2, 4
        if ( link_atoms(i)%mm(j) == 0 ) cycle
        do k = 1, mm_top%nangles
          tyangle = mm_top%angles(k)%type
          ty1 = tyangle(1:2)
          ty2 = tyangle(4:5)
          ty3 = tyangle(7:8)

          if ( (link_atoms(i)%qm_type(1) == ty1)     .and. &
               (mm_atoms(link_atoms(i)%mm(1))%attype == ty2) .and. &
               (mm_atoms(link_atoms(i)%mm(j))%attype == ty3) ) then
            link_atoms(i)%reson(resonance_id)%params(3+j,1) = k
          elseif ( (link_atoms(i)%qm_type(1) == ty3)     .and. &
                   (mm_atoms(link_atoms(i)%mm(1))%attype == ty2) .and. &
                   (mm_atoms(link_atoms(i)%mm(j))%attype == ty1) ) then
            link_atoms(i)%reson(resonance_id)%params(3+j,1) = k
          endif
        enddo
      enddo
      enddo

      !	dihedral X--Cq -- Cmm--Y parameters 5 to 13
      do i = 1, numlink
      do jq = 2, 4
        if ( link_atoms(i)%qm(jq) == 0 ) cycle ! no more neighbours

        do jm = 2, 4
          if ( link_atoms(i)%qm(jm) == 0 ) cycle ! only two MM neighbours

          m = 0
          do k = 1, mm_top%ndihe
            tydihe = mm_top%dihedrals(k)%type
            ty1 = tydihe(1:2)
            ty2 = tydihe(4:5)
            ty3 = tydihe(7:8)
            ty4 = tydihe(10:11)

            if ( ty1 == 'X ' ) then ! if X the same constant
              if ( (link_atoms(i)%qm_type(1) == ty2) .and. &
                   (mm_atoms(link_atoms(i)%mm(1))%attype == ty3) ) then
                link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,1) = k
              elseif ( (link_atoms(i)%qm_type(1) == ty3) .and. &
                       (mm_atoms(link_atoms(i)%mm(1))%attype == ty2) )  then
                link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,1) = k
              endif

            else ! if not X different constants
              if ( (link_atoms(i)%qm_type(jq)     == ty1) .and. &
                   (link_atoms(i)%qm_type(1)      == ty2) .and. &
                   (mm_atoms(link_atoms(i)%mm(1))%attype  == ty3) .and. &
                   (mm_atoms(link_atoms(i)%mm(jm))%attype == ty4) ) then
                link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,1) = k

                if ( perdihe2(k) < 0 ) then
                  m = m +1
                  link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,m)   = k
                  link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,m+1) = k+1
                endif
              elseif ( (link_atoms(i)%qm_type(jq) == ty4) .and. &
                       (link_atoms(i)%qm_type(1)  == ty3) .and. &
                       (mm_atoms(link_atoms(i)%mm(1))%attype  == ty2).and.&
                       (mm_atoms(link_atoms(i)%mm(jm))%attype == ty1) ) then
                link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,1) = k

                if ( perdihe2(k) < 0 ) then
                  m = m+1
                  link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,m)   = k
                  link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,m+1) = k+1
                endif
              endif
            endif
          enddo ! k
        enddo   ! jm
      enddo     ! jq
      enddo     ! linkatoms

      !	dihedral Cq--Cmm--Y--Z parameters 14 to 22
      do i  = 1, numlink
      do jm = 2, 4
        if ( link_atoms(i)%mm(jm) == 0 ) cycle

        do j = 1, 3
          if ( link_atoms(i)%mm2(jm,j) == 0 ) cycle

          m = 0
          do k = 1, mm_top%ndihe
            tydihe = mm_top%dihedrals(k)%type
            ty1 = tydihe(1:2)
            ty2 = tydihe(4:5)
            ty3 = tydihe(7:8)
            ty4 = tydihe(10:11)

            if ( ty1 == 'X ' ) then ! if X the same constant
              if ( (mm_atoms(link_atoms(i)%mm(1))%attype  == ty2) .and. &
                   (mm_atoms(link_atoms(i)%mm(jm))%attype == ty3) )  then
                link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,1) = k
              elseif ( (mm_atoms(link_atoms(i)%mm(1))%attype  == ty3) .and. &
                       (mm_atoms(link_atoms(i)%mm(jm))%attype == ty2) ) then
                link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,1) = k
              endif

            else ! if not X different constants
              if ( (link_atoms(i)%qm_type(1)         == ty1) .and. &
                   (mm_atoms(link_atoms(i)%mm(1))%attype     == ty2) .and. &
                   (mm_atoms(link_atoms(i)%mm(jm))%attype    == ty3) .and. &
                   (mm_atoms(link_atoms(i)%mm2(jm,j))%attype == ty4) ) then

                if ( perdihe2(k) < 0 ) then
                  link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,1) = k

                  if ( perdihe2(k+1) < 0 ) then
                    link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,2) = k+1

                    if ( perdihe2(k+2) < 0 ) then
                      link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,3) = k+2
                      link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,4) = k+3
                    else
                      link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,3) = k+2
                      exit
                    endif
                  else
                    link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,2) = k+1
                    exit
                  endif

                else
                  link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,1) = k
                  exit
                endif

              elseif ( (link_atoms(i)%qm_type(1)         == ty4) .and. &
                       (mm_atoms(link_atoms(i)%mm(1))%attype     == ty3) .and. &
                       (mm_atoms(link_atoms(i)%mm(jm))%attype    == ty2) .and. &
                       (mm_atoms(link_atoms(i)%mm2(jm,j))%attype == ty1) ) then
                if ( perdihe2(k) < 0 ) then
                  link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,1) = k

                  if ( perdihe2(k+1) < 0 ) then
                    link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,2) = k+1

                    if ( perdihe2(k+2) < 0 ) then
                      link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,3) = k+2
                      link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,4) = k+3
                    else
                      link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,3) = k+2
                      exit
                    endif
                  else
                    link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,2) = k+1
                    exit
                  endif
                else
                  link_atoms(i)%reson(resonance_id)%params(13+3*(jm-2)+j,1) = k
                  exit
                endif
              endif ! one side or other
            endif ! with or without X
          enddo ! end dihedrals
        enddo ! for neigh 2, neigh 1, numlink
      enddo
      enddo ! numlink
    enddo  ! resonance_id

    call de_alloc( perdihe2, 'perdihe2', 'get_link_ff' )
  end subroutine get_link_ff_parameters

  ! calculates Energy and Forces of HLink.
  subroutine link2( r_mm, natot, n_mm, frc, stress, mm_top, &
                    Ener, lattice_type, cell )

    use alloc      , only : de_alloc, re_alloc
    use functions  , only : angle_v2, scalar_v2, norm_v2
    use mm_topology, only : fftopology
    use qmmm_pbc   , only : pbc_displ_vector, reccel
    use units      , only : Ang, pi

    type(fftopology), intent(in) :: mm_top
    integer  , intent(in)    :: natot
    integer  , intent(in)    :: n_mm
    real(dp) , intent(in)    :: r_mm(3,natot)
    real(dp) , intent(in)    :: cell(3,3)
    character, intent(in)    :: lattice_type
    real(dp) , intent(out)   :: Ener
    real(dp) , intent(inout) :: frc(3,natot)
    real(dp) , intent(inout) :: stress(3,3)

    integer  :: i, j, k, l, m, jm, jq, at_id(4), resonance_id
    real(dp) :: Elink_for_h_atom, fce(3), rij, scal, r12,     &
                r32, dr12r32, angulo, fpp(3), ip, dversor(3), &
                rhq, r_qm, dr(3), dang(3), dr12(3), dr32(3),   &
                stress_fact, h_kangle, amber_cell(3,3),       &
                amber_kcell(3,3), ang_fac, dr43(3), dr42(3), dr31(3)

    real(dp), pointer  :: perdihe2(:), flink(:,:)
    real(dp), external :: volcel

    nullify( perdihe2, flink )
    call re_alloc( perdihe2, 1, mm_top%ndihe, 'perdihe2', 'link2')
    call re_alloc( flink   , 1, 3, 1, natot , 'flink'   , 'link2')

    do i = 1, mm_top%ndihe
      perdihe2(i) = mm_top%dihedrals(i)%per
    enddo

    ! link mm2 asignment: linkmm atoms neighbours
    ! change units
    amber_cell(:,:) = cell(:,:) / Ang
    call reccel( 3, amber_cell, amber_kcell, 0 )

    ! reasignation of perdihe2
    do i = 1, mm_top%ndihe
      if ( perdihe2(i) < 0 ) perdihe2(i) = -perdihe2(i)
    enddo

    ang_fac     = pi / 180.0_dp
    stress_fact = 1.0_dp / (Ang * volcel(amber_cell))

    do resonance_id = 1, num_resonances
      Ener = 0.0_dp
      flink(1:3,1:natot) = 0.0_dp

      do i = 1, numlink ! calculation of E and F
        Elink_for_h_atom = 0.0_dp

        if ( h_link_scaling_type /= 'splam' ) then
          ! correction to bond Cqm -- X parameter 1
          at_id(1) = link_atoms(i)%qm(1)
          at_id(2) = natot - n_mm + link_atoms(i)%mm(1)
          k   = link_atoms(i)%reson(resonance_id)%params(1,1)

          dr = (r_mm(:,at_id(1)) - r_mm(:,at_id(2))) / Ang
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr )

          rij = norm_v2( dr )
          Elink_for_h_atom = Elink_for_h_atom + mm_top%bonds(k)%energy( dr )
          call mm_top%bonds(k)%grad( dr, dang )

          flink(:,at_id(1)) = -dang(:)
          flink(:,at_id(2)) =  dang(:)
          do m = 1, 3
            stress(m,:) = stress(m,:) + stress_fact * r_mm(m,at_id(1)) * dang(:)
          enddo
        endif

        ! angle  X -- Cqm -- Cmm parameters 2 to 4
        !	3 angles ( one for each x) j, x index
        do j = 2, 4
          if ( link_atoms(i)%qm(j) == 0 ) cycle

          at_id(1) = natot - n_mm + link_atoms(i)%mm(1)
          at_id(2) = link_atoms(i)%qm(1)
          at_id(3) = link_atoms(i)%qm(j)
          k = link_atoms(i)%reson(resonance_id)%params(j,1)

          dr12(:) = ( r_mm(:,at_id(1)) - r_mm(:,at_id(2)) ) / Ang
          dr32(:) = ( r_mm(:,at_id(3)) - r_mm(:,at_id(2)) ) / Ang
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr12 )
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr32 )
          angulo = angle_v2( dr12, dr32 )

          h_kangle = mm_top%angles(k)%k - link_atoms(i)%reson(resonance_id)%angle_k(j-1)

          Elink_for_h_atom = Elink_for_h_atom + h_kangle * &
                             ( (angulo - mm_top%angles(k)%r_eq) * ang_fac ) **2

          scal = scalar_v2( dr12, dr32 )
          r12  = norm_v2( dr12 )
          r32  = norm_v2( dr32 )
          dr12r32 = r12 * r32

          dang(:) = ( dr32(:) - scal * dr12(:) / (r12*r12) ) / dr12r32
          dang(:) = - dang(:) / ( sqrt( 1.0_dp - (scal / dr12r32)**2 ) )
          dang(:) = 2.0_dp * h_kangle * ang_fac * dang(:) * &
                    ( angulo - mm_top%angles(k)%r_eq )

          ! atom qm cq force, sum of each particular x
          flink(:,at_id(1)) = flink(:,at_id(1)) - dang(:)

          ! force over x, change at_id(3) by at_id(1)
          dang(:) = ( dr12(:) - scal * dr32(:) / (r32*r32) ) / dr12r32
          dang(:) = - dang(:) / ( sqrt( 1.0_dp - (scal / dr12r32)**2 ) )
          dang(:) = 2.0_dp * h_kangle * ang_fac * dang(:) * &
                    ( angulo - mm_top%angles(k)%r_eq )

          ! sum each x
          flink(:,at_id(3)) = flink(:,at_id(3)) - dang(:)

          ! middle atom force (change in derivate)
          dang(:) = scal * ( r32 * dr12(:) / r12 + r12 * dr32(:) / r32 ) &
                    - (dr12(:) + dr32(:)) * r12 * r32
          dang(:) = dang(:) / dr12r32 * dr12r32
          dang(:) = - dang(:) / ( sqrt( 1.0_dp - (scal / dr12r32)**2 ) )
          dang(:) = 2.0_dp * h_kangle * ang_fac * dang(:) * &
                    ( angulo - mm_top%angles(k)%r_eq )

          ! central atom force
          flink(:,at_id(2)) = flink(:,at_id(2)) - dang(:)
        enddo ! enddo for each possible atom in angle

        ! angle  Cqm -- Cmm -- X parameters 5 to 7
        !	3 angles ( one for each x) j, x index
        do j = 2, 4
          if ( link_atoms(i)%mm(j) == 0 ) cycle

          at_id(1) = link_atoms(i)%qm(1)
          at_id(2) = natot - n_mm + link_atoms(i)%mm(1)
          at_id(3) = natot - n_mm + link_atoms(i)%mm(j)
          k = link_atoms(i)%reson(resonance_id)%params(3+j,1)

          dr12(:) = ( r_mm(:,at_id(1)) - r_mm(:,at_id(2)) ) / Ang
          dr32(:) = ( r_mm(:,at_id(3)) - r_mm(:,at_id(2)) ) / Ang
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr12 )
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr32 )
          angulo = angle_v2( dr12, dr32 )

          Elink_for_h_atom = Elink_for_h_atom + mm_top%angles(k)%k * &
                        ((angulo - mm_top%angles(k)%r_eq)* ang_fac)**2

          scal = scalar_v2( dr12, dr32 )

          r12 = norm_v2( dr12 )
          r32 = norm_v2( dr32 )
          dr12r32  = r12 * r32
          h_kangle = mm_top%angles(k)%k

          dang(:) = ( dr32(:) - scal * dr12(:) / (r12*r12) ) / dr12r32
          dang(:) = - dang(:) / ( sqrt( 1.0_dp - (scal / dr12r32)**2 ) )
          dang(:) = 2.0_dp * h_kangle * ang_fac * dang(:) * &
                    ( angulo - mm_top%angles(k)%r_eq )

          ! atom qm cq force, sum of each particular x
          flink(:,at_id(1)) = flink(:,at_id(1)) - dang(:)


          ! force over x, change at_id(3) by at_id(1)
          dang(:) = ( dr12(:) - scal * dr32(:) / (r32*r32) ) / dr12r32
          dang(:) = - dang(:) / ( sqrt( 1.0_dp - (scal / dr12r32)**2 ) )
          dang(:) = 2.0_dp * h_kangle * ang_fac * dang(:) * &
                    ( angulo - mm_top%angles(k)%r_eq )

          ! sum each x
          flink(:,at_id(3)) = flink(:,at_id(3)) - dang(:)

          ! middle atom force (change in derivate)
          dang(:) = scal * ( r32 * dr12(:) / r12 + r12 * dr32(:) / r32 ) &
                    - (dr12(:) + dr32(:)) * r12 * r32
          dang(:) = dang(:) / dr12r32 * dr12r32
          dang(:) = - dang(:) / ( sqrt( 1.0_dp - (scal / dr12r32)**2 ) )
          dang(:) = 2.0_dp * h_kangle * ang_fac * dang(:) * &
                    ( angulo - mm_top%angles(k)%r_eq )

          ! central atom force
          flink(:,at_id(2)) = flink(:,at_id(2)) - dang(:)
        enddo

        !	dihedrals  X--Cq -- Cmm--Y parameters 7 to 16
        !	for each 3 X, for each 3 Y
        do jq = 2, 4
        do jm = 2, 4
          if ( link_atoms(i)%mm(jm) == 0 ) cycle
          at_id(1) = link_atoms(i)%qm(jq)
          at_id(2) = link_atoms(i)%qm(1)
          at_id(3) = natot - n_mm + link_atoms(i)%mm(1)
          at_id(4) = natot - n_mm + link_atoms(i)%mm(jm)

          if ( any(at_id(:) == 0 ) ) cycle
          dr12 = (r_mm(:,at_id(1)) - r_mm(:,at_id(2))) / Ang
          dr32 = (r_mm(:,at_id(3)) - r_mm(:,at_id(2))) / Ang
          dr31 = (r_mm(:,at_id(3)) - r_mm(:,at_id(1))) / Ang
          dr43 = (r_mm(:,at_id(4)) - r_mm(:,at_id(3))) / Ang
          dr42 = (r_mm(:,at_id(4)) - r_mm(:,at_id(2))) / Ang
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr12 )
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr32 )
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr31 )
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr43 )
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr42 )

          ! who many parameters for dihedrals 7 to 16
          do m = 1, 4
            k = link_atoms(i)%reson(resonance_id)%params(7+3*(jq-2)+jm-1,m)

            if ( k == 0 ) cycle
            Elink_for_h_atom = Elink_for_h_atom + &
                               mm_top%dihedrals(k)%energy( dr12, dr32, dr43 )

            ! force for each of the 4 atoms
            do j = 1, 4
              call mm_top%dihedrals(k)%grad( fce, j, dr12, dr32, dr43, &
                                                           dr31, dr42 )

              ! force assignation
              flink(:,at_id(j)) = flink(:,at_id(j)) - fce(:)
            enddo
          enddo ! enddo more than 1 parameter
        enddo ! enddo X and Y
        enddo

        ! dihedrals Cq -- Cmm--Y--Z  parameters 17 to 25
        ! for each 3 X, for each 3 Y
        do jm = 2, 4
          if ( link_atoms(i)%mm(jm) == 0 ) cycle

          do jq = 1, 3
            if ( link_atoms(i)%mm2(jm,jq) == 0) cycle ! is other Z atom?

            at_id(1) = link_atoms(i)%qm(1)
            at_id(2) = natot - n_mm + link_atoms(i)%mm(1)
            at_id(3) = natot - n_mm + link_atoms(i)%mm(jm)
            at_id(4) = natot - n_mm + link_atoms(i)%mm2(jm,jq)

            if ( any(at_id(:) == 0 ) ) cycle
            dr12 = (r_mm(:,at_id(1)) - r_mm(:,at_id(2))) / Ang
            dr32 = (r_mm(:,at_id(3)) - r_mm(:,at_id(2))) / Ang
            dr31 = (r_mm(:,at_id(3)) - r_mm(:,at_id(1))) / Ang
            dr43 = (r_mm(:,at_id(4)) - r_mm(:,at_id(3))) / Ang
            dr42 = (r_mm(:,at_id(4)) - r_mm(:,at_id(2))) / Ang
            call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr12 )
            call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr32 )
            call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr31 )
            call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr43 )
            call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dr42 )

            l = 16 + 3 *( jm - 2 ) + jq
            ! how many parameters for dihedrals 17 to 25
            do m = 1, 4 ! more than 1 paramter
              k = link_atoms(i)%reson(resonance_id)%params(l,m)
              if ( k == 0 ) cycle

              Elink_for_h_atom = Elink_for_h_atom + &
                                 mm_top%dihedrals(k)%energy( dr12, dr32, dr43 )

              do j = 1, 4 ! force over 4 atoms
                call mm_top%dihedrals(k)%grad( fce, j, dr12, dr32, dr43, &
                                               dr31, dr42 )

                ! force assignation
                flink(:,at_id(1)) = flink(:,at_id(1)) - fce(:)
              enddo
            enddo ! enddo more than 1 parameter
          enddo ! enddo to Y and Z
        enddo

        ! For each LA Force over HL scaling and sum where corresponds
        ! HL force division in parallel fpar and perpendicular fpp
        ! unitary versor dversor(3) at_id(3)=LA

        at_id(1) = link_atoms(i)%qm(1)
        at_id(2) = natot - n_mm + link_atoms(i)%mm(1)
        dversor(:) = (r_mm(:,at_id(1)) - &
                   link_atoms(i)%reson(resonance_id)%r(:)) / Ang

        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dversor )

        rhq = norm_v2( dversor )
        dversor(1:3) = dversor(1:3) / rhq

        ip = dversor(1) * link_atoms(i)%reson(resonance_id)%fa(1) + &
             dversor(2) * link_atoms(i)%reson(resonance_id)%fa(2) + &
             dversor(3) * link_atoms(i)%reson(resonance_id)%fa(3)

        if ( h_link_scaling_type == 'splam' ) then
          flink(1:3,at_id(2)) = flink(1:3,at_id(2)) + &
                                link_atoms(i)%reson(resonance_id)%bond_k * ip * dversor(1:3)
          flink(1:3,at_id(1)) = flink(1:3,at_id(1)) + ip * dversor(1:3) - &
                                link_atoms(i)%reson(resonance_id)%bond_k * ip * dversor(1:3)
        else
          flink(1:3,at_id(1)) = flink(1:3,at_id(1)) + ip * dversor(1:3)
        endif
        fpp(1:3) = link_atoms(i)%reson(resonance_id)%fa(1:3) - ip * dversor(1:3)
        dang(:)  = ( r_mm(:,at_id(1)) - r_mm(:,at_id(2)) ) / Ang

        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, dang )
        r_qm = norm_v2( dang )

        flink(:,at_id(2)) = flink(:,at_id(2)) + fpp(:) * rhq / r_qm
        flink(:,at_id(1)) = flink(:,at_id(1)) + fpp(:) * (1.0_dp - rhq / r_qm)

        Ener = Ener + resonance(resonance_id)%weight * Elink_for_h_atom
      enddo ! numlink

      frc(1:3,1:natot) = frc(1:3,1:natot) + &
          resonance(resonance_id)%weight * flink(1:3,1:natot)
    enddo

    call de_alloc( perdihe2, 'perdihe2', 'link2' )
    call de_alloc( flink   , 'flink'   , 'link2' )
  end subroutine link2

  ! calculates HL position
  subroutine link3( na_qm, na_mm, qm_atoms, mm_atoms, mm_top, cell, latt_type )
    use functions  , only : norm_v2
    use mm_topology, only : fftopology, mm_atom_t, qm_atom_t
    use qmmm_pbc   , only : pbc_displ_vector, reccel
    use units      , only : Ang

    integer  , intent(in) :: na_qm
    integer  , intent(in) :: na_mm
    type(qm_atom_t) , intent(in) :: qm_atoms(na_qm)
    type(mm_atom_t) , intent(in) :: mm_atoms(na_mm)
    type(fftopology), intent(in) :: mm_top
    real(dp) , intent(in) :: cell(3,3)
    character, intent(in) :: latt_type

    integer resonance_id, i, at_id(2), k
    real(dp) :: dist, rlink2qm, kcell(3,3), dr(3)

    call reccel( 3, cell, kcell, 0 )

    ! dist HL-CQM
    do resonance_id = 1, num_resonances
    do i = 1, numlink
      at_id(1) = link_atoms(i)%qm(1)
      at_id(2) = link_atoms(i)%mm(1)

      dr(:) = mm_atoms(at_id(2))%r(:) - qm_atoms(at_id(1))%r(:)
      call pbc_displ_vector( latt_type, cell, kcell, dr )
      dist = norm_v2( dr )
      dr(:) = dr(:) / dist

      k = link_atoms(i)%reson(resonance_id)%params(1,1)
      if ( h_link_scaling_type == 'splam' ) then
        rlink2qm = link_atoms(i)%reson(resonance_id)%bond_r * Ang + &
                   link_atoms(i)%reson(resonance_id)%bond_k * &
                   ( dist - mm_top%bonds(k)%r_eq * Ang )
      else
        rlink2qm = link_atoms(i)%reson(resonance_id)%bond_r * Ang
      endif

      link_atoms(i)%reson(resonance_id)%r(:) = &
        qm_atoms(at_id(1))%r(:) + rlink2qm * dr(:)
    enddo
    enddo
  end subroutine link3

  logical function check_out_of_range( atomi, atomj, chi, chj, rmin )
    use sys, only : die

    integer         , intent(in) :: atomi
    integer         , intent(in) :: atomj
    character(len=4), intent(in) :: chi
    character(len=4), intent(in) :: chj
    real(dp)        , intent(in) :: rmin

    real(dp), parameter :: bond_range=1.2_dp
    character :: ch1, ch2

    ! List of typical single bond between species.
    ! Bond lengths taken from:
    ! http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html

    check_out_of_range = .false.
    ch1 = chi(1:1)
    ch2 = chj(1:1)

    if ( (ch1 == 'C') .and. (ch2 == 'C') ) then
      if ( rmin > 1.54_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'C') .and. (ch2 == 'H')) .or. &
             ((ch1 == 'H') .and. (ch2 == 'C')) ) then
      if ( rmin > 1.09_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'C') .and. (ch2 == 'O')) .or. &
             ((ch1 == 'O') .and. (ch2 == 'C')) ) then
      if ( rmin > 1.43_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'C') .and. (ch2 == 'N')) .or.&
             ((ch1 == 'N') .and. (ch2 == 'C')) ) then
      if ( rmin > 1.47_dp * bond_range ) check_out_of_range = .true.
    elseif ( (ch1 == 'N') .and. (ch2 == 'N') ) then
      if ( rmin > 1.45_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'N') .and. (ch2 == 'H')) .or. &
             ((ch1 == 'H') .and. (ch2 == 'N')) ) then
      if ( rmin > 1.01_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'N') .and. (ch2 == 'O')) .or. &
             ((ch1 == 'O') .and. (ch2 == 'N')) ) then
      if ( rmin > 1.47_dp * bond_range ) check_out_of_range = .true.
    elseif ( (ch1 == 'H') .and. (ch2 == 'H') ) then
      if ( rmin > 0.74_dp * bond_range ) check_out_of_range = .true.
    elseif ( (ch1 == 'O') .and. (ch2 == 'O') ) then
      if ( rmin > 1.48_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'C') .and. (ch2 == 'S')) .or. &
             ((ch1 == 'S') .and. (ch2 == 'C')) ) then
      if ( rmin > 1.82_dp * bond_range ) check_out_of_range = .true.
    elseif ( ((ch1 == 'C') .and. (ch2 == 'F')) .or. &
             ((ch1 == 'F') .and. (ch2 == 'C')) ) then
      if ( rmin > 1.35_dp * bond_range ) check_out_of_range = .true.
    endif

    if ( check_out_of_range ) return

    write(*,*) "No suitable parameter for Cqm-LA bond. Please set a value "//&
               "in the bonds block of amber.parm. This value is needed to "//&
               "decide if the atom are close enough to make a bond."
    write(*,*) 'Pair distance(A): ', rmin
    write(*,*) 'QM atom: ', atomi
    write(*,*) 'MM atom: ', atomj
    write(*,*) 'Parameters for bond equilibrium distance: ', chi, '-', chj
    call die( 'No suitable parameters found for Cqm-LA bond in amber.parm.' )
  end function check_out_of_range

  subroutine set_default_H_link_bond_parameters( qmtype, H_bond_length, &
                                                 MM_kbond, H_kbond )
    use sys, only: die
    implicit none
    character(len=4), intent(in)  :: qmtype
    real(dp)        , intent(in)  :: MM_kbond
    real(dp)        , intent(out) :: H_bond_length
    real(dp)        , intent(out) :: H_kbond

    if ( qmtype(1:1) == 'C' ) then
      H_bond_length = 1.09_dp
      H_kbond = sqrt( MM_kbond / 413.0_dp )
    elseif ( qmtype(1:1) == 'N' ) then
      H_bond_length = 1.01_dp
      H_kbond = sqrt( MM_kbond / 391.0_dp )
    elseif ( qmtype(1:1) == 'O' ) then
      H_bond_length = 0.96_dp
      H_kbond = sqrt( MM_kbond / 366.0_dp )
    else
      call die('No default parameters fof H link bond')
    endif
  end subroutine set_default_H_link_bond_parameters
end module linkatoms
