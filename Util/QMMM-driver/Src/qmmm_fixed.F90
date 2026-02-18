module qmmm_constraints
  !! This module contains routines to apply and read atomic
  !! constraints.
  !! The first routine, fixed1, reads the input file and builds a list of
  !! "blocked" atoms, which are then used when fixed2 removes their
  !! forces and velocities.

  ! Input explanation:
  !   You can have up to 4 different constraint groups at the same time,
  ! chosen from the following menu:
  !   ( 1) Constraint all ATOMS selected (by atom index).
  !   (-1) Constraint all ATOMS, except those selected.
  !   ( 2) Constraint all RESIDUES selected (by residue index).
  !   (-2) Constraint all RESIDUES, except those selected.
  !   ( 3) Constraint all MM non-H atoms.
  !   ( 4) Constraint all MM alpha-C, C and N atoms.
  !   ( 5) Constraint all non-water residues.
  !   ( 6) Constraint all water molecules outside a given distance from the
  !        center of mass of a given residue.
  !   ( 7) Adds a restraint to water molecules that are far away from the center
  !        of mass. Presumably, this prevents the solvent box from exploding...
  !
  ! The input should be basically as follows:
  !  * First a character and number of constraint groups available. Then for
  !    each group:
  !     * A character followed by the type of each constraint group.
  !     * If type is <= 2, the masks indicating the selection.
  !
  ! These masks should:
  !    1) start by F/f/P/p
  !    2) a character (usually f, as in "From").
  !    3) the starting index of the constraint.
  !    4) a character (usually t, as in "To")
  !    5) the final index of the constraint.
  ! This might be space-sensitive, so I added "|" characters to
  ! indicate line start in a file.
  !
  ! |%block PositionConstraints
  ! | n 5
  ! | t 1
  ! |F f 5  t 6
  ! |F f 10 t 17
  ! |F f 23 t 25
  ! | t 2
  ! |F f 84 t 87
  ! |F f 90 t 91
  ! | t 4
  ! | t 6
  ! |  115
  ! |  13
  ! | t 7
  ! |%end block
  !
  ! In this case we are constraining:
  !  * atoms 5-6, 10-17, 23-25
  !  * residues 84-87, 90-91
  !  * all MM alpha-C, C and N atoms
  !  * all water molecules more than 13 Ang away from residue number 115
  !  * all water molecules that are within 2.5 Ang from the maximum distance
  !  * towards de center of mass.
  ! (Type -1 works as type 1, and type -2 works as type 2. Also, types 3, 4, 5
  !  and 7 do not need any additional input)
  !
  ! For all constraint groups of the type < 2, you can have as many
  ! atom selection ranges as you want, so long as they start with the
  ! f, F, p or P letters.

  implicit none
  public :: qmmm_fixed1
  public :: qmmm_fixed2
  public :: add_siesta_qmmm_ntcon

  private
  logical :: first_call = .true.
contains

  subroutine qmmm_fixed1( na_qm, na_mm, natot, nroaa, blocklist, mm_atoms, wat )
    !! Reads the constrained atom block.
    use alloc      , only : re_alloc, de_alloc
    use fdf        , only : block_fdf, parsed_line
    use fdf        , only : fdf_block, fdf_bline, fdf_bintegers, fdf_brewind, &
                            fdf_bvalues, fdf_bnames, fdf_bclose, fdf_bbackspace
    use functions  , only : dist2_v2
    use mm_topology, only : mm_atom_t
    use precision  , only : dp
    use sys        , only : die
    use units      , only : Ang

    implicit none
    integer         , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)  :: natot
      !! Total (QM+MM) mumber of atoms.
    integer         , intent(in)  :: nroaa
      !! Total number of MM residues.
    type(mm_atom_t) , intent(in)  :: mm_atoms(na_mm)
      !! Index of the residue containing a given atom.
    logical         , intent(out) :: wat
      !! Constrain water cap.
    integer         , intent(out) :: blocklist(natot)
      !! The list of "blocked" (i.e. constrained) atoms.

    integer          :: itype, icont, nt, aa, iat, ires, ncon(4), ctype(4)
    logical          :: reuse_line, last_line
    real(dp)         :: cut, cqm(3), mdist, dist, dist2
    character(len=1) :: ch1, exp
    character(len=4) :: ch4

    integer , pointer :: con1(:,:), con2(:,:)
    real(dp), pointer :: r_bak(:,:)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    ! These first are just to avoid warnings.
    cut = 0.0_dp; aa = 0


    wat          = .false.
    blocklist(:) = 0

    ! Reads the 'PositionConstraints' block. Returns if not found.
    if ( .not. fdf_block('PositionConstraints',bfdf) ) return
    write(6,'(/,a)') 'fixed: reading "PositionConstraints" block'

    if ( .not. fdf_bline(bfdf,pline) ) then
      write(6,'(A)') "fixed: No constraints found."
      call fdf_bclose( bfdf )
      return
    endif

    nt = fdf_bintegers(pline, 1)

    nullify( con1, con2 )
    call re_alloc( con1, 1, 4, 1, natot, 'con1', 'qmmm_fixed1' )
    call re_alloc( con2, 1, 4, 1, natot, 'con2', 'qmmm_fixed1' )

    reuse_line = .false.
    last_line  = .false.
    do itype = 1, nt
      if ( last_line ) exit
      if ( .not. reuse_line ) then
        if ( .not. fdf_bline(bfdf,pline) ) &
          call die('fixed: Error reading constraint type.')
      else
        reuse_line = .false.
      endif
      ctype(itype) = fdf_bintegers(pline, 1)

      if ( (na_mm == 0) .and. (ctype(itype) > 1) ) &
        call die('fixed: constrained type for only QM atoms must be only 1.')

      if ( ctype(itype) <= 2 ) then
        icont = 1

        do while (.true.)
          if ( icont > natot ) &
            call die('fixed: sets of constraints must not exceed natot.')

          if ( .not. fdf_bline(bfdf,pline) ) then
            last_line = .true.
            exit
          endif

          exp = trim(fdf_bnames(pline,1))

          if ( (exp == 'f') .or. (exp == 'F') .or. (exp == 'p') .or. &
                (exp == 'P') ) then
            con1(itype,icont) = fdf_bintegers(pline,1)
            con2(itype,icont) = fdf_bintegers(pline,2)
            icont = icont +1
          else
            reuse_line = .true.
            exit
          endif
        enddo
        ncon(itype) = icont -1

        if ( ncon(itype) == 0 ) &
          call die('fixed: sets of constraints must not be zero')

        do icont = 1, ncon(itype)
          if ( con2(itype,icont) < con1(itype,icont) ) &
            call die( 'fixed: sets of constraints must be defined in the '//&
                      ' correct order.')
        enddo
      elseif ( ctype(itype) == 6 ) then
        if ( .not. fdf_bline(bfdf,pline) ) &
          call die('fixed: Error reading constraint type 6.')
        aa = fdf_bintegers(pline,1)
        if ( .not. fdf_bline(bfdf,pline) ) &
          call die('fixed: Error reading constraint type 6.')
        cut = fdf_bvalues(pline,1)

      elseif ( (ctype(itype) < -2) .or. (ctype(itype) > 7) ) then
        call die('fixed: Wrong type of constraint')

      endif !type
    enddo   !nt
    call fdf_bclose( bfdf )

    ! Blocklist assignation
    do itype = 1, nt
      select case ( ctype(itype) )
      case (1)
        do icont = 1, ncon(itype)
          if ( con2(itype,icont) > natot ) &
            call die('fixed: atoms in constraint must not exeed natot.')

          do iat = con1(itype,icont), con2(itype,icont)
            blocklist(iat) = 1
          enddo
        enddo

      case (-1)
        blocklist(:) = 1
        do icont = 1, ncon(itype)
          if ( con2(itype,icont) > natot ) &
            call die('fixed: atoms in constraint must not exeed natot.')

          do iat = con1(itype,icont), con2(itype,icont)
            blocklist(iat) = 0
          enddo
        enddo

      case (2)
        do icont = 1, ncon(itype)
          if ( con2(itype,icont) > nroaa ) &
            call die('fixed: residues in constraint must not exeed nroaa.')

          do ires = con1(itype,icont), con2(itype,icont)
          do iat  = 1, na_mm
            if ( ires == mm_atoms(iat)%aanum ) blocklist(na_qm+iat) = 1
          enddo
          enddo
        enddo

      case (-2)
        blocklist(:) = 1
        do icont = 1, ncon(itype)
          if ( con2(itype,icont) > nroaa ) &
            call die('fixed: residues in constraint must not exeed nroaa.')

          do ires = con1(itype,icont), con2(itype,icont)
          do iat  = 1, na_mm
            if ( ires == mm_atoms(iat)%aanum ) blocklist(na_qm+iat) = 0
          enddo
          enddo
        enddo

      case (3)
        do iat = 1, na_mm
          ch4 = mm_atoms(iat)%atname
          ch1 = ch4(1:1)
          if ( ch1 /= 'H') blocklist(na_qm+iat) = 1
        enddo

      case (4)
        do iat = 1, na_mm
          if ( (mm_atoms(iat)%atname == 'CA') .or. &
               (mm_atoms(iat)%atname == 'C') .or. &
               (mm_atoms(iat)%atname == 'N') ) blocklist(na_qm+iat) = 1
        enddo

      case (5)
        do iat = 1, na_mm
          if ( mm_atoms(iat)%aaname /= 'HOH' ) blocklist(na_qm+iat) = 1
        enddo

      case (6)
        nullify( r_bak )
        call re_alloc( r_bak, 1, 3, 1, na_mm, 'r_bak', 'qmmm_fixed1' )

        do iat = 1, na_mm
          r_bak(1:3,iat) = mm_atoms(iat)%r(1:3) / Ang
        enddo
        cqm   = 0.0_dp
        icont = 0

        do iat = 1, na_mm
          if ( mm_atoms(iat)%aanum /= aa ) cycle
          icont  = icont  +1
          cqm(1:3) = cqm(1:3) + r_bak(1:3,iat)
        enddo
        cqm(1:3) = cqm(1:3) / icont

        mdist = 0.0_dp
        dist  = 0.0_dp
        do iat = 1, na_mm
          if ( mm_atoms(iat)%aanum /= aa ) cycle

          dist = dist2_v2( r_bak(:,iat), cqm(:) )
          if ( dist > mdist ) mdist = dist
        enddo
        mdist = sqrt(mdist)
        cut   = cut + mdist

        dist  = 0.0_dp
        dist2 = cut * cut
        do iat = 1, na_mm
          if ( (mm_atoms(iat)%aaname == 'HOH') .and. &
               (mm_atoms(iat)%atname == 'O') ) then
            dist = dist2_v2( r_bak(:,iat), cqm(:) )

            if ( dist > dist2 ) then
              blocklist(iat+na_qm)   = 1
              blocklist(iat+na_qm+1) = 1
              blocklist(iat+na_qm+2) = 1
            endif
          elseif ( mm_atoms(iat)%aaname /= 'HOH' ) then
            dist = dist2_v2( r_bak(:,iat), cqm(:) )

            if ( dist > dist2 ) blocklist(iat+na_qm) = 1
          endif
        enddo
        call de_alloc( r_bak, 'r_bak', 'qmmm_fixed1' )

      case (7)
        wat = .true.
        write(6,'(/,a)') 'fixed: Running restraining water cap'

      end select ! type
    enddo !nt

    call de_alloc( con1, 'con1', 'qmmm_fixed1' )
    call de_alloc( con2, 'con2', 'qmmm_fixed1' )
  end subroutine qmmm_fixed1

  subroutine  qmmm_fixed2( na_qm, na_mm, natot, nfree, blocklist, &
                           mm_atoms, forces, cforces, vat )
    !! Imposes force and velocity constraints.
    use mm_topology, only : mm_atom_t
    use precision  , only : dp

    implicit none
    integer , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer , intent(in)    :: na_mm
      !! Number of MM atoms.
    integer , intent(in)    :: natot
      !! Total (QM+MM) mumber of atoms.
    integer , intent(in)    :: blocklist(natot)
      !! The list of "blocked" (i.e. constrained) atoms.
    type(mm_atom_t), intent(in) :: mm_atoms(na_mm)
      !! The list of "blocked" (i.e. constrained) MM atoms.
    real(dp), intent(in)    :: forces(3,natot)
      !! Unconstrained atomic forces.
    real(dp), intent(out)   :: cforces(3,natot)
      !! Constrained atomic forces.
    real(dp), intent(out)   :: vat(3,natot)
      !! Atomic velocities.
    integer , intent(inout) :: nfree
      !! Number of free atoms.

    integer :: iat
    cforces(:,:) = forces(:,:)

    ! Nullify forces and velocities.
    do iat = 1, na_qm
      if ( blocklist(iat) == 1 ) then
        cforces(1:3,iat) = 0.0_dp
        vat(1:3,iat)     = 0.0_dp
      endif
    enddo
    do iat = 1, na_mm
      if ( (mm_atoms(iat)%is_blocked) .or. (blocklist(iat+na_qm) == 1) ) then
        cforces(1:3,iat+na_qm) = 0.0_dp
        vat(1:3,iat+na_qm)     = 0.0_dp
      endif
    enddo

    ! Set total number of free atoms
    if ( first_call ) then
      nfree = 0
      do iat = 1, na_qm
        if ( blocklist(iat) == 0 ) nfree = nfree +1
      enddo
      do iat = 1, na_mm
        if ( (blocklist(iat+na_qm) == 0) .and. &
             (.not. mm_atoms(iat)%is_blocked) ) &
          nfree = nfree +1
      enddo

      write(6,'(/a,2x,i5)') 'siesta-qmmm: Total Free (Unconstrained) Atoms:', &
                            nfree
      first_call = .false.
    endif
  end subroutine qmmm_fixed2

  subroutine add_siesta_qmmm_ntcon( cell, na_qm, na_mm, isa, amass, xa, &
                                    blocklist, mm_atoms, ntcon )
    ! Returns the total amount of atoms constrained.
    use alloc      , only : de_alloc, re_alloc
    use m_fixed    , only : fixed
    use mm_topology, only : mm_atom_t
    use precision  , only : dp

    implicit none
    real(dp), intent(in)  :: cell(3,3)
      !! The unit cell vectors.
    integer , intent(in)  :: na_qm
      !! Number of QM atoms.
    integer , intent(in)  :: na_mm
      !! Number of MM atoms.
    integer , intent(in)  :: isa(na_qm+na_mm)
      !! Species index belonging to each atom.
    real(dp), intent(in)  :: amass(na_qm+na_mm)
      !! Atomic masses.
    real(dp), intent(in)  :: xa(3,na_qm+na_mm)
      !! Atomic coordinates.
    integer , intent(in)  :: blocklist(na_qm+na_mm)
      !! The list of "blocked" (i.e. constrained) atoms.
    type(mm_atom_t), intent(in) :: mm_atoms(na_mm)
      !! The list of "blocked" (i.e. constrained) MM atoms.
    integer , intent(out) :: ntcon
      !! Total number of atoms constrained.

    integer           :: iat
    ! These are all dummy variables for "fixed()"
    real(dp)          :: stress(3,3), cstress(3,3)
    real(dp), pointer :: fa(:,:), cfa(:,:)

    nullify( fa, cfa )
    call re_alloc( fa , 1, 3, 1, na_qm+na_mm, 'fa' , 'add_siesta_qmmm_ntcon' )
    call re_alloc( cfa, 1, 3, 1, na_qm+na_mm, 'cfa', 'add_siesta_qmmm_ntcon' )

    stress(:,:) = 0.0_dp
    fa(:,:)     = 1.0_dp

    call fixed( cell, stress, na_qm+na_mm, isa, amass, xa, fa, &
                cstress, cfa, ntcon )
    call de_alloc( fa, 'fa', 'add_siesta_qmmm_ntcon' )

    ! set total number of free atoms
    do iat = 1, na_qm
      if ( blocklist(iat) == 1 ) then
        if ( cfa(1,iat) > 0.0_dp ) ntcon = ntcon +1
        if ( cfa(2,iat) > 0.0_dp ) ntcon = ntcon +1
        if ( cfa(3,iat) > 0.0_dp ) ntcon = ntcon +1
      endif
    enddo
    do iat = 1, na_mm
      if ( (blocklist(iat+na_qm) == 1 ) .or. (mm_atoms(iat)%is_blocked) ) then
        if ( cfa(1,iat) > 0.0_dp ) ntcon = ntcon +1
        if ( cfa(2,iat) > 0.0_dp ) ntcon = ntcon +1
        if ( cfa(3,iat) > 0.0_dp ) ntcon = ntcon +1
      endif
    enddo

    call de_alloc( cfa, 'cfa', 'add_siesta_qmmm_ntcon' )
  end subroutine add_siesta_qmmm_ntcon
end module qmmm_constraints