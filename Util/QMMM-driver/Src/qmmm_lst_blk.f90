module qmmm_list_block
  !! This module contains the declaration of both blocked (i.e non-moving or
  !! constrained ) MM atoms and those atoms that are close enough to participate
  !! in QM-MM interactions (only used afterwards when writing the PDB).
  !!
  !! Basically we first attempt to read a ".lst" file with closest MM atoms. If
  !! that fails, we consider all atoms within rcut_qmm as being "QM-MM active".
  !! Then, we do the same, trying to read a ".blk" file; MM atoms specified in
  !! the file will be considered "blocked" (immovable). If no such file exists,
  !! we consider any atom BEYOND blockmm_radius to be frozen.
  !!
  !! If rcut_qmmm is greater than 99A, all MM atoms are considered to be
  !! "QM-MM active". In the same manner, if blockmm_radius is greater than
  !! 99A, no MM atoms are considered to be blocked.
  implicit none
  public :: qmmm_lst_blk

  private

contains
  subroutine qmmm_lst_blk( na_qm, na_mm, nroaa, atxres, rclas, rcut_qmmm, &
                           blockmm_radius, mm_atoms, slabel )

    !! Calculates the rcut and block QM-MM only in the first step. Attempts
    !! to first read an .lst or .blk file in order to avoid calculations.
    use alloc      , only : re_alloc, de_alloc
    use functions  , only : dist2_v2
    use mm_topology, only : mm_atom_t
    use precision  , only : dp
    use units      , only : Ang

    implicit none
    integer         , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)    :: na_mm
      !! Number of MM atoms.
    integer         , intent(in)    :: nroaa
      !! Number of MM residues.
    integer         , intent(in)    :: atxres(nroaa)
      !! Number of atoms per residue.
    real(dp)        , intent(in)    :: rclas(3,na_qm+na_mm)
      !! Atomic positions.
    real(dp)        , intent(in)    :: rcut_qmmm
      !! Cut-off radius for QM-MM interactions.
    real(dp)        , intent(in)    :: blockmm_radius
      !! Cut-off radius for blocking (constraining) MM atoms.
    character(len=*), intent(in)    :: slabel
      !! System label (name).
    type(mm_atom_t) , intent(inout) :: mm_atoms(na_mm)
      !! Information for MM atoms. Will set up "is_blocked" and
      !! "is_qm_neighbour".

    ! Internal variables.
    character(len=255) :: fname
    integer            :: ios, iat, jat, qat, ires, iu, count, tmp_blk
    logical            :: bloqmmm, liqmmm, found
    real(dp)           :: dist, dist2
    real(dp), pointer  :: cm(:,:)

    ! External subroutines and functions.
    external           :: io_assign, io_close

    ios = 0

    ! Calculate center of masses of all residues.
    nullify( cm )
    call re_alloc( cm, 1, 3, 1, nroaa, 'cm', 'qmmm_lst_blk' )

    iat = na_qm +1
    do ires = 1, nroaa
      do jat = 1, atxres(ires)
        cm(1:3,ires) = cm(1:3,ires) + rclas(1:3,iat)
        iat = iat +1
      enddo

      cm(1:3,ires) = cm(1:3,ires) / atxres(ires)
    enddo

    ! rcut QM-MM
    found  = .false.
    liqmmm = .false.

    ! Find file name and verify existence.
    fname  = trim(slabel) // '.lst'
    inquire( file = fname, exist = found )

    if ( found ) then
      call io_assign( iu )
      open( iu, file = fname, status = 'old' )

      ! Read QMMM list.
      do iat = 1, na_mm
        read( iu, *, iostat = ios) tmp_blk
        if ( check_ios( ios, '.lst' ) ) then
          call io_close( iu )
          call de_alloc( cm, 'cm', 'qmmm_lst_blk' )

          do jat = 1, iat
            mm_atoms(jat)%is_qm_neighbour = .false.
          enddo
          return
        endif

        if (tmp_blk > 0) mm_atoms(iat)%is_qm_neighbour = .true.
      enddo

      call io_close( iu )
      write( 6, '(/a)' ) 'qm-mm: Reading QM-MM neighbours list from file.'
    else
      if ( .not. (rcut_qmmm > 99.0_dp) ) liqmmm = .true.
      if ( liqmmm ) then
        ! QM-MM neigh list
        write( 6, '(/a,f12.6)' ) 'qm-mm: cut off radius QM-MM (Ang):', rcut_qmmm

        dist     = 0.0_dp
        dist2    = rcut_qmmm * rcut_qmmm * Ang * Ang
        do iat = 1, na_mm
          mm_atoms(iat)%is_qm_neighbour = .true.
        enddo

        do qat = 1, na_qm
          iat = 1

          do ires = 1, nroaa
            dist = dist2_v2( rclas(:,qat), cm(:,ires) )
            do jat = 1, atxres(ires)
              if ( dist2 > dist ) mm_atoms(iat)%is_qm_neighbour = .false.
              iat = iat +1
            enddo
          enddo
        enddo

        count = 0
        do iat = 1, na_mm
          if ( .not. (mm_atoms(iat)%is_qm_neighbour)) count = count +1
        enddo

        write( 6, '(/a,2x,i5)' ) 'qm-mm: Number of QM-MM interacting atoms:', &
                                 count

        ! Open file and write listqmmm.
        call io_assign( iu )
        open( iu, file = fname, form = 'formatted', status = 'unknown' )

        do iat = 1, na_mm
          write(iu,*) mm_atoms(iat)%is_qm_neighbour
        enddo

        call io_close( iu )

      else ! liqmmm = .false.
        write(6,'(/a)') 'qm-mm: Warning QM-MM cutoff too large.'
      endif
    endif

    ! block QM-MM
    found   = .false.
    bloqmmm = .false.

    ! Find file name and check if input exists.
    fname = trim(slabel) // '.blk'
    inquire( file = fname, exist = found )
    if ( found ) then
      call io_assign( iu )
      open( iu, file = fname, status = 'old' )

      ! read blockqmmm
      do iat = 1, na_mm
        read( iu, *, iostat = ios ) tmp_blk
        if ( tmp_blk > 0 ) mm_atoms(iat)%is_blocked = .true.

        if ( check_ios( ios, '.blk' ) ) then
          call io_close( iu )
          call de_alloc( cm, 'cm', 'qmmm_lst_blk' )

          do jat = 1, iat
            mm_atoms(jat)%is_blocked = .false.
          enddo
          return
        endif
      enddo

      call io_close( iu )
      write( 6, '(/a)' ) 'qm-mm: Reading blocked QM-MM atoms from file.'
    else
      if ( .not. (blockmm_radius > 99.0_dp) ) bloqmmm = .true.
      if( bloqmmm ) then
        ! Fixing MM atoms beyond block cut off
        write( 6, '(/a,f12.6)' ) 'qm-mm: cut off radius Block (Ang):', &
                                 blockmm_radius
        dist         = 0.0_dp
        dist2        = blockmm_radius * blockmm_radius * Ang * Ang

        do iat = 1, na_mm
          mm_atoms(iat)%is_blocked = .true.
        enddo

        do qat = 1, na_qm
          iat = 1

          do ires = 1, nroaa
            dist = dist2_v2( rclas(:,qat), cm(:,ires) )

            do jat = 1, atxres(ires)
              if ( dist2 > dist ) mm_atoms(iat)%is_blocked = .false.
              iat = iat +1
            enddo
          enddo
        enddo

        count = 0
        do iat = 1, na_mm
          if ( .not. (abs(blockmm_radius) > 0.0_dp) ) &
            mm_atoms(iat)%is_blocked = .true.
          if ( .not. mm_atoms(iat)%is_blocked ) count = count +1
        enddo

        write( 6, '(/a,2x,i5)' ) 'qm-mm: Number of QM-MM free atoms:', count

        ! Write blockQMMM to file.
        call io_assign( iu )
        open( iu, file = fname, form = 'formatted', status = 'unknown' )

        do iat = 1, na_mm
          write( iu, * ) mm_atoms(iat)%is_blocked
        enddo

        call io_close( iu )

      else ! bloqmmm=.false.
        write( 6, '(/a)' ) 'qm-mm: Warning Block cutoff too large.'
      endif
    endif
    call de_alloc( cm, 'cm', 'qmmm_lst_blk' )

    call pxfflush( 6 )
  end subroutine qmmm_lst_blk

  function check_ios( ios, f_ext ) result( exit_routine )
    !! Checks read status when attempting to read files.
    use sys, only : message
    implicit none
    integer         , intent(in) :: ios
      !! Read status.
    character(len=4), intent(in) :: f_ext
      !! File extention involved in read.
    logical :: exit_routine

    exit_routine = .false.
    if ( ios /= 0 ) then
      call message( 'WARNING',&
                    'qm-mm: Problem reading from '//f_ext//' file.' )
      exit_routine = .true.
    endif
  end function check_ios

end module qmmm_list_block
