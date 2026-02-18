module reinsert_atoms_m
  !! Subroutines used to deal with the PBC
  implicit none
  public :: reinserting_atoms_in_box

contains

  subroutine reinserting_atoms_in_box( latt_type, natot, na_qm, nac, qm_atoms, &
                                       mm_atoms, is_atom_blocked, cell )
    use mm_topology, only : mm_atom_t, qm_atom_t
    use qmmm_pbc   , only : reccel
    use precision  , only : dp

    implicit none
    character, intent(in)    :: latt_type
    integer  , intent(in)    :: natot
    integer  , intent(in)    :: na_qm
    integer  , intent(in)    :: nac
    integer  , intent(in)    :: is_atom_blocked(natot)
    real(dp) , intent(in)    :: cell(3,3)

    type(qm_atom_t), intent(inout) :: qm_atoms(na_qm)
    type(mm_atom_t), intent(inout) :: mm_atoms(nac)

    integer  :: iCrd, iat, jWat, kCrd, nk, shift
    real(dp) :: rnk, kcell(3,3)

    iat = 0
    shift = 0

    if ( latt_type == 'D' ) then
      do
        iat = iat + shift +1
        shift = 0

        if ( iat > natot ) exit
        if ( is_atom_blocked(iat) == 1) cycle


        if ( iat <= na_qm ) then
          do iCrd = 1, 3
            rnk = qm_atoms(iat)%r(iCrd) / cell(iCrd,iCrd)

            nk = INT( rnk )
            if ( rnk < 0.0_dp ) nk = INT(rnk-1.0_dp)

            qm_atoms(iat)%r(iCrd) = qm_atoms(iat)%r(iCrd) - nk*cell(iCrd,iCrd)
          enddo

        else
          if ( mm_atoms(iat-na_qm)%aaname == 'HOH' ) shift = 2

          do iCrd = 1, 3
            rnk = mm_atoms(iat-na_qm)%r(iCrd) / cell(iCrd,iCrd)

            nk = INT( rnk )
            if ( rnk < 0.0_dp ) nk = INT(rnk-1.0_dp)

            do jWat = 0, shift
              mm_atoms(iat+jWat-na_qm)%r(iCrd) = &
              mm_atoms(iat+jWat-na_qm)%r(iCrd) - nk * cell(iCrd,iCrd)
            enddo
          enddo
        endif
      enddo ! atoms

    else  ! lattice type
      call reccel( 3, cell, kcell, 0 )

      do
        iat = iat + shift +1
        if ( iat > natot ) exit
        shift = 0

        if ( iat <= na_qm ) then
          do kCrd = 1, 3
            rnk = 0.0_dp

            ! Here it is assumed that kcell is the matrix, whose rows are
            ! the reciprocal vectors without 2*pi factor.
            do iCrd = 1, 3
              rnk = rnk + qm_atoms(iat)%r(iCrd) * kcell(iCrd,kCrd)
            enddo

            nk = INT(rnk)
            if ( rnk < 0.0_dp ) nk = INT(rnk-1.0_dp)

            if ( nk == 0 ) cycle
            do iCrd = 1, 3
              qm_atoms(iat)%r(iCrd) = qm_atoms(iat)%r(iCrd) - cell(iCrd,kCrd)*nk
            enddo
          enddo

        else
          if ( mm_atoms(iat-na_qm)%aaname == 'HOH' ) shift = 2

          do kCrd = 1, 3
            rnk = 0.0_dp

            ! Here it is assumed that kcell is the matrix, whose rows are
            ! the reciprocal vectors without 2*pi factor.
            do iCrd = 1, 3
              rnk = rnk + mm_atoms(iat-na_qm)%r(iCrd) * kcell(iCrd,kCrd)
            enddo

            nk = INT(rnk)
            if ( rnk < 0.0_dp ) nk = INT(rnk-1.0_dp)

            if ( nk == 0 ) cycle
            do jWat = 0, shift
            do iCrd = 1, 3
              mm_atoms(iat+jWat-na_qm)%r(iCrd) = &
                mm_atoms(iat+jWat-na_qm)%r(iCrd) - cell(iCrd,kCrd) * nk
            enddo
            enddo
          enddo
        endif
      enddo ! atoms
    endif ! lattice type
  end subroutine reinserting_atoms_in_box

end module reinsert_atoms_m