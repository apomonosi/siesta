module center_all
  !! This module contains routines to center the QM system in the cell.
  use mm_topology, only : mm_atom_t, qm_atom_t

  implicit none
  public :: centermol
  public :: centerdyn
contains

  subroutine centermol( na_qm, na_mm, qm_atoms, mm_atoms, ucell )
    !! Centers quantum subsystem inside the cell.
    use precision, only: dp

    implicit none
    integer , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer , intent(in)    :: na_mm
      !! Number of MM atoms.
    real(dp), intent(in)    :: ucell(3,3)
      !! Cell vectors.
    type(qm_atom_t), intent(inout) :: qm_atoms(na_qm)
      !! QM atoms.
    type(mm_atom_t), intent(inout) :: mm_atoms(na_mm)
      !! MM atoms.

    integer  :: iAtom, iCrd
    real(dp) :: d0(3), min_r(3), max_r(3), delta(3)

    d0(:)      = 0.0_dp
    min_r(1:3) = qm_atoms(1)%r(1:3)
    max_r(1:3) = qm_atoms(1)%r(1:3)

    do iAtom = 1, na_qm
    do iCrd  = 1, 3
      if ( qm_atoms(iAtom)%r(iCrd) > max_r(iCrd) ) &
        max_r(iCrd) = qm_atoms(iAtom)%r(iCrd)
      if ( qm_atoms(iAtom)%r(iCrd) < min_r(iCrd) ) &
        min_r(iCrd) = qm_atoms(iAtom)%r(iCrd)
    enddo
    enddo

    do iCrd = 1, 3
      delta(iCrd) = max_r(iCrd) - min_r(iCrd)
      d0(iCrd)    = ( ucell(iCrd,iCrd) - delta(iCrd) ) * 0.5_dp
    enddo

    do iAtom = 1, na_qm
      qm_atoms(iAtom)%r(1:3) = qm_atoms(iAtom)%r(1:3) - min_r(1:3) + d0(1:3)
    enddo
    do iAtom = 1, na_mm
      mm_atoms(iAtom)%r(1:3) = mm_atoms(iAtom)%r(1:3) - min_r(1:3) + d0(1:3)
    enddo

    write(6,'(/a)') 'qmmm centermol: System centered in cell.'

  end subroutine centermol

  subroutine centerdyn( na_qm, na_mm, qm_atoms, mm_atoms, ucell )
    !! Subroutine called in each dynamic step to recenter the quantum subsystem
    !! inside the cell when the dynamics drives it near a border.
    use precision, only: dp

    implicit none
    integer , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer , intent(in)    :: na_mm
      !! Number of MM atoms.
    real(dp), intent(in)    :: ucell(3,3)
      !! Cell vectors.
    type(qm_atom_t), intent(inout) :: qm_atoms(na_qm)
      !! QM atoms.
    type(mm_atom_t), intent(inout) :: mm_atoms(na_mm)
      !! MM atoms.

    integer  :: iAtom, iCrd
    real(dp) :: d0(3), min_r(3), max_r(3), delta(3), lbord, dbord(3)
    logical  :: ctr

    d0(:)    = 0.0_dp
    dbord(:) = 0.0_dp
    lbord    = 2.0_dp
    ctr      = .false.

    min_r(1:3) = qm_atoms(1)%r(1:3)
    max_r(1:3) = qm_atoms(1)%r(1:3)

    do iAtom = 1, na_qm
    do iCrd  = 1, 3
      if ( qm_atoms(iAtom)%r(iCrd) > max_r(iCrd) ) &
        max_r(iCrd) = qm_atoms(iAtom)%r(iCrd)
      if ( qm_atoms(iAtom)%r(iCrd) < min_r(iCrd) ) &
        min_r(iCrd) = qm_atoms(iAtom)%r(iCrd)
    enddo
    enddo

    do iCrd = 1, 3
      dbord(iCrd) = ucell(iCrd,iCrd) - max_r(iCrd)
      if ( (min_r(iCrd) < lbord) .or. (dbord(iCrd) < lbord) ) ctr=.true.
    enddo

    if ( .not. ctr ) return

    write(6,'(/a)') 'centerdyn: System centered in cell.'
    do iCrd = 1, 3

      delta(iCrd) = max_r(iCrd) - min_r(iCrd)
      d0(iCrd)    = ( ucell(iCrd,iCrd) - delta(iCrd) ) * 0.5_dp

      if ( d0(iCrd) < lbord ) &
        write(6,'(/a)') 'centerdyn: QM atoms too close to cell border.'
    enddo

    do iAtom = 1, na_qm
      qm_atoms(iAtom)%r(1:3) = qm_atoms(iAtom)%r(1:3) - min_r(1:3) + d0(1:3)
    enddo
    do iAtom = 1, na_mm
      mm_atoms(iAtom)%r(1:3) = mm_atoms(iAtom)%r(1:3) - min_r(1:3) + d0(1:3)
    enddo
  end subroutine centerdyn
end module center_all