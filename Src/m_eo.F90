! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_eo

  use precision, only: dp
  implicit none

  public

  real(dp), pointer, save :: eo(:,:,:)  ! Hamiltonian eigenvalues
  real(dp), pointer, save :: qo(:,:,:)  ! Occupations of eigenstates

  logical, save :: scf_eigenvalues_available

contains

  subroutine eoqo_reset()
    use alloc, only : de_alloc
    implicit none

    call de_alloc( eo, 'eo', 'siesta_init' )
    call de_alloc( qo, 'qo', 'siesta_init' )
    nullify(eo, qo)
  end subroutine eoqo_reset

end module m_eo





