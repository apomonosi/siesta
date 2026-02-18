! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module fc_matrices_m

  use class_dSpData1D
  use class_dSpData2D
  use class_Sparsity

  implicit none

  private
  save

  !> Derivative of Hamiltonian with respect to current FC displacement
  type(dSpData2D), public :: dHdR_2D
  !> Hamiltonian for negative FC displacement
  type(dSpData2D), public :: fc_H_2D(2)
  !> Derivative of Overlap matrix with respect to current FC displacement
  type(dSpData1D), public :: dSdR_1D
  !> Overlap matrix for FC displacements
  type(dSpData1D), public :: fc_S_1D(2)
  !> Number of supercells in each direction for FC displacements
  integer, public :: fc_nsc(3,2)

end module fc_matrices_m
