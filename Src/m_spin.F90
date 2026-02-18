! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module t_spin

  use precision, only: dp
  
  implicit none
  private
  
  !> Type containing a simulations spin configuration
  !!
  !! Thus this type contains _all_ relevant information
  !! regarding the spin-configuration.
  type, public :: tSpin

     !> Dimensionality of the Hamiltonian
     integer :: H = 1
     !> Dimensionality of the density matrix
     integer :: DM = 1
     !> Dimensionality of the energy density matrix
     integer :: EDM = 1
     !> Dimensionality of the grid operations
     integer :: Grid = 1
     !> Dimensionality of the diagonal spin-components
     integer :: spinor = 1

     ! Storing the type of calculation
     
     !> Whether the simulation has no spin
     logical :: none = .true.
     !> Collinear spin
     logical :: Col = .false.
     !> Non-colinear spin
     logical :: NCol = .false.
     !> Spin-orbit coupling
     logical :: SO = .false.

     !> Flavors of off-site implementation
     !> We should use 'full' instead of 'offsite'
     logical :: SO_offsite = .false.
     logical :: SO_onsite  = .false.

     ! Perhaps one could argue that one may
     ! associate a symmetry to the spin which
     ! then denotes whether the spin-configuration
     ! is assumed time-reversal symmetric...
     ! Hmm... 

  end type tSpin

end module t_spin

module m_spin
  use precision, only: dp

  use t_spin, only: tSpin

  implicit none
  
  private

  !> Spin configuration for SIESTA
  type(tSpin), public, save :: spin

  ! Use plain integers instead of pointers, to avoid problems
  ! in the PEXSI-only nodes, which might not call the spin_init
  ! routine. The values are copied at the end of that routine.
  
  ! Create short-hands for the spin-configuration
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%Grid
    integer, save, public, pointer :: nspin => null() ! (Grid)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%spinor
    integer, save, public, pointer :: spinor_dim => null() ! (spinor)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%H, spin%DM
    integer, save, public, pointer :: h_spin_dim => null() ! (H and DM)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%EDM
    integer, save, public, pointer :: e_spin_dim => null() ! (EDM)
  

  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%none
  logical, save, public, pointer :: NoMagn ! (none)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%Col
  logical, save, public, pointer :: SPpol ! (Col)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%NCol
  logical, save, public, pointer :: NonCol ! (NCol)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%SO
  logical, save, public, pointer :: SpOrb ! (SO)

  ! TODO : this is unrelated to the spin-configuration...
  ! Consider moving this to some other place...
  logical, save, public :: TrSym   = .true.

  ! Whether we are performing spiral arrangement of spins
  logical, save, public :: Spiral  = .false.
  ! Pitch wave vector for spiral configuration
  real(dp), save, public :: qSpiral(3) = 0._dp


end module m_spin
