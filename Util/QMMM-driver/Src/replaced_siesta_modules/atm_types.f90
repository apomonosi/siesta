module atm_types
  !! Reduction of atm_types in siesta. Contains the data for QM atomic
  !! species.
  use precision, only : dp

  ! Species_info: Consolidate all the pieces of information in one place
  type, public :: species_info
    character(len=2)  ::  symbol       ! Atom symbol
    character(len=20) ::  label
    integer           ::  z            ! Atomic number
    real(dp)          ::  mass         ! Atomic mass
    real(dp)          ::  zval         ! Valence charge
    integer           ::  norbs = 0    ! Total number of orbitals
    real(dp)          ::  orb_pop(100) ! Orbital population
  endtype species_info

  integer, public :: nspecies
  type(species_info), allocatable, public :: species(:)
end module atm_types