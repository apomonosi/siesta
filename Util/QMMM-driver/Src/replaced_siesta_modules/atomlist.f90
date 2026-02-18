
module atomlist
  !! This module is intended as a simplification of SIESTA's atomlist module.
  !! This reduces the amount of dependencies needed for the driver.
  use precision, only : dp

  implicit none
  real(dp), public :: qtot    = 0.0_dp
    !! Total number of electrons.
  integer         , pointer, public :: iza(:)   => null()
    !! Atomic number of each atom.
  real(dp)        , pointer, public :: amass(:) => null()
    !! Atomic masses.
  character(len=2), pointer, public :: elem(:)  => null()
    !! Contains element names.

contains
  subroutine initatomlists()
    use alloc      , only : re_alloc, de_alloc
    use atm_types  , only : species
    use siesta_geom, only : na_u, na_s, isa, xa_last

    implicit none
    integer :: ia, is, io

    na_s = na_u ! We do not use a supercell within the QMMM driver.
    nullify( iza, amass, elem, xa_last )
    call re_alloc( iza    , 1, na_u, 'iza'  , 'atomlist' )
    call re_alloc( amass  , 1, na_u, 'amass', 'atomlist' )
    call re_alloc( elem   , 1, na_u, 'elem' , 'atomlist' )
    call re_alloc( xa_last, 1, 3, 1, na_u, 'xa_last' , 'atomlist' )

    qtot    = 0.0_dp
    do ia = 1, na_u
      is = isa(ia)

      amass(ia) = species(is)%mass
      iza(ia)   = species(is)%z
      elem(ia)  = species(is)%symbol

      do io = 1, species(is)%norbs
        qtot = qtot + species(is)%orb_pop(io)
      enddo
    enddo
  end subroutine initatomlists
end module atomlist