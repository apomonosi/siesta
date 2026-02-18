module atmfuncs
  !! This module is intended as a simplification of SIESTA's atmfuncs module.
  !! It contains only a few selected functions, working as an accessor for
  !! the data stored in atm_types.
  implicit none
  private
  public :: floating, izofis, labelfis, massfis

contains
  function floating(is)
    !! Returns whether this species is a floating (ghost) atom.
    logical floating
    integer, intent(in) :: is
      !! Species index.

    floating = izofis(is) < 0
  end function floating

  function izofis( is )
    !! Returns the species' atomic number.
    use atm_types, only : species
    integer :: izofis
    integer, intent(in) :: is
      !! Species index.

    izofis = species(is)%z
  end function izofis

  function massfis( is )
    !! Returns the species' atomic mass.
    use atm_types, only : species
    use precision, only : dp
    real(dp) :: massfis
    integer, intent(in) :: is
      !! Species index.

    massfis=species(is)%mass
  end function massfis

  function labelfis ( is )
    !! Returns the species' label.
    use atm_types, only : species
    character(len=20)   :: labelfis
    integer, intent(in) :: is
      !! Species index.

    labelfis= species(is)%label
  end function labelfis
end module atmfuncs