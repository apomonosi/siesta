!!
!! This file contains modules that are true dummies, used only
!! to avoid more complex dependencies within siesta. Most of these
!! are only used in read_options.
module siesta_master
  implicit none
  private
  character(len=132), public :: input_file = ' '
end module siesta_master

module option_charges_m
  implicit none

  type :: option_charge_t
    integer :: STEPS_ALLOWED
    integer :: step = 0
    integer :: format = 0
  end type
  public :: option_charge_t

  type :: option_charges_t
    type(option_charge_t) :: hirshfeld
    type(option_charge_t) :: voronoi
    type(option_charge_t) :: mulliken
    type(option_charge_t) :: spin
  end type
  public :: option_charges_t

end module option_charges_m