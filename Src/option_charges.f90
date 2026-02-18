! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module option_charges_m

   implicit none

   private

   !< Never do a charge output.
   integer, parameter, public :: OPT_CHARGE_NONE = 0

   !< Write out after every geometry initialization.
   !!
   !! This happens just after the gemoetry will start a new SCF
   !! cycle. It thus uses the just created density, either through
   !! interpolation, or through initialization.
   integer, parameter, public :: OPT_CHARGE_INIT = int(b'1')

   !< Write out after every geometry SCF convergence.
   integer, parameter, public :: OPT_CHARGE_GEOMETRY = int(b'10')

   !< Write charges during SCF.
   integer, parameter, public :: OPT_CHARGE_SCF = int(b'100')

   !< Write charges during SCF, but just after mixing.
   integer, parameter, public :: OPT_CHARGE_SCF_AFTER_MIX = int(b'1000')

   !< Write out charges after the full calculation has finished.
   integer, parameter, public :: OPT_CHARGE_END = int(b'10000')

   ! Alternatively, one could make
   ! format be an array with Options, if the array has zero length
   ! then nothing will be printed. In this way different charge
   ! populations could be shown differently. I.e. one could
   ! imagine doing hirshfeld only on pz orbitals (omitting Dscf for
   ! other orbitals).

   public :: option_charge_t

   !< The charge option for individual charge methods.
   !!
   !! It holds user-requested information such as:
   !! - when should it be calculated (see the above options)
   !! - the format it should be written in (charge method dependent)
   !!
   !! Likely more options will come.
   type :: option_charge_t

      !< A bit-set specifying which allowed steps it can do
      integer :: STEPS_ALLOWED

      !< Designate which steps the charge calculations will kick in.
      !< Currently only the above listed charge-steps are available.
      integer :: step = OPT_CHARGE_NONE

      !< Specification of the output format.
      !< The code may use this in any way it likes.
      !< Other values may be charge-method dependent.
      integer :: format = 0

   contains

      procedure, public :: run_string => charge_t_char
      procedure, public :: run => charge_t_run
      procedure, public :: can_run => charge_t_can_run
      procedure, private :: charge_t_enable_run, charge_t_enable_run_string
      generic, public :: enable_run => charge_t_enable_run, charge_t_enable_run_string
      procedure, public :: disable_run => charge_t_disable_run

   end type

   public :: option_charges_t

   !< Holder of all possible charge methods.
   !!
   !! A placeholder type, that defines the different methods
   !! avaible, and when they are available.
   type :: option_charges_t
      type(option_charge_t) :: hirshfeld = option_charge_t( &
         OPT_CHARGE_GEOMETRY + OPT_CHARGE_SCF + OPT_CHARGE_END)

      type(option_charge_t) :: voronoi = option_charge_t( &
         OPT_CHARGE_GEOMETRY + OPT_CHARGE_SCF + OPT_CHARGE_END)

      type(option_charge_t) :: mulliken = option_charge_t( &
         OPT_CHARGE_INIT + &
         OPT_CHARGE_GEOMETRY + &
         OPT_CHARGE_SCF + OPT_CHARGE_SCF_AFTER_MIX + &
         OPT_CHARGE_END)

      type(option_charge_t) :: spin = option_charge_t( &
         OPT_CHARGE_INIT + &
         OPT_CHARGE_GEOMETRY + &
         OPT_CHARGE_SCF + OPT_CHARGE_SCF_AFTER_MIX + &
         OPT_CHARGE_END)

      ! TODO: type(option_charge_t) :: moments_spin
   end type

contains

   !< Convert the steps it should print-out, to a string for readability.
   !!
   !! It checks all `%steps` and appends their string-equivalents
   !! to a single string.
   elemental function charge_t_char(this) result(chr)
      class(option_charge_t), intent(in) :: this
      character(len=56) :: chr

      chr = ' '
      if ( this%run(OPT_CHARGE_INIT) ) then
         chr = "init"
      end if
      if ( this%run(OPT_CHARGE_SCF) ) then
         chr = trim(chr) // "+scf"
      end if
      if ( this%run(OPT_CHARGE_SCF_AFTER_MIX) ) then
         chr = trim(chr) // "+scf-after-mix"
      end if
      if ( this%run(OPT_CHARGE_GEOMETRY) ) then
         chr = trim(chr) // "+geometry"
      end if
      if ( this%run(OPT_CHARGE_END) ) then
         chr = trim(chr) // "+end"
      end if
      if ( len_trim(chr) == 0 ) then
         chr = "none"
      else if ( chr(1:1) == "+" ) then
         chr = trim(chr(2:))
      end if

   end function

   !< Function used to check if a certain method can run a particular step.
   !!
   !! In only checks towards the `STEPS_ALLOWED` attribute of the type.
   !!
   !! Example:
   !! ```fortran
   !! type(option_charge_t) :: hello = option_charge_t(OPT_CHARGE_GEOMETRY)
   !! print *, hello%can_run(OPT_CHARGE_END) ! is False 
   !! ```
   elemental function charge_t_can_run(this, step) result(run)
      class(option_charge_t), intent(in) :: this
      integer, intent(in) :: step
      logical :: run
      run = iand(this%STEPS_ALLOWED, step) == step
   end function


   !< Check whether the charge-option has enabled a specific step.
   !!
   !! In only checks towards the `step` attribute of the type.
   !! It assumes that `step` has been populated without violating `can_run`.
   elemental function charge_t_run(this, step) result(run)
      class(option_charge_t), intent(in) :: this
      integer, intent(in), optional :: step
      logical :: run

      run = this%step /= OPT_CHARGE_NONE
      if ( present(step) ) then
         ! > 0 allows one to check for multiple options at the
         ! same time, potentially we could force a `strict`
         ! argument, which checks in a strict sense
         run = run .and. (iand(this%step, step) > 0)
      end if
   end function

   !< Turn on a charge-step calculation.
   !!
   !! It will only turn it on, if it is allowed to do so in the `STEPS_ALLOWED`
   !! variable.
   subroutine charge_t_enable_run(this, step)
      class(option_charge_t), intent(inout) :: this
      integer, intent(in) :: step

      if ( iand(this%STEPS_ALLOWED, step) == step ) then
         this%step = ior(this%step, step)
      end if

   end subroutine

   !< Turn on charge-step calculations based on a string.
   !!
   !! The step-options can be separated
   !! by `+`, `/`, `:` or `&`.
   !! Then it searches for substrings that matches the
   !! option names above.
   subroutine charge_t_enable_run_string(this, step)
      use m_char, only: lcase

      class(option_charge_t), intent(inout) :: this
      character(len=*), intent(in) :: step
      character(len=len_trim(step) + 1):: str

      ! Add + to end of string to ensure we can compare correctly
      str = lcase(trim(step)) // "+"
      if ( has(str, "never") .or. &
           has(str, "none") ) &
           call this%disable_run()
      if ( has(str, "init") ) &
         call this%enable_run(OPT_CHARGE_INIT)
      if ( has(str, "scf") ) &
         call this%enable_run(OPT_CHARGE_SCF)
      if ( has(str, "scf-after-mix") ) &
         call this%enable_run(OPT_CHARGE_SCF_AFTER_MIX)
      if ( has(str, "geometry") .or. has(str, "geom") ) &
         call this%enable_run(OPT_CHARGE_GEOMETRY)
      if ( has(str, "end") ) &
         call this%enable_run(OPT_CHARGE_END)

   contains

      function has(str, substr)
         character(len=*), intent(in) :: str, substr
         logical :: has

         has = index(str, substr // "+") > 0
         if ( has ) return
         has = index(str, substr // "/") > 0
         if ( has ) return
         has = index(str, substr // ":") > 0
         if ( has ) return
         has = index(str, substr // "&") > 0

      end function

   end subroutine

   !< Disables running charge calculations.
   !!
   !! The parameter `step` is optional, and if not supplied, it will
   !! disable all steps.
   subroutine charge_t_disable_run(this, step)
      class(option_charge_t), intent(inout) :: this
      integer, intent(in), optional :: step

      if ( .not. present(step) ) then
         this%step = OPT_CHARGE_NONE
         return
      end if

      if ( iand(this%step, step) == step ) then
         this%step = this%step - step
      end if

   end subroutine

end module
