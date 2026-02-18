! 
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units
  
#ifdef SIESTA__UNITS_ORIGINAL
  use units_legacy_m
#else
  use units_codata2018_m
#endif

  implicit none

  public

  ! Also public : unit_convfac

  integer, parameter, private :: dp = selected_real_kind(14,100)

  !--------------------------------
! Internally, siesta works with length: Bohr.
!                               energy: Rydberg.
!                                 time: femtosecond
!  real(dp), parameter :: Bohr   = 1.0_dp
!  real(dp), parameter :: Rydberg = 1.0_dp
!  real(dp), parameter :: Femtosecond = 1.0_dp
!
!  Ang = Bohr / 0.529177
!   eV = Rydberg / 13.60580
!  Joule = eV / 1.6e-19_dp
!  Meter = Ang / 1.0e-10_dp
!  Pascal = Joule/Meter**3
!   kBar  = Pascal * 1.0e4
!  Ryd^-1 (time) = fs/0.04837769
!   .... and so on.
!
! amu is the conversion from standard atomic mass (Da) to the units
! that SIESTA uses when calculating the kinetic energy of atomic
! nuclei, and for the atomic equations of motion in MD.
!
! The amu is defined 1/12 of the mass of a C-12 isotope (about 1822 times
! the mass of the electron); so in principle:
!
! 1 amu = 1.66053906660e-27 kg
! Thus, amu = 1.66053906660e-27 * Joule * second^2 / meter^2

! For consistency, these should be better defined in terms of the units table.
! Parameter initialization cannot use non-intrinsic functions, though.
! These could be then made into 'protected variables' and initialized
! in a call to a "units initialization function".

  !! real(dp), public, protected :: Ang
  !! ...
  !! subroutine units_initialization()
  !!    Ang = unit_conversion_factor('ang','bohr',...)
  !!    eV  = unit_conversion_factor('eV','ryd',...)
  !!    ...
  !! end subroutine units_initialization

CONTAINS
  
  ! This function is a drop-in replacement for 'fdf_convfac' that can be
  ! used independently of fdf (in particular, in phases before fdf is set up).
  !
    function unit_convfac(from, to) result(factor)
       use sys, only: die
       character(len=*), intent(in)   :: from, to
       real(dp)                       :: factor

       integer :: stat
       character(len=256) :: msg

       factor = unit_convfac_work(from, to, stat, msg)

       if (stat /= 0) then
          call die('unit_convfac' // trim(msg))
       endif
     end function unit_convfac
! ------
     
    function unit_convfac_work(from,to,stat,msg) result (factor)

      use units_common_m, only: leqi
      
    implicit none

    character(len=*), intent(in)   :: from, to
    integer, intent(out)           :: stat
    character(len=*), intent(out)  :: msg
    real(dp)                       :: factor

    character(len=20)      :: phys_dim_to, phys_dim_from
    character(len=20)      :: unit_name_to, unit_name_from
    character(len=40)      :: new_from
    real(dp)               :: value_to, value_from
    
    call inquire_unit(to, stat, phys_dim_to, unit_name_to, value_to)
    if (stat == -1) then
       msg = 'Unknown unit = ' // trim(to)
       RETURN
    else if (stat == 1) then
       msg = 'Ambiguous unit (please fix the code to specify physical dimension) = ' // trim(to)
       RETURN
    endif
       
    call inquire_unit(from, stat, phys_dim_from, unit_name_from, value_from)
    if (stat == -1) then
       msg = 'Unknown unit = ' // trim(from)
       RETURN
    else if (stat == 1) then
       ! "from" unit is ambiguous. 
       if (len_trim(phys_dim_from) > 0) then
          msg = 'Unit name ' // trim(from) // &
               ' is ambiguous, even with qualification.'
          RETURN
       else
          ! Try casting the physical dimension to that of 'to'
          new_from = trim(phys_dim_to) // ":" // trim(from)
          call inquire_unit(new_from, stat, phys_dim_from, unit_name_from, value_from)
          if (stat == -1) then
             msg = 'Unit name ' // trim(from) // &
                  ' is ambiguous and cast to target (' // trim(new_from) // ') does not exist'
             RETURN
          else if (stat == 1) then
             msg = 'Ambiguous unit even after casting! (case sensitivity needed?) = ' &
                   // trim(new_from)
             RETURN
          else
          ! Do nothing. 
          endif
       endif
       
    else
       ! Do nothing.
          
    endif

    ! Final checks
    if (.not. leqi(phys_dim_to, phys_dim_from)) then
       msg = "Incompatible dimensions: " &
            // trim(phys_dim_to) // ":" // trim(unit_name_to) // " , " &
            // trim(phys_dim_from) // ":" // trim(unit_name_from)
       stat = -1
       RETURN
    endif

    factor = value_from / value_to

  END FUNCTION unit_convfac_work

end module
