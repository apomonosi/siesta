! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This is a lightweight wrapper layer to offer a set of handlers for termination
! and messaging (there could be more).
! The main entries are actually procedure pointers, with defaults appropriate to
! 'small, serial' clients. Larger programs might use the handler setters to change
! the functionality. (See examples in siesta_handlers_m.F90 and tbt_handlers_m.F90)

      module sys
!
!     Termination and messaging routines, MPI aware
!      
      implicit none

      public :: die      ! Prints an error message and terminates the program
      public :: message  ! Prints a message string if node==0
      public :: reset_messages_file  ! 
      public :: bye      ! Cleans up and exits the program
      
      interface
         subroutine message_interf(level,str)
         character(len=*), intent(in) :: level
         character(len=*), intent(in) :: str
         end subroutine message_interf
      end interface

      interface
         subroutine reset_messages_file_interf()
         end subroutine reset_messages_file_interf
      end interface

      interface
         subroutine die_interf(str)
         character(len=*), intent(in) :: str
         end subroutine die_interf
      end interface

      interface
         subroutine bye_interf(str)
         character(len=*), intent(in) :: str
         end subroutine bye_interf
      end interface

      procedure(die_interf), pointer :: die => simple_die_routine
      procedure(bye_interf), pointer :: bye => simple_bye_routine
      procedure(message_interf), pointer :: message => simple_message_routine
      procedure(reset_messages_file_interf), pointer :: reset_messages_file=> simple_reset_routine

      public :: set_die_handler
      public :: set_bye_handler
      public :: set_message_handler
      public :: set_reset_message_handler
      
      private

      CONTAINS

        subroutine set_die_handler(func)
          procedure(die_interf) :: func
          die => func
        end subroutine set_die_handler
        
        subroutine set_bye_handler(func)
          procedure(bye_interf) :: func
          bye => func
        end subroutine set_bye_handler
        
        subroutine set_message_handler(func)
          procedure(message_interf) :: func
          message => func
        end subroutine set_message_handler

        subroutine set_reset_message_handler(func)
          procedure(reset_messages_file_interf) :: func
          reset_messages_file => func
        end subroutine set_reset_message_handler

! auxiliary routine to provide a non-zero exit code
  subroutine exit(code)
    use iso_c_binding, only: C_INT

    integer(C_INT), intent(in) :: code

    interface
       subroutine c_exit(code) bind(C,name="exit")
         use iso_c_binding, only: c_int
         integer(c_int), intent(in) :: code
       end subroutine c_exit
    end interface

    call c_exit(code)
  end subroutine exit

! --------------------------------
  subroutine simple_die_routine(str)
    character(len=*), intent(in)   :: str
    write(0,'(a,a)') "[error]: " // trim(str)
    write(6,'(a,a)') "[error]: " // trim(str)
    !
    call exit(1)
    ! alternatively
    ! stop
    ! In F2018:
    ! error stop 1
    ! or 'C' abort to get a traceback
  end subroutine simple_die_routine
! --------------------------------
  subroutine simple_bye_routine(str)
    character(len=*), intent(in)   :: str
    write(0,'(a,a)') "[bye]: " // trim(str)
    write(6,'(a,a)') "[bye]: " // trim(str)
    !
    call exit(0)
    ! Alternatively:
    ! stop
  end subroutine simple_bye_routine
! --------------------------------
  subroutine simple_message_routine(level,str)
    character(len=*), intent(in)   :: level
    character(len=*), intent(in)   :: str
    write(0,'(a,a,a)') trim(level) // ": " // trim(str)
    write(6,'(a,a,a)') trim(level) // ": " // trim(str)
    !
  end subroutine simple_message_routine
! --------------------------------
  subroutine simple_reset_routine()
    !
  end subroutine simple_reset_routine
! --------------------------------

      end module sys
