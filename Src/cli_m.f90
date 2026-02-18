! ---
! Copyright (C) 2021     	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module cli_m

  implicit none

  private

  public :: get_command_arg

contains

  subroutine get_command_arg(pos, val, length)
    ! Wrapper for get_command_argument that checks the status
    ! in order to prevent otherwise undetected truncations.

    use sys, only: die
    implicit none

    ! Arguments
    integer, intent(in) :: pos               ! argument number (position)
    character(len=*), intent(out) :: val     ! argument retrieved
    integer, intent(out), optional :: length ! argument length

    ! Internal variables
    integer :: ierr          ! status return value
    integer :: length_loc    ! length returned by get_command_argument
    character(len=10) :: str_pos, str_len_val, str_len_arg, str_err

    ! Call the intrinsic subroutine.
    call get_command_argument(pos, val, length_loc, ierr)

    ! Check returned status.
    select case (ierr)
       case (0)
         ! Success.
         if (present(length)) length = length_loc
         return
       case (-1)
         ! Error retrieving the argument: argument too long.
         write(str_pos,'(I0)') pos
         write(str_len_val,'(I0)') len(val)
         write(str_len_arg,'(I0)') length_loc
         call die ('Command argument in position '//trim(str_pos)//' has &
              &length '//trim(str_len_arg)//'. It is too long to be correctly &
              &retrieved by the command line parser. Please use arguments of &
              &at most '//trim(str_len_val)//' characters.')
       case default
         ! Capture other processor-dependent errors (error should be positive).
         write(str_pos,'(I0)') pos
         write(str_err,'(I0)') ierr
         call die ('Error code ('//trim(str_err)//') received when trying to &
              &retrieve the command line argument in position '// &
              trim(str_pos)//'.')
    end select

  end subroutine get_command_arg

end module cli_m

