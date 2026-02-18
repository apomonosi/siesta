!
! Copyright (C) 1996-2021	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine tbt_init_output(IO)

  ! Subroutine to initialize the output file for SIESTA
  ! This should *only* be used to close(6) and open(6) for
  ! the new output file
  !
  ! Taken from reinit, and adapted by Nick R. Papior, 2017

  use files, only: stdout_file
  use tbt_reinit_m, only: tbt_parse_command_line

  implicit none

  ! Whether this node is allowed to perform IO
  logical, intent(in) :: IO

  ! Quick return for non-IO node
  ! Perhaps this should be abber
  if ( .not. IO ) return

  ! First we determine whether all output should be
  ! written to stdout or to a file.
  call tbt_parse_command_line(output_file=stdout_file)
  if ( stdout_file /= ' ' ) then
     ! Close default output to create a new handle
     close(unit=6)
     open(unit=6, file=trim(stdout_file), form="formatted", &
          position="rewind", action="write", status="unknown")
  end if

end subroutine tbt_init_output
