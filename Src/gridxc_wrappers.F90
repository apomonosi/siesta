! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------------------------------------
subroutine gridxc_timer_start(str)
  character(len=*), intent(in)  :: str
  call timer("gridxc@"//trim(str),1)
end subroutine
subroutine gridxc_timer_stop(str)
  character(len=*), intent(in)  :: str
  call timer("gridxc@"//trim(str),2)
end subroutine
