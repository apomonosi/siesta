!
! These are 'external' versions of the routines in module sys
!
subroutine die(str)
  use sys, only: die_sys => die
  character(len=*), intent(in) :: str
  call die_sys(str)
end subroutine die

subroutine bye(str)
  use sys, only: bye_sys => bye
  character(len=*), intent(in) :: str
  call bye_sys(str)
end subroutine bye

subroutine message(level,str)
  use sys, only: message_sys => message
  character(len=*), intent(in) :: level
  character(len=*), intent(in) :: str
  call message_sys(level, str)
end subroutine message

subroutine reset_messages_file()
  use sys, only: reset_messages_file_sys => reset_messages_file
  call reset_messages_file_sys()
end subroutine reset_messages_file

! Versions of the alloc module handlers for libgridxc
! The alloc_memory_event is dummy (suitable for dependencies, but
! inhibiting Siesta's reporting of gridxc allocations/deallocations.
! This is a temporary kludge to enable compilation when the libgridxc
! library does not use procedure pointers for its handlers.
! 
#ifndef SIESTA__GRIDXC_HAS_PP
subroutine alloc_error_report(str,code)

  character(len=*), intent(in)  :: str
  integer, intent(in)  :: code
  if (code == 0) call die(str)

end subroutine
subroutine alloc_memory_event(bytes,name)
  integer, intent(in) :: bytes
  character(len=*), intent(in)  :: name
  ! do nothing
end subroutine
#endif

