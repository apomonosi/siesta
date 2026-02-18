subroutine set_local_handlers()

  use local_handlers_m, only: die, bye, message
  use sys, only: set_die_handler, set_bye_handler
  use sys, only: set_message_handler

!  use alloc_handlers_m, only: alloc_memory_event
!  use alloc, only: set_alloc_event_handler

  call set_die_handler(die)
  call set_bye_handler(bye)
  call set_message_handler(message)

!  call set_alloc_event_handler(alloc_memory_event)
  
end subroutine set_local_handlers

