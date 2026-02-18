module qmmm_timers
  !! This module handles the timers that are present in the
  !! SIESTA QMMM driver. Its main interface with the rest
  !! of the program is qmmm_timer, the rest is to be kept
  !! encapsulated here.
  use precision, only: dp

  implicit none
  public :: qmmm_timer

  private
  integer, parameter :: NMAX = 1000
    !! Maximum number of timers.
  character(len=12) :: progs(NMAX)
    !! List of timer names.
  integer :: nprogs
    !! Current total number of timers.
  integer :: ncalls(NMAX)
    !! Number of calls for each timer.
  real(dp) :: count_rate
    !! Integer to seconds conversion.
  real(dp) :: time0
    !! Starting time offset for all timers.
  real(dp) :: time1(NMAX)
    !! Current time start for each timer.
  real(dp) :: timet(NMAX)
    !! Current accumulated total time for each timer.
  logical  :: first = .true.
    !! Whether this is the first call.

contains

  subroutine qmmm_timer( prog, iopt )
    !! Timer handle for external routines.
    implicit none
    character(len=*), intent(in) :: prog
      !! Timer name.
    integer         , intent(in) :: iopt
      !! Option: 0 resets, 1 begins, 2 ends, 3 prints.

    call qmmm_timetot( prog, iopt )
  end subroutine qmmm_timer

  subroutine qmmm_timetot(prog, iopt )
    !! Finds and prints the total time spent in a given routine or program
    !! section. It must be called with iopt = 1 at the beginning of each
    !! routine and 0 at the end of it.
    !
    ! The following options are available with iopt:
    !   IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
    !   IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
    !   IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
    !   IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
    use precision , only : dp

    implicit none
    character(len=*), intent(in) :: prog
      !! Timer name.
    integer         , intent(in) :: iopt
      !! Option: 0 resets, 1 begins, 2 ends, 3 prints.

    real(dp), parameter :: ZERO = 0.0_dp, HUNDRED = 100.0_dp, &
                           TIMMIN = 1.0e-6_dp

    integer  :: time_int, count_rate_int, iprog
    logical  :: new_timer
    real(dp) :: time, timetl, timtot, avgtme, fractn

    if ( first ) then
      call system_clock ( count_rate = count_rate_int )
      count_rate = real( count_rate_int, kind=dp )
      first = .false.
    endif
    call system_clock (time_int)

    time = time_int / count_rate

    if ( iopt == 0 ) then
      nprogs = 0
      time0  = time
    elseif ( iopt == 1 .or. iopt == 2 ) then
      new_timer = .true.

      do iprog = 1, nprogs
        if ( progs(iprog) == prog ) then
          new_timer = .false.
          exit
        endif
      enddo

      if ( new_timer ) then
        nprogs = nprogs +1
        if ( nprogs > NMAX ) then
          write( 6, * ) 'timer: Maximum number of timers reached. Timer: ', prog
          return
        endif

        iprog = nprogs
        progs(iprog)  = prog
        ncalls(iprog) = 0
        timet(iprog)  = ZERO
      endif

      if ( iopt == 1 ) then
        ncalls(iprog) = ncalls(iprog) +1
        time1(iprog)  = time
      else
        timet(iprog)  = timet(iprog) + time - time1(iprog)
      endif

    elseif ( iopt == 3 ) then
      timtot = time - time0
      if ( timtot < TIMMIN) return
      if ( (prog == 'ALL') .or. (prog == 'all') ) then
        write( 6, '(/,A)' ) 'timer: CPU execution times (min):'
        write( 6, '(A,2X,A10,A9,2A12,A9)') 'timer:', 'Routine   ', 'Calls',&
                                           'Time/call', 'Tot.time', '%'
        do iprog = 1, nprogs
          timetl = timet(iprog)
          avgtme = timet(iprog) / ncalls(iprog)
          fractn = timet(iprog) / timtot * HUNDRED

          write( 6, '(A,2X,A10,I9,2F12.3,F9.2)' )  &
            'timer:', progs(iprog), ncalls(iprog), &
            avgtme / 60.0_dp, timetl / 60.0_dp, fractn
        enddo
        write( 6, * ) ' '
      else
        do iprog = 1, nprogs
          if ( progs(iprog) == prog ) then
            timetl = timet(iprog)
            fractn = timet(iprog) / timtot * HUNDRED
            write( 6, '(A,A10,I6,F12.3,F7.2)' ) &
              'timer: Routine,Calls,Time (min),% = ', &
              progs(iprog), ncalls(iprog), timetl / 60.0_dp, fractn
          endif
        enddo
      endif
    else
      write( 6, * ) 'timer: INVALID OPTION IOPT =', iopt, prog
    endif
  end subroutine qmmm_timetot
end module qmmm_timers
