module qmmm_ioxv_m
  !! This module handles XV restart input/output in the form of
  !! XV files.
  implicit none
  public :: qmmm_ioxv

  private
  logical :: first_call = .true.
    !! Indicates whether is the first time we are calling this routine.
  character(len=255) :: fname
    !! Name for the XV restart file.
contains

  subroutine qmmm_ioxv( slabel, task, istep, rststep, nat, ucell, vucell, &
                        r, v, fxv , fv )
    !! Saves positions and velocities.
    use precision, only : dp
    use sys      , only : die

    implicit none
    character(len=*), intent(in)    :: slabel
      !! The system label or name.
    character(len=1), intent(in)    :: task
      !! R for read, W for write.
    integer         , intent(in)    :: nat
      !! Number of atoms in unit cell.
    integer         , intent(in)    :: istep
      !! Current MD step.
    real(dp)        , intent(inout) :: ucell(3,3)
      !! Unit cell vectors.
    real(dp)        , intent(inout) :: vucell(3,3)
      !! Unit cell "velocity", i.e. the derivative of the
      !! unit cell vectors with respect to time.
    real(dp)        , intent(inout) :: r(3,nat)
      !! Atomic positions.
    real(dp)        , intent(inout) :: v(3,nat)
      !! Atomic velocities.
    integer         , intent(inout) :: rststep
      !! Step to read from restart file.
    logical         , intent(inout) :: fxv
      !! On output, indicates whether the XV file has been found.
    logical         , intent(inout) :: fv
      !! On output, indicates whether atomic velocities have been found.

    ! Externals
    external :: io_assign, io_close

    ! Other internal variables.
    integer  :: ia, iu, iv, na, ios

    ! Find name of file
    if ( first_call ) then
      fname      = trim(slabel)//'.qmmm.XV'
      first_call = .false.
    endif
    ios = 0

    ! Choose between read or write
    if ( (task == 'r') .or. (task == 'R') ) then
      ! Check if input file exists
      fxv = .false.
      inquire( file = fname, exist = fxv )

      if ( fxv ) then
        call io_assign( iu )
        open( iu, file = fname, status = 'old' )

        write( 6, '(/,a)' ) &
          'qmmm_ioxv: Reading coordinates and velocities from file.'

        do iv = 1, 3
          read( iu, *, iostat = ios ) ucell(1,iv) , ucell(2,iv) , ucell(3,iv),&
                                      vucell(1,iv), vucell(2,iv), vucell(3,iv)
          call check_ios( ios, -1 )
        enddo
        read( iu, *, iostat = ios ) na, rststep
        call check_ios( ios, -2 )

        if ( na /= nat ) &
          call die('qmmm_ioxv: Wrong number of atoms on XV restart!')

        do ia = 1, nat
          read( iu, *, iostat = ios ) r(1,ia), r(2,ia), r(3,ia), &
                                      v(1,ia), v(2,ia), v(3,ia)
          call check_ios( ios, ia )
        enddo

        call io_close( iu )

        ! Check if velocities are present.
        do iv = 1, 3
        do ia = 1, nat
          if ( abs( v(iv,ia) ) > 0.0_dp ) fv = .true.
        enddo
        enddo

      else
        write( 6, '(/,a)' ) 'qmmm_ioxv: WARNING: XV file not found.'
        return
      endif

    elseif ( (task == 'w') .or. (task == 'W') ) then
      ! Write coordinate and velocity restart.
      call io_assign( iu )
      open( iu, file = fname, form = 'formatted', status = 'unknown' )

      do iv = 1, 3
        write( iu, '(2(3x,3f18.9))' ) ucell(1,iv) , ucell(2,iv) , ucell(3,iv) ,&
                                      vucell(1,iv), vucell(2,iv), vucell(3,iv)
      enddo
      write( iu, * ) nat, istep

      do ia = 1, nat
        write( iu, '(3f18.9,3x,3f18.9)' ) r(1,ia), r(2,ia), r(3,ia), &
                                          v(1,ia), v(2,ia), v(3,ia)
      enddo
      call io_close( iu )
    endif

  end subroutine qmmm_ioxv

  subroutine check_ios( ios, ierr )
    !! Checks file IO status.
    use sys, only : die
    implicit none
    integer, intent(in) :: ios
      !! R/W io status.
    integer, intent(in) :: ierr

    character(len=255) :: msg
    character(len=6)   :: iat

    if ( ios == 0 ) return

    msg = 'qmmm_ioxv: problem reading from file.'
    if ( ierr == -1 ) msg = trim(msg) // ' Could not read cell vectors.'
    if ( ierr == -2 ) msg = trim(msg) // &
                      ' Could not read number of atoms and step.'
    if ( ierr > 1 ) then
      write( iat, '(I6)' ) ierr
      msg  = trim(msg) // ' Error while reading atom ' // trim(iat) // '.'
    endif
    call die( msg )
  end subroutine check_ios
end module qmmm_ioxv_m
