module ioxvrestr_m
  !! This module contains a single subroutine to save a position and
  !! velocity restart at each restrained MD step.
  implicit none
  public :: ioxvrestr

contains
  subroutine ioxvrestr( nAtm, ucell, rAtm, vAtm, step )
    !! Saves positions and velocities at each restrained step.
    use fdf      , only : fdf_string
    use precision, only : dp

    implicit none
    integer , intent(in) :: nAtm
      !! Total number of atoms.
    real(dp), intent(in) :: ucell(3,3)
      !! Unit cell vectors.
    real(dp), intent(in) :: ratm(3,nAtm)
      !! Atomic positions.
    real(dp), intent(in) :: vatm(3,nAtm)
      !! Atomic velocities.
    integer , intent(in) :: step
      !! Current step.

    character(len=164) :: sname
    character(len=255) :: fname
    character(len=20)  :: step_str
    integer            :: ia, iu

    ! Called from siesta.
    external           :: io_assign, io_close

    ! Converts current step number to string and then
    ! build the output filename.
    write( step_str, '(I20)' ) step

    sname = fdf_string( 'SystemLabel', 'siesta' )
    fname = trim(sname) // '.XV.'
    fname = trim(fname) // trim(adjustl(step_str))

    ! Open file
    call io_assign( iu )
    open( iu, file = fname, form = 'formatted', status = 'unknown' )

    ! Write data on file
    write( iu, '(3x,3f18.9)' ) ucell(1,1), ucell(2,1), ucell(3,1)
    write( iu, '(3x,3f18.9)' ) ucell(1,2), ucell(2,2), ucell(3,2)
    write( iu, '(3x,3f18.9)' ) ucell(1,3), ucell(2,3), ucell(3,3)

    write( iu, '(I8)') nAtm
    do ia = 1, nAtm
      write( iu, '(3f18.9,3x,3f18.9)' ) ratm(1,ia), ratm(2,ia), ratm(3,ia), &
                                        vatm(1,ia), vatm(2,ia), vatm(3,ia)
    enddo

    ! Close file
    call io_close( iu )
  end subroutine ioxvrestr
end module ioxvrestr_m
