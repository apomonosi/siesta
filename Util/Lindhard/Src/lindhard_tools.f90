module lindhard_tools_m
  !! This module contains the tools to compute the Lindhard function,
  !! including file I/O and fdf options.

  implicit none
  private

  public :: readlindhard
  public :: readeig
  public :: inver3
  public :: fermif

contains

  subroutine readlindhard( temp, q0, q1, q2, ngridx, ngridy, ngridz, nq1, nq2, &
                           firstband, lastband )
    !! Reads the input variables to compute the Lindhard function.
    !! P. Ordejon, September 2014

    use precision, only : dp
    use fdf      , only : fdf_get

    implicit none
    integer, intent(out)  :: nq1
      !! Number of interpolation points in first direction.
    integer, intent(out)  :: nq2
      !! Number of interpolation points in second direction.
    integer, intent(out)  :: ngridx
      !! Number of interpolation grid points in first direction.
    integer, intent(out)  :: ngridy
      !! Number of interpolation grid points in second direction.
    integer, intent(out)  :: ngridz
      !! Number of interpolation grid points in third direction.
    integer, intent(out)  :: firstband
      !! Lowest band to consider in computing the Lindhard Function.
    integer, intent(out)  :: lastband
      !! Highest band to consider in computing the Lindhard Function.
    real(dp), intent(out) :: q0(3)
      !! Point in k-space defining the origin of the 2D grid.
    real(dp), intent(out) :: q1(3)
      !! The vector from q0 to q1 is the first vector of the grid.
    real(dp), intent(out) :: q2(3)
      !! The vector from q0 to q2 is the second vector of the grid.
    real(dp), intent(out) :: temp
      !! Temperature to compute the Lindhard Function (in Kelvin)

    temp = fdf_get( 'Lindhard.Temperature', 0.0d0, 'K' )

    ngridx = fdf_get( 'Lindhard.ngridx', 0 )
    ngridy = fdf_get( 'Lindhard.ngridy', 0 )
    ngridz = fdf_get( 'Lindhard.ngridz', 0 )

    nq1 = fdf_get( 'Lindhard.nq1', 1 )
    nq2 = fdf_get( 'Lindhard.nq2', 1 )

    firstband = fdf_get( 'Lindhard.firstband', 0 )
    lastband  = fdf_get( 'Lindhard.lastband', 0 )

  end subroutine readlindhard

  subroutine readeig(eo, ef, no, ns, nk, get_no)
    !! Reads the eigenvalues from the file with extension .EIG.
    use fdf
    use files    , only : slabel, label_length
    use m_io     , only : io_close, io_assign
    use precision, only : dp
    use sys      , only : die

    implicit          none
    real(kind=dp), intent(inout) :: eo(no,ns,nk)
      !! Eigenvalues.
    real(kind=dp), intent(inout) :: ef
      !! Fermi energy.
    integer      , intent(inout) :: no
      !! Number of orbitals.
    integer      , intent(inout) :: ns
      !! Number of spins components.
    integer      , intent(inout) :: nk
      !! Number of k-points.
    logical      , intent(in)    :: get_no
      !! Flag to read the file for the first time.

    character(len=label_length+4), save :: fname
    integer :: iu, ik, is, jk, io, n, i1, i2, j, ndone

    fname = slabel
    fname = trim(fname) // '.EIG'

    call io_assign( iu )
    open( iu, file=fname, form='formatted', status='unknown' )

    read(iu,*) ef
    read(iu,*)   no, ns, nk
    if ( ns > 1 ) call die('Error: Lindhard not ready for spin components.')

    if ( get_no ) then
      call io_close( iu )
      return
    endif

    is = 1
    do jk = 1,nk
      n = min(no,10)
      ndone = 1
      read(iu,*) ik, (eo(io,is,jk),io=1,n)
      j= no / 10
      n=mod(no,10)
      do i1=1, no/10 -1
        ndone=ndone+1
        read(iu,*) (eo(i2+i1*10,is,jk), i2=1,10)
      enddo
      if (10*ndone < no) then
       read(iu,*) (eo(10*j+i2,is,jk),i2=1,n)
      endif
    enddo

    call io_close( iu )
  end subroutine readeig

  subroutine inver3( a , b )
    !! Inverts a 3x3 matrix
    use precision, only : dp

    implicit none

    real(kind=dp), intent(in)  :: a(3,3)
      !! Input matrix.
    real(kind=dp), intent(out) :: b(3,3)
      !! Inverted matrix.
    integer       :: i
    real(kind=dp) :: c

    b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
    b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
    b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
    b(2,1) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
    b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
    b(2,3) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
    b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
    b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
    b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
    do i = 1, 3
      c = 1.0_dp / (a(1,i)*b(i,1) + a(2,i)*b(i,2) + a(3,i)*b(i,3) )
      b(i,:) = b(i,:) * c
    enddo
  end subroutine inver3

  function fermif( x, temp ) result(f)
    !! Computes the Fermi-Dirac distribution function.
    use precision, only : dp
    use units    , only : Kelvin, eV

    implicit none
    real(kind=dp), intent(in)  :: x
      !! Energy.
    real(kind=dp), intent(in)  :: temp
      !! Temperature in K.
    real(kind=dp) f, conv

    conv = eV / Kelvin
    if ( abs(temp) < 1e-16_dp ) then
      f = 0.5_dp
      if ( x < 0.0_dp ) then
        f = 1.0_dp
      elseif ( x > 0.0_dp ) then
        f = 0.0_dp
      endif
    else
      f = 1.0_dp / ( 1.0_dp + exp((x*conv)/temp) )
    endif
  end function fermif

end module lindhard_tools_m