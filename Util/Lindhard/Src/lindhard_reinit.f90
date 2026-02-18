module lindhard_reinit_m
   !! This module contains the subroutine to initialize the reading of
   !! input data for the Lindhard utility. It reads the system name and
   !! the system label from the FDF file.

   implicit none
   private
   public :: lindhard_reinit

contains

  subroutine lindhard_reinit( sname )
    !! Initializes the reading of input data for the Lindhard utility. This
    !! must be called from the node 0.
    use fdf      , only : fdf_init, fdf_get
    use files    , only : slabel, label_length
    use m_io     , only : io_close, io_assign
    use precision, only : dp

    implicit none
    character(len=*), intent(inout) :: sname
      !! Name of the system.

    character(len=20)  :: string, fileout, sname_default
    character(len=150) :: line, filein
    character(len=59)  :: slabel_default
    integer            :: count, length, lun, lun_tmp, ios
    logical            :: debug_input, file_exists

    write(6,'(a)') '                  * 2D Lindhard Function Calculation  *    '

    inquire( file = 'INPUT_DEBUG', exist = debug_input )
    if ( debug_input ) then
      write(6,'(a)') 'reinit: Reading input from file INPUT_DEBUG.'

      call io_assign( lun )
      filein = 'INPUT_DEBUG'
      open( lun, file = 'INPUT_DEBUG', form = 'formatted', status = 'old' )
      rewind( lun )

    else ! Read from standard input
     write(6,'(/a)') 'reinit: Reading options from standard input'

      do
        call system_clock( count )
        write( string, * ) count
        filein = 'INPUT_TMP.'//adjustl( string )
        inquire( file = filein, exist = file_exists )
        if ( .not. file_exists ) exit
      enddo

      call io_assign(lun_tmp)
      open( lun_tmp, file = filein, form = 'formatted', status = 'replace' )
      rewind( lun_tmp )
    endif
    write( 6,'(a,23(1h*),a,28(1h*))' ) '***', ' Dump of input data file '

    lun = 5
    ios = 0
    do
      read( lun, iostat = ios, fmt = '(a)' ) line
      if ( ios /= 0 ) exit
      call chrlen( line, 0, length )
      if ( length /= 0 ) then
        write(6,'(a)') line(1:length)
        if ( .not. debug_input ) write(lun_tmp,'(a)') line(1:length)
      endif
    enddo
    write(6,'(a,23(1h*),a,29(1h*))') '***', ' End of input data file '

    ! Choose proper file for fdf processing
    if ( debug_input ) then
      call io_close(lun)
    else
      call io_close(lun_tmp)
    endif

    fileout = 'fdf.log'
    call fdf_init(filein,fileout)
    write(6,*) 'reinit: FDF initialized.'

    ! Get system name and label from FDF
    sname_default = ' '
    sname = fdf_get( 'SystemName', sname_default )
    write(6,'(a,a)') 'reinit: System Name: ', trim(sname)

    slabel_default  = 'siesta'
    slabel = fdf_get( 'SystemLabel', slabel_default )
    write(6,'(a,a)') 'reinit: System Label: ',slabel
  end subroutine lindhard_reinit

  subroutine chrlen( string, nchar, lchar )
    !! Returns the length of a string, ignoring trailing blanks and null.
    implicit none

    character(*), intent(in)  :: string
      !! String to be considered.
    integer     , intent(in)  :: nchar
      !! Length of the string to be considered.
    integer     , intent(out) :: lchar
      !! Length of the string without trailing blanks and null.

    lchar = nchar
    if ( lchar <= 0 ) lchar = len( string )

    if (lchar < 1) return
    do while( ( (string(lchar:lchar) == ' ') .or. &
                (string(lchar:lchar) == CHAR(0)) ) )
      lchar = lchar - 1

      if (lchar < 1) exit
    enddo
  end subroutine chrlen

end module lindhard_reinit_m
