module qmmm_reinit_m
  !! This module contains a single subroutine that initializes the
  !! fdf reading and writing for the QMMM driver.
  public :: qmmm_reinit

contains
  subroutine qmmm_reinit( sname )
    !! Initializes the data reading for SIESTA, setting the system label for
    !! the QMMM driver.
    use fdf         , only : fdf_init, fdf_set_unit_handler, fdf_get
    use files       , only : slabel, label_length
    use precision   , only : dp
    use parallel    , only : Node
    use qmmm_files_m, only : qmmm_files
    use units       , only : inquire_unit

    implicit none
    character(len=*), intent(out) :: sname
      !! System name.

    character(len=20)  :: tstring, fdfout, sname_default, slabel_default
    character(len=150) :: line
    integer            :: count, length, lun, lun_tmp, lun_block_tmp, ios
    logical            :: debug_input, file_exists

    if ( Node == 0 ) then
      write(6,'(/a)') '                           ****************************  '
      write(6,'(a)')  '                           *  WELCOME TO SIESTA-QMMM  *  '
      write(6,'(a)')  '                           ****************************  '

      ! Dump data file to output file and generate scratch file for FDF to read
      ! from (except if QMMM_INPUT_DEBUG exists)

      inquire( file = 'QMMM_INPUT_DEBUG', exist = debug_input )
      if ( debug_input ) then ! Read from debug input
        write(6,'(a)') 'WARNING: Reading input from file QMMM_INPUT_DEBUG'

        call io_assign( lun )
        qmmm_files%input = 'QMMM_INPUT_DEBUG'
        open( lun, file = 'QMMM_INPUT_DEBUG', form = 'formatted', status = 'old' )
        rewind( lun )
      else ! Read from standard input.
        write(6,'(/a)') 'qmmm_reinit: Reading from standard input'
        lun = 5
        call io_assign( lun_tmp )

        do while (.true.)
          call system_clock( count )
          write( tstring, * ) count

          qmmm_files%input = 'INPUT_TMP.' // adjustl( tstring )
          inquire( file = qmmm_files%input, exist = file_exists )

          if ( .not. file_exists ) exit
        end do

        open( lun_tmp, file = qmmm_files%input, form = 'formatted', &
              status = 'replace' )
        rewind( lun_tmp )
      endif

      write( 6, '(a,23(1h*),a,28(1h*))' ) '***', ' Dump of input data file.'
      call io_assign( lun_block_tmp )
      open( lun_block_tmp, file = 'qmmm_block_input', form = 'formatted', &
            status = 'replace')

      ios = 0
      do while ( .true. )
        read( lun, iostat = ios, fmt = '(a)' ) line
        if ( ios /= 0 ) exit

        length = len_trim( line )
        if ( length > 0 ) then
          write( 6, '(a)' ) line(1:length)

          if ( .not. debug_input ) write( lun_tmp, '(a)' ) line(1:length)
          write( lun_block_tmp, '(a)' ) line(1:length)
        endif
      enddo
      write( 6, '(a,23(1h*),a,29(1h*))' ) '***', ' End of input data file.'

      ! Choose proper file for fdf processing
      if ( debug_input ) then
        call io_close( lun )
        call io_close( lun_block_tmp )
      else
        call io_close( lun_tmp )
        call io_close( lun_block_tmp )
      endif

      fdfout = 'fdf.log'
      call fdf_init( trim(qmmm_files%input), fdfout )
      call fdf_set_unit_handler( inquire_unit )

      ! Defile Name of the system ...
      sname_default = 'SiestaSystem'
      sname = fdf_get( 'SystemName', sname_default )

      write( 6, '(/a,71(1h-))' ) 'qmmm_reinit: '
      write( 6, '(a,a)' ) 'qmmm_reinit: System Name: ', trim(sname)
      write( 6, '(a,71(1h-))' ) 'qmmm_reinit: '

      ! Defile System Label (short name to label files) ...
      slabel_default = 'siesta'
      slabel = fdf_get( 'SystemLabel', slabel_default )
      write( 6, '(a,a)' ) 'qmmm_reinit: System Label: ', slabel
      write( 6, '(a,71(1h-))' ) 'qmmm_reinit: '
    endif
  end subroutine qmmm_reinit
end module qmmm_reinit_m
