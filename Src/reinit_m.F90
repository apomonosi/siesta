!
! Copyright (C) 1996-2021       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module reinit_m

  implicit none

  private

  public :: reinit, parse_command_line

contains

  subroutine reinit(sname, skip_cmdargs)

    ! Subroutine to initialise the reading of the data for SIESTA
    !
    !     It uses the FDF (Flexible Data Format) package
    !     of J.M.Soler and A.Garcia
    !
    ! Taken from redata. Writen by P.Ordejon, December 1996
    ! **************************** OUTPUT *********************************
    ! character    slabel      : System Label (to name output files)
    ! character(len=*) sname       : System Name
    ! **********************************************************************

    ! Modules
    use parallel,    only : Node
    use fdf
    use files,       only : slabel, label_length
    use files,       only : stdin_file, stdout_file
    use siesta_master, only: input_file  ! fdf data file
    use units,         only: inquire_unit

#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
#endif
    implicit none

    character(len=*), intent(out) :: sname
    ! Prevent command-line processing by Siesta if the API is in use
    logical, intent(in), optional :: skip_cmdargs

    ! Internal variables .................................................
    character(len=50) :: fileout

    integer :: narg, count
    integer :: length, lun, lun_tmp, iostat, ierr
    character(len=256) :: line
    character(len=label_length) :: aux_str
    logical :: stdin_file_found

    logical :: debug_input, file_exists, skip
    character(len=8) :: mydate
    character(len=10) :: mytime

    external :: die
    
    ! Default input-file
    stdin_file = "fdf_input"

    skip = .false.
    if (present(skip_cmdargs)) then
       skip=skip_cmdargs
    endif


    if (Node.eq.0) then
       ! Print Welcome ........................................................
       write(6,'(/a)') '                           *********************** '
       write(6,'(a)')  '                           *  WELCOME TO SIESTA  * '
       write(6,'(a)')  '                           *********************** '
       write(6,'(/a)') ' NOTE: Siesta V5 introduces some new defaults for basis-set generation'
       write(6,'(a)')  ' NOTE: which might lead to slightly different numerical results.'
       write(6,'(a/)') ' NOTE: Please see the manual for full details and compatibility measures.'

       if( skip ) then
          ! skip the command line arguments
          narg = 0
       else
          ! Number of arguments provided in the command line.
          narg = command_argument_count()
       endif

       ! Initialisation: stdin file yet to be determined.
       stdin_file_found = .false.

       ! Set name of file to read from.
       ! Choose proper file for fdf processing
       ! (INPUT_DEBUG if it exists or "standard input",
       ! processed and dumped to a temporary file).
       ! Note this is only done in the master node.

       inquire(file='INPUT_DEBUG',exist=debug_input)
       if ( debug_input ) then
         write(6,'(a)') &
              'WARNING: Siesta is reading its input from file INPUT_DEBUG'
         stdin_file = 'INPUT_DEBUG'
         stdin_file_found = .true.

       else if ( narg > 0 ) then
         ! If the last argument may be the input file name,
         ! call the command line parser.
         call parse_command_line(input_file=aux_str)

         ! If the parser found an input file, this should not be
         ! an empty string.
         if (aux_str /= ' ') then

            ! Check that the file exists (the parser already checked that it
            ! is different from the output file and that it is short enough
            ! to be safely stored in stdin_file).
            inquire(file=aux_str, exist=file_exists)
            if (file_exists) then
               stdin_file = trim(aux_str)
               stdin_file_found = .true.
               write(*,'(/,2a)') 'reinit: Reading from file ' // &
                    trim(stdin_file)
            else
               ! The parser said this would be the input file.
               call die ('Cannot find requested input file "' // &
                    trim(aux_str)//'".')
            end if
          end if
       end if

       ! If stdin_file has not been found yet:
       if (.not. stdin_file_found) then
          if ( input_file/=' ' ) then
             ! Read data from given file
             ! This is the case for running siesta under SIESTA-master
             write(6,'(a)') 'reinit: Reading from file '//trim(input_file)
             call io_assign(lun)
             open(lun,file=input_file,form='formatted',status='old')
             rewind(lun)
          else
             ! Read from standard input
             write(6,'(/a)') 'reinit: Reading from standard input'
             lun = 5
          end if                 ! (input_file/=' ')

          ! In any of these two cases, we will copy the contents of
          ! input_file / standard input to a newly created stdin_file.

          ! Make sure we get a new file
          call io_assign(lun_tmp)
          file_exists = .true.
          do while ( file_exists )
             call system_clock( count )
             write(stdin_file,"(a,i5.5)") 'INPUT_TMP.', mod(count,100000)
             inquire(file=stdin_file, exist=file_exists )
          end do

          ! Open this file to dump the input data.
          open(lun_tmp, file=stdin_file, form='formatted', status='new')
          write(6,'(a)') 'reinit: Dumping input in ' // trim(stdin_file)
          write(6,"(a,23('*'),a,28('*'))") '***', ' Dump of input data file '

          ! Line by line, dump the data...
          do
             ! ... read from input_file/standard input (lun)...
             read(lun,iostat=iostat,fmt='(a)') line
             if ( iostat /= 0 ) exit
             length = len_trim(line)
             ! Skip empty lines
             if ( length /= 0 ) then
                write(6,'(a)') line(1:length)
                ! ... in stdin_file (lun_tmp).
                if (.not. debug_input) write(lun_tmp,'(a)') line(1:length)
             end if
          end do

          ! End data dump.
          write(6,"(a,23('*'),a,29('*'))") '***', ' End of input data file '
          call io_close(lun_tmp)

          ! Be sure to close the input master file,
          ! otherwise we may end up opening too many files.
          if ( lun /= 5 ) call io_close(lun)

          ! "stdin_file" for fdf is now the temporary file.
          ! This was necessary historically to allow
          ! the rewinds involved in fdf operation.

       end if ! block for copying input_file/standard input to stdin_file
    end if ! Node .eq. 0


    !!! Set up fdf. !!!

    ! Choose a 'unique' prefix for the log (and possible debug) fdf files.
    ! The time string may be slightly different in different processors,
    ! depending on the system time.
    call date_and_time(mydate,mytime)
    ! mydate has form ccyymmdd and mytime has form hhmmss.sss, so
    ! the string used below conforms to ISO 8601 with millisecond precision.
    write(fileout,"(a)") 'fdf.' // mydate // 'T' // mytime // ".log"

    if (Node .eq. 0) then
       call fdf_init(stdin_file, trim(fileout))
    else
       ! the other nodes will get the data structure
    endif

    ! Parse the command line (only if not using the Siesta API)
    ! Note that there can be modifications to the fdf structure
    if (.not. skip) then
       call parse_command_line(info=.false.)
    endif

#ifdef MPI
      call broadcast_fdf_struct(0,mpi_comm_world)
#endif

    call fdf_set_unit_handler(inquire_unit)
    
    ! Define Name of the system ...
    sname = fdf_string('SystemName', ' ')
    if (Node.eq.0) then
       write(6,'(/a,71("-"))') 'reinit: '
       write(6,'(a,a)') 'reinit: System Name: ',trim(sname)
       write(6,'(a,71("-"))') 'reinit: '
    end if

    ! Define System Label (short name to label files) ...
    slabel = fdf_string('SystemLabel', 'siesta')
    ! Check that the SystemLabel is not empty.
    count = len_trim(slabel)
    if ( count == 0 ) call die('SystemLabel must be at least 1 character!')
    ! Check that there are no spaces in the SystemLabel
    length = index(slabel, ' ')
    if ( length > 0 .and. length < count ) then
       call die('SystemLabel must *NOT* contain any spaces!')
    end if
    if (Node.eq.0) then
      write(6,'(a,a)') 'reinit: System Label: ',trim(slabel)
      write(6,'(a,71("-"))') 'reinit: '
    end if

  end subroutine reinit

!==============================================================================

  subroutine parse_command_line(info, input_file, output_file)

    ! This subroutine requires one and only one optional argument.
    !
    ! This subroutine may be called with up to 4 purposes:
    ! To retrieve the input file, this
    !    requires present(infile).
    ! To retrieve the output file (-o/-out), this
    !    requires present(outfile).
    ! To process help and version options, this
    !    requires info present and true.
    ! To process all the other options, this
    !    requires info present and false.
    !
    ! All 4 calls should parse the arguments in the same way,
    ! hence the single subroutine.
    !
    ! It is the responsibility of the programmer to ensure that,
    ! when info is present, all MPI processes call this subroutine;
    ! otherwise the call will hang.
    !
    ! Design choice: info is the first argument in order to enforce that
    ! the input_file or output_file keywords are always used
    ! when requesting a file name from the parser.

    use fdf
    use files,        only: label_length
    use cli_m,        only: get_command_arg
#ifdef MPI
    use mpi_siesta
#endif
    use parallel,     only: Node
    use siesta_version_info, only: siesta_print_version
    implicit none

    ! Arguments
    logical, intent(in), optional :: info
    character(len=*), intent(out), optional :: input_file
    character(len=*), intent(out), optional :: output_file

    ! Internal variables
    character(len=*), parameter :: myself = 'parse_command_line'
    character(len=label_length) :: infile
    character(len=label_length) :: outfile
    character(len=10) :: str
    integer :: narg ! number of arguments
    character(len=2048) :: line, line_orig, line2
    integer :: nopts ! number of options
    logical :: process_fdf, process_info
    integer :: ia
#ifdef MPI
    integer :: MPIerror
#endif

    ! Initialise files
    outfile = ' '
    infile = ' '

    ! Ensure that the number of optional arguments is valid.
    nopts = count([present(input_file), present(output_file), present(info)])
    if (nopts .NE. 1) call die('ERROR: badly formed call to '//myself//'.')

    ! Ensure that, if present, file names are long enough.
    if (present(input_file)) then
       if (len(input_file) .GT. len(infile)) &
          call die('ERROR: input_file variable too short.')
    else if (present(output_file)) then
       if (len(output_file) .GT. len(outfile)) &
          call die('ERROR: output_file variable too short.')
    end if

    ! Determine whether info or fdf will be processed in this call.
    if (present(info)) then
       process_info = info
       process_fdf = .not. info
    else
       process_info = .false.
       process_fdf = .false.
    end if

    ! Number of arguments provided in the command line.
    narg = command_argument_count()

    ! Read special variables from the command line
    ia = 0
    do while ( ia < narg )

      ia = ia + 1
      call get_command_arg(ia, line_orig)

      if ( line_orig(1:1) /= '-' ) then
        ! It is not an option it must be the input file
        ! With this the input file may be in between options
        if (len_trim(line_orig) > len(infile)) then
          ! Prevent truncation.
          write(str,'(I0)') len(infile)
          call die ('The argument ('//trim(line_orig)//') is too &
              &long to be used as the input file name, please use &
              &a file name of at most '//str//' characters.')
        else if ( len_trim(infile) > 0 ) then
          call die('There are two arguments thought to be input files: &
              &"'//trim(infile)//'" and "'//trim(line_orig)//'". &
              &Please only supply one input file.')
        else
          infile = trim(line_orig)
          cycle
        end if
      else
        line = line_orig
      end if

      ! Truncate '-' to no '-'
      do while ( line(1:1) == '-' )
        line = line(2:)
      end do

      ! We allow these line
      select case (line)

        !! Options that require a second argument.

      case ('out', 'o', 'fdf', 'L', 'V')
        if ( ia >= narg ) call die('Missing argument on command line, ' &
            // trim(line))
        ia = ia + 1
        call get_command_arg(ia,line2)

        if ( (line == 'out' .or. line == 'o') ) then
          if (len_trim(line2) .GT. len(outfile)) then
            ! Prevent truncation
            write(str,'(I0)') len(outfile)
            call die ('The "'// trim(line_orig) //'" argument (' // &
                trim(line2) // ') is too long to be used as the output &
                &file name, please use a file name of at most ' // &
                str // ' characters.')
          else
            outfile = trim(line2)
          end if
        end if

        if (process_fdf) then
          ! We allow these variations:
          !  FDFLabel=0.1:eV
          !  FDFLabel:0.1:eV
          !  FDFLabel=0.1=eV
          !  "FDFLabel 0.1 eV"
          line2 = cmd_tokenize(line2)

          select case (line)
          case ('L')
            line2 = 'SystemLabel '//trim(line2)
          case ('V')
            if ( index(trim(line2), ' ') == 0 ) then
              ! Default to eV argument; users expect this unit for
              ! applied bias.
              line2 = 'TS.Voltage '//trim(line2)//' eV'
            else
              line2 = 'TS.Voltage '//trim(line2)
            end if
          end select

          call fdf_overwrite(line2)
        end if

        !! Single-argument options.

      case ('electrode', 'elec')
        if (process_fdf) then
          line = 'SaveHS true'
          call fdf_overwrite(line)
          ! This is deprecated, so lets not force write *another* file
          !line = 'TS.HS.Save true'
          !call fdf_overwrite(line)
          line = 'TS.DE.Save true'
          call fdf_overwrite(line)
        end if

      case ('version', 'v')
        ! If a version option is found,
        ! print version information...
        if (process_info) then
          if (Node == 0) call siesta_print_version

          ! ... and quietly end execution
          ! (do not call die to avoid the bye message)
          call pxfflush(6)
#ifdef MPI
          call MPI_Finalize(MPIerror)
#endif
          stop
        end if

      case ('help', 'h')
        ! If a help option is found,
        ! print help information...
        if (process_info) then
          if (Node == 0) call print_help

          ! ... and quietly end execution
          ! (do not call die to avoid the bye message)
          call pxfflush(6)
#ifdef MPI
          call MPI_Finalize(MPIerror)
#endif
          stop
        end if

      case default
        call die('Error: Unknown command line option: "' // &
            trim(line_orig)//'".')

      end select

    end do

    ! If any files were requested, check the parsing of files
    ! and do the assignment if everything is OK.
    if ( present(output_file) .OR. present(input_file) ) then
       if ( outfile /= ' ' .AND. outfile == infile ) then
          call die('Error: Requested output file matches input file, "' // &
               trim(outfile) // '".')
       else if (present(output_file)) then
          output_file = outfile
       else ! present(input_file)
          input_file = infile
       end if
    end if

!------------------------------------------------------------------------------

  contains

    function cmd_tokenize(line) result(tline)
      implicit none

      character(len=*), intent(in) :: line
      character(len=len(line)) :: tline

      integer :: i, n
      n = len(tline)
      tline = line
      do i = 1 , n
         if ( tline(i:i) == ':' .or. &
              tline(i:i) == '=' ) then
            tline(i:i) = ' '
         end if
      end do
    end function cmd_tokenize

  end subroutine parse_command_line

!==============================================================================

  subroutine print_help

    use, intrinsic :: iso_fortran_env, only: stderr => ERROR_UNIT
!   use sys, only : bye
    implicit none

    write(stderr,'(a)')'Help for calling SIESTA'
    write(stderr,'(a)')''
    write(stderr,'(a)')'Usage:'
    write(stderr,'(a)')'  siesta [OPTIONAL ARGUMENTS...]'
    write(stderr,'(a)')''
    write(stderr,'(a)')'OPTIONAL ARGUMENTS:'
    write(stderr,'(a)')'  -help|-h'
    write(stderr,'(a)')'      Only print this help.'
    write(stderr,'(a)')'  -version|-v'
    write(stderr,'(a)')'      Only print version and compilation information.'
    write(stderr,'(a)')'  -out|-o <file>'
    write(stderr,'(a)')'      Instead of writing to stdout, write to file.'
    write(stderr,'(a)')'  -L <name>'
    write(stderr,'(a)')'      Short-hand for setting SystemLabel.'
    write(stderr,'(a)')'  -fdf <label>=<value>[:<unit>]'
    write(stderr,'(a)')'      Set the label to the corresponding value.'
    write(stderr,'(a)')'  -V <value>:<unit>'
    write(stderr,'(a)')'      Short-hand for setting TS.Voltage.'
    write(stderr,'(a)')'  -electrode|-elec'
    write(stderr,'(a)')'      Force Save.HS and TS.DE.Save to true.'
    write(stderr,'(a)')'  <fdf-file>'
    write(stderr,'(a)')'      Use file as fdf-input, you do not need to pipe it in.'
    write(stderr,'(a)')'      If not provided must be piped in: siesta < RUN.fdf'
    write(stderr,'(a)')''
    write(stderr,'(a)')'For further help, please see the SIESTA manual (pdf).'

  end subroutine print_help

end module reinit_m
