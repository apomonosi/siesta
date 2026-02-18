module siestaqmmm_m
  !! This module contains all the helper routines needed for
  !! communication with SIESTA.
  implicit none

  public :: send_coordinates
  public :: recv_forces
  public :: close_communicators
  public :: set_input_file

  private
  type pipes_t
    !! Variables for pipe interface.
    character(len=224)   :: qm_name = 'siesta.qmatoms'
      !! Pipe for QM atom information (coordinates, cell).
    character(len=224)   :: pc_name = 'siesta.mmatoms'
      !! Pipe for MM atom information (coordinates, charges).
    character(len=224)   :: fr_name = 'siesta.allforces'
      !! Pipe for forces and stress.
    integer              :: qm_uid  = 111
      !! File unit for QM information.
    integer              :: mm_uid  = 222
      !! File unit for MM information.
    integer              :: fr_uid  = 333
      !! File unit for QM forces and stress.
  end type pipes_t

  type socket_t
    !! Variables for socket interface.
    character(len=1024) :: host = 'localhost'
      !! Socket host address.
    integer             :: inet = 0
      !! Type of socket, 0 for UNIX socket, 1 for TCP/IP.
    integer             :: port = 10002
      !! Port number, only used for TCP/IP sockets.
    integer             :: id
      !! Socket identifier.
  end type socket_t

  logical :: siesta_launched = .false.
    !! Indicates whether SIESTA has already been launched.
  logical :: using_sockets   = .false.
    !! Indicates whether we are using the socket interface
    !! (if not, pipes).
  character(len=256) :: input_fname = "siestaqmmm.fdf"
    !! The name of the input file for this library.

  type(pipes_t)  :: qmmm_pipes
    !! Variables for pipe interface.
  type(socket_t) :: qmmm_socket
    !! Variables for socket interface.

contains

  subroutine set_input_file( input_file_name )
    !! Sets the name for this library's input file.
    implicit none
    character(len=*), intent(in) :: input_file_name
      !! Input file name.

    input_fname = trim( input_file_name )
  end subroutine set_input_file

  subroutine open_communicators( )
    !! Opens pipes and/or sockets used.
    use f90sockets, only : create_socket
    implicit none

    if ( using_sockets ) then
      call create_socket( qmmm_socket%id, qmmm_socket%inet, &
                          qmmm_socket%port, qmmm_socket%host )

    else
      call system( 'rm -f '//trim(qmmm_pipes%qm_name) )
      call system( 'rm -f '//trim(qmmm_pipes%pc_name) )
      call system( 'rm -f '//trim(qmmm_pipes%fr_name) )

      call system( 'mkfifo '//trim(qmmm_pipes%qm_name) )
      call system( 'mkfifo '//trim(qmmm_pipes%pc_name) )
      call system( 'mkfifo '//trim(qmmm_pipes%fr_name) )

      open( unit = qmmm_pipes%qm_uid, file = qmmm_pipes%qm_name, &
            form = "formatted", status = "old", position = "asis" )
      open( unit = qmmm_pipes%mm_uid, file = qmmm_pipes%pc_name, &
            form = "formatted", status = "old", position = "asis" )
      open( unit = qmmm_pipes%fr_uid, file = qmmm_pipes%fr_name, &
            form = "formatted", status = "old", position = "asis" )
    endif

    write(*,'(A)') "Communicators opened."
  end subroutine open_communicators

  subroutine close_communicators( )
    !! Closes pipes and/or sockets used.
    use f90sockets, only : close_socket

    implicit none
    if ( using_sockets ) then
      call close_socket( qmmm_socket%id )
      call system( 'unlink /tmp/ipi_'//trim(qmmm_socket%host) )

    else
      close( qmmm_pipes%qm_uid )
      close( qmmm_pipes%mm_uid )
      close( qmmm_pipes%fr_uid )

    endif

    siesta_launched = .false.
    write(*,'(A)') "Communicators closed."
  end subroutine close_communicators

  subroutine init_siesta( )
    !! Reads the inputfile and launches siesta.
    use fdf, only : fdf_init, fdf_shutdown, fdf_get, leqi

    implicit none
    character(len=100) :: run_command
    character(len=50)  :: inp_str
    character(len=256) :: siesta_input
    integer            :: n_tasks
    logical            :: is_serial


    write(*,'(A)') "Using "//trim(input_fname)//&
                  " as input file for the QMMM library."
    call fdf_init( trim(input_fname) )

    run_command  = fdf_get( 'ParallelCommand', 'mpirun' )
    siesta_input = fdf_get( 'SiestaInput'    , 'INPUT_DEBUG' )
    n_tasks      = fdf_get( 'nCpus'          , 1)
    is_serial    = fdf_get( 'runSerial'      , .false.)

    call fdf_shutdown( )

    write(*,'(A)') "Using "//trim(siesta_input)//&
                   " as input file for SIESTA."
    !! Checks the original siesta file for pipe/sockets information.
    call fdf_init( trim(siesta_input) )

    inp_str = fdf_get( 'QMMM.Driver.Type', 'socket' )
    if ( leqi( inp_str, 'socket' ) .or. leqi( inp_str, 'sockets' ) ) then
      using_sockets = .true.
    elseif  ( leqi( inp_str, 'pipe' ) .or. leqi( inp_str, 'pipes' ) ) then
      using_sockets = .false.
    else
      stop 'ERROR - SIESTAQMMM: Wrong Driver interface type. '//&
           'Must be "pipe" or "socket".'
    endif

    if ( using_sockets ) then
      qmmm_socket%host = fdf_get( 'QMMM.Driver.Address'   , qmmm_socket%host )
      qmmm_socket%port = fdf_get( 'QMMM.Driver.Port'      , qmmm_socket%port )
      inp_str          = fdf_get( 'QMMM.Driver.SocketType', 'inet' )

      qmmm_socket%host = trim(qmmm_socket%host)//achar(0)
      if ( leqi( inp_str, 'unix' ) ) then
        qmmm_socket%inet = 0
      else if ( leqi( inp_str, 'inet' ) ) then
        qmmm_socket%inet = 1
      else
        stop 'ERROR - SIESTAQMMM: Unrecognized socket type. '//&
             'Must be either "unix" or "inet". '
      endif

    else
      qmmm_pipes%qm_name = fdf_get( 'QMMM.Driver.QMRegionFile', &
                                    qmmm_pipes%qm_name )
      qmmm_pipes%pc_name = fdf_get( 'QMMM.Driver.MMChargeFile', &
                                    qmmm_pipes%pc_name )
      qmmm_pipes%fr_name = fdf_get( 'QMMM.Driver.ForceOutFile', &
                                    qmmm_pipes%fr_name )

    endif
    call fdf_shutdown( )

    !! For pipes, pipes must be created before the call to SIESTA.
    !! However, when using sockets, we must be sure SIESTA is already running.
    if ( .not. using_sockets ) call open_communicators( )

    if ( is_serial ) then
      call system( 'siesta < ' // trim(siesta_input) // ' > siestaqmmm.out &' )

    else
      write( inp_str, * ) n_tasks
      call system( trim(run_command) // ' -n ' // inp_str //&
                   ' siesta < ' // trim(siesta_input) // ' > siestaqmmm.out &' )

    endif

    if ( using_sockets ) call open_communicators( )

    siesta_launched = .true.
  end subroutine init_siesta

  subroutine send_coordinates( n_qm, n_mm, r_qm, r_mm, pc, cell, last )
    !! Sends the coordinates for both MM and QM regions (and related variables)
    !! to SIESTA.
    use f90sockets, only : writebuffer
    use precision , only : dp

    implicit none
    integer , intent(in)  :: n_qm
      !! Total number of atoms in the QM region.
    integer , intent(in)  :: n_mm
      !! Total number of atoms in the MM region.
    real(dp), intent(in)  :: r_qm(3,1:n_qm)
      !! Positions of atoms in the QM region.
    real(dp), intent(in)  :: r_mm(3,1:n_mm)
      !! Positions of atoms in the MM region.
    real(dp), intent(in)  :: pc(1:n_mm)
      !! Classical partial charges of atoms in the MM region.
    real(dp), intent(in)  :: cell(3,3)
      !! Size of the unit cell.
    logical , intent(in)  :: last
      !! Indicates whether this is the last MD step.

    integer :: max_size, iat
    real(dp), allocatable :: dpbuffer(:)

    if ( .not. siesta_launched ) call init_siesta( )

    max_size = 3
    max_size = max( max_size, n_qm + n_mm )
    allocate( dpbuffer(3*max_size) )
    dpbuffer = 0.0_dp

    if ( using_sockets ) then
      dpbuffer(1:9) = reshape( cell, (/9/) )
      call writebuffer( qmmm_socket%id, dpbuffer(1:9), 9 )
      call writebuffer( qmmm_socket%id, n_qm )

      dpbuffer(1:3*n_qm) = reshape( r_qm, (/3*n_qm/) )
      call writebuffer( qmmm_socket%id, dpbuffer(1:3*n_qm), 3*n_qm )
      if ( last ) then
        call writebuffer( qmmm_socket%id, 'last', 4 )
      else
        call writebuffer( qmmm_socket%id, 'goon', 4 )
      endif

      call writebuffer( qmmm_socket%id, n_mm )
      dpbuffer(1:3*n_mm) = reshape( r_mm, (/3*n_mm/) )
      call writebuffer( qmmm_socket%id, dpbuffer(1:3*n_mm), 3*n_mm )
      call writebuffer( qmmm_socket%id, pc, n_mm )

    else
      do iat = 1, 3
        write( qmmm_pipes%qm_uid, * ) cell(:,iat)
      enddo
      write( qmmm_pipes%qm_uid, * ) n_qm
      do iat = 1, n_qm
        write( qmmm_pipes%qm_uid, * ) r_qm(:,iat)
      enddo

      if ( last ) then
        write( qmmm_pipes%qm_uid, * ) 'last'
      else
        write( qmmm_pipes%qm_uid, * ) 'goon'
      endif
      call flush( qmmm_pipes%qm_uid )

      write( qmmm_pipes%mm_uid, * ) n_mm
      do iat = 1, n_mm
        write( qmmm_pipes%mm_uid, * ) r_mm(:,iat), pc(iat)
      enddo
      call flush( qmmm_pipes%mm_uid )
    endif

    deallocate( dpbuffer )
  end subroutine send_coordinates

  subroutine recv_forces( n_qm, n_mm, f_qm, f_mm, stress, energy )
    !! Receives QMMM forces and stress from SIESTA.
    use f90sockets, only : readbuffer
    use precision , only : dp

    implicit none
    integer , intent(in)  :: n_qm
      !! Total number of atoms in the QM region.
    integer , intent(in)  :: n_mm
      !! Total number of atoms in the MM region.
    real(dp), intent(out) :: f_qm(3,1:n_qm)
      !! Forces over QM atoms.
    real(dp), intent(out) :: f_mm(3,1:n_mm)
      !! Forces over MM atoms.
    real(dp), intent(out) :: stress(3,3)
      !! Forces over MM atoms.
    real(dp), intent(out) :: energy
      !! QM/MM energy.

    integer :: max_size, iat, n_tmp
    real(dp), allocatable :: dpbuffer(:)

    max_size = 3
    max_size = max( max_size, n_qm + n_mm )
    allocate( dpbuffer(3*max_size) )
    dpbuffer = 0.0_dp

    if ( using_sockets ) then
      call readbuffer( qmmm_socket%id, energy )
      call readbuffer( qmmm_socket%id, dpbuffer, 9 )
      stress = reshape( dpbuffer, (/3,3/) )

      call readbuffer( qmmm_socket%id, n_tmp ) ! Optional check for QM atoms
      call readbuffer( qmmm_socket%id, dpbuffer, 3*n_qm )
      f_qm = reshape( dpbuffer, (/3,n_qm/) )

      call readbuffer( qmmm_socket%id, n_tmp ) ! Optional check for MM atoms
      call readbuffer( qmmm_socket%id, dpbuffer, 3*n_mm )
      f_mm = reshape( dpbuffer, (/3,n_mm/) )

    else
      read( qmmm_pipes%fr_uid, * ) energy
      do iat = 1, 3
        read( qmmm_pipes%fr_uid, * ) stress(:,iat)
      enddo

      read( qmmm_pipes%fr_uid, * ) n_tmp ! Optional check for QM atoms
      do iat = 1, n_qm
        read( qmmm_pipes%fr_uid, * ) f_qm(:,iat)
      enddo

      read( qmmm_pipes%fr_uid, * ) n_tmp ! Optional check for MM atoms
      do iat = 1, n_mm
        read( qmmm_pipes%fr_uid, * ) f_mm(:,iat)
      enddo

    endif
    deallocate( dpbuffer )
  end subroutine recv_forces
end module siestaqmmm_m
