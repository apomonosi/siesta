!
! Copyright (C) 1996-2021	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module QMMM_interface
  !! This module handles information passage from the MM driver to SIESTA
  !! and vice-versa. Essentially:
  !!   - Siesta receives coordinates for QM atoms, coordinates for partial
  !!     charges and the values of said charges.
  !!   - After SCF, SIESTA sends forces over QM atoms and partial charges to
  !!     the MD driver (master) code.
  !!
  !! Possible interfaces include:
  !!   1 - pipes
  !!   2 - sockets

  implicit none
  public :: QMMM_communication_setup
  public :: QMMM_communication_close
  public :: QMMM_get_qm_coords
  public :: QMMM_get_partialcharges
  public :: QMMM_send_forces

  private
  type QMMM_comm_t
    integer            :: ctyp = 1
      !! Communication type for MM program-SIESTA interactions.
      !! 1 for pipes, 2 for sockets.
    character(len=224) :: driver = 'unknown'
      !! Placeholder in case specifying the MM driver program is needed.

    ! These are used for pipe or file communication.
    character(len=255) :: qmfile = 'siesta.qmatoms'
      !! File or pipe containing the QM atom coordinates.
    integer            :: qm_uid = 0
      !! Unit for the QM atoms file/pipe.
    character(len=255) :: pcfile = 'siesta.mmatoms'
      !! File or pipe containing the partial charges and their coordinates.
    integer            :: pc_uid = 0
      !! Unit for the partial charges file/pipe.
    character(len=255) :: frfile = 'siesta.allforces'
      !! File or pipe where siesta outputs forces over QM atoms and partial
      !! charges.
    integer            :: fr_uid = 0
      !! Unit for the forces file/pipe.

    ! The following options are only used in the socket implementation.
    integer             :: socket = 0
      !! Socket identificator for socket communications.
    character(len=1024) :: shost  = 'localhost'
      !! Socket host address.
    integer             :: sport  = 10002
      !! Socket port.
    integer             :: stype  = 1
      !! Socket type, usually 'inet' (=1) but might be 'unix' (=0).
  end type QMMM_comm_t

  type(QMMM_comm_t) :: QMMM_comm
    !! Structure containing communication information for QMMM.

  logical :: is_last_step = .false.
    !! Whether we are in the last step.
contains

  subroutine QMMM_communication_setup( )
    !! Sets up all the data needed for QM-MM communicators.
    use parallel  , only : IONode
    use f90sockets, only : open_socket

    implicit none
    logical :: pipe_exists

    ! We first read the inputs from FDF.
    call QMMM_communication_inputs( )

    select case ( QMMM_comm%ctyp )
    case (1)
      if ( IONode ) then
        ! We wait for the driver to create the pipes.
        pipe_exists = .false.
        do while (.not. pipe_exists)
          inquire( file = QMMM_comm%qmfile, exist = pipe_exists )
        enddo

        ! Open pipe for QM coordinates and cell vectors.
        call io_assign( QMMM_comm%qm_uid )
        open( unit = QMMM_comm%qm_uid, file = QMMM_comm%qmfile, &
              form = 'formatted', status = 'old', position = 'asis' )

        pipe_exists = .false.
        do while (.not. pipe_exists)
          inquire( file = QMMM_comm%pcfile, exist = pipe_exists )
        enddo

        ! Open pipe for MM coordinates and charges.
        call io_assign( QMMM_comm%pc_uid )
        open( unit = QMMM_comm%pc_uid, file = QMMM_comm%pcfile, &
              form = 'formatted', status = 'old', position = 'asis' )

        pipe_exists = .false.
        do while (.not. pipe_exists)
          inquire( file = QMMM_comm%frfile, exist = pipe_exists )
        enddo

        ! Open pipe for forces and stress output.
        call io_assign( QMMM_comm%fr_uid )
        open( unit = QMMM_comm%fr_uid, file = QMMM_comm%frfile, &
              form = 'formatted', status = 'old', position = 'asis' )
      endif

    case (2)
      if ( IONode ) then
        write(*,'(/,a,i4,i8,2x,a)') &
          ' QMMM - opening socket: inet, port, host = ', &
          QMMM_comm%stype, QMMM_comm%sport, trim( QMMM_comm%shost )

        call open_socket( QMMM_comm%socket, QMMM_comm%stype, &
                          QMMM_comm%sport , QMMM_comm%shost )
      endif

    case default
      call die( 'ERROR: QMMM_communication_setup '//&
                '- wrong type of communication.' )

    end select

  end subroutine QMMM_communication_setup

  subroutine QMMM_communication_close( )
    !! Closes all sockets or pipes used for QMMM communication.
    use parallel  , only : IONode
    use f90sockets, only : close_socket

    implicit none

    select case ( QMMM_comm%ctyp )
    case (1)
      if ( IONode ) then
        write(*,'(/,a,i4,i8,2x,a)') ' QMMM - closing pipes.'
        close( QMMM_comm%qm_uid )
        close( QMMM_comm%pc_uid )
        close( QMMM_comm%fr_uid )

      endif

    case (2)
      if ( IONode ) then
        write(*,'(/,a,i4,i8,2x,a)') ' QMMM - closing socket.'
        call close_socket( QMMM_comm%socket )
      endif

    case default
      call die( 'ERROR: QMMM_communication_setup '//&
                '- wrong type of communication.' )

    end select

  end subroutine QMMM_communication_close

  subroutine QMMM_communication_inputs( )
    !! Reads input variables and stores them in QMMM_comm
    use fdf, only : fdf_get, leqi
    use sys, only : die

    implicit none
    character(len=224) :: inp_str

    QMMM_comm%driver = fdf_get( 'QMMM.Driver'     , QMMM_comm%driver )
    inp_str          = fdf_get( 'QMMM.Driver.Type', 'socket'   )

    if ( leqi( inp_str, 'socket' ) .or. leqi( inp_str, 'sockets' ) ) then
      QMMM_comm%ctyp = 2
    else if  ( leqi( inp_str, 'pipe' ) .or. leqi( inp_str, 'pipes' ) ) then
      QMMM_comm%ctyp = 1
    endif

    select case ( QMMM_comm%ctyp )
    case (1)
      QMMM_comm%qmfile = fdf_get( 'QMMM.Driver.QMRegionFile', QMMM_comm%qmfile )
      QMMM_comm%pcfile = fdf_get( 'QMMM.Driver.MMChargeFile', QMMM_comm%pcfile )
      QMMM_comm%frfile = fdf_get( 'QMMM.Driver.ForceOutFile', QMMM_comm%frfile )

    case (2)
      QMMM_comm%shost = fdf_get( 'QMMM.Driver.Address'   , QMMM_comm%shost )
      QMMM_comm%sport = fdf_get( 'QMMM.Driver.Port'      , QMMM_comm%sport )
      inp_str         = fdf_get( 'QMMM.Driver.SocketType', 'inet' )

      if ( leqi( inp_str, 'unix' ) ) then
        QMMM_comm%stype = 0
      else if ( leqi( inp_str, 'inet' ) ) then
        QMMM_comm%stype = 1
      else
        call die( 'ERROR: QMMM_communication_inputs - unrecognized socket type')
      endif

    case default
      call die( 'ERROR: QMMM_communication_inputs - wrong QMMM.DriverType.' )

    end select

  end subroutine QMMM_communication_inputs

  subroutine QMMM_get_qm_coords( n_qm, qm_crd, ucell )
    !! Gets the coordinates for QM atoms in a QMMM simulation.
    use f90sockets, only : readbuffer
    use fdf       , only : leqi
    use parallel  , only : IONode
    use precision , only : dp
    use sys       , only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Double_Precision, MPI_Logical, MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast
#endif

    implicit none
    integer , intent(in)    :: n_qm
      !! Total number of QM atoms.
    real(dp), intent(inout) :: qm_crd(3,n_qm)
      !! QM atom coordinates.
    real(dp), intent(inout) :: ucell(3,3)
      !! Simulation box.

    character(len=4)      :: last_buffer
    integer               :: nat_in
    real(dp)              :: cell_in(9)
    real(dp), allocatable :: crd_in(:)
#ifdef MPI
    integer :: MPIerror
#endif

    select case( QMMM_comm%ctyp )
    case (1) ! Pipes
      if ( IONode ) then
        read( QMMM_comm%qm_uid, * ) ucell

        read( QMMM_comm%qm_uid, * ) nat_in
        if ( nat_in /= n_qm ) &
          call die( 'ERROR: QMMM_get_qm_coords - unexpected number of atoms.' )
        read( QMMM_comm%qm_uid, * ) qm_crd

        read( QMMM_comm%qm_uid, * ) last_buffer
        if ( leqi(last_buffer, 'last')) is_last_step = .true.
      endif

#ifdef MPI
      call MPI_Bcast( ucell(1,1)  ,      9, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
      call MPI_Bcast( qm_crd(1,1) , 3*n_qm, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
      call MPI_Bcast( is_last_step,      1, MPI_Logical         , 0, &
                      MPI_Comm_World, MPIerror )
#endif

    case (2) ! Sockets
      if ( IONode ) call readbuffer( QMMM_comm%socket, cell_in, 9 )
#ifdef MPI
      call MPI_Bcast( cell_in, 9, MPI_Double_Precision, 0, MPI_Comm_World, &
                      MPIerror )
#endif
      ucell = reshape( cell_in, (/3,3/) )

      ! Read and check number of atoms
      if ( IONode ) then
        call readbuffer( QMMM_comm%socket, nat_in )
        if ( nat_in /= n_qm ) &
          call die( 'ERROR: QMMM_get_qm_coords - unexpected number of atoms.' )
      endif

      ! Read and broadcast atomic coordinates
      allocate( crd_in(3*n_qm) )
      if ( IONode ) call readbuffer( QMMM_comm%socket, crd_in     , 3*n_qm )
      if ( IONode ) call readbuffer( QMMM_comm%socket, last_buffer,      4 )

      if ( leqi(last_buffer, 'last')) is_last_step = .true.
#ifdef MPI
      call MPI_Bcast( crd_in,  3*n_qm, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
      call MPI_Bcast( is_last_step, 1, MPI_Logical, 0, &
                      MPI_Comm_World, MPIerror )
#endif
      qm_crd = reshape( crd_in, (/3, n_qm/) )
      deallocate( crd_in )

    case default
      call die( 'ERROR: QMMM_get_qm_coords - wrong type of communication.' )

    end select

    ! Print for verification.
    if ( IONode ) then
      write(*,'(A)') 'QM/MM cell and QM coordinates received (Bohr).'
      write(*,'(A)') 'Cell (Bohr) ='
      write(*,'(3f12.6)') ucell
      write(*,'(A)') 'QM Coordinates (Bohr) ='
      write(*,'((3f12.6))') qm_crd
    end if
  end subroutine QMMM_get_qm_coords

  subroutine QMMM_get_partialcharges( n_pc, pc_crd, pc_chrg )
    !! Gets the coordinates and charge values for classical partial
    !! charges in a QMMM simulation. Pointers to pc_crd and pc_chrg
    !! should be already deallocated.
    use alloc     , only : re_alloc
    use f90sockets, only : readbuffer
    use parallel  , only : IONode
    use precision , only : dp
    use sys       , only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Integer, MPI_Double_Precision, MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast
#endif

    implicit none
    integer          , intent(out) :: n_pc
      !! Total number of classical partial charges.
    real(dp), pointer, intent(out) :: pc_crd(:,:)
      !! Partial charges' coordinates.
    real(dp), pointer, intent(out) :: pc_chrg(:)
      !! Partial charges' values.

    integer :: icrg
    real(dp), allocatable :: pc_data(:)
#ifdef MPI
    integer :: MPIerror
#endif

    select case( QMMM_comm%ctyp )
    case (1) ! Pipes
      if ( IONode ) then
        read( QMMM_comm%pc_uid, * ) n_pc

        if ( n_pc < 0 ) call die( ' ERROR: QMMM_get_partialcharges - number '//&
                                  'of partial charges cannot be less than 0.' )
      endif

#ifdef MPI
      call MPI_Bcast( n_pc, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror )
#endif
      nullify( pc_crd, pc_chrg )
      call re_alloc( pc_crd, 1, 3, 1, n_pc, 'pc_crd', &
                     'QMMM_get_partialcharges' )
      call re_alloc( pc_chrg,      1, n_pc, 'pc_chrg', &
                     'QMMM_get_partialcharges' )

      if ( n_pc < 1 ) return

      if ( IONode ) then
        do icrg = 1, n_pc
          read( QMMM_comm%pc_uid, * ) pc_crd(1:3,icrg), pc_chrg(icrg)
        enddo
      endif

#ifdef MPI
      call MPI_Bcast( pc_crd(1,1), 3*n_pc, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
      call MPI_Bcast( pc_chrg(1) ,   n_pc, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
#endif

    case (2) ! Sockets
      if ( IONode ) then
        call readbuffer( QMMM_comm%socket, n_pc )

        if ( n_pc < 0 ) call die( ' ERROR: QMMM_get_partialcharges - number '//&
                                  'of partial charges cannot be less than 0.' )
      endif

#ifdef MPI
      call MPI_Bcast( n_pc, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror )
#endif
      nullify( pc_crd, pc_chrg )
      call re_alloc( pc_crd, 1, 3, 1, n_pc, 'pc_crd', &
                     'QMMM_get_partialcharges' )
      call re_alloc( pc_chrg,      1, n_pc, 'pc_chrg', &
                     'QMMM_get_partialcharges' )

      if ( n_pc < 1 ) return

      if ( IONode ) then
        allocate( pc_data(3*n_pc) )
        call readbuffer( QMMM_comm%socket, pc_data, 3*n_pc)

        pc_crd = reshape( pc_data, (/3, n_pc/) )

        call readbuffer( QMMM_comm%socket, pc_data, n_pc)

        pc_chrg(1:n_pc) = pc_data(1:n_pc)

        deallocate( pc_data )
      endif

#ifdef MPI
      call MPI_Bcast( pc_crd(1,1), 3*n_pc, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
      call MPI_Bcast( pc_chrg(1) ,   n_pc, MPI_Double_Precision, 0, &
                      MPI_Comm_World, MPIerror )
#endif

    case default
      call die( 'ERROR: QMMM_get_partialcharges - '//&
                'wrong type of communication.' )

    end select

  end subroutine QMMM_get_partialcharges

  subroutine QMMM_send_forces( n_qm, n_pc, qm_frc, pc_frc, Etot, stress, &
                               maxstep )
    !! Sends forces over QM atoms and over classical partial charges
    !! to the MM program, in addition to cell stress and total energy.
    use f90sockets, only : writebuffer
    use parallel  , only : IONode
    use precision , only : dp

    implicit none
    integer , intent(in)    :: n_qm
      !! Total number of QM atoms.
    integer , intent(in)    :: n_pc
      !! Total number of classical partial charges.
    real(dp), intent(in)    :: qm_frc(3,n_qm)
      !! Forces over QM atoms.
    real(dp), intent(in)    :: pc_frc(3,n_pc)
      !! Forces over classical partial charges.
    real(dp), intent(in)    :: Etot
      !! Total energy.
    real(dp), intent(in)    :: stress(3,3)
      !! Cell stress.
    integer , intent(inout) :: maxstep
      !! Maximum possible number of MD steps. Used to end the simulation on
      !! SIESTA's side.

    integer :: iat
    real(dp), allocatable :: str_v(:), qmf_v(:), pcf_v(:)

    ! Here we guarantee that this will be the last step in SIESTA too.
    if ( is_last_step ) maxstep = 0

    select case( QMMM_comm%ctyp )
    case (1) ! Pipes
      if ( IONode ) then
        write( QMMM_comm%fr_uid, * ) Etot

        do iat = 1,3
          write( QMMM_comm%fr_uid, * ) stress(:,iat)
        end do

        write( QMMM_comm%fr_uid, * ) n_qm
        do iat = 1, n_qm
          write( QMMM_comm%fr_uid, * ) qm_frc(:,iat)
        end do

        if ( n_pc > 0 ) then
          write( QMMM_comm%fr_uid, * ) n_pc
          do iat = 1, n_pc
            write( QMMM_comm%fr_uid, * ) pc_frc(:,iat)
          enddo
        endif

        call pxfflush( QMMM_comm%fr_uid )
      endif

    case (2) ! Sockets
      if ( IONode ) then
        allocate( str_v(9), qmf_v(3*n_qm) )

        qmf_v(:) = reshape( qm_frc, (/3*n_qm/) )
        str_v(:) = reshape( stress, (/9/) )

        call writebuffer( QMMM_comm%socket, Etot )
        call writebuffer( QMMM_comm%socket, str_v, 9 )
        call writebuffer( QMMM_comm%socket, n_qm )
        call writebuffer( QMMM_comm%socket, qmf_v, 3*n_qm )

        deallocate( str_v, qmf_v )

        if ( n_pc > 0 ) then
          allocate( pcf_v(3*n_pc) )
          pcf_v = reshape( pc_frc, (/3*n_pc/) )

          call writebuffer( QMMM_comm%socket, n_pc )
          call writebuffer( QMMM_comm%socket, pcf_v, 3*n_pc )

          deallocate( pcf_v )
        endif
      endif

    case default
      call die( 'ERROR: QMMM_send_forces - wrong type of communication.' )

    end select
  end subroutine QMMM_send_forces

end module QMMM_interface
