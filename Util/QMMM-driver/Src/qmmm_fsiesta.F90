module qmmm_fsiesta
  use precision          , only : dp

! Support routines for siesta-as-a-subroutine in Unix/Linux.
! The routines that handle the other side of the communication are
! in module iopipes of siesta program.
! Usage:
!   call siesta_launch( label, nnodes )
!     character(len=*),intent(in) :: label  : Name of siesta process
!                                             (prefix of its .fdf file)
!     integer,optional,intent(in) :: nnodes : Number of MPI nodes
!
!   call siesta_units( length, energy )
!     character(len=*),intent(in) :: length : Physical unit of length
!     character(len=*),intent(in) :: energy : Physical unit of energy
!
!   call siesta_forces( label, na, xa, cell, energy, fa, stress )
!     character(len=*), intent(in) :: label      : Name of siesta process
!     integer,          intent(in) :: na         : Number of atoms
!     real(dp),         intent(in) :: xa(3,na)   : Cartesian coords
!     real(dp),optional,intent(in) :: cell(3,3)  : Unit cell vectors
!     real(dp),optional,intent(out):: energy     : Total energy
!     real(dp),optional,intent(out):: fa(3,na)   : Atomic forces
!     real(dp),optional,intent(out):: stress(3,3): Stress tensor
!   call siesta_quit( label )
!     character(len=*),intent(in) :: label  : Name of one siesta process,
!                                             or 'all' to stop all procs.
! Behaviour:
! - If nnodes is not present among siesta_launch arguments, or nnodes<2,
!   a serial siesta process will be launched. Otherwise, a parallel
!   mpirun process will be launched.
! - If siesta_units is not called, length='Ang', energy='eV' are
!   used by default. If it is called more than once, the units in the
!   last call become in effect.
! - The physical units set by siesta_units are used for all the siesta
!   processes launched
! - If siesta_forces is called without a previous call to siesta_launch
!   for that label, it assumes that the siesta process has been launched
!   (and the communication pipes created) externally in the shell.
!   In this case, siesta_forces only opens its end of the pipes and begins
!   communication through them.
! - If argument cell is not present in the call to siesta_forces, or if
!   the cell has zero volume, it is assumed that the system is a molecule,
!   and a supercell is generated automatically by siesta so that the
!   different images do not overlap. In this case the stress returned
!   has no physical meaning.
! - The stress is defined as dE/d(strain)/Volume, with a positive sign
!   when the system tends to contract (negative pressure)
! - The following events result in a stopping error message:
!   - siesta_launch is called twice with the same label
!   - siesta_forces finds a communication error trough the pipes
!   - siesta_quit is called without a prior call to siesta_launch or
!     siesta_forces for that label
! - If siesta_quit is not called for a launched siesta process, that
!   process will stay listening indefinitedly to the pipe and will need
!   to be killed in the shell.
! - siesta_units may be called either before or after siesta_launch
! J.M.Soler and A.Garcia. Nov.2003

   use fdf

  implicit none

PUBLIC :: siesta_launch, siesta_units, siesta_send_coordinates, &
               siesta_receive_forces, siesta_quit, create_fifos

PRIVATE ! Nothing is declared public beyond this point

! Holds data on siesta processes and their communication pipes
  type :: proc
    private
    character(len=80) :: label ! Name of process
    integer           :: iuc   ! I/O units for coords/forces commun.
    integer           :: iuf   ! I/O units for coords/forces commun.
    integer           :: iuq   ! I/O units for coords/forces commun.
  end type proc

! Global module variables
  integer, parameter :: max_procs = 100
  type(proc),   save :: p(max_procs)
  integer,      save :: np=0
  character(len=32), save :: xunit = 'Ang'
  character(len=32), save :: eunit = 'eV'
  character(len=32), save :: funit = 'eV/Ang'
  character(len=32), save :: sunit = 'eV/Ang**3'

CONTAINS

  subroutine create_fifos( label )
    !! Creates IO pipes.
    use linkatoms, only : num_resonances, resonance

    implicit none
    character(len=*),  intent(in) :: label

    character(len=255) :: cpipe, fpipe, qpipe
    character(len=120) :: task
    integer            :: resonance_id

    do resonance_id = 1, num_resonances
      ! Check that pipe does not exist already
      if ( idx(trim(resonance(resonance_id)%path)//label ) /= 0) &
        print*, 'siesta_launch: ERROR: process for label ', trim(label), &
          ' already launched'

      ! Create pipes
      cpipe = trim(resonance(resonance_id)%path)//trim(label)//'.coords'
      fpipe = trim(resonance(resonance_id)%path)//trim(label)//'.forces'

      call system( 'rm -f '//cpipe )
      call system( 'rm -f '//fpipe )
      task = 'mkfifo '//trim(cpipe)//' '//trim(fpipe)
      qpipe = trim(resonance(resonance_id)%path)//trim(label)//'.pc'
      task  = trim(task)//' '//qpipe

      call system( task )

      ! Open this side of pipes
      call open_pipes( trim(resonance(resonance_id)%path)//label )
    enddo
  end subroutine create_fifos

  subroutine siesta_launch( label )
    !! Launches a SIESTA process for QM calculations.
    use linkatoms          , only : num_resonances, resonance
    use siesta_qmmm_options, only : launch_siesta_flag, parallel_command

    implicit none
    character(len=*),  intent(in) :: label

    character(len=120) :: override_command
    character(len=512) :: task
    character(len=512) :: siesta_bin, siesta_bin_rel
    integer            :: resonance_id, siesta_nodes, envstat
    logical            :: lex

    ! set processor number for the QM siesta calculation
    siesta_nodes = fdf_integer( 'NumberMPInodes', 1 )

    if ( launch_siesta_flag ) then

      call get_environment_variable( "SIESTA", value = siesta_bin, &
                                     status = envstat )

      if ( envstat > 0 ) then
        write(*,*) "SIESTA environment variable not found when attempting "//&
                   "to launch siesta."
        write(*,*) "Defaulting to ./siesta in the working directory."
        siesta_bin     = "./siesta"
        siesta_bin_rel = "../siesta"
      elseif ( envstat < -1 ) then
        write(*,*) "Error while reading SIESTA environment variable."
        write(*,*) "Defaulting to ./siesta in the working directory."
        siesta_bin     = "./siesta"
        siesta_bin_rel = "../siesta"
      else
        write(*,*) "SIESTA environment variable set to "//trim(siesta_bin)
        siesta_bin_rel = siesta_bin
      endif

      do resonance_id = 1, num_resonances
        if ( resonance_id > 1 ) then
          inquire( file = trim(resonance(resonance_id)%path)//'siesta', &
                   exist = lex )
          if (lex) &
            call system( 'rm -f '//trim(resonance(resonance_id)%path)//'siesta')

          call system( 'ln -s '//trim(siesta_bin_rel)//' '//&
                       trim(resonance(resonance_id)%path) )
        endif

        ! Start siesta process
        if ( siesta_nodes > 1 ) then
          call get_environment_variable( "SIESTA_PARALLEL_COMMAND", &
                                         value  = override_command, &
                                         status = envstat )

          if ( (envstat > 0) .or. (envstat < -1) ) then
             write( override_command, * ) &
               trim(parallel_command)//' -n ', siesta_nodes
          endif
          write(*,*) "Using '"//trim(override_command)//&
                     "' for parallel excecution."

          write( task, * ) 'cd '//trim(resonance(resonance_id)%path)//' ; '//&
                           trim(override_command)//' '//&
                           trim(siesta_bin)//' < '//trim(label)//'.fdf > '//&
                           trim(label)//'.out &'
        else
          write( task, * ) 'cd '//trim(resonance(resonance_id)%path)//' ; '//&
                           trim(siesta_bin)//' < '//trim(label)//'.fdf > '//&
                           trim(label)//'.out &'
        endif

        call system( task )
        write( * , fmt = '(A,I4,A)' ) &
          'SIESTA has been launched on ', siesta_nodes, ' procs'
        call pxfflush( 6 )
      enddo

    else
      do resonance_id = 1, num_resonances
        call system( 'echo 0 > '//trim(resonance(resonance_id)%path)//&
                     'start_siesta' )
      enddo
    endif

  end subroutine siesta_launch

  subroutine siesta_units( length, energy )
    implicit none
    character(len=*), intent(in) :: length, energy
    xunit = length
    eunit = energy
    funit = trim(eunit)//'/'//trim(xunit)
    sunit = trim(eunit)//'/'//trim(xunit)//'**3'
  end subroutine siesta_units

  subroutine siesta_send_coordinates( label, na_qm, r_all, cell, &
                                      na_mm, mm_atoms, is_last )
    use linkatoms  , only : num_resonances, resonance, numlink, link_atoms
    use mm_topology, only : mm_atom_t

    implicit none
    character(len=*), intent(in) :: label
    integer         , intent(in) :: na_qm
    integer         , intent(in) :: na_mm
    real(dp)        , intent(in) :: r_all(3,na_qm+na_mm)
    type(mm_atom_t) , intent(in) :: mm_atoms(na_mm)
    real(dp)        , intent(in) :: cell(3,3)
    logical         , intent(in) :: is_last

    integer           :: i, ia, ip, iu, resonance_id
    real(dp)          :: c(3,3)
    logical, save     :: firstTime = .true.

    do resonance_id = 1, num_resonances

      ! Find system index
      ip = idx( trim(resonance(resonance_id)%path)//label )
      if ( ip == 0 ) then
        call open_pipes( label )
        ip = idx( label )
        stop 'Error while opening pipes.'
      endif

      ! Copy unit cell
      c = cell

      ! Print coords for debugging
      print'(/,2a)', 'siesta_send_coordinates: label = ', &
          trim(resonance(resonance_id)%path)//trim(label)
      print'(3a,/,(3f12.6))', 'siesta_send_coordinates: cell (', &
          trim(xunit), ') =', cell
      print'(3a,/,(3f12.6))', 'siesta_send_coordinates: rqm (', &
          trim(xunit), ') =', r_all(:,1:na_qm)
      print'(3a,/,(3f12.6))', 'siesta_send_coordinates: rlink (',&
          trim(xunit), ') =', &
            ( link_atoms(i)%reson(resonance_id)%r(1:3), i = 1, numlink )
      print*,'##################################################'

      if ( firstTime ) then
        call siesta_launch( label )
        firstTime = .false.
      endif

      ! Write coordinates to pipe
      iu = p(ip)%iuc

      do i = 1, 3
        write(iu,*) c(:,i)
      enddo
      write(iu,*) na_qm + numlink
      do ia = 1, na_qm
        write(iu,*) r_all(:,ia)
      enddo
      do ia = 1, numlink
        write(iu,*) link_atoms(ia)%reson(resonance_id)%r(:)
      enddo

      if ( is_last ) then
        write(iu,*) 'last'
      else
        write(iu,*) 'endq'
      endif
      call pxfflush( iu )

      iu = p(ip)%iuq
      write(iu,*) na_mm
      do ia = 1, na_mm
        write(iu,*) r_all(:,ia+na_qm), mm_atoms(ia)%pc
      enddo
      call pxfflush( iu )
    enddo

  end subroutine siesta_send_coordinates

  subroutine siesta_receive_forces( label, na_qm, na_mm, energy, fa, stress )

    use alloc    , only : re_alloc, de_alloc
    use linkatoms, only : num_resonances, resonance, numlink, link_atoms

    implicit none
    character(len=*),   intent(in) :: label
    integer,            intent(in) :: na_qm
    integer,            intent(in) :: na_mm
    real(dp), optional, intent(out):: energy
    real(dp), optional, intent(out):: fa(3,na_qm+na_mm)
    real(dp), optional, intent(out):: stress(3,3)

    real(dp)          :: e, s(3,3), wgt
    real(dp), pointer :: f(:,:)
    logical, save     :: firstTime = .true.

    integer :: i, ia, ip, iu, n, resonance_id

    nullify( f )
    call re_alloc( f, 1, 3, 1, na_qm+na_mm, 'f', 'siesta_receive_forces' )

    if (present(energy)) energy = 0.0_dp
    if (present(fa))     fa     = 0.0_dp
    if (present(stress)) stress = 0.0_dp

    do resonance_id = 1, num_resonances
      ! Find system index.
      ip = idx( trim(resonance(resonance_id)%path)//label )

      ! Print current label.
      print'(/,2a)', 'siesta_receive_forces: label = ', &
                   trim(resonance(resonance_id)%path)//trim(label)
      print*,'##################################################'

      ! Read forces from pipe
      iu = p(ip)%iuf

      read( iu, * ) e
      do i = 1, 3
        read(iu,*) s(:,i)
      enddo
      read(iu,*) n

      if ( n /= na_qm+numlink ) then
        print*, 'siesta_forces: ERROR: na_u mismatch: na_u, n =', &
                na_qm+numlink, n
        stop
      endif

      do ia = 1, na_qm
        read(iu,*) f(:,ia)
      enddo
      do ia = 1, numlink
        read(iu,*) link_atoms(ia)%reson(resonance_id)%fa(1:3)
      enddo

      if ( na_mm > 0 ) then
        read(iu,*) n
        if ( n /= na_mm ) then
          print*, 'siesta_forces: ERROR: na_mm mismatch: na_mm, n =', &
                  na_mm, n
          stop
        endif
        do ia = 1, na_mm
          read(iu,*) f(:,na_qm+ia)
        enddo
      endif

      ! Print forces for debugging
      print'(3a,f12.6)', 'siesta_receive_forces: energy (',&
                          trim(eunit),') =', e
      print'(3a,/,(3f12.6))', 'siesta_receive_forces: stress (',&
                              trim(sunit),') =', s
      print'(3a,/,(3f12.6))', 'siesta_receive_forces: forces (',&
                              trim(funit),') =', f

      ! Copy results to output arguments
      wgt = resonance(resonance_id)%weight

      if ( present(energy) ) energy  = energy  + wgt * e
      if ( present(stress) ) stress  = stress  + wgt * s
      if ( present(fa)     ) fa(:,:) = fa(:,:) + wgt * f(:,:)
    enddo

    firstTime = .false.
    if ( num_resonances > 1 ) then
      ! Print forces for debugging
      print'(3a,f12.6)', 'siesta_receive_forces: total QM energy (',&
                         trim(eunit),') =', energy
      print'(3a,/,(3f12.6))', 'siesta_receive_forces: total QM stress (',&
                              trim(sunit),') =', stress
      print'(3a,/,(3f12.6))', 'siesta_receive_forces: total QM forces (',&
                              trim(funit),') =', fa(:,1:na_qm)
    endif

    call de_alloc( f, 'f', 'siesta_receive_forces' )
  end subroutine siesta_receive_forces

  recursive subroutine siesta_quit( label )
    implicit none
    character(len=*), intent(in) :: label

    integer :: ip

    if ( label == 'all' ) then
      ! Stop all siesta processes
      do ip = 1, np
        call siesta_quit( p(ip)%label )
      enddo

      return
    endif

    ! Stop one siesta process
    ip = idx( label ) ! Find process index
    if ( ip == 0 ) then
      print*, 'siesta_quit: ERROR: unknown label: ', trim(label)
      stop
    endif

    ! We just close the pipes.
    close( p(ip)%iuc, status = "delete" )
    close( p(ip)%iuf, status = "delete" )
    close( p(ip)%iuq, status = "delete" )

    ! Move last process to this slot and remove process.
    if ( ip < np ) then
      p(ip)%label = p(np)%label
      p(ip)%iuc   = p(np)%iuc
      p(ip)%iuf   = p(np)%iuf
      p(ip)%iuq   = p(np)%iuq
    endif
     np = np - 1
  end subroutine siesta_quit

  subroutine open_pipes( label )
    implicit none
    character(len=*), intent(in) :: label

    integer            :: iuc, iuf, iuq
    character(len=255) :: cpipe, fpipe, qpipe

    ! Check that pipe does not exist already
    if ( idx(label) /= 0 ) &
      print*, 'open_pipes: ERROR: pipes for ', trim(label), &
              ' already opened'

    ! Get io units for pipes
    call get_io_units( iuc, iuf, iuq )

    ! Open pipes
    cpipe = trim(label)//'.coords'
    fpipe = trim(label)//'.forces'
    qpipe = trim(label)//'.pc'
    open( unit=iuc, file=cpipe, form="formatted", &
          status="old", position="asis")
    open( unit=iuf, file=fpipe, form="formatted", &
          status="old", position="asis")
    open( unit=iuq, file=qpipe, form="formatted", &
          status="old", position="asis")

    ! Store data of process
    np = np + 1
    if ( np > max_procs ) then
      stop 'siesta_launch: ERROR: parameter max_procs too small'
    else
      p(np)%label = label
      p(np)%iuc   = iuc
      p(np)%iuf   = iuf
      p(np)%iuq   = iuq
    endif

  end subroutine open_pipes

  subroutine get_io_units( iuc, iuf, iuq )
    ! Finds three available I/O unit numbers
    implicit none
    integer, intent(out)  :: iuc, iuf, iuq

    integer :: i, ip
    logical :: unit_used

    iuc = 0
    iuf = 0
    iuq = 0
    do i = 10, 999
      !! This does not really work with pipes...
      inquire( unit = i, opened = unit_used )

      do ip = 1, np
        unit_used = ( unit_used ) .or. ( i == p(ip)%iuc ) .or. &
                    ( i == p(ip)%iuf .or. ( i == p(ip)%iuq ) )
      enddo
      if ( .not. unit_used ) then
        if ( iuc == 0 ) then
          iuc = i
        else if ( iuf == 0 ) then
          iuf = i
        else
          iuq = i
          return
        endif
      endif
    enddo

    stop 'fsiesta:get_io_units: ERROR: cannot find free I/O unit'
  end subroutine get_io_units

  integer function idx( label )
    !! Finds which of the stored labels is equal to a given input label.

    implicit none
    character(len=*), intent(in) :: label
    integer :: i
    do i = 1,np
      if ( label == p(i)%label ) then
        idx = i
        return
      end if
    end do
    idx = 0
  end function idx

end module qmmm_fsiesta
