program tscontour

  use units, only: Kelvin, eV
  use files, only: slabel, stdin_file
  use precision, only: dp
  use siesta_options, only : Temp
  use fdf, only : fdf_init, fdf_get, leqi

  use m_ts_options, only: read_ts_chem_pot, IsVolt, Volt, ts_kT, N_Elec, Elecs, N_mu, mus
  use m_ts_contour, only: read_contour_options, io_contour, ts_contour_reset
  use m_ts_chem_pot, only: delete
  use ts_electrode_m, only: electrode_t
  use m_ts_global_vars, only: TSmode, onlyS

  implicit none

  character(len=64) :: arg = ''
  integer :: narg, iarg, i
  logical :: exists = .false.

  narg = command_argument_count()
  iarg = 1
  do while( iarg <= narg )
     arg = ' '
     call get_command_argument(iarg,arg)
     select case ( arg )
     case ( '-h', '--help', '-help' )
       call help
     case default
       if ( arg(1:1) == '-' ) then
         write(0,'(a)') 'Either of the two errors has been encountered for the option "'//trim(arg)//'":'

         write(0,'(a)') ' 1) The option is not recognised'
         write(0,'(a)') ' 2) Input fdf cannot start with a hyphen "-"'
         call nl(0)
         call help
       end if
       stdin_file = arg
     end select
     iarg = iarg + 1
  end do

  if ( leqi(stdin_file, 'none') ) then
    write(0,'(a)') 'Could not find input file on the command line'
    call nl(0)
    call help
  end if

  ! check whether the file exists
  inquire(file=stdin_file, exist=exists)
  if (.not. exists ) then
    write(0,'(a)') 'Input file can not be piped into this program, &
        &please supply FDF file on command line...'
    call nl(0)
    call help
  end if

  ! Initialize the fdf
  call fdf_init(stdin_file, "tscontour.log")

  ! Fake transiesta mode (ensures ts-read methods does not crash)
  TSmode = .true.
  onlyS = .false.

  ! Retrieve default label (for consistent output)
  slabel = fdf_get('SystemLabel', 'siesta')

  ! Read (default) electronic temperature, to be used in read_ts_chem_pot
  Temp = fdf_get('ElectronicTemperature', 1.9e-3_dp, 'Ry')

  ! Read in chemical potentials (and electronic temperature)
  ! Reads:
  !   ts_kT
  !   Volt
  !   IsVolt
  call read_ts_chem_pot()

  ! Instead of reading in electrodes we allocate and *fake* what is needed
  N_elec = N_mu
  allocate(Elecs(N_elec))
  Elecs(:)%Eta = fdf_get('TS.Elecs.Eta',0.001_dp*eV,'Ry')

  do i = 1 , N_elec
    Elecs(i)%name = mus(i)%name
    Elecs(i)%ID = i
    Elecs(i)%mu => mus(i)

    mus(i)%N_el = 1
    allocate(mus(i)%el(1))
    mus(i)%el(1) = i
  end do

  ! Read in contour options
  call read_contour_options(N_Elec, Elecs, N_mu, mus, ts_kT, IsVolt, Volt)

  ! Write out everything
  call io_contour(IsVolt, mus, slabel)

  do i = 1 , N_elec
    call delete(mus(i))
    call Elecs(i)%delete(all=.true.)
  end do

  deallocate(Elecs, mus)

  call ts_contour_reset()

contains

  subroutine nl(u)
    integer, intent(in), optional :: u
    if ( present(u) ) then
      write(u,*)
    else
      write(*,*)
    end if
  end subroutine nl

  subroutine help()
    write(0,'(a)') 'Helps writing contour as output'
    write(0,'(a)') 'Options:'
    write(0,'(a)') ' <fdf>      : the input fdf file that needs conversion.'
    write(0,'(a)') ' -h         : this help.'
    stop
  end subroutine help

end program tscontour
