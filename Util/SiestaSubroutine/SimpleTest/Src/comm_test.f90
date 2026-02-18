program main
  !
  ! Test case in which not all ranks are involved in the
  ! Siesta launch.
  ! Contributed by Miha Gunde, 2025
  !
  use mpi
  use fsiesta

  implicit none
  integer, parameter :: rp=kind(1.d0)
  integer :: ierr, world_me, group_me, world_np, group_np, group_comm, mygroup
  character(:), allocatable :: fname
  real(rp), dimension(3,3) :: x
  real(rp) :: etot


  call mpi_init(ierr)

  call mpi_comm_size( MPI_COMM_WORLD, world_np, ierr )
  call mpi_comm_rank( MPI_COMM_WORLD, world_me, ierr )

  ! filename with input
  fname="test_input"

  ! atomic positions
  x(:,1) = [ 0.000_rp,  0.000_rp,   0.00_rp ]
  x(:,2) = [ 0.757_rp,  0.586_rp,   0.00_rp ]
  x(:,3) = [ -0.757_rp,  0.586_rp,   0.00_rp ]

  ! split the comm such that me=0 is in group=0, everyone else in group=1
  if( world_me==0 ) then
     mygroup=0
  else
     mygroup=1
  end if
  call mpi_comm_split( MPI_COMM_WORLD, mygroup, world_me, group_comm, ierr )

  ! get size and ranks of the new groups
  call mpi_comm_size( group_comm, group_np, ierr)
  call mpi_comm_rank( group_comm, group_me, ierr )

  write(*,"(5(a,1x,i0,2x))")&
       "world_size:",world_np, &
       "world_me:",world_me, &
       "mygroup:",mygroup, &
       "group_size:",group_np, &
       "group_me:",group_me


  ! launch siesta with group=1 comm
  call siesta_units( "Ang", "eV" )
  if( mygroup==1 ) then
     call siesta_launch( trim(fname), comm=group_comm )
     call siesta_forces( trim(fname), 3, x, energy=etot )
     if( group_me==0 ) write(*,*) "etot:",etot
  end if


  call mpi_comm_free(group_comm, ierr)
  call mpi_finalize(ierr)
end program main
