! ---
! Copyright (C) 1996-2021       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module runinfo_m
  implicit none

  private
  public :: runinfo

contains

  subroutine runinfo

     ! This subroutine prints runtime information.
     ! This includes information about the cwd, MPI, and OpenMP.

!$   use omp_lib, only : omp_get_num_threads
!$   use omp_lib, only : omp_get_nested, omp_set_nested
!$   use omp_lib, only : omp_get_max_active_levels
!$   use omp_lib, only : omp_get_schedule, omp_set_schedule
!$   use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
!$   use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
#if defined(_OPENMP) && ( _OPENMP >= 201307 )
!$   ! The above value should correspond to OpenMP 4.0.
!$   use omp_lib, only : omp_get_proc_bind
!$   use omp_lib, only : OMP_PROC_BIND_FALSE, OMP_PROC_BIND_TRUE
!$   use omp_lib, only : OMP_PROC_BIND_MASTER
!$   use omp_lib, only : OMP_PROC_BIND_CLOSE, OMP_PROC_BIND_SPREAD
#endif
     use parallel, only: Nodes
     use posix_calls, only: getcwd

     implicit none

     character(len=2048) :: cwd
     integer :: istat
     integer :: i, is
!$   integer :: ns    ! Number of nested levels

     write(6,'(/,a)') 'Runtime information:'

     ! Current working directory
     call getcwd(cwd, istat)
     if ( istat == 0 ) then
        write(6,'(2a)') '* Directory : ', trim(adjustl(cwd))
     end if

#ifdef MPI
     ! Sadly PEXSI workers can only be determined after
     ! fdf has been opened, so we can't print-out that information
     ! (for now)
     if ( Nodes > 1 ) then
        write(6,'(a,i0,a)') '* Running on ', Nodes, ' nodes in parallel.'
     else
        write(6,'(a)') '* Running in serial mode (only 1 MPI rank).'
     endif
#else
     write(6,'(a)') '* Running in serial mode.'
#endif

!$OMP parallel default(shared)
!$OMP master

     ! Number of threads and processes.
!$   i = omp_get_num_threads()
!$   write(*,'(a,i0,a)') '* Running ',i,' OpenMP threads.'
!$   write(*,'(a,i0,a)') '* Running ',Nodes*i,' processes.'

     ! Thread affinity information.
#if _OPENMP >= 201307
!$   i = omp_get_proc_bind()
!$   select case ( i )
!$     case ( OMP_PROC_BIND_FALSE )
!$       write(*,'(a)') '* OpenMP NOT bound (please bind threads!)'
!$     case ( OMP_PROC_BIND_TRUE )
!$       write(*,'(a)') '* OpenMP bound'
!$     case ( OMP_PROC_BIND_MASTER )
!$       write(*,'(a)') '* OpenMP bound (master)'
!$     case ( OMP_PROC_BIND_CLOSE )
!$       write(*,'(a)') '* OpenMP bound (close)'
!$     case ( OMP_PROC_BIND_SPREAD )
!$       write(*,'(a)') '* OpenMP bound (spread)'
!$     case default
!$       write(*,'(a)') '* OpenMP bound (unknown)'
!$   end select
#endif

     ! Scheduling information.
!$   call omp_get_schedule(i,is)
!$   select case ( i )
!$      case ( OMP_SCHED_STATIC )
!$         write(*,'(a,i0)') '* OpenMP runtime schedule STATIC, chunks ', is
!$      case ( OMP_SCHED_DYNAMIC )
!$         write(*,'(a,i0)') '* OpenMP runtime schedule DYNAMIC, chunks ', is
!$         if ( is == 1 ) then
!$           ! this is the default scheduling, probably the user
!$           ! has not set the value, predefine it to 32
!$           is = 32
!$         write(*,'(a,i0)') '** Changing chunk size:'
!$         write(*,'(a,i0)') '** OpenMP runtime schedule DYNAMIC, chunks ',is
!$      end if
!$      case ( OMP_SCHED_GUIDED )
!$         write(*,'(a,i0)') '* OpenMP runtime schedule GUIDED, chunks ',is
!$      case ( OMP_SCHED_AUTO )
!$         write(*,'(a,i0)') '* OpenMP runtime schedule AUTO, chunks ',is
!$      case default
!$         write(*,'(a,i0)') '* OpenMP runtime schedule UNKNOWN, chunks ',is
!$   end select

     ! Nested levels.
!$   if ( omp_get_nested() ) then
!$      ns = omp_get_max_active_levels()
!$   else
!$      write(*,'(a)') '** OpenMP (trying to FORCE nesting)'
!$      call omp_set_nested(.true.)
!$      ns = omp_get_max_active_levels()
!$   end if
!$   write(*,'(a,i0,a)') '* OpenMP allows ',ns,' nested levels.'

!$OMP end master
!$OMP end parallel
!$   call omp_set_schedule(i,is)

  end subroutine runinfo

end module runinfo_m

