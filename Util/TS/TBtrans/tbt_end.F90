! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication
subroutine tbt_end()

  use fdf, only: fdf_shutdown
#ifdef MPI
  use mpi_siesta, only : MPI_Finalize
#endif
  use parallel, only : Node
  use memory_log,      only : memory_report
  use m_timestamp, only : timestamp
  use m_wallclock, only : wallclock

  ! Cleaning memory
  use m_tbt_kpoint, only: tbt_kpoint_reset
  use m_tbt_tri_init, only: tbt_tri_reset
  use m_tbt_regions, only: tbt_regions_reset
  use m_tbt_hs, only: tbt_hs_reset
  use m_tbt_options, only: tbt_options_reset
#ifdef NCDF_4
  use m_tbt_delta, only: delete_delta
  use m_tbt_dH, only: dH
  use m_tbt_dSE, only: dSE
#endif

  use m_tbt_contour, only: tbt_contour_reset
  use m_ts_sparse, only: ts_sparse_reset
  use m_ts_method, only: ts_method_reset
  use m_ts_contour_eq, only: ts_contour_eq_reset
  use m_ts_contour_neq, only: ts_contour_neq_reset

#ifdef MPI
  integer :: MPIerror
#endif

  call tbt_kpoint_reset()
  call tbt_tri_reset()
  call tbt_regions_reset()
  call tbt_hs_reset()

#ifdef NCDF_4
  call delete_delta(dH)
  call delete_delta(dSE)
#endif

  call tbt_contour_reset()
  call ts_contour_eq_reset()
  call ts_contour_neq_reset()

  call tbt_options_reset()

  call ts_sparse_reset()
  call ts_method_reset()

  ! Clean the fdf-files
  call fdf_shutdown()

  ! Stop time counter
  call timer( 'tbtrans', 2 )
  call timer( 'all', 3 )

  ! Print allocation report
  call memory_report( printNow=.true. )

  if ( Node == 0 ) then
     call timestamp("End of run")
     call wallclock("End of run")
  end if

#ifdef MPI
  call MPI_Finalize(MPIerror)
#endif

end subroutine tbt_end
