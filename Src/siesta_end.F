! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end( )
      use atm_types,       only : atm_types_reset
      use fdf,             only : fdf_shutdown
      use memory_log,      only : memory_report
      use precision,       only : dp
      use siesta_cml,      only : cml_p, mainXML, cmlEndModule
      use siesta_cmlsubs,  only : siesta_cml_exit
      use siesta_master,   only : siesta_subroutine  ! Is siesta a subroutine?
#ifdef DEBUG_XC
      USE siestaXC,        only : closeDebugOutputFile
#endif /* DEBUG_XC */
      use m_timestamp,     only : timestamp
      use m_wallclock,     only : wallclock, wallclock_shutdown
      use parallel,        only : IOnode, SIESTA_worker
#ifdef MPI
      use mpi_siesta
#endif
      use meshsubs,        only : meshsubs_reset
      use moreMeshSubs,    only : resetMeshDistr
      use densematrix,     only : resetDenseMatrix
      use sparse_matrices, only : resetSparseMatrices,
     &                            resetOrbitalDistributions
      use m_fixed,         only : resetFixedPointers
      use band,            only : reset_bands
      use meshphi,         only : resetMeshPhi
      use meshdscf,        only : resetDscfPointers
      use m_dscfcomm,      only : resetdscfComm
      use m_rhog,          only : resetRhog
      ! NP, new mixing
      use m_mixing_scf,    only : mixers_scf_reset
#ifdef SIESTA__CHESS
      USE siesta_options,  only : isolve, SOLVE_CHESS
#endif

      use kpoint_scf_m, only: reset_kpoint_scf
      use kpoint_dos_m, only: reset_kpoint_dos
      use kpoint_pdos_m, only: reset_kpoint_pdos
      use kpoint_ldos_m, only: reset_kpoint_ldos
      use ts_kpoint_scf_m, only: reset_ts_kpoint_scf

      use m_ts_method, only: ts_method_reset
      use m_ts_contour_eq, only: ts_contour_eq_reset
      use m_ts_contour_neq, only: ts_contour_neq_reset
      use m_ts_options, only: ts_options_reset
      use m_ts_sparse, only: ts_sparse_reset
      use zmatrix    , only: zmatrix_reset
      use spinorbit  , only: spinorbit_reset

#ifdef SIESTA__CHESS
      use m_chess, only: CheSS_finalize
#endif

#ifdef SIESTA__FLOOK
      use flook_siesta,    only : slua_call, LUA_FINALIZE
      use flook_siesta,    only : slua_close
      use siesta_options,  only : LUA
      use siesta_dicts,  only : dict_clean
#endif

      use m_diag, only: diag_exit

#ifdef TRACING
      use extrae_eventllist, only : writeList, deleteList, eventlist
#endif
      use m_io,            only : io_assign, io_close
      use m_io_yaml,       only : siesta_write_yaml
      use m_cite,          only : announce_citations, reset_citations
      use siesta_version_info,    only : siesta_version_str
      use atomlist,        only : reset_atomlists
      use matel_m,         only : cleanup_matel_tables
      use m_matel_registry,only : cleanup_matel_registry
      use listsc_module,   only : LISTSC_RESET
      use m_ksvinit,       only : ksv_reset
      use m_dhscf,         only : dhscf_reset
      use basis_specs,     only : basis_specs_reset
      use basis_types,     only : basis_parameters_reset
      use iogrid_netcdf,   only : iogrid_netcdf_reset
      use ncdf_io_m,       only : ncdfio_reset
      use siesta_geom,     only : siesta_geom_reset
      use m_forces,        only : forces_reset
      use m_eo,            only : eoqo_reset
      use chemical,        only : chemical_reset
      use m_fft,           only : fft_trigs_reset
      use gpfa_fft,        only : fft_gpfa_trigs_reset
      use m_dynamics,      only : md_reset
      use writewave,       only : reset_writewave
      use molecularmechanics, only: reset_twobody
      use m_dftu,          only : reset_dftu
      use m_occ_proj,      only : reset_proj_dftu
      use m_dminim  ,      only : reset_dminim_c
      use m_grdsam  ,      only : reset_grdsam
#ifdef SIESTA__DFTD3
      use grimme_dispersion_m, only : dftd3_finalize
      use siesta_options, only : want_dftd3_dispersion
#endif

      use QMMM_core     , only : QMMM_end
      use QMMM_structure, only : mm_charges

      implicit none

      real(dp) :: R12(3), S12, DSDR(3)
      integer :: iu
#ifdef MPI
      integer MPIerror
#endif
!------------------------------------------------------------------------- BEGIN
      if (SIESTA_worker) then

      if (ionode) then
        call siesta_write_yaml()
        call announce_citations(siesta_version_str)
        call reset_citations()
      end if

#ifdef SIESTA__FLOOK
      ! Final call for the end-user
      call slua_call(LUA,LUA_FINALIZE)
      ! Clean up the lua environment
      call slua_close(LUA)
      call dict_clean()
#endif

!     Free allocated memory
      ! Reset matrix element tables and registry
      call cleanup_matel_tables()
      call cleanup_matel_registry()
      R12 = 0._dp
      call LISTSC_reset()
      call reset_atomlists()
      call resetSparseMatrices( )
      call resetMeshDistr( )
      call resetdscfComm( )
      call resetDenseMatrix( )
      call resetFixedPointers( )
      call reset_bands( )
      call resetMeshPhi( )
      call resetDscfPointers( )
      call resetRhog()
      call mixers_scf_reset()
      call ksv_reset()
      call dhscf_reset()
      call spinorbit_reset()
      call resetOrbitalDistributions( )
      call reset_writewave()
      call reset_twobody()
      call reset_dftu()
      call reset_proj_dftu()
      call reset_dminim_c()
      call reset_grdsam()

      call reset_kpoint_ldos()
      call reset_kpoint_pdos()
      ! This *must* be called after LDOS+PDOS since
      ! they may be associated with DOS
      call reset_kpoint_dos()
      ! This *must* be called before SCF since
      ! they may be associated with SCF
      call reset_ts_kpoint_scf()
      call reset_kpoint_scf()

      ! Clean transiesta stuff
      call ts_contour_eq_reset()
      call ts_contour_neq_reset()
      call ts_sparse_reset()
      call ts_method_reset()
      call ts_options_reset()

      ! Clean diagonalization
      call diag_exit()
      call eoqo_reset()

      ! Clean species and atomic data
      call zmatrix_reset()
      call atm_types_reset()
      call basis_specs_reset()
      call basis_parameters_reset()
      call siesta_geom_reset()
      call forces_reset()
      call chemical_reset()

      ! Others
      call meshsubs_reset()
      call fft_trigs_reset()
      call fft_gpfa_trigs_reset()
      call md_reset()
      call iogrid_netcdf_reset()
      call ncdfio_reset()


#ifdef SIESTA__DFTD3
      if (IOnode .and. want_dftd3_dispersion) 
     &    call dftd3_finalize() 
#endif

#ifdef SIESTA__CHESS
      if (isolve == SOLVE_CHESS) then
          call CheSS_finalize()
      end if
#endif
      call mm_charges%clean( )
      call QMMM_end( )

!     Output memory use up to the end of the program
!!      call printmemory( 6, 1 )

!     Print allocation report
      call memory_report( printNow=.true. , shutdown=.true.)

#ifdef DEBUG_XC
      call closeDebugOutputFile()
#endif /* DEBUG_XC */

!     Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

#ifdef TRACING
!     Write extrae user function list
      if (IOnode) then
        call writeList(eventlist)
      end if
      call deleteList(eventlist)
#endif /* TRACING */


!     Print final date and time
      if (IOnode) then
        write(*,*) ! Newline
        call timestamp('End of run')
        call wallclock('End of run')
        ! Ensure clock file is closed
        call wallclock_shutdown()
      endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif
      endif  ! SIESTA_worker

! Finalize MPI, unless siesta is running as a subroutine
      !
      ! Add stamp to MESSAGES file
      call message('INFO','Job completed')

      ! Clean the fdf-files
      call fdf_shutdown()

#ifdef MPI
      call MPI_Barrier(MPI_Comm_World, MPIerror)
      MPIerror = MPI_SUCCESS
      if (.not.siesta_subroutine) call MPI_Finalize( MPIerror )
#endif

      if (IOnode) then
         call io_assign(iu)
         open(iu, file='0_NORMAL_EXIT', form='FORMATTED')
#ifdef MPI
         write(iu,'(a,i0)') 'SIESTA completed, MPI exit: ',
     &        MPIerror
#else
         write(iu,'(a)') 'SIESTA completed successfully'
#endif
         call io_close(iu)
      endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
