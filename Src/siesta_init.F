! ---
! Copyright (C) 1996-2021	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
#include "mpi_macros.f"

      MODULE m_siesta_init
      private
      public :: siesta_init

      CONTAINS

      subroutine siesta_init( skip_cmdargs )
      use kpoint_scf_m,       only: setup_kpoint_scf, kpoint_scf
      use kpoint_scf_m,       only: gamma_scf
      use kpoint_pdos_m,      only: gamma_pdos
      use kpoint_ldos_m,      only: gamma_ldos
      use Band,               only: gamma_bands, setup_bands
      use m_ksvinit,          only: gamma_polarization,
     &                              estimate_pol_kpoints
      use m_struct_init,      only: struct_init
      use units,              only: Kelvin
      use siesta_options
      use timer_options,      only: use_tree_timer, use_parallel_timer
      use timer_options,      only: time_mpi_calls
      use timer_options,      only: use_walltime
      use timer_options,      only: use_inbalance
      use sparse_matrices,    only: dm_history
      use class_Fstack_Pair_Geometry_dSpData2D, only: new, max_size
      use siesta_geom
      use atomlist,           only: no_u, rmaxkb, amass, lasto, qtot,
     &                              qtots,
     &                              iza, rmaxo, zvaltot, superc,
     &                              initatomlists, no_l, rmaxdftu
      use fdf
      use sys,                only: die, bye
      use gridxc,             only: gridxc_init
      use molecularmechanics, only: inittwobody
      use metaforce,          only: initmeta
      use m_mpi_utils,        only: broadcast
      use alloc,              only: re_alloc
      use memory_log,         only: memory_report
      use parallelsubs,       only: getnodeorbs
      use m_iostruct,         only: write_struct, read_struct
      use zmatrix,            only: lUseZmatrix
      use zmatrix,            only: write_canonical_ucell_and_Zmatrix
      use m_supercell,        only: exact_sc_ag
      use files,              only: slabel
      use siesta_cmlsubs,     only: siesta_cml_init
      use m_timestamp,        only: timestamp
      use m_wallclock,        only: wallclock
      use parallel,           only: Node, Nodes, PEXSINodes, IOnode
      use parallel,           only: SIESTA_worker
      use parallel,           only: SIESTA_Group, SIESTA_Comm
      use m_energies
      use m_steps
      use m_spin,             only: spin
      use spin_subs_m,        only: init_spin, init_spiral
      use m_rmaxh
      use m_forces
      use m_eo
      use m_fixed,            only: init_fixed, print_fixed
      use m_ioxv,             only: xv_file_read
      use m_local_DOS,    only: init_local_DOS
      use m_projected_DOS,    only: init_projected_DOS
      use writewave,          only: gamma_wavefunctions,
     &                              setup_wf_kpoints
      use m_new_dm,           only: get_allowed_history_depth
      use m_object_debug,     only: set_object_debug_level
      use siesta_cml
      use gridXC,  only : getXC=>gridxc_getXC   ! Returns the XC functional to be used

#ifdef DEBUG_XC
      USE siestaXC, only: setDebugOutputUnit
#endif /* DEBUG_XC */
      USE m_timer, only: timer_report ! Sets options for report of CPU times
      use m_check_walltime
#ifdef MPI
      use mpi_siesta
#endif
      use m_cite, only: add_citation, init_citation

      use m_ts_init, only : ts_init

      use m_diag_option, only: read_diag, print_diag
      use m_diag_option, only: diag_recognizes_neigwanted
#ifdef SIESTA__DFTD3
      use grimme_dispersion_m, only : dftd3_initialize
#endif
#ifdef SIESTA__FLOOK
      use siesta_dicts, only : dict_populate
      use flook_siesta, only : slua_init, slua_call, LUA_INITIALIZE
#endif

#ifdef SIESTA__WANNIER90
      use m_w90_wrapper,   only: setup_w90_wrapper
#endif

#ifdef NCDF_4
      use netcdf_ncdf
#endif

      use reinit_m, only: reinit, parse_command_line
      use runinfo_m, only: runinfo
      use QMMM_core     , only: QMMM_initialize, doing_QMMM
      use QMMM_interface, only: QMMM_get_QM_coords
      use QMMM_structure, only: mm_charges
      use siesta_version_info, only: siesta_print_version

      implicit none
      logical, intent(in), optional :: skip_cmdargs

#ifdef MPI
      logical :: initialized
      integer :: MPIerror  ! Return error code in MPI routines
      integer :: MPI_compare_result
#endif
      real(dp):: veclen       ! Length of a unit-cell vector
      integer :: i, is, ispin, n_dm_items
      integer :: neigmin      ! Min. number of eigenstates (per k point)
      integer :: ns           ! Number of species
      logical :: not_using_auxcell
      real(dp) :: walltime_m, walltime_s
      integer :: nstates_spin(2)  ! dummy intermediate variable
      MPI_GROUP_TYPE :: World_Group
      integer :: npSIESTA
      external :: read_options
      ! external :: reset_messages_file
      external :: automatic_cell
      real(dp), external :: volcel

      integer :: nf, nXCfunc
      integer, parameter :: MaxFunc = 10
      character(len=20) :: XCauth(MaxFunc)
      character(len=20) :: XCfunc(MaxFunc)
      real(dp) :: XCweightX(MaxFunc)
      real(dp) :: XCweightC(MaxFunc)
      character(len=2) :: idx_string
      logical :: skip

!----------------------------------------------------------------------BEGIN
! Initialise MPI unless siesta is running as a subroutine
! of a driver program that may have initialized MPI already

      skip = .false.
      if (present(skip_cmdargs)) then
         skip=skip_cmdargs
      endif

#ifdef MPI
C     Initialise MPI and set processor number
      call MPI_Initialized( initialized, MPIerror )
      if (.not.initialized) then
#ifdef _OPENMP
        call MPI_Init_Thread(MPI_Thread_Serialized, i, MPIerror)
        if ( MPI_Thread_Serialized /= i ) then
          ! the requested threading level cannot be asserted
          ! Notify the user
          write(*,'(a)') '!!! Could not assert serialized threads'
        end if
#else
        call MPI_Init( MPIerror )
#endif
#ifdef _TRACE_
        call Extrae_init( )
        call Extrae_Restart( )
#endif
      endif ! (.not.initialized)

      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
#endif

      ! As soon as each node knows its rank,
      ! deal with info command options
      if (.not. skip) then
         call parse_command_line(info=.true.)
      endif

#ifdef MPI
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )

      ! Keeper for PEXSI-runs..
      PEXSINodes = Nodes
#endif

      !!! now done in next block...  call reset_messages_file()

      if (Node==0) then
         ! Make sure we delete this sentinel file
         open(unit=1,file="0_NORMAL_EXIT",status="unknown")
         close(unit=1,status="delete")
         ! reset_messages_file
         open(unit=1,file="MESSAGES",status="replace",
     $        position="rewind",action="write")
         close(unit=1)
      endif

      if (.not. skip) then
         ! Initialize the output file
         ! This will close and open unit=6 in case a command-line
         ! option is specified
         call init_output(Node == 0)
      endif
#ifdef DEBUG_XC
! Set output file unit for debug info
      call setDebugOutputUnit()
#endif /* DEBUG_XC */

      ! We have to immediately write out the options used
      ! to compile and run SIESTA
      if ( Node == 0 ) then

         ! Print version and other compilation information.
         call siesta_print_version

         ! Print parallelisation runtime information.
         call runinfo

         ! Write time-stamp at the top of the output
         call timestamp('Start of run')
      end if

! Initialise input
      call reinit( sname, skip )

#ifdef MPI

#ifdef SIESTA__PEXSI
      ! Optionally, restrict the number of processors that SIESTA uses, out
      ! of the total pool (but only for PEXSI, not yet for other solvers)

      npSIESTA  = fdf_get("MPI.Nprocs.SIESTA",Nodes)
      call MPI_Comm_Group(MPI_Comm_World, World_Group, MPIerror)
      call MPI_Group_incl(World_Group, npSIESTA,
     $                    (/ (i,i=0,npSIESTA-1) /),
     $                    SIESTA_Group, MPIerror)
      call MPI_Comm_create(MPI_Comm_World, SIESTA_Group,
     $                     SIESTA_Comm, MPIerror)

      SIESTA_worker = (Node < npSIESTA)
#else
      SIESTA_Comm = MPI_Comm_World
      SIESTA_worker = .true.
#endif
      ! Swap communicators
      ! This is needed in case that the whole program is run with
      ! a subset of the global processes. In that case, the DFT instance
      ! communicator (a better name than 'true world' communicator)
      ! is at this point MPI_COMM_WORLD
      MPI_COMM_DFT = MPI_COMM_WORLD
      ! This is needed since MPI_COMM_WORLD is used implicitly by
      ! Siesta for all operations.
      MPI_COMM_WORLD = SIESTA_Comm
      if (SIESTA_worker) then
         call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
         call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
      endif
#else /* not MPI */
      SIESTA_worker = .true.
      Node = 0
      Nodes = 1
      PEXSINodes = 1
#endif
      IOnode = SIESTA_worker .and. (Node .eq. 0)

      if ( IONode ) then
        ! Add citations
        call init_citation(trim(slabel))
        call add_citation("10.1088/0953-8984/14/11/302")
      end if

      ! Be sure to initialize the spin-configuration data
      ! for all processors.
      call init_spin( )

      if (.not. SIESTA_worker) RETURN

#ifdef DEBUG
!     Generates a debug file for every process to track the execution.
!     The file is called debug.$(PID) Where PID is the process number.
!     It also works in secuencial mode.
      call debugMpiOn( )
#endif
#ifdef NCDF_4
      call ncdf_IONode( IONode )
#endif

      if ( IOnode ) then
#ifdef MPI
         if ( PEXSINodes /= Nodes ) then
            write(6,'(/,a,2(i0,a))')
     &           '* Running on ', Nodes, ' SIESTA-nodes and ',
     &           PEXSINodes, ' PEXSI-nodes in parallel'
         end if
#endif

         call wallclock('Start of run')
      end if

!     Initialize some variables
      call init_Energies()


!     Set object debugging level
      call set_object_debug_level(0)
      i = fdf_get('DebugObjects.Node',0)
      if (Node==i) then
         if (fdf_get('DebugObjects',.false.)) then
            call set_object_debug_level(1)
         endif
      endif

!     Initialize CML (relies on reinit)
      call siesta_cml_init( )

! Set timer report file and threshold .................................
      call timer_report( file=trim(slabel)//'.times',
     .        threshold=fdf_get('timer_report_threshold', 0._dp))
      ! Note that the parallel timer might be inefficient for reports
      ! when large numbers of processors are used
      if ( PEXSINodes /= Nodes ) then
         use_parallel_timer = fdf_get('UseParallelTimer', .false.)
         use_tree_timer = fdf_get('UseTreeTimer', .true.)
      else
         use_parallel_timer = fdf_get('UseParallelTimer', .true.)
         use_tree_timer = fdf_get('UseTreeTimer', .false.)

      end if
      ! Enable tree-timer for this case
      if (fdf_get("timing-split-scf-steps",.false.)) then
         use_tree_timer = fdf_get('UseTreeTimer', .true.)
      endif
      use_walltime = fdf_get('use-walltime', use_tree_timer)
      ! Whether to print min/ave/max timings across MPI processes
      ! in certain sections
      use_inbalance = fdf_get('timer-use-inbalance', .false.)

!     Start time counter
!     Note new placement of this first use, so that
!     initialization is done after the setup of relevant variables

      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )

! Set allocation report level .........................................
      call memory_report( level=fdf_get('alloc_report_level',0),
     $           file=trim(slabel)//'.alloc',
     .           threshold=fdf_get('alloc_report_threshold', 0._dp),
     &           printNow=.false. )


!     Initialise exchange-correlation functional information
      call read_xc_info()
#ifdef MPI
      call gridxc_init(comm=MPI_COMM_ID(MPI_Comm_World))
#else
      call gridxc_init()
#endif

!     Write exchange-correlation functional to .cml file
      call getXC( nXCfunc, XCfunc, XCauth, XCweightX, XCweightC )

      if ( cml_p ) then
        call cmlStartPropertyList(xf=mainXML, title='XC.Mix',
     .                              dictRef='siesta:xc')
        do nf = 1,nXCfunc
          write(idx_string, "(I0.2)") nf
          call cmlStartPropertyList(mainXML, id=idx_string)
          call cmlAddProperty(xf=mainXML, dictRef='siesta:xc_func',
     .                    title='Functional', value=trim(XCfunc(nf)))
          call cmlAddProperty(xf=mainXML, dictRef='siesta:xc_auth',
     .                    title='Authors', value=trim(XCauth(nf)))
          call cmlAddProperty(xf=mainXML, dictRef='siesta:weight_ex',
     .                    title='Weight(Ex)', value=XCweightX(nf))
          call cmlAddProperty(xf=mainXML, dictRef='siesta:weight_ec',
     .                    title='Weight(Ec)', value=XCweightC(nf))
          call cmlEndPropertyList(xf=mainXML)
        enddo
        call cmlEndPropertyList(xf=mainXML)
      endif

!     Initialise force field component
      call inittwobody( )


!     Initialize pseudopotentials and atomic orbitals
      if (IOnode) call initatom( ns )
      call broadcast( ns )

      call broadcast_basis( )

      call register_rfs()

      atmonly = fdf_get( 'Atom-Setup-Only', .false. )
      if (atmonly) call bye( 'End of atom setup' )

!     Read geometry
      call struct_init( )      ! Sets na_u, isa, ucell

      ! Initialize the spin-spiral settings
      call init_spiral( ucell )

!     Initialize atom lists
      call initatomlists( )    ! Sets iza

!     early exit if only checking the structure
      struct_only = fdf_get( 'Output-Structure-Only', .false. )
      if (IONode) then
         call write_struct( ucell, na_u, isa, iza, xa )
         if (fdf_boolean('WriteCoorXmol',.false.)) then
            call coxmol(iza, xa, na_u)
         endif
         if (lUseZmatrix) then
            call write_canonical_ucell_and_Zmatrix(
     &                        filename="OUT.UCELL.ZMATRIX.INITIAL")
         endif
      endif
      if (struct_only) then
         call bye('End of structure processing')
      endif

!     Walltime control (numbers in seconds)
      if ( fdf_isphysical("MaxWalltime") ) then
         walltime_m = fdf_get("MaxWalltime", huge(1._dp), 's')
      else
         walltime_m = fdf_get("MaxWalltime", huge(1._dp))
      end if
      ! Period for clean-up operations
      if ( fdf_isphysical("MaxWalltime.Slack") ) then
         walltime_s = fdf_get("MaxWalltime.Slack", 5._dp, 's')
      else
         walltime_s = fdf_get("MaxWalltime.Slack", 5._dp)
      end if
      ! Note that the slack detracts from net available time
      walltime_warning = walltime_m - walltime_s

!-------------- Now we have the initial geometry

!     End of Initial Structure Processing
      if (Node.eq.0) then
        write(6,'(/,a,20("*"),a,28("*"))')
     &    'siesta: ', ' Simulation parameters '
        write(6,'(a)')  'siesta:'
        write(6,'(a)')  'siesta: The following are some of the '//
     &                           'parameters of the simulation.'
        write(6,'(a)')  'siesta: A complete list of the parameters '//
     &                           'used, including default values,'
        write(6,'(a,a)')'siesta: can be found in file out.fdf'
        write(6,'(a)')  'siesta:'
      endif

!     Allocate other arrays based on read sizes
      nullify(fa,cfa)
      call re_alloc( fa, 1, 3, 1, na_u, 'fa', 'siesta_init' )
      call re_alloc( cfa, 1, 3, 1, na_u, 'cfa', 'siesta_init' )

!     Read simulation data
      call read_options( na_u, ns, spin%H )

#ifdef SIESTA__DFTD3
      if (want_dftd3_dispersion) then
         !     Initialise D3 corrections module.
         if (IOnode) then
            call dftd3_initialize( )
         endif
      endif
#endif

      call get_allowed_history_depth(n_dm_items)
      call new(DM_history,n_dm_items,"(DM history stack)")
      if (ionode) print "(a,i0)", "Size of DM history Fstack: ",
     $                            max_size(DM_history)

      qtot = qtot - charnet     ! qtot set in initatomlists
                                ! charnet set in redata
      if (IOnode) then
        write(6,fmt="(a,f12.6)") 'Total number of electrons: ', qtot
        write(6,fmt="(a,f12.6)") 'Total ionic charge: ', zvaltot
      endif
!
!     Warn the user: if not doing a direct optimization, the Zmatrix
!     coordinates are no longer updated. Only coordinates are treated.
!
      if (lUseZmatrix) then
         if (idyn .ne. 0) then
            write(6,"(a)")
     &      " *** WARNING: Zmatrix form will be used only for input !!"
            write(0,"(a)")
     &      " *** WARNING: Zmatrix form will be used only for input !!"
         endif
      endif

! Calculate spin populations for fixed spin case...
      if (fixspin) then
        if (.not.spin%Col)
     &    call die( 'siesta: ERROR: ' //
     &              'You can only fix the spin of the system' //
     &              ' for collinear spin polarized calculations.' )
        qtots(1) = (qtot + total_spin) / 2.0_dp
        qtots(2) = (qtot - total_spin) / 2.0_dp
        if (IOnode) then
          write(6,"(a,f12.6)") 'Total number of UP electrons: ',
     &        qtots(1)
          write(6,"(a,f12.6)") 'Total number of DOWN electrons: ',
     &        qtots(2)
        end if
      else
        ! Distribute electrons equally per spin
        ! (For spin%spinor = 1, qtots(2) should not be referenced)
        qtots(:) = qtot / spin%spinor
      end if

      ! Get number of eigenstates that need to be calculated This
      ! generally needs to be set to half the number of charges plus a
      ! few extra states (such that the temperature tail can be found).
      ! When states are spinors (for collinear or non-collinear spin,
      ! i.e., for spin%spinor=2), to be occupied by a single electron
      ! each, this number refers to half the number of spinors.
      neigwanted = fdf_get('NumberOfEigenStates',no_u)

      ! This is the rough number of states to be occupied (per spin channel)
      ! For spin%spinor=1, this is half the number of electrons
      nstates_spin(1:2) = ceiling ( (qtots(1:2)/2.0_dp) *  spin%spinor )

      ! Calculate the number of eigenstates
      ! dependent on the # of spin-components.
      neigmin = 0
      do is = 1 , spin%spinor
        neigmin = max(neigmin, nstates_spin(is))
      end do

      ! Always at least have one additional state per 200 K.  This
      ! choice is rather arbitrary. If the system is large it is likely
      ! to have many states "close" to the Fermi level, and then the
      ! number of extra states should be higher. See below
      ! Note that we need at least one more state in order for the
      ! Fermi-level-determination routines to work

      neigmin = neigmin + max(1, ceiling((temp / Kelvin) / 200._dp))

      ! Correct the number of calculated eigenstates
      if ( neigwanted < 0 ) then

         ! If this number is negative it corresponds to the number of
         ! states ABOVE the charge neutrality point to be included.  For
         ! normal temperatures and small/medium systems, a number such
         ! as 10 might be enough, but for very high temperatures and/or
         ! large systems, this number might need to be increased significantly.

         ! Add the user requested states.
         neigwanted = neigmin + abs(neigwanted)

      end if

      ! Check number of eigenstates - cannot be larger than number of
      ! basis functions or smaller than number of occupied states + 1
      ! so that the Fermi level can be estimated
      neigwanted = max(neigwanted, neigmin)
      ! Truncate value to number of basis orbitals
      neigwanted = min(neigwanted,no_u)


!     Find maximum interaction range
      if ( negl ) then
        rmaxh = 2.0_dp*(rmaxo + rmaxdftu)
      else
        rmaxh = 2.0_dp*(rmaxo + max(rmaxkb,rmaxdftu))
      endif

!     Madelung correction for charged systems
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif

!     Parallel initialisation
      call initparallel( no_u, na_u, lasto, xa, ucell, rmaxh )
      if (IOnode) call show_distribution( )

!     Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

!     Find k-grid for Brillouin zone integration
!     NOTE: We need to know whether gamma is .true. or
!     not early, in order to decide whether to use an
!     auxiliary supercell for the calculation of matrix elements.
      call setup_kpoint_scf( ucell )
      not_using_auxcell = gamma_scf
#ifdef SIESTA__ELSI
      if (isolve == SOLVE_ELSI) then
         if (mod(Nodes,kpoint_scf%N*spin%spinor) /=0) then
            call die("The number of MPI processes " //
     $           "must be a multiple of nk*nspin")
         endif
      endif
#endif
!     Call initialisation of PDOS here since we need to check if
!     the auxiliary supercell is needed for a non-gamma calculation
      call init_projected_DOS( ucell )
      if ( do_pdos ) then
        not_using_auxcell = not_using_auxcell.and. gamma_pdos
      endif

!     Call initialisation of LDOS here since we need to check if
!     the auxiliary supercell is needed for a non-gamma calculation
      call init_local_DOS( ucell )
      if ( do_ldos ) then
        not_using_auxcell = not_using_auxcell .and. gamma_ldos
      endif

      ! Read in diagonalization options
      ! Note that only the sampled BZ determines the value of
      ! ParallelOverK. Oother routines performing diagonalization
      ! (i.e., pdos, bands...) might reset this option.
      call read_diag(Gamma_scf, spin%H)

      ! Do not print option values now unless we are using diagonalization
      ! as the ground-state solver
      if (isolve == SOLVE_DIAGON) call print_diag()

      if ( diag_recognizes_neigwanted() .or. isolve == SOLVE_ELSI ) then
         if ( IONode .and. neigwanted /= no_u )
     &      write(6,"(a,t53,'= ',i0,' / ',i0)")
     &      'diag: Number of calculated states [no_calc / no_u]',
     &      neigwanted, no_u
      else
         if ( IONode .and. isolve == SOLVE_DIAGON .and.
     $        neigwanted /= no_u) then
            write(6,"(a)")
     &      'diag: The diagonalization method chosen ' //
     $      'cannot compute a reduced number ' //
     $      'of eigenvectors.'
            write(6,"(a,i0,a,i0)")
     &           'diag: Resetting neigwanted ', neigwanted, ' to ', no_u
         endif
        neigwanted = no_u
      end if

      nullify(eo,qo)
      call re_alloc(eo, 1, no_u, 1, spin%spinor, 1, kpoint_scf%N, 'eo',
     &              'siesta_init')
      call re_alloc(qo, 1, no_u, 1, spin%spinor, 1, kpoint_scf%N, 'qo',
     &              'siesta_init')

      call setup_bands( )
      not_using_auxcell = not_using_auxcell .and. gamma_bands

      call setup_wf_kpoints( )
      not_using_auxcell = not_using_auxcell .and. gamma_wavefunctions

      call estimate_pol_kpoints( ucell )
      not_using_auxcell = not_using_auxcell .and. gamma_polarization
!
!     User can request that the calculation is done with an explicit
!     auxiliary supercell and hermitian version, even if using only
!     the gamma point
!
!     In case the user has requested the use of auxiliary
!     cell we simply use that option.
      if ( .not. not_using_auxcell ) use_aux_cell = .true.

!     Find required supercell
!     2*rmaxh is used to guarantee that two given orbitals in the
!     supercell can only overlap once
      if ( .not. use_aux_cell ) then
        nsc(1:3) = 1
      else
        ! Determine the tightest auxiliary supercell using
        ! also the atomic positions
        call exact_sc_ag(negl,ucell,na_u,isa,xa,nsc)
      end if

      mscell = 0.0_dp
      do i = 1, 3
        mscell(i,i) = nsc(i)
        nsc_old(i)  = nsc(i)
      enddo

!     Find auxiliary supercell (required only for k sampling)
      call superc( ucell, scell, nsc )

!     Initialise metadynamic forces if required
      call initmeta( )

      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else if (idyn .ge. 1 .and. idyn .le. 5) then
        inicoor = istart
        fincoor = ifinal
      else if (idyn .eq. 6) then
        inicoor = 0
        fincoor = (ia2-ia1+1)*3*2
      else if (idyn .eq. 7) then
        call die( "'PHONON' support is deprecated" )
      else if ( (idyn == 8) .or. (idyn == 11) ) then
        inicoor = 0
        fincoor = huge(1)
#ifdef NCDF_4
      else if (idyn == 9) then
         if ( IONode ) then
            write(*,'(/,a)')'expcoord: '//repeat('*',69)
            write(*,'(a,i0)')
     &           'expcoord: Currently reached coordinate step = ',
     &           inicoor
            write(*,'(a,i0)')
     &           'expcoord: Final coordinate step = ',
     &           fincoor
            write(*,'(a,/)')'expcoord: '//repeat('*',69)
         end if
#endif
#ifdef SIESTA__FLOOK
      else if (idyn == 10) then
         inicoor = 0
         ! Controlled in external LUA machine
         fincoor = huge(1)
#endif
      else
        call die( 'siesta: wrong idyn' )
      endif

      ! initialize the fixed positions (we need to check whether the
      ! fixed positions make sense)
      call init_fixed( ucell, na_u , isa, iza )
      call print_fixed( )

      ! Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. xv_file_read))
     &    call vmb( na_u, tempinit, amass, xa, isa, va )

      istp = 0

      call ts_init(spin%Grid,ucell,na_u,xa,lasto,no_u,inicoor,fincoor)

      ! Initialize QMMM module.
      call QMMM_initialize( idyn == 11 )
      if ( idyn == 11 ) then
        ! If we are using an external driver, we need to reset
        ! coordinates.
        call QMMM_get_QM_coords( na_u, xa, ucell )
        if ( volcel(ucell) < 1.0e-8_dp )
     &     call automatic_cell( ucell, scell, na_u, xa, isa, charnet )
      endif
      ! Updates partial charge distribution for QM/MM (if needed)
      if ( doing_QMMM() ) call mm_charges%update( na_u, xa, ucell )

#ifdef SIESTA__FLOOK

      ! Populate dictionaries
      call dict_populate()

      ! Initialize LUA handle, also die if the Lua relaxation
      ! is requested
      call slua_init(LUA, idyn == 10)

      ! Call a preprocess step right after initialization
      call slua_call(LUA,LUA_INITIALIZE)

#endif

!     Initialize WANNIER90 to be run directly in SIESTA
#ifdef SIESTA__WANNIER90
      if( w90_wrapper_processing ) then
!
!       Set up all the variables required to run WANNIER90 within SIESTA
!       (not as a postprocessing tool)
!       At this step, these variables are only known at the IONode
!
        call setup_w90_wrapper()
      endif

#endif

!     Output memory use before main loop
!!      call printmemory( 6, 0 )

!     Initialization now complete. Flush stdout.
      if (ionode) call pxfflush( 6 )

      call timer( 'Setup', 2 )

      END subroutine siesta_init

!--------------------------------------------------------------------------- END

      END MODULE m_siesta_init
