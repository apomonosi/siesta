! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine amn( ispin, index_manifold )
!
!     In this subroutine we compute the overlaps between Bloch states
!     onto trial localized orbitals
!
!     These matrices will depend on three indices (see paragraph after Eq. (16) 
!     of the review by
!     N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
!     or Eq. (1.8) of the Wannier Users Guide, Version 1.2
!     $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vert g_{n} \rangle =
!    \sum_{\mu} \sum_{lcell} c_{\mu m}^{\ast} (\vec{k}) \times
!    \exp^{i \vec{k} \cdot \left( \tau_{\mu} + \vec{R}_{lcell} \right)} \times
!    \langle \phi_{\mu} (\vec{r}-\tau_{\mu} - \vec{r}_{lcell} \vert g_{n}\rangle
!     where $m$ runs between 1 and the number of bands for wannierization,
!     $n$ runs between 1 and the number of projection functions, 
!     $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
!     will be computed. 
!
!     OUTPUT: 
!     File called seedname.amn, where the overlap matrices are written in 
!     the format required by Wannier90
!   
!     BEHAVIOUR:
!     The coefficients of the wave functions are introduced
!     in amn through the matrix coeffs, computed in diagonalizeHk
!     
!     PARALLELIZATION: 
!     This subroutine is parallelized on the projectors, so each node
!     takes care on a projection.     
!
!     Implemented by J. Junquera and R. Korytar, July 2013
!

  use precision,          only: dp                  ! Real double precision type
  use parallel,           only: Nodes               ! Total number of Nodes
  use parallel,           only: Node                ! Local Node
  use parallel,           only: IONode              ! Input/output node
  use parallelsubs,       only: WhichNodeOrb        ! Which node handles the
                                                    !   the information of a 
                                                    !   given orbital in the
                                                    !   unit cell
  use parallelsubs,       only: GlobalToLocalOrb    ! Transformation from a 
                                                    !   global index of a given
                                                    !   orbital in the unit cell
                                                    !   to the local index 
                                                    !   within a node
  use atomlist,           only: rmaxo               ! Max. cutoff atomic orbital
  use siesta_geom,        only: scell               ! Lattice vector of the
                                                    !   supercell in real space
  use siesta_geom,        only: ucell               ! Unit cell lattice vectors
  use siesta_geom,        only: nsc                 ! Diagonal element of 
                                                    !   the supercell
  use siesta_geom,        only: na_s                ! Number of atoms in the
                                                    !   supercell
  use siesta_geom,        only: xa                  ! Atomic positions
  use cellsubs,           only: reclat              ! Finds reciprocal unit cell
                                                    !   vector
  use m_switch_local_projection, only: numproj      ! Total number of projectors
  use m_switch_local_projection, only: numkpoints   ! Total number of k-points
                                                    !   for which the overlap of
                                                    !   the periodic part of the
                                                    !   wavefunct with a 
                                                    !   neighbour k-point will
                                                    !   be computed
  use m_switch_local_projection, only: kpointsfrac  ! List of k points relative
                                                    !   to the reciprocal 
                                                    !   lattice vectors.
                                                    !   First  index: component
                                                    !   Second index: k-point  
                                                    !      index in the list
  use m_switch_local_projection, only: projections  ! Trial projection functions
  use m_switch_local_projection, only: numincbands  ! Number of bands for 
                                                    !   wannierization
                                                    !   after excluding bands  
  use m_switch_local_projection, only: nincbands_loc! Number of bands for 
                                                    !   wannierization
                                                    !   after excluding bands  
                                                    !   in the local Node
  use m_switch_local_projection, only: coeffs       ! Coefficients of the
                                                    !   wavefunctions.
                                                    !   First  index: orbital
                                                    !   Second index: band
                                                    !   Third  index: k-point
  use m_switch_local_projection, only: Amnmat       ! Matrix of the overlaps of 
                                                    !   trial projector funtions
                                                    !   with Eigenstates of the
                                                    !   Hamiltonian
  use trialorbitalclass,  only: trialorbital        ! Derived type to define the
                                                    !    localized trial
                                                    !    orbitals
  use m_matel_registry,   only: register_in_tf_pool ! Subroutine that assigns a
                                                    !    global index to the 
                                                    !    trial projection
                                                    !    functions
  use m_matel_registry,   only: show_pool
  use matel_m,            only: matel               ! MATEL implementation 
                                                    !   with the global indices 
                                                    !   of the radial functions 
                                                    !   as inputs
  use matel_m,            only: init_matel_wannier  ! Table initializer
  use atmfuncs,           only: orb_gindex          ! Subroutine that gives
                                                    !   the global index of an
                                                    !   atomic orbital
  use w90_wrapper_types,only: compute_chempotwann_after_scf 
                                                    ! Are we going to use the 
                                                    !   unperturbed Wannier functions
                                                    !   as trial guess to compute 
                                                    !   the matrix elements of the 
                                                    !   Wanniers with chemical potentials?
  use w90_wrapper_types,only: first_chempotwann   ! First time the calculation 
                                                    !   of the matrix elements of the 
                                                    !   Hamiltonian with the chemical
                                                    !   chemical potential of the 
                                                    !   Wanniers is called?
  use atomlist,           only: no_s                ! Number of atomic orbital
                                                    !   in the supercell
  use atomlist,           only: iaorb               ! Atomic index of each orbital
  use atmfuncs,           only: rcut                ! Returns orbital 
                                                    !   cutoff radius
  use sparse_matrices,    only: S                   ! Overlap matrix in 
                                                    !   sparse format
  use sparse_matrices,    only: xijo                ! Relative position between
                                                    !   neighbour atoms
  use sparse_matrices,    only: numh                ! Number of non-zero element
                                                    !   per row in the matrices
                                                    !   (local to MPI Node)
  use sparse_matrices,    only: listhptr            ! Pointer to the start of
                                                    !   ach row (-1) of the
                                                    !   hamiltonian matrix
  use sparse_matrices,    only: listh               ! Nonzero hamiltonian-matrix
                                                    !   element column indices 
                                                    !   for each matrix row

  use w90_wrapper_types, only: w90_manifold_t
  use w90_wrapper_types, only: manifold_bands
  use w90_wrapper_types,   only: numh_man_proj
                                               ! Number of projections that will
                                               !  be  handled in the local node
  use w90_wrapper_types,   only: listhptr_man_proj
                                       ! Index pointer to listh_man_proj such
                                       ! listh_man_proj(listhptr_man_proj(1)+1)
                                       !   is the first projector of the first
                                       !   manifold handled by the local node
                                       ! listh_man_proj(listhptr_man_proj(io)+1)                                       !   is thus the first projector of
                                       !   of manifold 'io' while
                                       ! listh_man_proj(listhptr_man_proj(io) +                                        !                numh_man_proj(io))
                                       !   is the last projectors of manifold
                                       !   'io'.
                                       ! Dimension: number of manifolds
  use w90_wrapper_types,   only: listh_man_proj
                                               ! The column indices for 
                                               !   the projectors
                                               !   of all the manifolds handled
                                               !   by the local node
  use w90_wrapper_types, only: coeffs_wan_nao
                                               ! Coefficients of the
                                               !   Wannier functions in a basis
                                               !   of NAO
                                               !   First  index: Index of the
                                               !       manifold and Wannier func
                                               !       handled by numh_man_proj,
                                               !       listhptr_man_proj, and
                                               !       listh_man_proj, and
                                               !   Second index: NAO in the
                                               !       supercell
                                               !   Third index: Spin component

!
! Variables for the diagonalization
!
  use atomlist,           only: no_u         ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: indxuo       ! Index of equivalent orbital in 
                                             !   the unit cell
!
! Allocation/Deallocation routines
!
  use alloc,              only: re_alloc     ! Reallocation routines
  use alloc,              only: de_alloc     ! Deallocation routines

  use sys,                only: die          ! Termination routine

#ifdef MPI
  use m_orderbands,       only: which_band_in_node  ! Given a node and a 
                                                    !   local index,
                                                    !   this array gives the
                                                    !   global index of the band
                                                    !   stored there
  use m_orderbands,       only: sequential_index_included_bands 
                                                    ! Sequential number of the
                                                    !   bands included for
                                                    !   wannierization
                                                    !   (the bands are listed
                                                    !   in order of incremental
                                                    !   energy)

  use mpi_siesta
#endif

  implicit none

  integer, intent(in)     :: ispin               ! Spin component
  integer, intent(in)     :: index_manifold      ! Index of the manifold
                                                 !   that is wannierized

  type orbitallinkedlist
    real(dp),dimension(3)           :: center
    integer                         :: specie
    integer                         :: specieindex
    integer                         :: globalindex
    type(orbitallinkedlist),pointer :: nextitem
  end type


! Local variables

  type(w90_manifold_t), pointer :: w90man    ! For ease of reference

  integer  :: ik              ! Counter for loop on k-points
  integer  :: iproj           ! Counter for loop on projections
  integer  :: io              ! Counter for loop on orbital 
  integer  :: iorb            ! Counter for loop on orbital 
  integer  :: iat             ! Atomic index of each orbital
  integer  :: ind_proj        ! Counter for sequential indices
                              !    of projections
  integer  :: iband           ! Counter for loop on bands
  integer  :: nincbands       ! Number of bands for wannierization
  integer  :: iua             ! Equivalent atom in the
                              !    first unit cell
  integer  :: iuo             ! Equivalent orbital in the
                              !    first unit cell
  integer  :: ix              ! Counter for cartesian directions
  integer  :: icell           ! Counter for unit cell lattice
                              !    vectors
  integer  :: gindex          ! Global index of the trial projector function
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: globalindexorbital ! Global index of the neighbour atomic orbital
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: globalindexproj ! Global index of the neighbour atomic orbital
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: indexproj       ! Index of the projector 
                              !   This index runs from 1 to the total number
                              !   of projections
  integer  :: isc(3)
  integer  :: iproj_local     ! Local index within one node of the 
                              !   Wannier function
  integer  :: iproj_global    ! Global index of the Wannier function
  real(dp) :: trialcenter(3)  ! Position where the trial function is centered
                              !   (in Bohr)
  real(dp) :: trialrcut       ! Cutoff radius of the trial function
  real(dp) :: r12(3)          ! Relative position of the trial function with
                              !   respect to the neighbour orbital
  real(dp) :: overlap         ! Overlap between the trial function and a 
                              !   given atomic orbital
  real(dp) :: gradient(3)     ! Grandient of the overlap between 
                              !   the trial function and a given atomic orbital
  real(dp) :: phase           ! Product of the k-vector with the position
                              !   where the neighbour orbital is centered
  real(dp) :: kvector(3)      ! k-point vector in the Wannier90 grid
  real(dp) :: dxa(3)          ! Cell vector that translates a
                              !   given atom in the unit cell
                              !   to the equivalent in the supercell
  real(dp) :: rcell(3,3)      ! Reciprocal unit cell vectors
                              !   (without 2*pi factor)
  real(dp) :: xatorb(3)       ! Position of the atomic orbital
 
  integer  :: jneig, juo, ind, jo
  real(dp) :: xneig(3)        ! Position of the atomic orbital


  complex(dp), dimension(:,:), pointer :: psiloc => null() ! Coefficients of the wave
                                             !   function (in complex format)
  complex(dp) :: exponential                 ! Exponential of exp( i * phase )
  complex(dp) :: cstar                       ! Conjugate of the coefficient
                                             !   of the wave function      
  complex(dp), parameter :: iu = cmplx(0.0_dp,1.0_dp,kind=dp) ! Imaginary unit

  type(trialorbital), save :: tf             ! Projetion function

  type(orbitallinkedlist),pointer   :: item  ! Variable that contains all the
                                             !   information of the atomic
                                             !   orbitals neighbours of a given
                                             !   trial projection function
!
! Variable required to globalize the neighbour list
!
  integer  :: io_local        ! Counter for the local atomic orbitals in a
                              !    given node
  integer  :: io_global       ! Global index of an atomic orbital
  integer  :: BNode           ! Node that contains the information of a given
                              !    orbital in the unit cell
  integer  :: maxnhg          ! Compute the maximum number of
                              !    interacting atomic orbitals neighbours
  integer,  pointer :: numhg(:) => null()
                       ! Temporal array used to broadcast the array numh,
                       ! with number of non-zero elements
                       ! each orbital connects to
  integer,  pointer :: listhptrg(:) => null()
                       ! Temporal array used to broadcast the array listhptr
                       ! index pointer to listh such that listh(listhptr(1) + 1)
                       ! is the first non-zero element of orbital 1.
                       ! listh(listhptr(io) + 1) is thus the
                       ! first non-zero element
                       ! of orbital 'io' while listh(listhptr(io) + numh(io))
                       ! is the last non-zero element of orbital 'io'.
  integer,  pointer :: listhg(:) => null()
                       ! Temporal array used to broadcast the array listh
                       !      the column indices for the non-zero elements
  real(dp), pointer :: xijloc(:,:) => null()
                       ! Temporal array used to broadcast the relative
                       !    position between neighbours
  real(dp), pointer :: Sloc(:) => null()
                       ! Temporal array used to broadcast the overlap matrix
#ifdef MPI
  integer     :: iband_global                ! Global index for a band
  integer     :: iband_sequential            ! Sequential index of the band
  integer     :: MPIerror
  complex(dp), dimension(:,:), pointer :: auxloc => null()! Temporal array for the
                                             !   the global reduction of Amnmat
  integer, external :: numroc
#endif 

  external :: timer

! Start time counter
  call timer( 'Amn', 1 )

  w90man => manifold_bands(index_manifold)  !------------ Shortcut
  
!! For debugging
!  write(6,'(a,l5,2i5)') 'amn: compute_chempotwann = ', compute_chempotwann_after_scf, Node, Nodes
!  write(6,'(a,l5,2i5)') 'amn: first_chempotwann = ',   first_chempotwann, Node, Nodes
!! End debugging

! Allocate memory related with the overlap matrix between the trial projection
! function and the Hamiltonian eigenstate.
! These matrices will depend on three indices (see paragraph after Eq. (16) 
! of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
! or Eq. (1.8) of the Wannier Users Guide, Version 1.2
! $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vert g_{n} \rangle =
!    \sum_{\mu} \sum_{lcell} c_{\mu m}^{\ast} (\vec{k}) \times
!    \exp^{i \vec{k} \cdot \left( \tau_{\mu} + \vec{R}_{lcell} \right)} \times
!    \langle \phi_{\mu} (\vec{r}-\tau_{\mu} - \vec{r}_{lcell} \vert g_{n}\rangle
! where $m$ runs between 1 and the number of bands considered for wannierization
! $n$ runs between 1 and the number of projection functions, 
! $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
! will be computed 

  nincbands = numincbands( ispin )

  call re_alloc( Amnmat,         &
 &               1, nincbands,   &
 &               1, numproj,     &
 &               1, numkpoints,  &
 &               'Amnmat',       &
 &               'Amn'           )

! Initialize Amnmat
  Amnmat = cmplx( 0.0_dp, 0.0_dp, kind =dp )

! Allocate memory related with a local variable where the coefficients 
! of the eigenvector at the k-point will be stored
! Only nincbands are retained for wannierization, that is why the
! second argument is made equal to nincbands
  call re_alloc( psiloc, 1, no_u, 1, nincbands, 'psiloc', 'Amn' )

! Assign a global index to the trial functions
! The same global index is assigned to a given trial function in all the nodes

  ! The projectors are registered manifold by manifold,
  ! and it is not clear whether, once the wannierization has been
  ! carried out, the previous manifolds's matel-related data are still needed.
  ! Just in case, now the wannier projectors are added to the registry cumulatively,
  ! and the wannier-related matel tables are re-generated afresh with progressively
  ! larger sizes. If later on it is decided that the previous manifold(s)' data can
  ! be removed, calls to the 'reset_wannier_tables' and
  ! 'remove_wannier_projectors_from_registry' routines can be added.
  
  if (.not. w90man%projectors_are_registered ) then 
    if (allocated(w90man%projector_gindex)) deallocate(w90man%projector_gindex)
    allocate ( w90man%projector_gindex(numproj) )
    do iproj = 1, numproj
       tf = projections( iproj )
       ! Check whether it is a "Siesta orbital", In that case, generate a
       ! new trial orbital with the right orientation, and register it properly:
       if (tf%from_basis_orbital) then
          block
            use radial, only: radial_rescale_radfunc, rad_func
            use m_matel_registry, only: peek_at_registered_radfunc
            use m_matel_registry, only: register_in_rf_pool
            
            integer :: orb_gindex, l, m
            type(rad_func), pointer :: rfunc, func_new => null()
            
            orb_gindex = tf%iorb_gindex
            call peek_at_registered_radfunc(orb_gindex,rfunc,l,m)
            if ( modulo(m,2) == 1 ) then
               allocate(func_new)
               call radial_rescale_radfunc(rfunc,func_new,scale=-1.0_dp)
               call register_in_rf_pool(func_new, l, m, "wann_twin_orb", [iproj], gindex)
            else
               call register_in_rf_pool(rfunc, l, m, "wann_orb", [iproj], gindex)
            endif
         end block
      else
         call register_in_tf_pool( tf, gindex )
      endif
      w90man%projector_gindex( iproj ) = gindex
    enddo
    if (ionode) then
       write(*, "(a,i0,a)") "Registered Wannier projectors for manifold ", &
                        index_manifold,":"
       call show_pool(first = w90man%projector_gindex(1), &
            last = w90man%projector_gindex(numproj))
    endif

    w90man%projectors_are_registered = .true.
    ! (re)initialize the matel tables for orbital-wannier matrix elements
    call init_matel_wannier( numproj )
  endif
  
! Find reciprocal unit cell vectors (without 2*pi factor)
  call reclat( ucell, rcell, 0 )

  if( compute_chempotwann_after_scf .and. (.not. first_chempotwann) ) then
!   Copy of the procedure implemented in subroutine diagkp
    nullify(numhg, listhptrg)

!   Allocate local memory for global list arrays
    call re_alloc( numhg, 1, no_u, name='numhg',                               &
 &                 routine= 'amn' )
    call re_alloc( listhptrg, 1, no_u, name='listhptrg',                       &
 &                 routine= 'amn' )

!
!   Globalise numh
!   Loop over all the orbitals in the unit cell
    do io_global = 1, no_u
!     Localize which node handles the information related with the orbital io
      call WhichNodeOrb(io_global,Nodes,BNode)
      if ( Node .eq. BNode ) then
!       Identify the local index for the orbital in the unit cell in the
!       node that handles its information
        call GlobalToLocalOrb(io_global,Node,Nodes,io_local)
!       Assign the value of the number of neighbours of that particular orbital
        numhg(io_global) = numh(io_local)
      endif
!     Transfer the information from the node that contains the information
!     on the orbital io_global to all the other nodes
#ifdef MPI
      call MPI_Bcast( numhg(io_global),1,MPI_integer,BNode,                    &
 &                    MPI_Comm_World,MPIerror )
#endif
    enddo

!   Build global listhptr
    listhptrg(1) = 0
    do io_global = 2, no_u
     listhptrg(io_global) = listhptrg(io_global-1) + numhg(io_global-1)
    enddo

!   Globalise listh
!   Compute the maximum number of interacting atomic orbitals neighbours
!   considering all the orbitals in the unit cell
    maxnhg = listhptrg(no_u) + numhg(no_u)
    nullify(listhg)
    call re_alloc( listhg, 1, maxnhg, name='listhg',                           &
 &                 routine= 'amn' )
    do io_global = 1, no_u
      call WhichNodeOrb(io_global,Nodes,BNode)
      if (Node.eq.BNode) then
        call GlobalToLocalOrb(io_global,Node,Nodes,io_local)
        do jo = 1, numhg(io_global)
          listhg(listhptrg(io_global)+1:listhptrg(io_global)+numhg(io_global))=&
 &          listh(listhptr(io_local)+1:listhptr(io_local)+numh(io_local))
        enddo
      endif
#ifdef MPI
      call MPI_Bcast( listhg(listhptrg(io_global)+1),numhg(io_global),         &
 &                    MPI_integer,BNode,MPI_Comm_World,MPIerror )
#endif
    enddo

! Globalize the relative position of neighbours
    call re_alloc( xijloc, 1, 3, 1, maxnhg, name='xijloc', routine= 'amn' ) 
    call re_alloc( Sloc, 1, maxnhg, name='Sloc', routine= 'amn' ) 
    do io_global = 1, no_u
      call WhichNodeOrb(io_global,Nodes,BNode)
      if (Node.eq.BNode) then
        call GlobalToLocalOrb(io_global,Node,Nodes,io_local)
        do jo = 1, numh(io_local)
          Sloc(listhptrg(io_global)+jo) = S(listhptr(io_local)+jo)
        enddo 
        do jo = 1, numh(io_local)
          xijloc(1:3,listhptrg(io_global)+jo) = xijo(1:3,listhptr(io_local)+jo)
        enddo 
      endif 
#ifdef MPI
      call MPI_Bcast(Sloc(listhptrg(io_global)+1),numhg(io_global),        &
 &       MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
      call MPI_Bcast(xijloc(1,listhptrg(io_global)+1),3*numhg(io_global),  &
 &       MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
#endif
    enddo 
  endif

kpoints:                 &
  do ik = 1, numkpoints
!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units, 
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )

!   Initialize the local coefficient matrix for every k-point
    psiloc(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

#ifdef MPI
!   Store the local bands in this node on a complex variable
    do iband = 1, nincbands_loc
      iband_global = which_band_in_node(Node,iband)
      iband_sequential = sequential_index_included_bands(iband_global)

      do io = 1, no_u
        psiloc(io,iband_sequential) = coeffs(io,iband,ik)
      enddo

    enddo
!   Allocate workspace array for global reduction
    call re_alloc( auxloc, 1, no_u, 1, nincbands,   &
 &                 name='auxloc', routine='Amn' )
!   Global reduction of auxloc matrix
    auxloc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    call MPI_AllReduce( psiloc(1,1), auxloc(1,1),   &
 &                      no_u*nincbands,             &
 &                      MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
!   After this reduction, all the nodes know the coefficients of the
!   wave function for the point ik, for all the bands and for all atomic
!   orbitals
    psiloc(:,:) = auxloc(:,:)
#else
    do iband = 1, nincbands
      do io = 1, no_u
        psiloc(io,iband) = coeffs(io,iband,ik)
      enddo
    enddo
#endif

!   If the calculation has arrived to self-consistency after including the
!   chemical potential for the Wanniers, and we want to compute the 
!   Hamiltonian tight-binding matrix elements,
!   we are going to take as the localized guess functions 
!   the Wannier functions obtained before including the chemical potential
    if( compute_chempotwann_after_scf .and. (.not. first_chempotwann) ) then
      do iproj_local = 1, numh_man_proj(index_manifold)
!       Identify the global index of the Wannier
        ind_proj     = listhptr_man_proj(index_manifold) + iproj_local
        iproj_global = listh_man_proj(ind_proj)
        do iorb = 1, no_s
          if( abs(real(coeffs_wan_nao(ind_proj,iorb,ispin))) .gt. &
 &            w90man%threshold ) then 
            iuo = indxuo (iorb)                ! Equivalent orbital in 
                                               !   first unit cell
            iua = iaorb(iuo)                   ! Equivalent atom in 
                                               !   first unit cell
            iat = iaorb(iorb)                  ! Atom to which orbital belongs
            dxa(:) = xa(:,iat) - xa(:,iua)     ! Cell vector of atom iat
            isc(:) = nint( matmul(dxa,rcell) ) ! Cell index of atom iat
!           Find the index of the unit cell within the supercell where this
!           atom is located, centered on isc = 0
            do ix = 1,3
               ! Same centered in isc=0
               if (isc(ix)>nsc(ix)/2) isc(ix) = isc(ix) - nsc(ix) 
            enddo
!           Find the translated position of the atom in the supercell that
!           really takes a non-zero value in the unit cell
            xatorb(:) = xa(:,iua)
            do icell = 1, 3
              do ix = 1, 3
                xatorb(ix) = xatorb(ix) + isc(icell) * ucell(ix,icell)
              enddo
            enddo

            do jneig = 1, numhg(iuo)
              ind  = listhptrg(iuo) + jneig
              jo   = listhg(ind)
              juo  = indxuo(jo)
              xneig(:) = xijloc(:,ind) + xatorb(:)
              phase = -1.0_dp * dot_product( kvector, xneig )
              exponential = exp( iu * phase )
Band_loop2:                                                             &
              do iband = 1, nincbands
                cstar = conjg( psiloc(juo,iband) )
                  Amnmat(iband,iproj_global,ik) =         & 
 &                  Amnmat(iband,iproj_global,ik)   +     &
 &                  exponential * cstar * Sloc(ind) *     &
 &                  coeffs_wan_nao(ind_proj,iorb,ispin)
              enddo Band_loop2  ! End loop on the bands
            enddo 
          endif  ! End if the coefficients are larger than a threshold
        enddo    ! End loop on orbitals in the supercell
      enddo      ! End loop on projectors
      goto 999
    endif


!   Loop on the projections that will be computed in the local node
!   It would be better if each node computes all projections for
!   all its locally stored bands, instead of some projections for
!   all the bands...  In this way we will save the globalization of
!   band data.
!
#ifdef MPI
    do iproj = 1+Node, numproj, Nodes
#else
    do iproj = 1, numproj
#endif
      indexproj = iproj


!     Find the global index of the projector in the list of radial functions
!     that will be evaluated by MATEL
      globalindexproj = w90man%projector_gindex(indexproj) 

!     Find where the trial function is centered
      trialcenter = projections(indexproj)%center
  
!     Find the cutoff radius of the trial function
      trialrcut   = projections(indexproj)%rcut

!     For debugging
!      write(6,'(a,6i5,4f12.5,2i5,f12.5)')' ik, Node, nincbands, nincbands_loc, iproj, indexproj = ',  &
! &                ik, Node, nincbands, nincbands_loc, iproj,      &
! &                indexproj, trialcenter, trialrcut, globalindexproj, na_s, rmaxo
!     End debugging

!     Find the atomic orbitals that ovelap with our radial orbital
!     centered at trialcenter and with range trialrcut
      item => get_overlapping_orbitals( scell, rmaxo, trialrcut, na_s, &
 &                                      xa, trialcenter )

OrbitalQueue:                                                        &
      do while (associated(item))
        r12 = trialcenter - item%center
        globalindexorbital = orb_gindex( item%specie, item%specieindex )
        ! Compute the overlap between orbital with globalindexorbital
        ! And projector with globalindexproj
        call matel( 'S', globalindexorbital, globalindexproj, r12, &
     &                    overlap, gradient )
!!     For debugging
!      write(6,'(a,4i5,4f12.5)')' ik, Node, globalindices, r12, overlap = ',  &
! &                ik, Node, globalindexorbital, globalindexproj,      &
! &                r12, overlap
!!     End debugging
        phase = -1.0_dp * dot_product( kvector, item%center )
        exponential = exp( iu * phase )

!       Loop over occupied bands
Band_loop:                                                           &
        do iband = 1, nincbands
          cstar = conjg( psiloc(item%globalindex,iband) )
          Amnmat(iband,indexproj,ik) =      & 
 &          Amnmat(iband,indexproj,ik) +    &
 &          exponential * cstar * overlap
        enddo Band_loop  ! End loop on the bands

        item => item%nextitem
      enddo OrbitalQueue

    enddo   ! Loop on projections on the local node

999 continue

!   We need all the matrix Amnmat in IOnode
!   (to dump it into a file),but the results for some of the bands might
!   be computed in other nodes.
#ifdef MPI
!   Allocate workspace array for reduction
    if (IONode) then
       call re_alloc( auxloc, 1, nincbands, 1, numproj,    &
            &                 name='auxloc', routine='Amn' )
       auxloc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    endif

    call MPI_Reduce( Amnmat(1,1,ik), auxloc(1,1),                       &
 &                      nincbands*numproj,                              &
 &                      MPI_double_complex,MPI_sum,0,MPI_Comm_World,MPIerror )
    if (IONode) then
       Amnmat(:,:,ik) = auxloc(:,:)
    endif
#endif
  enddo kpoints

! Write the Amn overlap matrices in a file, in the format required
! by Wannier90
  if( IONode ) call writeamn( ispin )

! Deallocate some of the arrays
  call de_alloc( psiloc,  'psiloc',  'Amn' )
#ifdef MPI
   call de_alloc( auxloc,  'auxloc',  'Amn' )
#endif
  if( compute_chempotwann_after_scf .and. (.not. first_chempotwann) ) then
    call de_alloc( listhg,    name='listhg',    routine= 'amn' )
    call de_alloc( listhptrg, name='listhptrg', routine= 'amn' )
    call de_alloc( numhg,     name='numhg',     routine= 'amn' )
    call de_alloc( xijloc,    name='xijloc',    routine= 'amn' )
    call de_alloc( Sloc,      name='Sloc',      routine= 'amn' )
  endif

!! For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIError)
!#endif
!    call die()
!! End debugging


! End time counter
  call timer( 'Amn', 2 )

  return

  contains

  function get_overlapping_orbitals( latvec, atomrcut, trialrcut, numatoms, &
 &                                   atomcoords, trialcenter )              &
 &  result(firstitem)

    use precision,          only: dp            ! Real double precision type
    use sys,                only: die                 
    use neighbour,          only: maxnna        ! Maximum number of neighbours
    use neighbour,          only: jan           ! Atom-index of neighbours
    use neighbour,          only: r2ij          ! Squared distances to neighbors
    use neighbour,          only: xij           ! Vector from a given atom
                                                !   to its neighbours
    use neighbour,          only: mneighb       ! Subroutine to compute the
                                                !   number of neighbours
    use siesta_geom,        only: isa           ! Species index of each atom
    use atomlist,           only: lasto         ! Position of last orbital 
                                                !   of each atom
    use atomlist,           only: iphorb        ! Orbital index of each  orbital
                                                !   in its atom
    use atomlist,           only: indxuo        ! Index of equivalent orbital  
                                                !   in "u" cell
    use atmfuncs,           only: rcut          ! Function that determines the
                                                !   cutoff radius of a given
                                                !   orbital of a given specie

!   This subroutine yields a list of basis orbitals that overlap 
!   with a given trial orb.
!   This subroutine bridges siesta's mneighb() with our purposes.

    implicit none
  
!   Passing variables
    real(dp),dimension(3,3),intent(in) :: latvec      ! Lattice vectors of the 
                                                      !   supercell in 
                                                      !   real space.
    real(dp)               ,intent(in) :: atomrcut    ! Maximum cutoff radius in
                                                      !   the atomic orbital 
                                                      !   basis
    integer                ,intent(in) :: numatoms    ! Number of atoms in 
                                                      !   the supercell
    real(dp),dimension(3,numatoms),intent(in) :: atomcoords  ! Atomic positions 
    real(dp),dimension(3)  ,intent(in) :: trialcenter ! Center of the 
                                                      !  trial function 
    real(dp)               ,intent(in) :: trialrcut   ! Cutoff of the 
                                                      !  trial function 

!   Passing variables
    real(dp), dimension(:,:), allocatable,save :: coords 
!   A new array containing the coordinates of the na_s atoms in the unit cell
!   plus the position of the trial function will be required to call mneighb.
!   This new array, coords, will have (3,na_s+1) dimensions

    integer          :: jneig   ! Counter for the loop on neighbours
    integer          :: nneig   ! Number of neighbors of a given trial projector
    integer          :: atom    ! Atomic index of the neighbour
    integer          :: specie  ! Atomic species of the neighbour
    integer          :: iorb    ! Counter for loop on neighbour orbitals
    integer          :: joa     ! Index for the atomic orbital within 
                                !    a given atom
    real(dp)         :: radius  ! Radius that determine the scope of the search
    real(dp)         :: rij     ! Squared distance to neighbours

    type(orbitallinkedlist), pointer     :: firstitem
    type(orbitallinkedlist), pointer     :: newitem

! Initialize mneighb
    if( .not. allocated(coords) ) then !set up static "save" variables
      allocate( coords(3,numatoms+1) )
      coords(1:3,1:numatoms) = atomcoords(1:3,1:numatoms)
    endif

    coords(:,numatoms+1) = trialcenter(:)
    radius = trialrcut + atomrcut
    call mneighb( latvec, radius, numatoms+1, coords, 0 , 0, nneig )

!
! Look for atoms that overlap with the trial orbital
!
    call mneighb( latvec, radius, numatoms+1, coords, numatoms+1, 0, nneig )
    if (nneig.gt.maxnna)                                             &
 &    call die("amn: insufficient array shapes; see mneighb(..)")

!
! Prepare output list 
!

    firstitem => null()
    do jneig = 1, nneig
      atom = jan(jneig)
      if ( atom .gt. numatoms ) cycle !it is a trial projection function
      specie = isa(atom)              !it is an atom
      rij = dsqrt( r2ij(jneig) )
      do iorb = lasto( atom-1 )+1, lasto( atom )
        joa = iphorb( iorb )
        if ( (rcut(specie,joa) + trialrcut) .gt. rij ) then
          allocate(newitem)
          newitem%center(1:3) = xij(1:3,jneig) + trialcenter(1:3)
          newitem%specie      = specie
          newitem%specieindex = joa
          newitem%globalindex = indxuo( iorb )
          newitem%nextitem    => firstitem
          firstitem => newitem
        endif
      enddo
    enddo    ! End loop on number of neighbours
  end function get_overlapping_orbitals
end subroutine amn

