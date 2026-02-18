!
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_wannier_in_nao module:
!! Here we write the expansion of a Wannier function, 
!! output of the WANNIER-90 code, in the numerical atomic orbital (NAO)
!! basis set of SIESTA.
!!
!! The coefficients of the expansion are stored in the public variable
!! coeffs_wan_nao
!!
!! All the atomic orbitals *in the supercell* participate in the expansion
!! although the coefficients will decay very fast with the distance from the
!! center of the Wannier to the center of the atomic orbital
!!
!! It is parallelized in such a way that a given node knows all the coefficients
!! for a subset of numh_man_proj Wanniers
!!
!! Since both the Wannier functions and the localized atomic orbitals are real,
!! the coefficients of this expansion are also expected to be real
!!
!! The coefficients are written in files with the .WANNX format that can be
!! laterly read by DENCHAR for plotting
!! 
!! The analytical expressions can be found in
!!
!! <https://personales.unican.es/junqueraj/JavierJunquera_files/Notes/Wannier/wannier_in_nao.pdf>
!!

module m_wannier_in_nao

  use precision,      only: dp                ! Real double precision type

  implicit none

CONTAINS

#ifdef SIESTA__WANNIER90

  subroutine wannier_in_nao( ispin, index_manifold )
    use parallel,       only: Nodes          ! Total number of Nodes
    use parallel,       only: Node           ! Local Node
    use parallel,       only: BlockSize      ! BlockSize
    use parallel,       only: IONode         ! Input/output node
    use siesta_geom,    only: xa             ! Atomic positions for all the
                                             !   atoms in the supercell
                                             !   In Bohrs
    use siesta_geom,    only: isa            ! Species index of each atom
    use siesta_geom,    only: ucell          ! Unit cell lattice vectors
    use siesta_geom,    only: nsc            ! Diagonal element of the supercell
    use cellsubs,       only: reclat         ! Finds reciprocal unit cell vector
    use atomlist,       only: no_s           ! Number of atomic orbitals in the
                                             !   supercell
    use atomlist,       only: no_u           ! Number of atomic orbitals in the
                                             !   unit cell
    use atomlist,       only: indxuo         ! Index of equivalent orbital in 
                                             !   the unit cell
    use atomlist,       only: iaorb          ! Atomic index of each orbital
    use atomlist,       only: iphorb         ! Orbital index of each orbital 
                                             !   in its atom
    use atmfuncs,       only: labelfis       ! Atomic label
    use atmfuncs,       only: cnfigfio       ! Principal quantum number of the
                                             !   shell to what the orbital 
                                             !   belongs (for polarization 
                                             !   orbitals the quantum number 
                                             !   corresponds to the shell which
                                             !   is polarized by the orbital io)
    use atmfuncs,       only: symfio         ! Symmetry of the orbital

    use units,          only: Ang
    
    ! data to map to info formerly inside wannier modules ----
    use w90_wrapper_types,       only: w90_manifold_t
    use w90_wrapper_types,       only: manifold_bands
    use w90_wrapper_types,       only: reclatvec_w90_in, &
                                         numkpoints_w90_in, &
                                         kpointsfrac_w90_in
    !-------------
    use w90_wrapper_types,       only: numh_man_proj
                                             ! Number of projectors that will be
                                             !   handled in the local node
    use w90_wrapper_types,       only: listhptr_man_proj
                                       ! Index pointer to listh_man_proj such
                                       ! listh_man_proj(listhptr_man_proj(1)+1)
                                       !   is the first projector of the first
                                       !   manifold handled by the local node
                                       ! listh_man_proj(listhptr_man_proj(io)+1)                                                                   !   is thus the first projector of
                                       !   of manifold 'io' while
                                       ! listh_man_proj(listhptr_man_proj(io) +                                                                    !   numh_man_proj(io))
                                       !   is the last projectors of manifold
                                       !   'io'.
                                       ! Dimension: number of manifolds
    use w90_wrapper_types,       only: listh_man_proj
                                       ! The column indices for the projectors
                                       !   of all the manifolds handled by
                                       !   the local node
    use w90_wrapper_types,       only: coeffs_wan_nao
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
    use m_switch_local_projection, only: coeffs
                                             ! Coefficients of the wavefunctions
                                             !   First  index: orbital
                                             !   Second index: band
                                             !   Third  index: k-point
    use m_switch_local_projection, only: seedname
                                             ! Name of the file where the 
                                             !   WANNIER90 code
                                             !   reads or dumps the
                                             !   information.
    use writewave,      only: writew         ! Subroutine to dump the 
                                             !   coefficients of the 
                                             !   wavefunction
                                             !   in a .WANNX file
    use writewave,      only: setup_wfs_list ! Subroutine to setup the
                                             !   list of wave functions 
                                             !   (in this case, Wannier func.)
                                             !   to be plotted
    use writewave,      only: wfs_filename   ! Name of the file where the 
                                             !   coefficients of the Wanniers
                                             !   will be dumped
    use iso_c_binding, only : c_f_pointer, c_loc

!
! Allocation/Deallocation routines
!
    use alloc,          only: re_alloc       ! Reallocation routines
    use alloc,          only: de_alloc       ! Deallocation routines

!
! Termination routines
!
    use sys,            only: die            ! Termination routine

#ifdef MPI
    use mpi_siesta
    use m_orderbands,   only: which_band_in_node  
                                             ! Given a node and a local index,
                                             !   this array gives the
                                             !   global index of the band
                                             !   stored there
    use m_orderbands,   only: sequential_index_included_bands
                                             ! Sequential number of the
                                             !   bands included for
                                             !   wannierization
                                             !   (the bands are listed
                                             !   in order of incremental
                                             !   energy)
    use m_switch_local_projection, only: nincbands_loc
                                             ! Number of bands for
                                             !   wannierization
                                             !   after excluding bands
                                             !   in the local Node
#endif

!
! Input variables
!
    integer, intent(in) :: ispin            ! Spin component
    integer, intent(in) :: index_manifold   ! Index of the manifold
                                            !   that is wannierized

! 
! Internal variables 
! 
    integer  :: iu                           ! Logical unit
    integer  :: ik                           ! Counter for loop on k-points
    integer  :: iband                        ! Counter for loop on bands
    integer  :: iproj                        ! Counter for loop on projectors
    integer  :: iproj_global                 ! Counter for loop on projectors
    integer  :: iprojn                       ! Counter for loop on projectors
    integer  :: iprojm                       ! Counter for loop on projectors
    integer  :: ind                          ! Counter for sequential indices
                                             !    of projections
    integer  :: iorb                         ! Counter for loop on atomic 
                                             !   orbitals in the supercell
    integer  :: ia                           ! Atom to which orbital belongs
    integer  :: iua                          ! Equivalent atom in the
                                             !    first unit cell
    integer  :: iuo                          ! Equivalent orbital in the
                                             !    first unit cell
    integer  :: ix                           ! Counter for cartesian directions
    integer  :: icell                        ! Counter for unit cell lattice
                                             !    vectors
    integer  :: nk                           ! Number of k-points written in
                                             !   the .WANNX file
                                             !   Since the number of Wanniers 
                                             !   functions to be plotted do not 
                                             !   depend in k-points,
                                             !   it is set to 1
    integer  :: maxspn                       ! Maximum number of spin components
    logical  :: gamma                        ! Whether only the Gamma-point is 
                                             !    sampled.
                                             !    Since the dependency of the 
                                             !    phase is already included in 
                                             !    the coefficients, we set it
                                             !    up to .true.
    real(dp) :: kpoint(3)                    ! Coordinates of the k-point
    real(dp) :: kdummy(3)                    ! Dummy variable for the k-points
    real(dp) :: kxmu                         ! Dot product of the k-point and
                                             !   the position of the atom in the
                                             !   supercell
    real(dp) :: ckxmu                        ! Cosine of kxmu
    real(dp) :: skxmu                        ! Sine of kxmu
    real(dp) :: xatorb(3)                    ! Position of the atomic orbital
    real(dp), pointer :: psi(:,:,:)          ! Dummy variable to store the 
                                             !   coefficients of the Wanniers
                                             !   in a basis of atomic orbitals
                                             !   to call writew
    real(dp), target :: aux(2,no_s*5)        ! Dummy variable that will play
                                             !   the role of the eigenvalues
                                             !   in the call to writew
    complex(dp), dimension(:,:), pointer :: psiloc => null() 
                                             ! Coefficients of the wave
                                             !   function (in complex format)
    integer  :: isc(3)
    real(dp) :: dxa(3)                       ! Cell vector that translates a
                                             !   given atom in the unit cell
                                             !   to the equivalent in the 
    real(dp) :: rcell(3,3)                   ! Reciprocal unit cell vectors
                                             !   (without 2*pi factor)
    complex(dp), pointer :: coeffs_opt(:,:,:) => null()
                                             ! eigenvalues at the Wannier90 
                                             !    k-point mesh. 
                                             !    They are the ones than comes
                                             !    out for the diagonalization
                                             !    but if a disentanglement 
                                             !    is required, a farther 
                                             !    transformation is required 
                                             !    to generate the optimal 
                                             !    coefficients as explained 
                                             !    by Souza et al. 
                                             !    as explained below
    real(dp), pointer :: psi_real_1d(:)      ! Temporal array to store the
                                             !    coefficients of the Wanniers
                                             !    in a basis of NAO in a
                                             !    1D array in order to print
                                             !    them


#ifdef MPI
    integer     :: iband_global              ! Global index for a band
    integer     :: iband_sequential          ! Sequential index of the band
    integer     :: MPIerror
    complex(dp), dimension(:,:), pointer :: auxloc => null()
                                             ! Temporal array for the
                                             !   the global reduction of Amnmat
    integer, external :: numroc
#endif

    type(w90_manifold_t), pointer :: mnf

    ! mapping variables
#ifdef SIESTA__UNITS_ORIGINAL
    real(kind=dp), parameter :: bohr_angstrom_internal = 0.52917720859_dp
#endif
    integer :: num_proj, num_bands, num_kpts
    logical :: disentanglement
    real(dp) :: recip_lattice(3,3)
    real(dp), allocatable :: kpt_latt(:,:)

    ! Convenience device to point to manifold data
    mnf => manifold_bands(index_manifold)
    
    num_proj = mnf%numbands
    num_bands = mnf%number_of_bands
    num_kpts = numkpoints_w90_in
    disentanglement = mnf%disentanglement
    recip_lattice(:,:) = reclatvec_w90_in(:,:)
    allocate(kpt_latt,source=kpointsfrac_w90_in)
    
!   Find reciprocal unit cell vectors (without 2*pi factor)
    call reclat( ucell, rcell, 0 )

!   Allocate the array where the coefficients of the 
!   bands that will be wannierized
!   If there is no-disentanglement required, the coefficients are
!   the same that come out from the diagonalization.
!   If a disentanglement procedure is required 
!   (more bands that Wannier functions),
!   then we follow the recipe described in Sec. III A of the paper by
!   I. Souza et al. Phys. Rev. B 65 035109 (2001)
!   To find the N-dimensional subspace if the number of bands for
!   a k-point, N_k, is larger than N, we have to multiply the
!   coefficients by a matrix that comes from Wannier90
    nullify( coeffs_opt )
    call re_alloc( coeffs_opt,                                      &
 &                 1, no_u,                                         &
 &                 1, num_proj,                                     &
 &                 1, num_kpts,                                     &
 &                 name='coeffs_opt', routine='wannier_in_nao')
    coeffs_opt = cmplx(0.0_dp,0.0_dp,kind=dp)

!   Allocate memory related with a local variable where the coefficients
!   of the eigenvector at the k-point will be stored
!   Only num_bands are retained for wannierization, that is why the
!   second argument is made equal to num_bands
    call re_alloc( psiloc, 1, no_u, 1, num_bands,     &
 &                'psiloc', 'wannier_in_nao' )

#ifdef MPI
!   Store the local bands in this node on a complex variable
    do ik = 1, num_kpts

!     Initialize the local coefficient matrix for every k-point
      psiloc(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

      do iband = 1, nincbands_loc
        iband_global     = which_band_in_node(Node,iband)
        iband_sequential = sequential_index_included_bands(iband_global)

!!       For debugging
!        write(6,'(a,7i5)')                        &
! &       'Nodes, Node,ik, nbands_loc, iband, iband_global, iband_sequential=',&
! &        Nodes, Node, ik, nincbands_loc, iband, iband_global, iband_sequential
!!       End debugging

        do iorb = 1, no_u
          psiloc(iorb,iband_sequential) = coeffs(iorb,iband,ik)
        enddo 

      enddo 
!     Allocate workspace array for global reduction
      call re_alloc( auxloc, 1, no_u, 1, num_bands,                  &
 &                   name='auxloc', routine='wannier_in_nao' )
!     Global reduction of auxloc matrix
      auxloc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
      call MPI_AllReduce( psiloc(1,1), auxloc(1,1),        &
 &                        no_u*num_bands,                  &
 &                        MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
!     After this reduction, all the nodes know the coefficients of the
!     wave function for the point ik, for all the bands and for all atomic
!     orbitals
      psiloc(:,:) = auxloc(:,:)

      if( disentanglement ) then
        do iorb = 1, no_u
          do iprojn = 1, num_proj
            do iband = 1, num_bands
              coeffs_opt(iorb,iprojn,ik) = coeffs_opt(iorb,iprojn,ik) +  &
 &              mnf%u_matrix_opt(iband,iprojn,ik) * psiloc(iorb,iband)
            enddo 
          enddo
        enddo
      else
        coeffs_opt(:,:,ik) = psiloc(:,:) 
      endif

!!     For debugging
!      write(6,'(a,i5,3f12.5)')      &
! &      ' ik, kpt_latt(:,ik) = ',   &
! &        ik, kpt_latt(:,ik) 
!      if( disentanglement ) then
!        do iprojn = 1, num_proj
!          do iband = 1, num_bands
!            write(6,'(a,3i5,2f12.5)')                & 
! &           'Node, iband, iprojn, mnf%u_matrix_opt = ', &
! &            Node, iband, iprojn, mnf%u_matrix_opt(iband,iprojn,ik)
!          enddo 
!        enddo 
!      endif
!      do iband = 1, num_bands
!        do iorb = 1, no_u
!          write(6,'(3i5,2f12.5)') Node, iband, iorb, psiloc(iorb,iband)
!        enddo 
!      enddo 
!      do iprojn = 1, num_proj
!        do iorb = 1, no_u
!          write(6,'(a,3i5,2f12.5)')                                     &
! &          'Node, iprojn, iorb, coeffs_opt(iorb,iprojn,ik) = ',   &
! &          Node, iprojn, iorb, coeffs_opt(iorb,iprojn,ik)
!        enddo 
!      enddo 
!!     End debugging

    enddo    ! End loop on k-points
#else
    do ik = 1, num_kpts
      do iband = 1, num_bands
        do iorb = 1, no_u
          psiloc(iorb,iband) = coeffs(iorb,iband,ik)
        enddo 
      enddo
      if( disentanglement ) then
        do iorb = 1, no_u
          do iprojn = 1, num_proj
            do iband = 1, num_bands
              coeffs_opt(iorb,iprojn,ik) = coeffs_opt(iorb,iprojn,ik) +  &
 &              mnf%u_matrix_opt(iband,iprojn,ik) * psiloc(iorb,iband)
            enddo 
          enddo
        enddo
      else
        coeffs_opt(:,:,ik) = psiloc(:,:) 
      endif
    enddo
#endif

!   Loop over all the local Wannier functions handled locally in this node
    do iprojn = 1, numh_man_proj(index_manifold)
      ind          = listhptr_man_proj(index_manifold) + iprojn
      iproj_global = listh_man_proj(ind)

!!     For debugging
!      write(6,'(a,4i5)')                                                 &
! &      'wannier_in_nao: Node, index_manifold, iprojn, iproj_global = ', &
! &                       Node, index_manifold, iprojn, iproj_global 
!!     End debugging

!     Loop on all the orbital of the supercell to compute the corresponding
!     coefficient of the Wannier function on it.
      do iorb = 1, no_s
        iuo = indxuo(iorb)            ! Equivalent orbital in first unit cell
        ia  = iaorb(iorb)             ! Atom to which orbital belongs
        iua = iaorb(iuo)              ! Equivalent atom in first unit cell
        dxa(:) = xa(:,ia) - xa(:,iua) ! Cell vector of atom ia
        isc(:) = nint( matmul(dxa,rcell) )  ! Cell index of atom ia
!       Find the index of the unit cell within the supercell where this
!       atom is located, centered on isc = 0
        do ix = 1,3
          if (isc(ix)>nsc(ix)/2) isc(ix) = isc(ix) - nsc(ix) ! Same centered 
                                                             ! in isc=0
        enddo
!       Find the translated position of the atom in the supercell that 
!       really takes a non-zero value in the unit cell 
        xatorb(:) = xa(:,iua)
        do icell = 1, 3
          do ix = 1, 3
            xatorb(ix) = xatorb(ix) + isc(icell) * ucell(ix,icell)
          enddo 
        enddo 

!!       For debugging          
!        write(6,'(a,10i5,3f12.5)')                                    &
! &        'wannier_in_nao: Node, iproj_local, iproj_global, iorb, iuo, ia, iua, isc = ', &
! &                         Node, iprojn, iproj_global, iorb, iuo, ia, iua, isc(:), xatorb(:)
!!       End debugging          

!       Up to here, we know:
!       - the projector that will be expressed in a basis of NAO
!       - the atomic orbital for which the coefficient will be found
!       - the position of the atom where the atomic orbital is centered
!       Now, we perform the sum on k-point in the Equation written
!       in the header of the subroutine

        do ik = 1, num_kpts
!         Compute the coordinates of the k-point (in Ang^-1)
          kpoint(:) = kpt_latt(1,ik) * recip_lattice(1,:) +      &
 &                    kpt_latt(2,ik) * recip_lattice(2,:) +      &
 &                    kpt_latt(3,ik) * recip_lattice(3,:) 
!         Transform the coordinates of the k-point to Bohr^-1

#ifdef SIESTA__UNITS_ORIGINAL
          kpoint(:) = kpoint(:) * bohr_angstrom_internal
#else
          kpoint(:) = kpoint(:) / Ang
#endif

!         Compute the dot product between the k-point and the
!         atomic position
          kxmu = kpoint(1) * xatorb(1) +                        &
 &               kpoint(2) * xatorb(2) +                        &
 &               kpoint(3) * xatorb(3) 
          ckxmu = dcos(kxmu)
          skxmu = dsin(kxmu)

          do iprojm = 1, num_proj
            coeffs_wan_nao(ind,iorb,ispin) =                   &
 &            coeffs_wan_nao(ind,iorb,ispin) +                 &
 &            mnf%u_matrix(iprojm,iproj_global,ik)           *     &
 &            coeffs_opt(iuo,iprojm,ik)                  *     &
 &            cmplx(ckxmu,skxmu,kind=dp) 
          enddo 

        enddo ! End loop on k-points

      enddo   ! End loop on atomic orbitals
!     Divide by the number of k-points
      coeffs_wan_nao(ind,:,ispin) =                                  &
 &      coeffs_wan_nao(ind,:,ispin) / num_kpts
     
    enddo 

!!   For debugging
!    do iproj = 1, numh_man_proj(index_manifold) 
!      ind          = listhptr_man_proj(index_manifold) + iproj
!      iproj_global = listh_man_proj(ind)
!      do iorb = 1, no_s
!        write(6,'(a,8i5,2f12.5)')   &
! &        'Node, Nodes, ispin, iproj, ind, iproj_global, iorb, indxuo, coeffs_wan=',&
! &         Node, Nodes, ispin, iproj, ind, iproj_global, iorb, indxuo(iorb),        &
! &         coeffs_wan_nao(ind,iorb,ispin)
!      enddo 
!    enddo 
!!   End debugging

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( psi )
    call re_alloc( psi,                                             &
 &                 1, 2,                                            &
 &                 1, no_s,                                         &
 &                 1, numh_man_proj(index_manifold),                &
 &                 name='psi', routine='wannier_in_nao')
    psi = 0.0_dp

!!    For debugging
!    do ik = 1, num_kpts
!      do iproj = 1, num_proj
!        do iorb = 1, no_u
!          write(6,'(a,3i5,4f12.5)')' ik, iproj, iorb, coeff = ', &
! &          ik, iproj, iorb, coeffs(iorb,iproj,ik),              &
! &          mnf%u_matrix(iproj,iproj,ik)
!        enddo 
!      enddo 
!    enddo 
!
!     do iproj = 1, numh_man_proj(index_manifold)
!       do iorb = 1, no_s
!!!         if( indxuo(iorb) .eq. 15) then
!!         if( aimag(coeffs_wan_nao(index_manifold,iorb,ispin)) .gt. 1.d-5 ) then
!         if( Node .eq. 1) then
!         write(6,'(a,4i5,5f12.5)') ' Node, Nodes, iproj, iorb, coeffs_wan=', &
! &         Node, Nodes, iproj, iorb,                                         &
! &         coeffs_wan_nao(index_manifold,iorb,ispin), xa(:,iaorb(iorb))
!         endif
!!         endif
!       enddo 
!     enddo 
!!   End debugging

!   Set up the variables to call writew
    do iproj = 1, numh_man_proj( index_manifold )
      ind          = listhptr_man_proj(index_manifold) + iproj
      iproj_global = listh_man_proj(ind)
      do iorb = 1, no_s
        psi(1,iorb,iproj) = real(coeffs_wan_nao(ind,iorb,ispin))
        psi(2,iorb,iproj) = aimag(coeffs_wan_nao(ind,iorb,ispin))
!!       For debugging
!        if( abs(psi(1,iorb,iproj)) .gt. 1.d-5) then
!        if (iproj .eq. 4 ) then
!        write(6,'(a,2i5,8f12.5)')' iproj, iorb, psi = ', &
! &        iproj, iorb, psi(1,iorb,iproj), psi(2,iorb,iproj), xa(:,iaorb(iorb)), xafold(:,iaorb(iorb)) - xa(:,1)
!        endif 
!        if (iproj .eq. 5 ) then
!        write(6,'(a,2i5,8f12.5)')' iproj, iorb, psi = ', &
! &        iproj, iorb, psi(1,iorb,iproj), psi(2,iorb,iproj), xa(:,iaorb(iorb)), xafold(:,iaorb(iorb)) - xa(:,2)
!        endif 
!        endif
!!       End debugging
      enddo 
    enddo

    aux          = 0.0_dp
    kdummy       = 0.0_dp
    wfs_filename = trim(seedname)//".WANNX"
    nk           = 1
    gamma        = .false.
    maxspn       = 1
    

! 
!   Open the WANNX file and print the header of the file
!   with information of the atomic orbitals
!   This was done in the subroutine wwave when the coefficients of the 
!   wave functions at particular k-points are required
!
    if( IONode ) then
      call io_assign( iu )
      open(iu, file=wfs_filename,form="unformatted",status='unknown')

      write(iu) nk, gamma
      write(iu) maxspn
      write(iu) no_s
      write(iu) (iaorb(iorb),labelfis(isa(iaorb(iorb))),                    &
 &              iphorb(iorb), cnfigfio(isa(iaorb(iorb)),iphorb(iorb)),      &
 &              symfio(isa(iaorb(iorb)),iphorb(iorb)), iorb=1,no_s)
      call io_close( iu )
    endif

!!     For debugging
!      write(6,'(a,i5,l5)')' nk, gamma = ', nk, gamma
!      write(6,'(a,i5)')   ' maxspn = ', maxspn
!      write(6,'(a,i5)')   ' no_s   = ', no_s
!      write(6,'(a,i5)')   ' Node, BlockSize   = ', Node, BlockSize
!      do iorb = 1, no_s
!        write(6,*) iorb, iaorb(iorb),labelfis(isa(iaorb(iorb))),            &
! &                 iphorb(iorb), cnfigfio(isa(iaorb(iorb)),iphorb(iorb)),   &
! &                 symfio(isa(iaorb(iorb)),iphorb(iorb))
!      enddo 
!!     End debugging

    call setup_wfs_list(nk,no_s,1,num_proj,.false.,.false.)

    call c_f_pointer(c_loc(psi(1,1,1)),psi_real_1d,                        &
 &                   [size(psi,1)*size(psi,2)*size(psi,3)])

    call writew( no_s, num_proj, 1, kdummy, 1,                                &
 &               aux,psi_real_1d, &
 &               gamma, non_coll=.false., BlockSize=BlockSize  )

!   Deallocate some pointers
    call de_alloc( psiloc     )
    call de_alloc( psi        )
    call de_alloc( coeffs_opt )
#ifdef MPI
    call de_alloc( auxloc )
#endif

!!   For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!    call die()
!!   End debugging

  end subroutine wannier_in_nao

#else

  subroutine wannier_in_nao

  end subroutine wannier_in_nao
#endif

end module m_wannier_in_nao
