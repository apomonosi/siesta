!
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the switch_local_projection module:
!! The subroutines that compute the overlap of the periodic part of the 
!! wave functions at neighbour k-points, and the overlap between the wave 
!! functions and the initial guess functions might be called coming from
!! different inputs:
!!
!! 1) using WANNIER90 as a postprocessing tool (the "legacy" interface)
!! 2) calling WANNIER90 directly from SIESTA ("w90 wrapper" interface).
!!
!! Depending of the procedence, the input information is read in a
!! different way.  Here, depending on the kind of projection to
!! localized orbitals that will be carried out, the relevant matrices
!! and arrays that will be used in the subroutines mmn and amn are
!! populated.

module m_switch_local_projection

  use precision,      only: dp                ! Real double precision type
  use siesta_options, only: w90_processing    ! Will we call the interface with
                                              !    WANNIER90 as a postprocessing
                                              !    tool?
  use siesta_options, only: w90_wrapper_processing
                                              ! Will we call the WANNIER90 
                                              !   subroutines directly from 
                                              !   SIESTA (not as a 
                                              !   postprocessing tool)
                                              !   the Wannier transformation
  use atomlist,       only: no_u              ! Number of orbitals in unit cell
  use siesta_geom,    only: ucell             ! Unit cell lattice vectors
  use files,          only: label_length      ! Number of characters in slabel
  use trialorbitalclass

  implicit none

  integer           :: numkpoints        ! Total number of k-points
                                         !   in the chosen method of projection
  real(dp), pointer :: kpointsfrac(:,:)  ! List of k points relative
                                         !   to the reciprocal lattice vectors.
                                         !   First  index: components.
                                         !   Second index: k-point index in 
                                         !   the list
  integer           :: numbands(2)       ! Number of bands for wannierization
                                         !   before excluding bands
  integer           :: numincbands(2)    ! Number of included bands in the 
                                         !   calculation of the chosen 
                                         !   projection method after excluding
                                         !   bands
  logical, pointer :: isexcluded(:)      ! Masks excluded bands
  integer          :: blocksizeincbands  ! Maximum number of bands
                                         !   considered for wannierization 
                                         !   per node
  integer          :: nincbands_loc      ! Number of included bands in the 
                                         !   calculation of the overlap and 
                                         !   projection matrices.
  real(dp)         :: latvec(3,3)        ! Real space lattice vectors
  real(dp)         :: reclatvec(3,3)     ! Reciprocal lattice vectors
                                         !  Cartesian coordinates in Bohr^-1 
                                         !  First  index: component
                                         !  Second index: vector
!
! Variables related with the k-point list for which the overlap
! matrices Mmn between a k-point and its neighbor will be computed
!
  integer          :: nncount            ! The number of nearest
                                         !   neighbours belonging to
                                         !   each k-point of the 
                                         !   Monkhorst-Pack mesh
  integer, pointer :: nnlist_neig(:,:)   ! nnlist_neig(ikp,inn) is the index 
                                         !   of the inn-neighbour of ikp-point
                                         !   in the Monkhorst-Pack grid 
                                         !   folded to the first Brillouin zone
  integer, pointer :: nnfolding(:,:,:)   ! nnfolding(i,ikp,inn) is the
                                         !   i-component of the reciprocal 
                                         !   lattice vector
                                         !   (in reduced units) that brings
                                         !   the inn-neighbour specified in 
                                         !   nnlist_neig (which is in 
                                         !   the first BZ)
                                         !   to the actual \vec{k} + \vec{b} 
                                         !   that we need.
                                         !   In reciprocal lattice units.
  real(dp), pointer :: bvectorsfrac(:,:) ! The vectors b that connect
                                         !   each mesh-point k
                                         !   to its nearest neighbours


!
! Variables related with the coefficients of the wavefunctions and
! eigenvalues at the Wannier90 k-point mesh
!
  complex(dp), pointer :: coeffs(:,:,:) => null() 
                                         ! Coefficients of the wavefunctions.
                                         !   First  index: orbital
                                         !   Second index: band
                                         !   Third  index: k-point
  real(dp),    pointer :: eo(:,:) => null()  
                                         ! Eigenvalues of the Hamiltonian
                                         !   at the numkpoints introduced in
                                         !   kpointsfrac
                                         !   First  index: band index
                                         !   Second index: k-point index

!
! Variables related with the projections with trial functions,
! initial approximations to the MLWF
!
  integer  :: numproj                    ! Total number of projection centers,
                                         !   equal to the number of MLWF

  type(trialorbital), target, allocatable  :: projections(:)


!
! Output matrices
!

  complex(dp), pointer :: Mmnkb(:,:,:,:) => null()  
                                         ! Matrix of the overlaps of
                                         !   periodic parts of Bloch waves.
                                         !   <u_{ik}|u_{jk+b}>
                                         !   The first two indices refer to
                                         !   the number of occupied bands
                                         !   (indices m and n in standard
                                         !   notation, see for instance,
                                         !   Eq. (27) of the paper by
                                         !   Marzari et al., RMP 84, 1419 (2012)
                                         !   The third index refer to the kpoint
                                         !   The fourth index refer to the neig
  complex(dp), pointer :: Amnmat(:,:,:) => null() 
                                         ! Projections of a trial function
                                         !   with a Bloch orbital
                                         !   <\psi_{m k}|g_n>

!
! Variables related with the input/output files
!
  character(label_length+3)  :: seedname ! Name of the file where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, reads or dumps the
                                         !   information.



! Routines
  public :: switch_local_projection

  CONTAINS

  subroutine switch_local_projection( index_manifold )
    use w90_wrapper_types, only: manifold_bands
                                                  ! Derived type where all 
                                                  !   the details of the 
                                                  !   band manifolds for 
                                                  !   Wannier transformation
                                                  !   are stored
    use w90_wrapper_types, only: numkpoints_w90_in  
                                                  ! Number of k-points in 
                                                  !   the Monkhorst-Pack grid
                                                  !   that will be used in
                                                  !   Wannier90
    use w90_wrapper_types, only: kpointsfrac_w90_in  
                                                  ! List of k-points in 
                                                  !   the Monkhorst-Pack grid
                                                  !   that will be used in
                                                  !   the Wannier
                                                  !   transformation
                                                  !   (in fractional units, 
                                                  !   i.e. relative to the
                                                  !   reciprocal space 
                                                  !   lattice vectors)
    use w90_wrapper_types, only: nncount_w90_in    
                                                  ! The number of nearest
                                                  !   neighbours belonging to
                                                  !   each k-point of the 
                                                  !   Monkhorst-Pack mesh
    use w90_wrapper_types, only: nnlist_w90_in
                                                  ! nnlist_w90_in(ikp,inn)
                                                  !   is the index of the
                                                  !   inn-neighbour of ikp-point
                                                  !   in the Monkhorst-Pack grid
                                                  !   folded to the
                                                  !   first Brillouin zone
    use w90_wrapper_types, only: nnfolding_w90_in
                                                  ! nnfolding(i,ikp,inn) is the 
                                                  !   i-component
                                                  !   of the reciprocal lattice 
                                                  !   vector
                                                  !   (in reduced units) that 
                                                  !   brings the inn-neighbour 
                                                  !   specified in nnlist_w90_in
                                                  !   (which is in the first BZ)
                                                  !   to the actual
                                                  !   \vec{k} + \vec{b} that we
                                                  !   need.
                                                  !   In reciprocal 
                                                  !   lattice units.
    use w90_wrapper_types, only: bvectorsfrac_w90_in
    use w90_wrapper_types, only: latvec_w90_in
 
    use wannier90_types, only: numbands_wannier   ! Number of bands for 
    use wannier90_types, only: numkpoints_wannier ! Number of k-points in 
                                                  !   the Monkhorst-Pack grid
                                                  !   for which the overlap of
                                                  !   the periodic part of the
                                                  !   wavefunct with a 
                                                  !   k-point will be computed
    use wannier90_types, only: kpointsfrac_wannier! List of k-points in 
                                                  !   the Monkhorst-Pack grid
                                                  !   relative to the 
                                                  !   reciprocal lattice 
                                                  !   vectors.
                                                  !   First  index: components
                                                  !   Second index: k-point 
                                                  !   index in the list
    use wannier90_types, only: nncount_wannier    ! The number of nearest
                                                  !   neighbours belonging to
                                                  !   each k-point of the 
                                                  !   Monkhorst-Pack mesh
    use wannier90_types, only: nnlist_wannier     ! nnlist_wannier(ikp,inn)
                                                  !   is the index of the
                                                  !   inn-neighbour of ikp-point
                                                  !   in the Monkhorst-Pack grid
                                                  !   folded to the
                                                  !   first Brillouin zone
    use wannier90_types, only: nnfolding_wannier  ! nnfolding(i,ikp,inn) is the
                                                  !   i-component
                                                  !   of the reciprocal lattice 
                                                  !   vector (in reduced units) 
                                                  !   that brings the 
                                                  !   inn-neighbour specified 
                                                  !   in nnlist_wannier
                                                  !   (which is in the first BZ)
                                                  !   to the actual 
                                                  !   \vec{k} + \vec{b} that we
                                                  !    need.
                                                  !   In reciprocal lattice unit
    use wannier90_types, only: numbands_wannier   ! Number of bands for 
                                                  !   wannierization before
                                                  !   excluding bands
    use wannier90_types, only: numincbands_wannier! Number of included bands
                                                  !   in the calculation of
                                                  !   the overlap and projec.
                                                  !   matrices after excluding
                                                  !   bands
    use wannier90_types, only: nincbands_loc_wannier
                                                  ! Number of bands in the 
                                                  !   local node for 
                                                  !   wannierization after
                                                  !   after excluding bands
    use wannier90_types, only: blocksizeincbands_wannier
                                                  ! Maximum number of bands per
                                                  !   node considered for 
                                                  !   for wannierization
    use wannier90_types, only: isexcluded_wannier ! Masks excluded bands for
                                                  !   Wannier
    use wannier90_types, only: latvec_wannier     ! Reciprocal lattice vectors
    use wannier90_types, only: reclatvec_wannier  ! Reciprocal lattice vectors
    use wannier90_types, only: bvectorsfrac_wannier
    use wannier90_types, only: numproj_wannier    ! Number of projectors
    use wannier90_types, only: seedname_wannier
    use wannier90_types, only: projections_wannier

    use parallel,        only: IONode             ! Input/Output node
    use parallel,        only: Node               ! Local processor number
    use m_mpi_utils,     only: broadcast          ! Broadcasting routines
    use alloc,           only: re_alloc           ! Reallocation routines

    integer, intent(in) :: index_manifold         ! Index of the band manifold
                                                  !   in the Wannier
                                                  !   transformation

!
! Internal variables
!
    integer :: ik             ! Counter for loop on ik points
    integer :: iw             ! Counter for loop on projections
    integer :: nn             ! 
    integer :: i              ! 
    integer :: iband          ! Counter for loop on bands
    integer :: iorb           ! Counter for loop on atomic orbitals
    integer :: ivec           ! Counter for loop on vectors

    if( w90_wrapper_processing ) then
      if(IONode) then
        write(6,'(/a)')                                                     & 
 &        'switch_local_projection: Populating the relevant matrices for '  
        write(6,'(a)')                                                      & 
 &        'switch_local_projection: calling WANNIER90 directly from SIESTA'
        write(6,'(2a)') 'switch_local_projection: band manifold = ', &
            trim(manifold_bands(index_manifold)%name)
        numkpoints = numkpoints_w90_in
!       Initialize the list of k-points
        nullify( kpointsfrac )
        call re_alloc( kpointsfrac, 1, 3, 1, numkpoints,  &
 &                     name='kpointsfrac', routine='switch_local_projection')
        do ik = 1, numkpoints
          kpointsfrac(:,ik) = kpointsfrac_w90_in(:,ik)
        enddo

!       Initialize the list of neighbour k-points
        nullify( nnlist_neig   )
        nullify( nnfolding     )

!       Broadcast information regarding the number of k-points neighbours
!       and allocate in all nodes the corresponding arrays containing 
!       information about k-point neighbours

        nncount = nncount_w90_in

        call re_alloc( nnlist_neig, 1, numkpoints, 1, nncount,      &
 &                     name = "nnlist_neig", routine = "read_nnkp" )
        call re_alloc( nnfolding, 1, 3, 1, numkpoints, 1, nncount,  &
 &                     name = "nnfolding", routine = "read_nnkp" )
        nnlist_neig = nnlist_w90_in
        nnfolding   = nnfolding_w90_in

        numbands(:)   = manifold_bands(index_manifold)%final_band
        numincbands(:)= manifold_bands(index_manifold)%number_of_bands
        blocksizeincbands =  &
 &         manifold_bands(index_manifold)%blocksizeincbands

!       Initialize the list of excluded bands
        nullify( isexcluded )
        call re_alloc( isexcluded, 1, no_u, name='isexcluded', &
 &                     routine='switch_local_projection' )
        isexcluded = manifold_bands(index_manifold)%isexcluded

        latvec = ucell

!       Initialize number of projectors
        numproj = manifold_bands(index_manifold)%numbands
        if( allocated(projections) ) deallocate( projections )
        allocate(projections(numproj))
        projections = manifold_bands(index_manifold)%proj

!       Reciprocal lattice vectors
        call reclat( ucell, reclatvec, 1 )

        nullify( bvectorsfrac )
        call re_alloc( bvectorsfrac, 1, 3, 1, nncount,    &
                       name="bvectorsfrac", routine = "chosing_b_vectors")
        bvectorsfrac = bvectorsfrac_w90_in

        seedname = manifold_bands(index_manifold)%seedname
      endif   ! if (IONode .eq. 0)
      nincbands_loc = manifold_bands(index_manifold)%nincbands_loc
      goto 100
    endif

    if( w90_processing ) then
      ! We need this to reuse the registry index mechanism
      ! introduced in the new interface.
      ! If using the legacy interface, the manifold_bands structure
      ! will not be allocated. We need just a couple of fields
      ! in the "manifold 1" element of the structure       
      if (.not. allocated(manifold_bands)) then
         allocate(manifold_bands(1))
      else
         ! The user is trying to use *both* interfaces. While it might
         ! be possible, since the legacy interface is called at the
         ! end of 'siesta_analysis', it would be confusing.
         call message("WARNING","Manifold_bands array found allocated by legacy w90 interface")
         call die("Trying to use both w90 interfaces at the same time")
      endif
      if(IONode) then
        write(6,'(/a)')                                                     & 
 &        'switch_local_projection: Populating the relevant matrices for '
        write(6,'(a)')                                                      & 
 &        'switch_local_projection: the post-processing interface with Wannier 90'
        numkpoints = numkpoints_wannier

!       Initialize the list of k-points
        nullify( kpointsfrac )
        call re_alloc( kpointsfrac, 1, 3, 1, numkpoints,  &
 &                     name='kpointsfrac', routine='switch_local_projection')
        do ik = 1, numkpoints
          kpointsfrac(:,ik) = kpointsfrac_wannier(:,ik)
        enddo

!       Initialize the list of neighbour k-points
        nullify( nnlist_neig   )
        nullify( nnfolding     )

!       Broadcast information regarding the number of k-points neighbours
!       and allocate in all nodes the corresponding arrays containing 
!       information about k-point neighbours

        nncount = nncount_wannier

        call re_alloc( nnlist_neig, 1, numkpoints, 1, nncount,      &
 &                     name = "nnlist_neig", routine = "read_nnkp" )
        call re_alloc( nnfolding, 1, 3, 1, numkpoints, 1, nncount,  &
 &                     name = "nnfolding", routine = "read_nnkp" )
        nnlist_neig = nnlist_wannier
        nnfolding   = nnfolding_wannier

        numbands          = numbands_wannier
        numincbands       = numincbands_wannier
        blocksizeincbands = blocksizeincbands_wannier

!       Initialize the list of excluded bands
        nullify( isexcluded )
        call re_alloc( isexcluded, 1, no_u, name='isexcluded', &
 &                     routine='switch_local_projection' )
        isexcluded = isexcluded_wannier

        latvec    = latvec_wannier
!       Reciprocal lattice vector
        reclatvec = reclatvec_wannier

        nullify( bvectorsfrac )
        call re_alloc( bvectorsfrac, 1, 3, 1, nncount,    &
                       name="bvectorsfrac", routine = "chosing_b_vectors")
        bvectorsfrac = bvectorsfrac_wannier

        numproj = numproj_wannier
        if( allocated(projections) ) deallocate( projections )
        allocate(projections(numproj))
        projections = projections_wannier

        seedname = seedname_wannier
      endif  ! endif (IONode)
      nincbands_loc     = nincbands_loc_wannier
    endif

100 continue

#ifdef MPI
    call broadcast( seedname    )

    call broadcast( numbands          )
    call broadcast( numincbands       )
    call broadcast( blocksizeincbands )

    call broadcast( latvec      )
    call broadcast( reclatvec   )

    call broadcast( numkpoints  )
    if( Node .ne. 0) then
      call re_alloc( kpointsfrac, 1, 3, 1, numkpoints,                       &
 &                   name='kpointsfrac', routine='switch_local_projection')
    endif
    call broadcast( kpointsfrac )
    call broadcast( nncount     )
    if( Node .ne. 0) then
      call re_alloc( nnlist_neig, 1, numkpoints, 1, nncount,      &
 &                   name = "nnlist_neig", routine = 'switch_local_projection' )
      call re_alloc( nnfolding, 1, 3, 1, numkpoints, 1, nncount,  &
 &                   name = "nnfolding", routine = 'switch_local_projection' )
      call re_alloc( bvectorsfrac, 1, 3, 1, nncount,    &
                     name="bvectorsfrac", routine = 'switch_local_projection' )
    endif
    call broadcast( nnlist_neig  )
    call broadcast( nnfolding    )
    call broadcast( bvectorsfrac )

    if( Node .ne. 0) then
      call re_alloc( isexcluded, 1, no_u, name='isexcluded', &
 &                   routine='switch_local_projection' )
    endif
    call broadcast( isexcluded )

    call broadcast( numproj      )
    call broadcast_projections

#endif

!!   For debugging
!    write(6,'(/a,2i5)')                                         & 
! &    'switch_local_projection: index_manifold    = ', Node, index_manifold
!    write(6,'(a,i5,2x,a)')                                      & 
! &    'switch_local_projection: numkpoints        = ', Node, seedname
!    write(6,'(a,2i5)')                                          & 
! &    'switch_local_projection: numkpoints        = ', Node, numkpoints
!    do ik = 1, numkpoints
!      write(6,'(a,2i5,3f12.5)')                                &
! &      'switch_local_projection: ik, kpointsfrac = ',         &
! &      Node, ik, kpointsfrac(:,ik) 
!    enddo
!    write(6,'(a,3i5)')                                         &
! &    'switch_local_projection: numbands          = ',         &
! &      Node, numbands
!    write(6,'(a,3i5)')                                         &
! &    'switch_local_projection: numincbands       = ',         &
! &     Node, numincbands
!    write(6,'(a,2i5)')                                         &
! &    'switch_local_projection: nincbands_loc     = ',         &
! &     Node, nincbands_loc
!    write(6,'(a,2i5)')                                         &
! &    'switch_local_projection: blocksizeincbands = ',         &
! &     Node, blocksizeincbands
!    do iband = 1, no_u
!      write(6,'(a,2i5,l5)')                                    &
! &      'switch_local_projection: iband, isexcluded = ',       &
! &      Node, iband, isexcluded(iband)
!    enddo
!    do ivec = 1, 3
!      write(6,'(a,2i5,3f12.5)')                                 &
! &      'switch_local_projection: lattice vectors   = ',        &
! &      Node, ivec, latvec(:,ivec)
!    enddo
!    do ivec = 1, 3
!      write(6,'(a,2i5,3f12.5)')                                 &
! &      'switch_local_projection: reciprocal lattice= ',        &
! &      Node, ivec, reclatvec(:,ivec)
!    enddo
!    write(6,'(a,i5)') 'switch_local_projection: projections in Node: ', Node
!    do iw = 1, numproj
!       write(6,'(i5,3f12.6,3i3,f12.6)')   &
! &       Node,                            &
! &       projections(iw)%center(1:3),     &
! &       projections(iw)%l,               &
! &       projections(iw)%mr,              &
! &       projections(iw)%r,               &
! &       projections(iw)%zovera
!    enddo
!
!    write(6,'(a,i5)') 'switch_local_projection: nnkpts in Node: ', Node
!    write(6,'(2i4)') Node, nncount
!    do ik = 1,numkpoints
!      do nn = 1, nncount
!        write(6,'(3i6,3x,3i4)') &
! &        Node,ik,nnlist_neig(ik,nn),(nnfolding(i,ik,nn),i=1,3)
!      end do
!    end do
!    write(6,'(a,i5)') 'switch_local_projection: bvectorsfrac in Node: ', Node
!    do nn = 1, nncount
!      write(6,'(2i6,3x,3f12.5)') &
! &      Node, nn,(bvectorsfrac(i,nn),i=1,3)
!    enddo
!    write(6,'(a,i5/)')                                                        &
! &    'switch_local_projection: end bvectorsfrac in Node: ', Node
!      
!!   End debugging

  end subroutine switch_local_projection

end module m_switch_local_projection
