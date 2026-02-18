! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_siesta2wannier90

! OUTPUT: 
! File called seedname.mmn, where the overlap matrices are written in 
! the format required by Wannier90
!
! File called seedname.eig, where the eigenvalues of the Hamiltonian
! at the required k-points are written according to the format required
! by Wannier90.

  use precision, only : dp                   ! Real double precision type
  use parallel,  only : IOnode               ! Input/output node
  use parallel,  only : Node                 ! This process node
  use parallel,  only : Nodes                ! Total number of processor nodes
  use parallel,  only : BlockSize            ! Total number of processor nodes
  use sys,       only : die                  ! Termination routine
  use atomlist,  only: no_l                  ! Number of orbitals in local node
                                             ! NOTE: When running in parallel,
                                             !   this is core dependent
                                             !   Sum_{cores} no_l = no_u
  use atomlist,  only: no_u                  ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use files,     only : label_length         ! Number of characters in slabel
  use siesta_options, only: w90_write_mmn    ! Write the Mmn matrix for the
                                             !   interface with Wannier
  use siesta_options, only: w90_write_amn    ! Write the Amn matrix for the
                                             !   interface with Wannier
  use siesta_options, only: w90_write_eig    ! Write the eigenvalues for the
                                             !   interface with Wannier
  use siesta_options, only: w90_write_unk    ! Write the unks for the
                                             !   interface with Wannier
  use wannier90_types, only: latvec_wannier  ! Lattice vectors
  use wannier90_types, only: reclatvec_wannier  ! Reciprocal lattice vectors.
                                                !   Cartesian coordinates.
                                                !   Readed in Angstroms^-1 
                                                !   and transformed
                                                !   to Bohr^-1 internally
                                                !   First  index: component
                                                !   Second index: vector
                                                !   This is consistent with the 
                                                !   reciprocal lattice read in 
                                                !   siesta but it has changed 
                                                !   with respect to the 
                                                !   first version implemented 
                                                !   by R. Korytar
  use wannier90_types, only: numkpoints_wannier ! Total number of k-points
                                                !  for which the overlap of the
                                                !  periodic part of the wavefun
                                                !  with a neighbour k-point will
                                                !  be computed
  use wannier90_types, only: kpointsfrac_wannier! List of k points relative
                                                !   to the reciprocal lattice 
                                                !vectors.
                                                !   First  index: component
                                                !   Second index: k-point index 
                                                !                 in the list
  use wannier90_types, only: nncount_wannier
                             ! The number of nearest neighbours belonging to
                             !   each k-point of the Monkhorst-Pack mesh
  use wannier90_types, only: nnlist_wannier
                             ! nnlist(ikp,inn) is the index of the
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the
                             !   first Brillouin zone
  use wannier90_types, only: nnfolding_wannier
                             ! nnfolding(i,ikp,inn) is the i-component
                             !   of the reciprocal lattice vector
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units.
  use wannier90_types, only: nincbands_loc_wannier  
                             ! Number of bands in the local
                             !   node for wannierization 
                             !   after excluding bands
  use wannier90_types, only: blocksizeincbands_wannier
                                                ! Maximum number of bands
                                                !   considered for 
                                                !   wannierization per node
  use wannier90_types, only: bvectorsfrac_wannier
  use wannier90_types, only: numproj_wannier    ! Total number of projection centers,
                                                !   equal to the number of MLWF



  use wannier90_types, only: numbands_wannier
  use wannier90_types, only: numincbands_wannier
  use wannier90_types, only: isexcluded_wannier
  use wannier90_types, only: seedname_wannier
  use wannier90_types, only: projections_wannier
  use m_switch_local_projection, only: switch_local_projection
  use m_switch_local_projection, only: numbands
  use m_switch_local_projection, only: numincbands
  use m_switch_local_projection, only: isexcluded

#ifdef MPI
  use parallelsubs,       only : set_blocksizedefault
#endif

!
! Variables related with the number of bands considered for wannierization
!
  integer          :: numexcluded 
                             ! Number of bands to exclude from the calculation
                             !   of the overlap and projection matrices.
                             !   This variable is read from the .nnkp file

  integer, pointer :: excludedbands(:)
                             ! Bands to be excluded
                             !   This variable is read from the .nnkp file
  integer, pointer :: isincluded(:) ! Masks included bands




CONTAINS

subroutine siesta2wannier90

  use m_spin,        only: nspin         ! Number of spin components
  use files,         only: slabel        ! Short system label, 
                                         !   used to generate file names
  use m_digest_nnkp, only: read_nnkp     ! Subroutine that reads the .nnkp file
  use m_digest_nnkp, only: chosing_b_vectors ! Subroutine that computes the b
                                         ! vectors that connect each mesh 
                                         ! k-point to its nearest neighbours.
  use m_digest_nnkp, only: set_excluded_bands   ! Subroutine that chooses the 
                                                !   bands that are excluded from
                                                !   the wannierization procedure
  use m_digest_nnkp, only: number_bands_wannier ! Subroutine that computes the
                                         ! number of bands for wannierization

  implicit none

  integer                    :: ispin    ! Spin counter
  integer                    :: ik
  integer                    :: numbandswan(2)
                                         ! Number of bands for wannierization

#ifdef MPI
  integer, external :: numroc
  integer           :: nincbands
#endif

  external  :: timer, mmn

  call timer("siesta2wannier90",1)

  do ispin = 1, nspin
!   Append _up or _dn to the seedname if spin polarized
    seedname_wannier = trim(getFileNameRoot(ispin,nspin,slabel))
!   A priori, the .nnkp file generated by Wannier90 used as a 
!   postprocessing tool is independent of spin.
!   However, they provide examples (example08, bulk Fe) where
!   two different .nnkp files are generated, 
!   one for each component of the spin.
!   This inspired us to read the .nnkp file for each component of the spin
!   (note that the seedname will be different in every case).
!   The information is read only by the master node,
!   and broadcast to the rest of the nodes later. 
    if (IOnode) then
       write(6,'(/,a)')  &
 &       'siesta2wannier90: Reading the ' // trim(seedname_wannier) // '.nnkp file'
    endif
    call read_nnkp( seedname_wannier, latvec_wannier, reclatvec_wannier, numkpoints_wannier,  &
 &                  kpointsfrac_wannier, nncount_wannier, nnlist_wannier,                     &
 &                  nnfolding_wannier, numproj_wannier, projections_wannier, numexcluded,             &
 &                  excludedbands, w90_write_amn )

!!   For debugging
!    write(6,'(a,i5)') ' numkpoints_wannier = ', numkpoints_wannier
!    do ik = 1, numkpoints_wannier
!      write(6,'(a,i5,3f12.5)') ' ik, kpointsfrac_wannier = ',  &
! &                               ik, kpointsfrac_wannier(:,ik)
!    enddo
!!   End debugging

!   Compute the vectors that connect each mesh k-point to its nearest neighbours
    call chosing_b_vectors( kpointsfrac_wannier, nncount_wannier,  &
 &                          nnlist_wannier, nnfolding_wannier,     &
 &                          bvectorsfrac_wannier )

!   Compute the number of bands for wannierization
    call number_bands_wannier( numbandswan )
    numbands_wannier(ispin) = numbandswan(ispin)

!   Chose which bands are excluded from the wannierization procedure
    call set_excluded_bands(  ispin, numexcluded, excludedbands,      &
 &                            numbands_wannier, isexcluded_wannier,   &
 &                            isincluded, numincbands_wannier )

!   Allocate memory related with the coefficients of the wavefunctions
#ifdef MPI
!   Find the number of included bands for Wannierization that will be stored
!   per node. Use a block-cyclic distribution of nincbands over Nodes.
!
    call set_blocksizedefault( Nodes, numincbands_wannier(ispin),     &
 &                             blocksizeincbands_wannier )

!!    For debugging
!     write(6,'(a,3i5)')' diagonalizeHk: Node, Blocksize = ', &
! &     Node, numincbands_wannier(ispin), blocksizeincbands_wannier
!!    End debugging

     nincbands_loc_wannier = numroc( numincbands_wannier(ispin),      &
  &                          blocksizeincbands_wannier, Node, 0, Nodes )
!    For debugging
     write(6,'(a,3i5)')' diagonalizeHk: Node, nincbands_loc_wannier = ', &
 &     Node, nincbands_loc_wannier
!    End debugging
#else
     nincbands_loc_wannier = numincbands_wannier(ispin)
#endif

!
!   Populate the right arrays and matrices to call the generic routines
!   to compute the projections
!
    call switch_local_projection( 0 )

!   Compute the matrix elements of the plane wave,
!   for all the wave vectors that connect a given k-point to its nearest
!   neighbours
    call compute_pw_matrix( nncount_wannier, bvectorsfrac_wannier )

!   Compute the coefficients of the wavefunctions at the 
!   k-points of the Wannier90 mesh
    call diagonalizeHk( ispin )

!   Compute the overlap between periodic parts of the wave functions
    if( w90_write_mmn ) call mmn( ispin )

!   Compute the overlaps between Bloch states onto trial localized orbitals
    if( w90_write_amn ) call amn( ispin, 1 )

!   Write the eigenvalues in a file, in the format required by Wannier90
    if( IOnode .and. w90_write_eig ) call writeeig( ispin )

!   Write the values of the Bloch states in a box
    if( w90_write_unk ) call writeunk( ispin )

  enddo

  if (IOnode) then
    write(6,'(/,a)')  &
 &    'siesta2wannier90: All the information dumped in the corresponding files'
    write(6,'(a)')  &
 &    'siesta2wannier90: End of the interface between Siesta and Wannier90'
  endif

  call timer("siesta2wannier90",2)

end subroutine siesta2wannier90

function getFileNameRoot(ispin,nspin,root)
!
! Constructs filenames for various cases of spin:
! Nonpolarized, down and up by adding an extension to
! "root". Behavior:
! (a) unpolarized: no change
! (b) spin up: _up suffix
! (c) spin down: _dn suffix
!
  character(*),intent(in)      :: root
  integer,intent(in)           :: ispin, nspin
  character(len_trim(root)+3)  :: getFileNameRoot

  select case (nspin)
!   If only one spin component, do nothing
    case (1) 
      getFileNameRoot = trim(root)
!   If two spin components, 
!   append "_up" (for the first one) or "_dn" (for the second spin component)
!   to the seed name
    case (2) 
      select case (ispin)
         case (1)
           getFileNameRoot = trim(root)//"_up"
         case (2)
           getFileNameRoot = trim(root)//"_dn"
      end select 
!   If more than two spin components, stop the program 
!   and print and error message
    case default 
      call die("getFileNameRoot: non-collinear spin not implemented yet")
  end select 
end function getFileNameRoot

endmodule m_siesta2wannier90
