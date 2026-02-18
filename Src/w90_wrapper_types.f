! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors

!> \brief General purpose of the w90_wrapper types module:
!! defines data structures to handle the specification
!!     of the Wannier90 transformation within SIESTA. 
!!
!! This specification
!! is read from the fdf file by routine 'read_w90_in_specs'.
!!
!! Here we define the derived type w90_in_manifold_t, where all the
!! relevant information for a manifold of bands is stored.
!!
!! Then, we also define some relevant information regarding the k-point
!! sampling, the atomic structure, and eventually,
!! the chemical potential that can be applied to a particular
!! Wannier function

      module w90_wrapper_types
!
!

      use precision, only: dp          ! Double precision
      use files,     only: label_length! Number of characters in slabel
      use trialorbitalclass            ! Determines all the variables
                                       !   that define the localized
                                       !   functions used as first guess
                                       !   during the minimization of the
                                       !   spreading

      implicit none

      type, public ::  w90_manifold_t
          character(64) :: name
          character(label_length+32) :: seedname
                                       ! Seed of the name of the file 
                                       !   that will be used to dump the 
                                       !   information of the matrix elements
                                       !   produced by the WANNIER90 subroutines
          integer                :: initial_band
                                       !   Initial band of the manifold
                                       !     to be wannierized
          integer                :: final_band
                                       !   Final band of the manifold
                                       !     to be wannierized
          integer                :: number_of_bands
                                       ! Number of bands within the manifold
                                       !     passed to WANNIER90
          integer                :: numbands
                                       ! Number of bands that will be 
                                       !     orthonormalized.
                                       !     It might be different from the 
                                       !     number of bands in the manifold.
                                       !     In this case, a disentanglement
                                       !     procedure is required
          integer                :: num_iter = 500
                                       ! Number of iterations for th
                                       !     minimization of \Omega
          integer, pointer       :: orbital_indices(:) => null()
                                       ! Indices of the orbitals that will be
                                       !     used as localized functions in the
                                       !     initial guess of the minimization
                                       !     of the spreading
          integer, allocatable   :: projector_gindex(:)   
          logical                :: projectors_are_registered = .false.
                                       ! Registry indexes for localized functions
          logical, pointer       :: isexcluded(:) => null()
                                       ! List of bands excluded for
                                       !     Wannier transformation
          integer, pointer       :: isincluded(:) => null()
                                       ! List of bands included for
                                       !     Wannier transformation
          logical, pointer       :: orbexcluded(:) => null()
                                       ! List of orbitals excluded 
                                       !     for Wannier transformation
          integer, pointer       :: orb_in_manifold(:) => null()
                                       ! Sequential index of the orbitals that
                                       !     will be used in the Wannier 
                                       !     transformation
          integer                :: blocksizeincbands
                                       ! Blocking size for orbital
                                       !     distribution across the nodes
          integer                :: nincbands_loc
                                       ! Number of included bands 
                                       !     in the calculation 
                                       !     of the overlap and 
                                       !     projection matrices
                                       !     in the local node
          real(dp)               :: dis_win(2) = (/-1._dp, 0._dp/)
                                       ! Bottom/Top of the outer energy
                                       !     window for band disentanglement
                                       !     Units: energy
                                       !     Read from the fdf input
                                       !     in the 
                                       !     %block WannierProjections
          real(dp)               :: dis_froz(2) = (/-1._dp, 0._dp/)
                                       ! Bottom/Top of the inner (frozen) energy
                                       !     window for band disentanglement
                                       !     Units: energy
                                       !     Read from the fdf input
                                       !     in the 
                                       !     %block WannierProjections
          logical                :: frozen_states
                                       ! Flag that determines whether an inner
                                       !     energy window has been defined
          logical                :: disentanglement  
                                       ! Is the disentanglement 
                                       !   procedure required for this manifold?
          logical                :: dis_win_siesta
                                       ! Flag that determines whether the
                                       !    outer window will be read from the
                                       !    fdf file
          logical                :: dis_win_froz_siesta
                                       ! Flag that determines whether the
                                       !    inner frozen window will be
                                       !    read from the fdf file
          logical                :: wannier_plot = .false.
                                       ! Plot the Wannier functions?
          integer              :: wannier_plot_supercell(3) = (/1,1,1/)
                                       ! Size of the supercell for 
                                       !   plotting the WF
          logical                :: fermi_surface_plot = .false.
                                       ! Plot the Fermi surface?
          logical                :: write_hr = .false.
                                       ! Write the Hamiltonian in 
                                       !   the WF basis?
          logical                :: write_tb = .false.
                                       ! Write lattice vectors, 
                                       !   Hamiltonian, and position
                                       !   operator in WF basis
          logical                :: write_unk = .false.
                                       ! Write the UNK matrices
          type(trialorbital), allocatable :: proj(:)
                                       ! Initial guess for the localized 
                                       !   functions
          real(dp)               :: threshold = 1.d-6
                                 ! Threshold for the real part of the 
                                 !   coefficients of a Wannier in a basis 
                                 !   of NAO to compute the contribution to 
                                 !   the tight-binding matrix elements

          ! Unitary matrices coming out of wannier90
          complex(dp), allocatable :: u_matrix(:,:,:)
          complex(dp), allocatable :: u_matrix_opt(:,:,:)
          
      end type w90_manifold_t

      type(w90_manifold_t), public,
     .     allocatable, save, target     :: manifold_bands(:)

!
! Variables related with the neighbours of the k-points
!

      integer :: numkpoints_w90_in     ! Total number of k-points
                                       !   used in Wannier90

      real(dp), pointer :: kpointsfrac_w90_in(:,:)
                                       ! List of k points relative
                                       !   to the reciprocal lattice vectors.
                                       !   First  index: component
                                       !   Second index: k-point index in 
                                       !      the list
      integer           :: nncount_w90_in
                                       ! The number of nearest neighbours
                                       !      belonging to each k-point of the 
                                       !      Monkhorst-Pack mesh
      integer, allocatable  :: nnlist_w90_in(:,:)
                                       ! nnlist(ikp,inn) is the index of the
                                       !   inn-neighbour of ikp-point
                                       !   in the Monkhorst-Pack grid folded 
                                       !   to the first Brillouin zone
      integer, allocatable  :: nnfolding_w90_in(:,:,:)
                                       ! nnfolding(i,ikp,inn) is the i-component
                                       !   of the reciprocal lattice vector
                                       !   (in reduced units) that brings
                                       !   the inn-neighbour specified in nnlist
                                       !   (which is in the first BZ) to the
                                       !   actual \vec{k} + \vec{b} that we need
      real(dp), pointer :: bvectorsfrac_w90_in(:,:)
                                       ! The vectors b that connect
                                       !   each mesh-point k
                                       !   to its nearest neighbours
                                       !   In reciprocal lattice units. 
      integer           :: kmesh_w90_in(3)    
                                       ! Number of divisions along the three
                                       !   reciprocal lattice vectors that
                                       !   will be used in Wannier90

!
! Variables related with the atomic structure
!
      real(dp) :: latvec_w90_in(3,3)   ! Lattice vectors (in Ang???)
      real(dp) :: reclatvec_w90_in(3,3)! Reciprocal lattice vectors
                                       !  Cartesian coordinates in Ang^-1 ???
                                       !  First  index: component
                                       !  Second index: vector
!
! Variables related with the expansion of a Wannier in a basis of NAO
!
      complex(dp), pointer :: coeffs_wan_nao(:,:,:) 
                                       ! Coefficients of the
                                       !   Wannier functions in a basis
                                       !   of NAO
                                       !   First  index: Index of the 
                                       !       manifold and Wannier function,
                                       !       handled by numh_man_proj, 
                                       !       listhptr_man_proj, and
                                       !       listh_man_proj, and
                                       !   Second index: NAO in the
                                       !       supercell
                                       !   Third index: Spin

!
! Variables related with the computation of the matrix elements
! including the chemical potential for some Wannier functions
!
      logical, public :: wannier_coefficients_needed = .false.

      logical, public :: compute_chempotwann = .false.
                                       ! Compute the Hamiltonian matrix
                                       !   elements between NAO
                                       !   if a chemical potential in 
                                       !   a Wannier is applied
      logical, public :: compute_chempotwann_after_scf = .false.
                                       ! If compute_chempotwann = .true.
                                       !   this switch determines when
                                       !   the new matrix elements are
                                       !   added (only after regular scf)
      logical, public :: first_chempotwann   = .false.
                                       !  First time the calculation of the 
                                       !   matrix elements of the Hamiltonian
                                       !   with the chemical potential of the
                                       !   Wanniers is called?
      real(dp), pointer :: chempotwann_val(:,:)
                                       ! Chemical potential 
                                       !   applied to shift the energy of a 
                                       !   the matrix elements in real space
                                       !   associated with a given 
                                       !   Wannier function
                                       !   First index: Wannier function
                                       !   Second index: Spin component 
      integer, pointer :: numh_man_proj(:)     => null()
                                       ! Number of projections that will be
                                       !   handled in the local node
                                       !   for a given manifolds
                                       !   Dimension: number of manifolds
      integer, pointer :: listhptr_man_proj(:) => null()
                                       ! Index pointer to listh_man_proj such 
                                       ! listh_man_proj(listhptr_man_proj(1)+1)
                                       !   is the first projector of the first
                                       !   manifold handled by the local node 
                                       ! listh_man_proj(listhptr_man_proj(io)+1)
                                       !   is thus the first projector of 
                                       !   of manifold 'io' while 
                                       ! listh_man_proj(listhptr_man_proj(io) + 
                                       !                numh_man_proj(io)) 
                                       !   is the last projectors of manifold 
                                       !   'io'.
                                       ! Dimension: number of manifolds
      integer, pointer :: listh_man_proj(:)    => null()
                                       ! The column indices for the projectors
                                       !   of all the manifolds handled by 
                                       !   the local node


      end module w90_wrapper_types

