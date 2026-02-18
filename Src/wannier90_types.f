! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors
!
      module wannier90_types
!
!=======================================================================
!
!     This module defines data structures to handle the specification
!     of the Wannier transformation. This specification
!     is read from the fdf file by routine 'read_w90_in_specs'.

      use precision, only: dp
      use files,     only : label_length         ! Number of characters in slabel
      use trialorbitalclass

      implicit none


!
! Variables related with the atomic structure of the system
!

      real(dp) :: latvec_wannier(3,3)    ! Real lattice vectors.
                                 !   Cartesian coordinates.
                                 !   Readed in Angstroms and transformed to Bohr
                                 !   internally
                                 !   First  index: component
                                 !   Second index: vector
                                 !   This is consistent with the unit cell read
                                 !   in Siesta, but it has changed with respect
                                 !   the first version implemented by R. Korytar
      real(dp) :: reclatvec_wannier(3,3) ! Reciprocal lattice vectors.
                                         !   Cartesian coordinates.
                                         !   Readed in Angstroms^-1 
                                         !   and transformed
                                         !   to Bohr^-1 internally
                                         !   First  index: component
                                         !   Second index: vector
                                         !   This is consistent with the 
                                         !     reciprocal lattice read in Siesta
                                         !     but it has changed with respect
                                         !     the first version implemented 
                                         !     by R. Korytar


!
! Variables related with the k-point list for which the overlap
! matrices Mmn between a k-point and its neighbor will be computed
!
      integer  :: numkpoints_wannier ! Total number of k-points
                                     !   for which the overlap of the
                                     !   periodic part of the wavefunct
                                     !   with a neighbour k-point will
                                     !   be computed
      real(dp), pointer :: kpointsfrac_wannier(:,:)
                                     ! List of k points relative
                                     !   to the reciprocal lattice vectors.
                                     !   First  index: component
                                     !   Second index: k-point index in the list

!
! Variables related with the neighbours of the k-points
!
      integer           :: nncount_wannier
                             ! The number of nearest neighbours belonging to
                             !   each k-point of the Monkhorst-Pack mesh
      integer, pointer  :: nnlist_wannier(:,:)
                             ! nnlist(ikp,inn) is the index of the
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the
                             !   first Brillouin zone
      integer, pointer  :: nnfolding_wannier(:,:,:)
                             ! nnfolding(i,ikp,inn) is the i-component
                             !   of the reciprocal lattice vector
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units.
!
! Variables related with the number of bands considered for wannierization
!
      integer :: numbands_wannier(2)    ! Number of bands for wannierization
                                        !    before excluding bands
      integer :: numincbands_wannier(2) ! Number of included bands in the calc.
                                        !   of the overlap and projection 
                                        !   matrices after excluding the bands
      integer :: nincbands_loc_wannier  ! Number of bands for wannierizatio
                                        !   after excluding bands
                                        !   in the local node
      integer :: blocksizeincbands_wannier
                                        ! Maximum number of bands
                                        !   considered for wannierization
                                        !   per node
      logical, pointer :: isexcluded_wannier(:) 
                                        ! Masks excluded bands

      real(dp), pointer :: bvectorsfrac_wannier(:,:)
                                         !! The vectors b that connect 
                                         !!   each mesh-point k
                                         !!   to its nearest neighbours

!
! Variables related with the projections with trial functions,
! initial approximations to the MLWF
!
      integer  :: numproj_wannier        !! Total number of projection centers,
                                         !!   equal to the number of MLWF
      type(trialorbital), target, allocatable  :: projections_wannier(:)


!
! Variables related with the input/output files
!
      character(label_length+3)  :: seedname_wannier  ! Name of the file where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, reads or dumps the
                                         !   information.




      end module wannier90_types
