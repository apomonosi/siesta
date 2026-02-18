
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_w90_wrapper module:
!! interface between SIESTA and WANNIER90 code.
!!
!! In this module we perform a Wannier transformation of the Bloch functions
!! corresponding to a given manifold of bands,
!! following the recipe given in \cite Marzari-97.
!! In order to compute the transformation, SIESTA uses the low-level modules
!! of the WANNIER90 code(\cite Wannier90), version 3.0.0, as building blocks
!! A complete description of the theory behind the maximally localized
!! Wannier functions and its applications can be found in Ref. \cite Marzari-12
!!
!! The tutorial with a full explanation of the new labels in the input file
!! can be found in
!!
!! <https://personales.unican.es/junqueraj/JavierJunquera_files/Metodos/Wannier/Exercise-Wannier90-within-siesta.pdf>
!!
!! Example: one of the Wannier functions obtained after a transformation of the
!! top of the valence band manifold of SrTiO3 (O-2p in character)
!! \image html SrTiO3.manifold.1_00001.png
#include "mpi_macros.f"

module m_w90_wrapper
#ifdef SIESTA__WANNIER90


  use precision,      only: dp            ! Real double precision type
  use siesta_options, only: n_wannier_manifolds
                                          ! Number of bands manifolds
                                          !  that will be considered
                                          !  for Wannier transformation
  use siesta_options, only: w90_index_perturbed_manifold
                                          ! Index of the manifold that will
                                          !   be perturbed with an external
                                          !   chemical potential
  use siesta_options, only: w90_r_between_manifolds
                                          ! Will the position operator matrix
                                          !  elements be computed between bands
                                          !  of different manifolds
  use siesta_options, only: w90_mmn_diagonal
                                          ! Will the Mmn matrix be forced to
                                          !  be diagonal between bands of different
                                          !  manifolds?
  use parallel,       only: Node          ! Local processor number
  use parallel,       only: Nodes         ! Total number of processors in a
                                          !  parallel run
  use parallel,       only: IONode        ! Node for Input/output
  use w90_wrapper_types,   only: w90_manifold_t
  use w90_wrapper_types,   only: manifold_bands
                                       ! Variable where the initial
                                       !   and final band of each
                                       !   manifold are stored

  use w90_wrapper_types,   only: wannier_coefficients_needed

  use w90_wrapper_types,   only: compute_chempotwann
                                       ! Compute the Hamiltonian matrix
                                       !   elements between NAO
                                       !   if a chemical potential in
                                       !   a Wannier is applied
  use w90_wrapper_types,   only: first_chempotwann
                                       !  First time the calculation of the
                                       !   matrix elements of the Hamiltonian
                                       !   with the chemical potential of the
                                       !   Wanniers is called?
  use w90_wrapper_types,   only: chempotwann_val
                                       ! Chemical potential
                                       !   applied to shift the energy of a
                                       !   the matrix elements in real space
                                       !   associated with a given
                                       !   Wannier function
  use w90_wrapper_types,   only: numh_man_proj
                                       ! Number of projections that will be
                                       !   handled in the local node
                                       !   for a given manifolds
                                       !   Dimension: number of manifolds
  use w90_wrapper_types,   only: listhptr_man_proj
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
  use w90_wrapper_types,   only: listh_man_proj
                                       ! The column indices for the projectors
                                       !   of all the manifolds handled by
                                       !   the local node
  use w90_wrapper_types,       only: coeffs_wan_nao
                                       ! Coefficients of the
                                       !   Wannier functions in a basis
                                       !   of NAO
                                       !   First  index: Index of the
                                       !       manifold and projector,
                                       !       handled by numh_man_proj,
                                       !       listhptr_man_proj, and
                                       !       listh_man_proj lists
                                       !   Second index: NAO in the
                                       !       supercell
                                       !   Third index: Spin component
  use atomlist,           only: no_u   ! Number of orbitals in unit cell
  use atomlist,           only: no_u   ! Number of orbitals in unit cell
                                       ! NOTE: When running in parallel,
                                       !   this is core independent
  use atomlist,           only: no_s   ! Number of atomic orbitals in the
                                       !   supercell
  use sys,                only: die    ! Termination routine
  use fdf

!
! Allocation/Deallocation routines
!
  use alloc,              only: re_alloc  ! Reallocation routines
  use alloc,              only: de_alloc  ! Deallocation routines

#ifdef MPI
  use mpi_siesta
  use parallelsubs,       only: LocalToGlobalOrb
#endif

  implicit none


! Routines
  public :: setup_w90_wrapper
  public :: compute_matrices
  public :: compute_wannier
  public :: deallocate_wannier

  private

  MPI_COMM_TYPE :: w90_comm    ! Communicator for wannier90 wrapper

  CONTAINS

!> \brief General purpose of the subroutine setup_w90_wrapper:
!! set up the interface between SIESTA and WANNIER90
!!
!! Within this subroutine:
!!
!! 1. Include in the SystemName.bib file a citation to the papers where
!!    the implementation is based. Here:
!!    the original paper by Marzari and Vanderbilt,
!!    the last paper on the implementation of WANNIER90,
!!    and the pre-print by Íñiguez and Yildirim in the arxiv where the
!!    Lowdin orthogonalization is explained.
!!
!! 2. Read the WannierManifolds block in the fdf
!!    and set up the general information for every manifold,
!!    including initial (initial_band) and final (final_band) bands,
!!    number of bands that will be orthogonalized (numbands_w90_in),
!!    the orbitals that will be used as initial guess in the projections
!!    during the wannierization (orbital_indices),
!!    number of iterations that will be performed in the Wannierization
!!    procedure (num_iter),
!!    the energy windows for the disentanglement,
!!    and other output options (if the wannier functions will be plotted,
!!    or the Fermi surface computed).
!!    This is done in the subroutine read_w90_wrapper_specs.
!!
!! 3. Set up the indexing of the bands that will be excluded from the
!!    wannierization (isexcluded),
!!    the orbitals that will be used as the localized projector functions
!!    (orbexcluded),
!!    and the new indexing of the orbitals within a manifold (orb_in_manifold).
!!    This is done in the subroutine set_excluded_bands_w90_in.
!!
!! 4. Generate the trial localized functions from the atomic orbitals in the
!!    basis set of SIESTA.
!!    This is done in the subroutine define_trial_orbitals.
!!
!! 5. Determine the number of bands that will be treated per node in a
!!    parallel run.
!!
!! 6. Populate the lists required to handle the matrix of coefficients of
!!    a Wannier functions in terms of the numerical atomic orbitals of the
!!    supercells. These matrices are very similar to numh, listhptr, and
!!    listh to handle the Hamiltonian and Overlap sparse matrices
!!
!! 7. Read all the information regarding the
!!    k-point sampling (block Wannier.k),
!!    and the neighbours for all the k-points in the BZ,
!!    required for the wannierization.
!!    This is done in the subroutine read_kpoints_wannier
!!    by all the nodes simultaneously, so no need for broadcasting
!!    these variables.
!!
  subroutine setup_w90_wrapper()

    use m_cite,        only: add_citation          ! Subroutine used to cite
                                                   !   the proper papers where
                                                   !   the implementation is
                                                   !   based on
    use m_spin,        only: spin                  ! Spin configuration
                                                   !   for SIESTA
#ifdef MPI
    use parallelsubs,  only: set_blocksizedefault  ! Subroutine to find
                                                   !   a sensible default value
                                                   !   for the blocksize default
    use parallel,      only: BlockSize             ! Blocking factor used to
                                                   !   divide up the arrays over
                                                   !   the processes for the
                                                   !   Scalapack routines.
#endif

!
! Internal parameters
!

    integer :: index_manifold
    integer :: iproj_local
    integer :: iproj_global
    integer :: index
    integer :: numproj
    integer :: maxnh_man_proj

#ifdef MPI
    integer, external :: numroc                    ! Scalapack routine for
                                                   !  block-cyclic distributions
    integer, allocatable :: blocksizeprojectors(:)
    integer :: blocksize_tmp
    integer :: ierr
#endif

#ifdef MPI
    call MPI_Comm_Dup(mpi_comm_world, w90_comm, ierr)
#endif

    if( IONode ) then

!     Add a citation ...
!     ... to the Marzari and Vanderbilt paper on
!     Maximally Localized Wannier Functions
      call add_citation("10.1103/PhysRevB.56.12847")
!     ... to the Review of Modern Physics on
!     Maximally Localized Wannier Functions
      call add_citation("10.1103/RevModPhys.84.1419")
!!     ... to the paper describing the WANNIER90 code, version 2.x
!      call add_citation("10.1016/j.cpc.2014.05.003")
!!     ... to the paper describing the Löwdin orthogonalization
      call add_citation("10.1103/PhysRevB.71.180415")
!!     ... to the paper describing the WANNIER90 code, version 3.x
      call add_citation("10.1088/1361-648x/ab51ff")

    end if   ! close if( IONode )

!
!   Read and set up the general information for every manifold,
!   including initial (initial_band) and final (final_band) bands,
!   number of bands that will be orthogonalized (numbands_w90_in),
!   number of iterations that will be performed in the Wannierization
!   procedure (num_iter),
!   the orbitals that will be used in the wannierization (orbital_indices),
!   the energy windows for the disentanglement,
!   and other output options (if the wannier functions will be plotted,
!   or the Fermi surface computed).
!
    call read_w90_wrapper_specs

!   NOTE: This should be done only if the coefficients are going to be
!         required later on.

  if (wannier_coefficients_needed) then

!   Compute the lists that will be required to handle the coefficients of
!   the Wannier functions as an expansion of Numerical Atomic Orbitals

    nullify( numh_man_proj )
    call re_alloc( numh_man_proj,  1, n_wannier_manifolds, &
 &                 'numh_man_proj', 'setup_w90_wrapper' )
    numh_man_proj(:) = 0

    nullify( listhptr_man_proj )
    call re_alloc( listhptr_man_proj,  1, n_wannier_manifolds, &
 &                'listhptr_man_proj', 'setup_w90_wrapper' )
    listhptr_man_proj(:) = 0

#ifdef MPI
    allocate(blocksizeprojectors(n_wannier_manifolds))
#endif

    maxnh_man_proj = 0
    do index_manifold = 1, n_wannier_manifolds
      numproj = manifold_bands(index_manifold)%numbands
#ifdef MPI
!     Find the number of projectors that will be stored
!     per node. Use a block-cyclic distribution of numproj over Nodes.
      call set_blocksizedefault( Nodes, numproj,                  &
 &                               blocksizeprojectors(index_manifold) )
      numh_man_proj(index_manifold) = numroc( numproj,           &
 &               blocksizeprojectors(index_manifold), Node, 0, Nodes )
#else
      numh_man_proj(index_manifold) = numproj
#endif
      maxnh_man_proj = maxnh_man_proj + numh_man_proj(index_manifold)
    enddo

    listhptr_man_proj(1) = 0
    do index_manifold = 2, n_wannier_manifolds
      listhptr_man_proj(index_manifold)=listhptr_man_proj(index_manifold-1) +&
 &       numh_man_proj(index_manifold-1)
    enddo

    call re_alloc( listh_man_proj,  1, maxnh_man_proj, &
 &                'listh_man_proj', 'setup_w90_wrapper' )

    index = 0
#ifdef MPI
    blocksize_tmp = BlockSize
#endif
    do index_manifold = 1, n_wannier_manifolds
      do iproj_local = 1, numh_man_proj( index_manifold )
        index = index + 1
#ifdef MPI
        BlockSize = blocksizeprojectors(index_manifold)
        call LocalToGlobalOrb(iproj_local, Node, Nodes, iproj_global)
#else
        iproj_global = iproj_local
#endif
        listh_man_proj(index) = iproj_global
      enddo
    enddo
#ifdef MPI
    BlockSize = blocksize_tmp
#endif

!!   For debugging
!    do index_manifold = 1, n_wannier_manifolds
!       write(6,'(a,4i5)') &
! &       'setup_w90_wrapper: index_manifold, Node,numproj,numh_man_proj=',&
! &       index_manifold, Node, numproj, numh_man_proj(index_manifold)
!       write(6,'(a,4i5)') &
! &       'setup_w90_wrapper: index_manifold, Node, numproj, listhptr_man_proj =',&
! &       index_manifold, Node, numproj, listhptr_man_proj(index_manifold)
!    enddo
!    write(6,'(a,3i5)') &
! &    'setup_w90_wrapper: Node, maxnh_man_proj =',               &
! &                          Node, maxnh_man_proj
!    do index = 1, maxnh_man_proj
!      write(6,'(a,3i5)')                                             &
!        'setup_w90_wrapper: Node, index, listh_man_proj  =',       &
! &                            Node, index, listh_man_proj(index)
!    enddo
!!   End debugging

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( coeffs_wan_nao )
    call re_alloc( coeffs_wan_nao,                                  &
 &                 1, maxnh_man_proj,                               &
 &                 1, no_s,                                         &
 &                 1, spin%H,                                       &
 &                 name='coeffs_wan_nao', routine='wannier_in_nao')
    coeffs_wan_nao = cmplx(0.0_dp,0.0_dp,kind=dp)

#ifdef MPI
    deallocate(blocksizeprojectors)
#endif

 endif   ! wannier_coefficients_needed

!! For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIError)
!#endif
!    call die()
!! End debugging

!
!   Read all the information regarding the k-point sampling,
!   and the neighbours for all the k-points in the BZ
!   required for the wannierization.
!   This is done by all the nodes simultaneously, so no need for broadcasting
!   these variables
!
    call read_kpoints_wannier


  end subroutine setup_w90_wrapper

!
! ------------------------------------------------------------------------------
!


!> \brief  General purpouse of the subroutine read_w90_wrapper_specs:
!!         read all the info in the fdf file related with the
!!         band manifolds that will be treated in the
!!         Wannier transformation.
!!
!!         This information is contained in the
!!         block WannierManifolds.
!!         The derived variable manifold_bands
!!         is populated.
!!         The type of this variable is defined in the
!!         module w90_wrapper_types,
!!         in the file w90_in_types.f.
!!
!!         This subroutine is used only by the Node responsible for the
!!         input/output.

   subroutine read_w90_wrapper_specs()
!
!   Processes the information in an fdf file
!   regarding the bands that will enter into the orthonormalization procedure
!   and the atomic orbitals that will be orthonormalized.
!
    use fdf
    use files,    only: slabel          ! Short system label,
                                        !   used to generate file names
    use m_spin,   only: spin            ! Spin configuration
                                        !   for SIESTA
#ifdef MPI
    use parallelsubs,  only: set_blocksizedefault  ! Subroutine to find
                                                   !   a sensible default value
                                                   !   for the blocksize default
#endif

!   Internal variables
    integer :: index_manifold           ! Counter for the number of manifolds
    integer :: index_wannier            ! Variable to identify a given Wannier
    real(dp):: factor                   ! Conversion factor for energy units
    character(len=64) :: string_manifold

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

!   Allocate the pointer where all the data required for every manifold
!   will be stored

    allocate(manifold_bands(n_wannier_manifolds))

    ! Global default: threshold for the real part of the
    ! coefficients of a Wannier in a basis
    ! of NAO to compute the contribution to
    ! the tight-binding matrix elements
    manifold_bands(:)%threshold = &
        fdf_get('Wannier.Manifolds.Threshold', 1.d-6)

    ! Global default: Shall SIESTA compute the files which contain
    ! the periodic part of a Bloch function
    ! in the unit cell on a grid?
    manifold_bands(:)%write_unk = &
        fdf_get('Wannier.Manifolds.Unk', .false.)


!   Read the WannierManifolds block
!   If it is not present, do nothing
    if (.not. fdf_block('Wannier.Manifolds',bfdf)) RETURN
    !
    index_manifold = 0
    do while(fdf_bline(bfdf, pline))
      index_manifold = index_manifold + 1
      ! Read in the manifold name ...
      string_manifold = fdf_bnames(pline, 1)
      ! ... and process the corresponding named block
      call read_manifold(trim(string_manifold), &
          manifold_bands(index_manifold))
    end do

    ! Pending specific info in the manifold blocks
    wannier_coefficients_needed = fdf_get('Wannier.Compute.Coeffs',.false.)

!   Read the chemical potential associated with a given Wannier function
!   First check whether the block is present in the fdf file.
!   If it is not present, do nothing
    if (.not. fdf_block('Wannier.ChemicalPotential',bfdf)) RETURN

!   --------------------------- Wannier 'Chemical Potential' shifts
!   Allocate the pointer where the chemical potential associated with
!   a given Wannier function will be stored.
!   The correction will be applied only on the bands associated with
!   the first manifold, that is why we are using this number of bands
!   to allocate the variable.
    nullify( chempotwann_val )
    call re_alloc( chempotwann_val,                              &
 &                 1, manifold_bands(w90_index_perturbed_manifold)%numbands,  &
 &                 1, spin%H,                                    &
 &                 name='chempotwann_val',                       &
 &                 routine='read_w90_wrapper_specs')
    chempotwann_val = 0.0_dp


!   If the block with the chemical potentials is present,
!   set to true the flag to compute the shifts of the matrix elements
    compute_chempotwann = .true.
    first_chempotwann   = .true.
    wannier_coefficients_needed = .true.

!   Read the content of the block, line by line
    do while(fdf_bline(bfdf, pline))
      if( spin%H .eq. 1) then
        if (.not. fdf_bmatch(pline,'ivn')) &   ! We expect that the first line
 &        call die('Wrong format in Wannier.ChemicalPotential')
      elseif( spin%H .eq. 2) then
        if (.not. fdf_bmatch(pline,'ivvn')) &  ! We expect that the first line
 &        call die('Wrong format in ChemicalPotentialWannier')
      endif
      factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
      index_wannier = fdf_bintegers(pline,1)
      if( spin%H .eq. 1) then
        chempotwann_val(index_wannier,1) = fdf_breals(pline,1) * factor
      else if( spin%H .eq. 2) then
        chempotwann_val(index_wannier,1) = fdf_breals(pline,1) * factor
        chempotwann_val(index_wannier,2) = fdf_breals(pline,2) * factor
      endif
    enddo ! end loop over all the lines in the block ChemicalPotentialWannier

!!   For debugging
!    do index_wannier = 1, manifold_bands(1)%numbands
!      write(6,'(a,i5,a,2f12.5)')                                              &
! &      'read_w90_wrapper_specs: Chemical potential for Wannier number ',  &
! &      index_wannier, ' = ' , chempotwann_val(index_wannier,:)
!    enddo
!!   End debugging

  contains

    subroutine read_manifold(name, w90man)
      character(len=*), intent(in) :: name
      type(w90_manifold_t), intent(inout) :: w90man

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer :: number_unit_cells, number_wannier_unit_cell
      character(len=64) :: key, val
      logical :: key_found(100)
      integer :: ni, nl, i, j
      integer, allocatable :: orbs(:), tmp(:)
      real(dp) :: conv
      integer :: blocksizeincbands_tmp
      integer :: numincbands_tmp
      integer :: nincbands_loc_tmp
#ifdef MPI
      integer, external :: numroc                    ! Scalapack routine for
                                                     !  block-cyclic distributions
#endif

      if ( .not. fdf_block('Wannier.Manifold.'//trim(name), bfdf) ) then
        call die("WannierManifold."//trim(name)//" could not be found &
            &, this name was found in WannierManifolds and should exist!")
      end if

      w90man%name = trim(name)

      number_unit_cells = 1
      number_wannier_unit_cell = 1
      key_found(:) = .false.

      do while( fdf_bline(bfdf, pline) )

        ! Retrieve the key
        key = fdf_bnames(pline,1)

        if ( leqi(key, "bands") ) then
          key_found(1) = .true.

          ! Read inital and final bands
          w90man%initial_band = fdf_bintegers(pline, 1)
          w90man%final_band = fdf_bintegers(pline, 2)
          w90man%number_of_bands = w90man%final_band - w90man%initial_band + 1

          if ( w90man%number_of_bands < 0 ) then
            call die('W90[bands]: The final band should be larger than the initial one')
          end if

        else if ( leqi(key, "trial-orbitals") ) then
          key_found(2) = .true.

          ! Now we have to look for every unit-cell
          ni = fdf_bnintegers(pline)
          ! First read integers
          if ( ni > 0 ) then
            allocate(tmp(ni))
            do i = 1, ni
              tmp(i) = fdf_bintegers(pline, i)
            end do
            call int_array_extend(orbs, tmp)
            deallocate(tmp)
          end if

          nl = fdf_bnlists(pline)
          do i = 1, nl
            j = -1
            allocate(tmp(1))
            ! Query number of elements in the list
            call fdf_bilists(pline, i, j, tmp)
            deallocate(tmp)
            ! In case there are no entries, simply skip
            if ( j <= 0 ) cycle

            allocate(tmp(j))
            call fdf_bilists(pline, i, j, tmp)
            call int_array_extend(orbs, tmp)
            deallocate(tmp)
          end do

        else if ( leqi(key, "spreading.nitt") ) then
          key_found(3) = .true.

          w90man%num_iter = fdf_bintegers(pline, 1)

        else if ( leqi(key, "wannier_plot") .or. leqi(key, "wannier-plot") ) then
          key_found(4) = .true.

          w90man%wannier_plot = .true.
          if ( fdf_bnintegers(pline) == 3 ) then
            do i = 1, 3
              w90man%wannier_plot_supercell(i) = &
                  fdf_bintegers(pline, i)
            end do
          else
            w90man%wannier_plot_supercell(:) = fdf_bintegers(pline, 1)
          end if

        else if ( leqi(key, "fermi_surface_plot") .or. leqi(key, "fermi-surface-plot") ) then
          key_found(5) = .true.

          if ( fdf_bnnames(pline) == 2 ) then
            w90man%fermi_surface_plot = fdf_bboolean(pline, 1, after=1)
          else
            w90man%fermi_surface_plot = .true.
          end if

        else if ( leqi(key, "write_hr") .or. leqi(key,"write-hr") ) then
          key_found(6) = .true.

          if ( fdf_bnnames(pline) == 2 ) then
            w90man%write_hr = fdf_bboolean(pline, 1, after=1)
          else
            w90man%write_hr = .true.
          end if

        else if ( leqi(key, "write_tb") .or. leqi(key,"write-tb") ) then
          key_found(7) = .true.

          if ( fdf_bnnames(pline) == 2 ) then
            w90man%write_tb = fdf_bboolean(pline, 1, after=1)
          else
            w90man%write_tb = .true.
          end if

        else if ( leqi(key, "window") ) then
          key_found(8) = .true.

          val = fdf_bnames(pline, 2)
          conv = fdf_convfac(val, "eV")
          w90man%dis_win(1) = fdf_bvalues(pline, 1) * conv
          w90man%dis_win(2) = fdf_bvalues(pline, 2) * conv
          w90man%dis_win_siesta = .true.

        else if ( leqi(key, "window.frozen") ) then
          key_found(9) = .true.

          val = fdf_bnames(pline, 2)
          conv = fdf_convfac(val, "eV")
          w90man%dis_froz(1) = fdf_bvalues(pline, 1) * conv
          w90man%dis_froz(2) = fdf_bvalues(pline, 2) * conv
          w90man%frozen_states = .true.
          w90man%dis_win_froz_siesta = .true.

        else if ( leqi(key, "threshold") ) then
          key_found(10) = .true.

          w90man%threshold = fdf_bvalues(pline, 1)

        else if ( leqi(key, "write_unk") .or. leqi(key,"write-unk") ) then
          key_found(11) = .true.

          if ( fdf_bnnames(pline) == 2 ) then
            w90man%write_unk = fdf_bboolean(pline, 1, after=1)
          else
            w90man%write_unk = .true.
          end if

        else

          call die("W90[...]: Unknown key: "//trim(key)//" perhaps you mispelled?")

        end if

      end do

      if ( .not. w90man%write_unk ) then
        w90man%wannier_plot = .false.
      end if

      ! Calculate the true size of the number of w90 functions
      if ( .not. allocated(orbs) ) then
        call die("W90[trial-orbitals]: No trial-orbitals were found, please at least have 1 in the manifold!")
      end if
      if ( size(orbs) <= 0 ) then
        call die("W90[trial-orbitals]: No trial-orbitals were found, please at least have 1 in the manifold!")
      end if
      w90man%numbands = size(orbs)
      call re_alloc(w90man%orbital_indices, 1, w90man%numbands)
      w90man%orbital_indices(:) = orbs(:)
      deallocate(orbs)

      if ( w90man%numbands > w90man%number_of_bands ) then
        call die("W90[w90-functions/bands]: More Wannier functions requested than bands in the manifold")
      end if

      ! Determine whether disentanglement would be required
      w90man%disentanglement = w90man%number_of_bands /= w90man%numbands

      ! Update seedname for wannier90
      write(w90man%seedname, '(3a)') trim(slabel), ".manifold.", trim(name)

      ! Check option logic
      if ( w90man%disentanglement ) then
        if ( .not. (key_found(8) .and. key_found(9)) ) then
          call die("W90[window+window.frozen]: Disentanglement windows *must* be present.")
        end if
      end if

      if ( .not. key_found(1) ) then
        call die("W90[bands]: No bands are defined, must be present")
      end if

      if ( .not. key_found(2) ) then
        call die("W90[trial-orbitals]: No trial-orbitals are defined, must be present")
      end if

!    For debugging
!      write(6,'(/,a,i2)')                                                      &
! &      'read_w90_wrapper_specs: Dumping the information for the manifold number ', &
! &       index_manifold
!
!!     Write the seed of the name of the file to dump the output
!      write(6,'(a,i1,1x,a)')                                                   &
! &      'read_w90_wrapper_specs: Seed of the file to dump the output for manifold:',&
! &       index_manifold,                                                       &
! &       w90man%seedname
!
!!     Write the initial and the final band of the manifold
!      write(6,'(a,2i5)')                                                      &
! &      'read_w90_wrapper_specs: Initial band, Final band              = ', &
! &      w90man%initial_band,                   &
! &      w90man%final_band
!
!!     Write the total number of bands in the manifold
!      write(6,'(a,i5)')                                                       &
! &      'read_w90_wrapper_specs: Total number of bands in the manifold = ', &
! &      w90man%number_of_bands
!
!!     Write the total number of bands to be orthogonalized
!      write(6,'(a,i5)')                                                        &
! &'read_w90_wrapper_specs: Total number of bands in the Wannier transformation =', &
! &      w90man%numbands
!
!      if(w90man%disentanglement) then
!         write(6,'(a)')                                                  &
! &        'read_w90_wrapper_specs: Number of bands in the manifold is different of'
!         write(6,'(a)')                                                  &
! &        'read_w90_wrapper_specs: the number of bands to be orthogonalized.'
!         write(6,'(a)')                                                  &
! &        'read_w90_wrapper_specs: Disentanglement procedure required '
!      endif
!
!     Write the indices of the orbitals used as initial localized guess
!      do iorb = 1, w90man%numbands
!        write(6,'(a,i5)')                                                      &
! & 'read_w90_wrapper_specs: Index of the orbital used as initial localized guess =',&
! &      w90man%orbital_indices(iorb)
!      enddo
!
!!     Write the values of the outer energy window for band disentanglement
!      write(6,'(a,l5)')                                             &
! &      'read_w90_wrapper_specs: Disentanglement                    = ', &
! &      w90man%disentanglement
!      write(6,'(a,i5,2f12.5)')                                      &
! &      'read_w90_wrapper_specs: Outer energy window  (eV)          = ', &
! &      index_manifold,                                             &
! &      w90man%dis_win
!
!!     Write the values of the inner (frozen) energy window for
!!     band disentanglement
!      write(6,'(a,i5,2f12.5,l5)')                                   &
! &      'read_w90_wrapper_specs: Inner (frozen) energy window  (eV) = ', &
! &      index_manifold,                                             &
! &      w90man%dis_froz_min,         &
! &      w90man%dis_froz_max,         &
! &      w90man%frozen_states         &
!
!!     Write the values of the inner (frozen) energy window for
!!     band disentanglement
!      write(6,'(a,i5)')                                                    &
! &      'read_w90_wrapper_specs: Number of iterations for the minimization = ', &
! &      w90man%num_iter
!
!!     Write the values of the inner (frozen) energy window for
!!     band disentanglement
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_wrapper_specs: Plot the Wannier functions?               = ', &
! &      w90man%wannier_plot
!      write(6,'(a,i5)')                                                    &
! &      'read_w90_wrapper_specs: Plot the Wannier functions?               = ', &
! &      w90man%wannier_plot_supercell(1)
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_wrapper_specs: Plot the Fermi surface?              ', &
! &      w90man%fermi_surface_plot
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_wrapper_specs: Write the Hamiltonian in the basis of WF? = ', &
! &      w90man%write_hr
!      write(6,'(a,l5)')                                                    &
! &      'read_w90_wrapper_specs: Write the tight-binding elements in the basis of WF? = ', &
! &      w90man%write_tb
!!     End debugging


      ! Set up the indexing of the bands that will be excluded from the
      ! wannierization (isexcluded),
      ! the orbitals that will be used as the localized projector functions
      ! (orbexcluded) ,
      ! and the new indexing of the orbitals within a manifold (orb_in_manifold)
      call set_excluded_bands_w90_in(w90man)

      ! Generate the trial localized functions from the atomic orbitals in the
      ! basis set of SIESTA
      call define_trial_orbitals(name, w90man )

      ! Figure out the block-sizes for the manifold
      numincbands_tmp = w90man%number_of_bands
      ! Allocate memory related with the coefficients of the wavefunctions
#ifdef MPI
      ! Find the number of included bands for Wannierization that will be stored
      ! per node. Use a block-cyclic distribution of nincbands over Nodes.

      call set_blocksizedefault( Nodes, numincbands_tmp, blocksizeincbands_tmp )

      nincbands_loc_tmp = numroc( numincbands_tmp, blocksizeincbands_tmp, Node, 0, Nodes )
!!      For debugging
!       write(6,*)
!       write(6,'(a,l5)')                                                   &
! &       'setup_w90_wrapper: r_between_manifolds = ',                    &
! &       w90_r_between_manifolds
!       write(6,'(a,6i5)')                                                  &
! &       'setup_w90_wrapper: index_manifold, Node, Nodes, nincbands_loc_tmp, numincbands_tmp, Blocksize = ',  &
! &       index_manifold, Node, Nodes, nincbands_loc_tmp, numincbands_tmp,                  &
! &       blocksizeincbands_tmp
!!      End debugging
#else
      nincbands_loc_tmp = numincbands_tmp
      blocksizeincbands_tmp = 1
#endif

      w90man%nincbands_loc = nincbands_loc_tmp
      w90man%blocksizeincbands = blocksizeincbands_tmp

    end subroutine read_manifold

    !> \brief General purpose of the subroutine set_excluded_bands_w90_in:
    !! identify the bands that will be considered for wannierization
    !! and the atomic orbitals that will be used as initial guess for the
    !! projections.
    !!
    !! Within this subroutine we:
    !! 1. Identify the bands that will be considered for the wannierization of
    !! a given manifold. These can be found in the variable:
    !! manifold_bands(index_manifold)%isexcluded.
    !! 2. Identify the numerical atomic orbitals of the basis set that will
    !! be considered as local functions for the initial guess.
    !! These can be found in the variable
    !! manifold_bands(index_manifold)%orbexcluded.
    !! 3. Set a sequential index for the orbitals that will be used
    !! as initial guess, stored in the variable
    !! manifold_bands(index_manifold)%orb_in_manifold.
    !! This routine is called only by the IOnode.
    subroutine set_excluded_bands_w90_in(w90man)
      type(w90_manifold_t), intent(inout) :: w90man

      integer :: iband
      integer :: iorb
      integer :: index
      integer :: index_orb_included

      call re_alloc( w90man%isexcluded, 1, no_u, &
          name='isexcluded', routine='set_excluded_bands_w90_in' )

      call re_alloc( w90man%orbexcluded, 1, no_u, &
          name='orbexcluded', routine='set_excluded_bands_w90_in' )

      call re_alloc( w90man%orb_in_manifold, 1, no_u, &
          name='orb_in_manifold', routine='set_excluded_bands_w90_in' )

      ! By default, all the bands are excluded from the calculation
      w90man%isexcluded(:) = .true.

      do iband = w90man%initial_band, w90man%final_band
        ! Exclude the corresponding bands from the computation
        w90man%isexcluded( iband ) = .false.
      end do

!!   For debugging
!    write(6,*)
!    do iterex = 1, no_u
!     write(6,'(a,2i5,l5)')                                                     &
! &  'set_excluded_bands_w90_in: Node, excluded_bands_w90_in (iband, excluded)',&
! &     Node, iterex, w90man%isexcluded(iterex)
!    enddo
!!   End debugging

      ! By default, all the orbitals are excluded from the calculation
      w90man%orbexcluded(:) = .true.
      do iorb = 1, w90man%numbands
        index_orb_included = w90man%orbital_indices(iorb)
        if( index_orb_included > 0 ) &
            w90man%orbexcluded(index_orb_included) = .false.
      end do

      index = 0
      w90man%orb_in_manifold(:) = 0
      do iorb = 1, w90man%numbands
        index = index + 1
        index_orb_included = w90man%orbital_indices(iorb)
        if( index_orb_included > 0 ) &
            w90man%orb_in_manifold(index_orb_included) = index
      end do

!!   For debugging
!    write(6,*)
!    do iterex = 1, no_u
!     write(6,'(a,2i5,l5,i5)')                                            &
! &     'set_excluded_bands_w90_in: orbitals excluded: manifold, orb_unit cell, excluded, orb_mani', &
! &     index_manifold, iterex,                                           &
! &     w90man%orbexcluded(iterex),        &
! &     w90man%orb_in_manifold(iterex)
!    enddo
!!   End debugging

    end subroutine set_excluded_bands_w90_in

    !> \brief General purpose of the define_trial_orbitals subroutine:
    !! populate the derived variable where all the
    !! information regarding the trial localized functions will be stored.
    !!
    !! In this subroutine we populate the derived variable where all the
    !! information regarding the trial localized functions will be stored
    !! This variable, called proj_w90_in, is defined within the derived variable
    !! manifold_bands.
    !! The dimension of proj_w90_in should be the same as the number
    !! of bands to be wannierized
    !! Some default values (x-axis, z-axis, zovera, r) to define the
    !! hydrogenoid functions are taken directly from WANNIER90.
    !! Indeed, this is not a major issue since they will not be used
    !! in SIESTA, where we will take directly the shape of the numerical
    !! atomic orbitals in the basis set.
    !! The center of the localized projection functions are directly given in
    !! cartesian Bohrs, and not in fractional units as it is done in WANNIER90.
    !! The m-angular quantum number index is different in SIESTA and WANNIER90.
    !! In WANNIER90 they are defined as in Tables 3.1 and 3.2 of the User Guide.
    !! In SIESTA, they are defined as in the SystemName.ORB_INDX file.
    subroutine define_trial_orbitals(name, w90man)
      use trialorbitalclass,  only: trialorbital
                                        ! Derived type to define the
                                        !    localized trial
                                        !    orbitals
      use trialorbitalclass, only: cutoffs
                                        ! Cut-off radii in units of \alpha^-1
      use units,       only: Ang        ! Conversion factor
      use atmfuncs,    only: lofio      ! Returns angular momentum number
      use atmfuncs,    only: mofio      ! Returns magnetic quantum number
      use atmfuncs,    only: rcut       ! Returns orbital cutoff radius
      use atmfuncs,    only: orb_gindex ! Returns the global index of a
                                        !    basis orbital
      use atomlist,    only: iaorb      ! Atomic index of each orbital
      use atomlist,    only: iphorb     ! Orbital index of each orbital
                                        !    in its atom
      use siesta_geom, only: xa         ! Atomic coordinates
      use siesta_geom, only: isa        ! Species index of each atom
      use siesta_geom, only: ucell      ! Unit cell lattice vectors
                                        !   in real space as used inside
                                        !   SIESTA
                                        !   First  index: component
                                        !   Second index: vector
                                        !   In Bohrs

      character(len=*), intent(in) :: name
      type(w90_manifold_t), intent(inout) :: w90man

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer  :: number_projections
      integer  :: iproj        ! Counter for the number of projections
      integer  :: iorb         ! Index of the atomic orbital
      integer  :: ia           ! Pointer to atom to which orbital belongs
      integer  :: iao          ! Orbital index within atom
      integer  :: is           ! Species index
      real(dp) :: rc           ! Orbital's cutoff radius
      integer  :: l            ! Angular momentum of the localized trial function
      integer  :: m            ! Magnetic quantum number of the localized trial
                               !   function.
                               ! The corresponding to the real spherical harmonics
                               !   can be found in Table 3.1 and 3.2 of the
                               !   WANNIER90 Users Guide
      integer  :: r            ! radial (optional):
                               ! If r=2 – use a radial function with one node
                               !   (ie second highest pseudostate with that
                               !   angular momentum).
                               !   Default is r=1.
                               !   Radial functions associated with different
                               !   values of r should be orthogonal to each other.
      real(dp) :: xaxis(3)     ! Sets the x-axis direction. Default is x=(1,0,0)
      real(dp) :: xaxis_tmp(3) ! Sets the x-axis direction, read from WannierProjectors
      real(dp) :: yaxis(3)     ! Sets the y-axis direction
      real(dp) :: yaxis_tmp(3) ! Sets the y-axis direction, read from WannierProjectors
      real(dp) :: zaxis(3)     ! Sets the z-axis direction. Default is z=(0,0,1)
      real(dp) :: zaxis_tmp(3) ! Sets the z-axis direction, read from WannierProjectors
      real(dp) :: center_tmp(3)! Sets the center of the Wannier function,
                               !   read from WannierProjectors
      real(dp) :: zovera       ! zovera (optional):
                               !   the value of Z for the radial part
                               !   of the atomic orbital
                               !   (controls the diffusivity of the radial
                               !   function).
                               !   Units always in reciprocal Angstrom.
                               !   Default is zovera = 1.0.
      integer  :: nprojs_wannier90_style
                               ! Number of projectors whose information
                               !   must be included in the WannierProjectors."name" block
      integer  :: nprojs_in_block
                               ! Number of projectors actually
                               !   included in the WannierProjectors."name" block

      integer  :: idx
      integer, allocatable :: trial_orbs_index(:)

      integer  :: i, j
      real(dp) :: xnorm_tmp    ! Norm of the x-vector read from WannierProjectors
      real(dp) :: znorm_tmp    ! Norm of the z-vector read from WannierProjectors
      real(dp) :: cosine       ! Cosine between the x and z vectors

      real(dp), parameter :: eps4  = 1.0e-4_dp
      real(dp), parameter :: eps6  = 1.0e-6_dp

      ! Sets the default for the x-axis
      xaxis(1) = 1.0_dp
      xaxis(2) = 0.0_dp
      xaxis(3) = 0.0_dp

      ! Sets the default for the y-axis
      yaxis(1) = 0.0_dp
      yaxis(2) = 1.0_dp
      yaxis(3) = 0.0_dp

      ! Sets the default for the z-axis
      zaxis(1) = 0.0_dp
      zaxis(2) = 0.0_dp
      zaxis(3) = 1.0_dp

      ! Sets the default for zovera
      zovera = 1.0_dp
      ! zovera = 1.0_dp/0.529177_dp

      ! Set up the default value for r (use radial functions without node)
      r = 1

      ! There should be as many localized trial orbitals as the number of bands
      ! to be projected
      number_projections = w90man%numbands

      ! Allocate the derived variable where all the data for the
      ! localized trial orbitals will be stored
      if (allocated(w90man%proj)) deallocate( w90man%proj )
      allocate(w90man%proj(number_projections))

      nprojs_wannier90_style = 0
      ! trial_orbs are "classic Wannier trial functions" (as opposed to Siesta PAOs)
      allocate(trial_orbs_index(number_projections))  ! for bookeeping below

      do iproj = 1, number_projections
        ! Identify the atomic orbitals on which we are going to project
        iorb = w90man%orbital_indices(iproj)
        if( iorb > 0 ) then      ! Siesta orbital
          ia = iaorb(iorb)               ! Atom to which orbital belongs
          is = isa(ia)                   ! Atomic species of the atom where the
                                         !   orbital is centered
          iao = iphorb( iorb )           ! Orbital index within atom
          l  = lofio( is, iao )          ! Orbital's angular mumentum number
          m  = mofio( is, iao )          ! (Real) orbital's magnetic quantum number
          rc = rcut(  is, iao )          ! Orbital's cutoff radius
          w90man%proj(iproj)%center = xa(:,ia)
          w90man%proj(iproj)%zaxis  = zaxis
          w90man%proj(iproj)%xaxis  = xaxis
          w90man%proj(iproj)%yaxis  = yaxis
          w90man%proj(iproj)%zovera = zovera
          w90man%proj(iproj)%r      = r
          w90man%proj(iproj)%l      = l
          w90man%proj(iproj)%mr     = m
          w90man%proj(iproj)%rcut   = rc
          w90man%proj(iproj)%lmax   = l
          w90man%proj(iproj)%from_basis_orbital = .true.
          w90man%proj(iproj)%iorb   = iorb
          w90man%proj(iproj)%iorb_gindex = orb_gindex(is,iao)

        else   ! Wannier90-style trial orbital, to be specified in block

           nprojs_wannier90_style = nprojs_wannier90_style + 1
           trial_orbs_index(nprojs_wannier90_style) = iproj
        end if
      end do

    ! Read the information for the projectors from the Block,
    ! in case some hybrid trial functions are required
    if( nprojs_wannier90_style > 0 ) then
      ! Count that there is the correct number of lines present
      nprojs_in_block = fdf_block_linecount('WannierProjectors.'//name, 'rrriiirrrrrrr')
      if ( nprojs_wannier90_style /= nprojs_in_block ) then
        write(6,'(3a)') 'Wrong number of projectors in the WannierProjectors.',name,' block'
        write(6,'(a,i5)')' nprojs_wannier90_style = ', nprojs_wannier90_style
        write(6,'(a,i5)')' projectors_in_block   = ', iproj
        call die('Wrong number of projectors in the WannierProjectors.'//name//' block')
      end if

      if (.not. fdf_block('WannierProjectors.'//name,bfdf)) then
        call die('Block WannierProjectors to define the trial guess functions required')
      end if

      ! Read the content of the block, line by line
      ! Using the trial_orbs_index array we can relax the requirement
      ! that the 'negative indices' appear last in the manifold
      iproj = 0
      do while ( fdf_bline(bfdf, pline) )
        iproj = iproj + 1
        if ( iproj > number_projections ) &
            call die("Error in logic, should never be reached (check above)")

        center_tmp(1) = fdf_breals(pline,1)
        center_tmp(2) = fdf_breals(pline,2)
        center_tmp(3) = fdf_breals(pline,3)

        idx = trial_orbs_index(iproj)
        ASSOCIATE ( proj => w90man%proj(idx) ) ! Compact notation

        proj%zaxis(1)  = fdf_breals(pline,4)
        proj%zaxis(2)  = fdf_breals(pline,5)
        proj%zaxis(3)  = fdf_breals(pline,6)
        proj%xaxis(1)  = fdf_breals(pline,7)
        proj%xaxis(2)  = fdf_breals(pline,8)
        proj%xaxis(3)  = fdf_breals(pline,9)
        proj%zovera    = fdf_breals(pline,10)
        proj%l         = fdf_bintegers(pline,1)
        proj%mr        = fdf_bintegers(pline,2)
        proj%r         = fdf_bintegers(pline,3)
        proj%from_basis_orbital = .false.

        xaxis_tmp = proj%xaxis
        zaxis_tmp = proj%zaxis

        ! Check that the xaxis and the zaxis are orthogonal to each other
        xnorm_tmp = sqrt( dot_product(xaxis_tmp, xaxis_tmp) )
        if (xnorm_tmp < eps6) call die('define_trial_orbitals: |xaxis_tmp| < eps ')
        znorm_tmp = sqrt( dot_product(zaxis_tmp, zaxis_tmp) )
        if (znorm_tmp < eps6) call die('define_trial_orbitals: |zaxis_tmp| < eps ')

        cosine = dot_product(xaxis_tmp, zaxis_tmp) / (xnorm_tmp * znorm_tmp)
        if ( abs(cosine) > eps6 ) &
            call die('define_trial_orbitals: xaxis_tmp and zaxis_tmp are not orthogonal')
        if ( proj%zovera < eps6 ) &
            call die('define_trial_orbitals: zovera value must be positive')

        ! Compute the y-axis
        call cross(zaxis, xaxis, yaxis_tmp)
        proj%yaxis = yaxis_tmp / sqrt(dot_product(yaxis_tmp,yaxis_tmp))

        ! Now convert "center" from ScaledByLatticeVectors to Bohrs
        do i = 1 , 3
          proj%center(i) = 0.0_dp
          do j = 1 , 3
            proj%center(i) = center_tmp(j)*ucell(i,j) + proj%center(i)
          end do
        end do

        proj%zovera = proj%zovera * Ang  !To Bohr^-1

        ! Select the maximum angular momentum
        select case(proj%l)
        case(0:3)
          proj%lmax = proj%l
        case(-3:-1)
          proj%lmax = 1 ! spN hybrids
        case(-5:-4)
          proj%lmax = 2 ! spdN hybrids
        case default
          call die("Invalid l in define_trial_orbitals")
        end select

        ! Further checks
        if ( proj%mr < 1 ) call die("Invalid mr in define_trial_orbitals ")
        if ( proj%r > 3 .or. proj%r < 1 ) then
          call die("Invalid r in define_trial_orbitals")
        else if ( proj%l  >= 0 .and. proj%mr > 2*proj%l+1 ) then
          call die("Invalid mr_w in define_trial_orbitals")
        else if ( proj%l < 0 .and. proj%mr > 1-proj%l ) then
          call die("Invalid mr_w in define_trial_orbitals")
        end if

        ! Cut-off initialization (in Bohr)
        proj%rcut = cutoffs(proj%r) / proj%zovera

        END ASSOCIATE   ! --- cosmetic device to simplify code
      end do !  do while lines in the fdf block
    end if

  end subroutine define_trial_orbitals

  subroutine int_array_extend(orig, add)
      integer, intent(inout), allocatable :: orig(:)
      integer, intent(in) :: add(:)

      integer, allocatable :: orig2(:)
      integer :: norig, nadd

      if ( allocated(orig) ) then
        norig = size(orig)
      else
        norig = 0
      end if
      nadd = size(add)

      if ( nadd <= 0 ) return

      if ( norig == 0 ) then
        allocate(orig(nadd))
        orig(:) = add(:)
      else
        allocate(orig2(norig))
        orig2(:) = orig(:)
        deallocate(orig)
        allocate(orig(norig + nadd))
        orig(1:norig) = orig2(:)
        deallocate(orig2)
        orig(norig+1:) = add(:)
      end if

    end subroutine int_array_extend

  end subroutine read_w90_wrapper_specs

!> \brief General purpose of the subroutine read_kpoints_wannier:
!! process the information in the fdf file
!! required to generate the k-point sampling for Wannierization.
!!
!! In this subroutine, we process the information in the fdf file
!! required to generate the k-point sampling that will be used inside
!! WANNIER90 for the Wannier transformation.
!! The block that is readed and digested is Wannier.k.
!! Example:
!!
!! %block Wannier.k
!!
!!   20  20  1
!!
!! %endblock
!! Alternatively it may be provided as an integer list.
!! Wannier.k [20 20 1]
!!
!! where the three integers are the number of k-points along the corresponding
!! reciprocal lattice vector
!! The algorithm to generate the k-points in reciprocal units is borrowed from
!! the utility kmesh.pl of WANNIER90.
!!
!! Then, and following the recipes given in the Appendix of \cite Marzari-97,
!! we compute the distance to nearest-neighbour shells of k-points
!! the vectors connecting neighbour k-points in the Monkhorst-Pack mesh,
!! and check whether the completeness relation is fully satisfied
!! [Eq. (B1) of Ref. \cite Marzari-97]
!!
!! Finally, we populate the variables related with the neighbour k-points in
!! the Monkhorst-Pack mesh in the module where all the Wannier90 parameters are
!! stored.

  subroutine read_kpoints_wannier

    use fdf
    use siesta_geom,    only: ucell            ! Unit cell lattice vectors
                                               !   in real space as used inside
                                               !   SIESTA
                                               !   First  index: component
                                               !   Second index: vector
                                               !   In Bohrs
    use w90_wrapper_types,   only: latvec_w90_in   ! Lattice vectors
                                               !   in real space as used inside
                                               !   wannier90
                                               !   First  index: vector
                                               !   Second index: component
                                               !   In Ang.
    use w90_wrapper_types,   only: reclatvec_w90_in
                                               ! Reciprocal lattice vectors
                                               !   computed from real_lattice
                                               !   The factor 2.0 * pi
                                               !   is included
                                               !   First  index: component
                                               !   Second index: vector
                                               !   In Angstroms^-1
    use w90_wrapper_types,   only: kmesh_w90_in
                                               ! Number of divisions along the
                                               !   reciprocal lattice vectors
    use w90_wrapper_types,   only: numkpoints_w90_in
                                               ! Total number of k-points
    use w90_wrapper_types,   only: kpointsfrac_w90_in
                                               ! Coordinates of the k-points
                                               !   In fractional units, i. e.
                                               !   in units of the reciprocal
                                               !   lattice vectors
                                               !   First  index: component
                                               !   Second index: vector
    use w90_wrapper_types,   only: nncount_w90_in
                                               ! Same as nntot
                                               !    but in the module where the
                                               !    variables that controls the
                                               !    Wannier90 in SIESTA
                                               !    are stored
    use w90_wrapper_types,   only: nnlist_w90_in
                                               ! Same as nnlist
                                               !    but in the module where the
                                               !    variables that controls the
                                               !    Wannier90 in SIESTA
                                               !    are stored
    use w90_wrapper_types,   only: nnfolding_w90_in
                                               ! Same as nncell
                                               !    but in the module where the
                                               !    variables that controls the
                                               !    Wannier90 in SIESTA
                                               !    are stored
    use w90_wrapper_types,   only: bvectorsfrac_w90_in
                                               ! Vectors that connect each mesh
                                               !    k-point
                                               !    to its nearest neighbours
    use m_digest_nnkp,  only: chosing_b_vectors! Subroutine that computes the b
                                               !    vectors that connect each
                                               !    mesh k-point to its
                                               !    nearest neighbours.

    use units, only : Ang

!   Internal variables
    integer :: ik
    integer :: ikx, iky, ikz
    integer :: i
    integer :: nkp
    integer :: nn


    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    ! Read the data to generate the grid in reciprocal space that will be used
    ! for the Wannier Projections
    if ( fdf_block("Wannier.k", bfdf) ) then

      do while( fdf_bline(bfdf, pline) )
        ! We expect 3 integers for the line
        ! Each integer corresponds to the number of k-points
        ! along the given lattice vector.
        if (.not. fdf_bmatch(pline, 'iii') ) &
            call die('Wrong format in Wannier.k')
        kmesh_w90_in(1) = fdf_bintegers(pline,1)
        kmesh_w90_in(2) = fdf_bintegers(pline,2)
        kmesh_w90_in(3) = fdf_bintegers(pline,3)
      end do

    else if ( fdf_islist("Wannier.k") ) then

      call fdf_list("Wannier.k", 3, kmesh_w90_in)

    else

       ! Just gamma...
       ! return (do not...)
       kmesh_w90_in(:) = 1
    end if

!   Define the total number of k-points used in the Wannier projection
    numkpoints_w90_in = kmesh_w90_in(1) * kmesh_w90_in(2) * kmesh_w90_in(3)

!   Compute and store the components of the k-points in fractional units
    nullify( kpointsfrac_w90_in )
    call re_alloc( kpointsfrac_w90_in, 1, 3, 1, numkpoints_w90_in,  &
 &                 name='kpointsfrac_w90_in', routine='read_kpoints_wannier')

!   The algorithm to generate the k-points in fractional units
!   is borrowed from the utility kmesh.pl of Wannier90
    ik = 0
    do ikx = 0, kmesh_w90_in(1) - 1
      do iky = 0, kmesh_w90_in(2) - 1
        do ikz = 0, kmesh_w90_in(3) - 1
          ik = ik + 1
          kpointsfrac_w90_in(1,ik) = (ikx*1.0_dp)/kmesh_w90_in(1)
          kpointsfrac_w90_in(2,ik) = (iky*1.0_dp)/kmesh_w90_in(2)
          kpointsfrac_w90_in(3,ik) = (ikz*1.0_dp)/kmesh_w90_in(3)
        enddo
      enddo
    enddo

!   Transform the units of the unit cell lattice vector from Bohrs to Ang,
!   as required by Wannier90
!   In WANNIER90 the real space lattice vector matrix is the
!   transpose one used in SIESTA
!   In WANNIER90 the first component is the vector and the second the componen
!   In SIESTA the first component is the component and the second the vector
!   We make the transposition here to transfer the info to WANNIER90

    latvec_w90_in = transpose(ucell) / Ang

!   Compute the reciprocal lattice vectors as required by WANNIER90
!   The factor 2.0 * pi is included
    call reclat( latvec_w90_in, reclatvec_w90_in, 1 )

!   Now, we generate the vectors that connect nearest-neighbour shells
!   in reciprocal space, together with the completeness relation,
!   Eq. (B1) in N. Marzari et al. Physical Review B 56, 12847 (1997)

    !
    block
      use wannier90_m, only: wannier90_wrapper
      integer:: w90_lun, ik, w90_comm_f

      if (IOnode) then
         open(newunit=w90_lun,file="_nnkp.win",status="replace",action="write")
         ! This is not used, but needed in the input. The number is arbitrary
         write(w90_lun,*) "num_wann  =  120"
         write(w90_lun,*) "begin unit_cell_cart"
         do ik = 1, 3
            write(w90_lun,'(3f14.6)') latvec_w90_in(ik,:)
         enddo
         write(w90_lun,*) "end unit_cell_cart"
         write(w90_lun,"(a,1x,3i4)") "mp_grid", kmesh_w90_in(1:3)
         write(w90_lun,*) "begin kpoints"
         do ik = 1, numkpoints_w90_in
            write(w90_lun,'(3f18.12)') kpointsfrac_w90_in(:,ik)
         enddo
         write(w90_lun,*) "end kpoints"
         close(w90_lun)

         write(6,"(a)") "... Calling wannier90 to generate " // &
                          "k-point neighbor information"
         write(6,"(a)") "... See file _nnkp.wout for information"
      end if

      !   Call wannier90 in "nnkp" mode, and
      !   store:
      !    - the number of nearest neighbours belonging to each k-point of the
      !      Monkhorst-Pack mesh
      !    - the list of nearest neighoburs for each k-point
      !    - the vector, in fractional reciprocal lattice coordinates,
      !      that brings the nnth nearest neighbour of k-point nkp to
      !      its periodic image that is needed for computing the
      !      overlap matrices M_mn(k,b)

#ifdef MPI
      w90_comm_f = MPI_COMM_ID(w90_comm)
#endif
      call wannier90_wrapper("_nnkp",        &
#ifdef MPI
                       mpi_comm=w90_comm_f, &
#endif
                       nnkp_mode=.true., &
                       nntot_out=nncount_w90_in, &
                       nnlist_out=nnlist_w90_in, &
                       nncell_out=nnfolding_w90_in)


    end block

!   Compute the vectors that connect each mesh k-point
!   to its nearest neighbours
    call chosing_b_vectors( kpointsfrac_w90_in, nncount_w90_in,  &
 &                          nnlist_w90_in, nnfolding_w90_in,     &
 &                          bvectorsfrac_w90_in )

  end subroutine read_kpoints_wannier


! ------------------------------------------------------------------------------
!
!> \brief General purpose of the subroutine compute_matrices:
!! Compute the overlap between periodic parts of wave functions at neighbour
!! k-points, and the projections between the wave functions and the initial
!! guesses of the localized functions.
!!
!! In this subroutine we compute and dump in the corresponding files:
!! 1. The matrix elements of the plane waves connecting neighbour
!! points in the first-Brillouin zone for Wannierization in a basis
!! of atomic orbitals
!! (done in the subroutine compute_pw_matrix)
!! 2. The overlap matrices between the periodic part of wave functions
!! at neighbour k-points,
!! \f{eqnarray*}{
!! M_{m n}^{(\vec{k}, \vec{b})} =
!!        \langle u_{m \vec{k}} \vert u_{n \vec{k} + \vec{b}} \rangle,
!! \f}
!! that is the Eq. (27) of the Ref. \cite Marzari-12, or
!! Eq. (25) of Ref. \cite Marzari-97
!! (done in the subroutine mmn)
!! 3. The overlap matrices between the eigenstates of the one-particle
!! Hamiltonian and the localized trial orbitals, taken as initial guess,
!! \f{eqnarray*}{
!! A_{mn} (\vec{k}) = \langle \psi_{m\vec{k}} \vert g_{n} \rangle,
!! \f}
!! as presented right before Eq. (17) of Ref. \cite Marzari-12
!! (done in the subroutine amn)
!! 4. The eigenvalues for the bands to be wannierized
!! (done in the subroutine writeeig)
!! 5. The values of the periodic part of the wavefunctions at the points
!! of a real space grid. This must be done if we want to plot the
!! Wanniers
!! (done in the subroutine writeunk)

  subroutine compute_matrices( ispin, index_manifold )

    use m_switch_local_projection, only: seedname
    use m_switch_local_projection, only: nncount
    use m_switch_local_projection, only: numkpoints
    use m_switch_local_projection, only: bvectorsfrac
    use m_switch_local_projection, only: Amnmat
    use m_switch_local_projection, only: Mmnkb
    use files,                     only: slabel        ! Short system label,
                                                       !   used to generate file
                                                       !   names
    use m_spin,                    only: spin          ! Spin configuration
                                                       !   for SIESTA
    use sys,                       only: die

    integer, intent(in) :: ispin           ! Spin component
    integer, intent(in) :: index_manifold  ! Index of the manifold

!   Internal variables
    integer :: ik                          ! Counter for the loops on k-points
    integer :: n                           ! Counter for the loops on bands
    integer :: m                           ! Counter for the loops on bands

    integer :: number_of_bands_in_manifold_local
    integer :: number_of_bands_to_project

!   Variables related with the reading of the amn files
    integer            :: iman          ! Counter for the manifolds
    character(len=100) :: filename      ! Name of the files where the amn
                                        !   matrices of the manifolds to be
                                        !   mixed are stored
    character(len=50)  :: dummy         ! Dummy variable
    integer            :: nb_tmp        ! Number of bands (read from the file)
    integer            :: nkp_tmp       ! Number of k-points (read from the file
    integer            :: nw_tmp        ! Number of Wanniers (read from the file
!    integer            :: m             ! Band index
!                                        !    (read from the loop in the file)
!    integer            :: n             ! Projector index
!                                        !    (read from the loop in the file)
    integer            :: nkp           ! k-point index
                                        !    (read from the loop in the file)
    integer            :: num_amn       ! Total number of entries in the Amn
                                        !   matrix
    integer            :: iproj         ! Counter for the loop on projectors
    integer            :: icount        ! Counter for the loop on k-neighbouts
    integer            :: mband         ! Counter for the loop on bands
    integer            :: nband         ! Counter for the loop on bands

    external     :: io_assign              ! Assign a logical unit
    external     :: io_close               ! Close a logical unit

#ifdef MPI
    integer            :: MPIError
#endif

    type(w90_manifold_t), pointer :: mnf

    mnf => manifold_bands(index_manifold)

    if( spin%H .eq. 1) then
      seedname = mnf%seedname
    else if( spin%H .gt. 1) then
      write(seedname, "(a,'.spin.',i1.1)") trim(mnf%seedname), ispin
    end if

    number_of_bands_in_manifold_local = mnf%nincbands_loc
    number_of_bands_to_project = mnf%numbands

!   Compute the matrix elements of the plane wave,
!   for all the wave vectors that connect a given k-point to its nearest
!   neighbours
!    if( first_chempotwann ) call compute_pw_matrix( nncount, bvectorsfrac )
    call compute_pw_matrix( nncount, bvectorsfrac )

!   Compute the overlaps between Bloch orbitals at neighboring k points
!    if( first_chempotwann ) call mmn( ispin )
    call mmn( ispin )

!   Compute the overlaps between Bloch states onto trial localized orbitals
!    if( first_chempotwann ) call amn( ispin, index_manifold )
    call amn( ispin, index_manifold )

!   Check if the computations of the position operator matrix elements
!   between manifolds has to be carried out
    if( w90_r_between_manifolds ) then
      if( index_manifold .eq. 3) then
        if( (.not. compute_chempotwann) .or. first_chempotwann ) then
#ifdef MPI
          call MPI_barrier(MPI_Comm_world,MPIError)
#endif
          if( IONode ) then
            do nkp = 1, numkpoints
              do iproj = 1, manifold_bands(1)%numbands
                do mband = manifold_bands(1)%number_of_bands+1,  &
                           manifold_bands(1)%number_of_bands +   &
                           manifold_bands(2)%number_of_bands
                  Amnmat(mband, iproj, nkp) = cmplx(0.0_dp,0.0_dp,kind=dp)
                enddo
              enddo
              do iproj = manifold_bands(1)%numbands+1,    &
                         manifold_bands(3)%numbands
                do mband = 1, manifold_bands(1)%number_of_bands
                  Amnmat(mband, iproj, nkp) = cmplx(0.0_dp,0.0_dp,kind=dp)
                enddo
              enddo
            enddo

            if( w90_mmn_diagonal ) then
              do nkp = 1, numkpoints
                do icount = 1, nncount
                  do mband = 1, manifold_bands(1)%number_of_bands
                    do nband = manifold_bands(1)%number_of_bands+1, &
                               manifold_bands(1)%number_of_bands +   &
                               manifold_bands(2)%number_of_bands
                      Mmnkb(mband,nband,nkp,icount) = cmplx(0.0_dp,0.0_dp,kind=dp)
                    enddo
                  enddo
                  do mband = manifold_bands(1)%numbands+1,    &
                             manifold_bands(3)%numbands
                    do nband = 1, manifold_bands(1)%number_of_bands
                      Mmnkb(mband,nband,nkp,icount) = cmplx(0.0_dp,0.0_dp,kind=dp)
                    enddo
                  enddo
                enddo
              enddo
            endif

!           Write the Mmn overlap matrices in a file, in the format required
!           by Wannier90
            call writemmn( ispin )

!           Write the Amn overlap matrices in a file, in the format required
!           by Wannier90
            call writeamn( ispin )

          endif
#ifdef MPI
          call MPI_barrier(MPI_Comm_world,MPIError)
#endif
        endif
      endif
    endif

!   Write the selected eigenvalues within the manifold in the file
!   SystemLabel.eigW
    if( IOnode ) call writeeig( ispin )

    if( mnf%write_unk ) call writeunk( ispin )

    if (IONode) then
      write(6,'(/,a)')  &
 &     'compute_matrices: All the information dumped in the corresponding files'
      write(6,'(a)')  &
 &     'compute_matrices: End of the interface between Siesta and Wannier90'
    endif

  end subroutine compute_matrices


!> \brief General purpose of the subroutine compute_wannier:
!!
!! Within this subroutine:
!! 1. We create a .win file required by wannier90 using
!!    variables transferred from different modules in SIESTA.
!! 2. We call the wannier90 wrapper in "full" mode.

  subroutine compute_wannier( ispin, index_manifold )

!
!   Variables related with the atomic structure coming from SIESTA
!
    use siesta_geom,    only : na_u            ! number of atoms in the
                                               !   unit cell
    use siesta_geom,    only : xa              ! atomic positions
                                               !   in cartesian coordinates
                                               !   units in Bohrs
    use siesta_geom,    only : isa             ! species index of each atom
    use units,          only : Ang             ! conversion factor from
                                               !   Ang to Bohrs
    use atm_types,      only : species         ! information about the different
                                               !   chemical species
!
!   Variables related with the k-points sampling coming from SIESTA
!
    use w90_wrapper_types,   only: latvec_w90_in   ! Lattice vectors
    use w90_wrapper_types,   only: kpointsfrac_w90_in
    use w90_wrapper_types,   only: numkpoints_w90_in
    use w90_wrapper_types,   only: kmesh_w90_in
                                               ! Number of divisions along the
                                               !   reciprocal lattice vectors
!
!   Variables related with the input/output coming from SIESTA
!
    use files,          only: slabel           ! Short system label,
                                               !   used to generate file names
    use m_spin,         only: spin             ! Spin configuration for SIESTA
!
!   Variables related with the post-processing coming from SIESTA
!
    use m_energies,     only: ef               ! Fermi energy
    use units,          only: eV               ! Conversion factor from Ry to eV

!
!   Variables related with the parallelization, coming from SIESTA
!
    use parallel,        only: IOnode          ! input/output node?
!
! Input variables
!
    integer, intent(in) :: ispin            ! Spin index
    integer, intent(in) :: index_manifold   ! Index of the manifold
                                            !   that is wannierized

    character(len=256)  :: seedname
    type(w90_manifold_t), pointer :: mnf
    integer :: w90_comm_f

    mnf => manifold_bands(index_manifold)

!   Set up the variables related with the writing of the Hamiltonian
    if( spin%H == 1) then
      seedname = mnf%seedname
    else if( spin%H > 1) then
      write(seedname, "(a,'.spin.',i1.1)") &
          trim(mnf%seedname), ispin
    end if

    block
      use wannier90_m, only: wannier90_wrapper
      integer:: w90_lun, w90_eig, ik, ia, iband
      character(len=256) :: filename

      if (IOnode) then
         filename = trim(seedname) // ".win"
         open(newunit=w90_lun,file=filename,status="replace",action="write")

         write(w90_lun,*) "num_wann =", mnf%numbands
         write(w90_lun,*) "num_bands =", mnf%number_of_bands
         write(w90_lun,*) "num_iter =", mnf%num_iter
         ! numproj set by w90
         if (mnf%write_hr) write(w90_lun,*) "write_hr = T"
         if (mnf%write_tb) write(w90_lun,*) "write_tb = T"

         if (numkpoints_w90_in==1) write(w90_lun,*) "gamma_only = T"

         if (mnf%dis_win_siesta) then
            write(w90_lun,"(a,f14.6)") "dis_win_min =", mnf%dis_win(1)
            write(w90_lun,"(a,f14.6)") "dis_win_max =", mnf%dis_win(2)
         endif
         if (mnf%dis_win_froz_siesta) then
            write(w90_lun,"(a,f14.6)") "dis_froz_min =", mnf%dis_froz(1)
            write(w90_lun,"(a,f14.6)") "dis_froz_max =", mnf%dis_froz(2)
         endif

         write(w90_lun,*) "begin unit_cell_cart"
         do ik = 1, 3
            write(w90_lun,'(3f14.6)') latvec_w90_in(ik,:)
         enddo
         write(w90_lun,*) "end unit_cell_cart"

         write(w90_lun,*) "begin atoms_cart"
         do ia = 1, na_u
            write(w90_lun,'(a,1x,3f12.6)') trim(species(isa(ia))%label), &
                                           xa(:,ia)/Ang
         enddo
         write(w90_lun,*) "end atoms_cart"


         write(w90_lun,"(a,1x,3i4)") "mp_grid", kmesh_w90_in(1:3)
         write(w90_lun,*) "begin kpoints"
         do ik = 1, numkpoints_w90_in
            write(w90_lun,'(3f18.12)') kpointsfrac_w90_in(:,ik)
         enddo
         write(w90_lun,*) "end kpoints"

         if (mnf%wannier_plot) write(w90_lun,*) "wannier_plot = T"
         write(w90_lun,"(a,3i4)") "wannier_plot_supercell =", &
                               mnf%wannier_plot_supercell(1:3)

         if (mnf%fermi_surface_plot) write(w90_lun,*) "fermi_surface_plot = T"
         write(w90_lun,'(a,f14.6)') "fermi_energy = ", ef/eV

         close(w90_lun)

         write(6,"(a)") "... Calling wannier90 for this manifold"
         write(6,"(a)") "... See file " // trim(seedname) // &
                        ".wout for information"

      end if

      !   Call wannier90 in "full" mode,
      !   store the U and U_opt matrices
      !   Files .amn, .mmn, and .eigW have by
      !   now been created in compute_matrices
#ifdef MPI
      w90_comm_f = MPI_COMM_ID(w90_comm)
#endif
      call wannier90_wrapper(seedname,          &
#ifdef MPI
                          mpi_comm=w90_comm_f, &
#endif
                          eigfile_ext =".eigW",               &
                          u_matrix_out=mnf%u_matrix,          &
                          u_matrix_opt_out=mnf%u_matrix_opt)

    end block

  end subroutine compute_wannier


!> \brief General purpose of the subroutine deallocate_wannier:
!! free data structures (wip) and communicator
!!

  subroutine deallocate_wannier()

#ifdef MPI
    integer :: ierr

    call MPI_Comm_Free(w90_comm, ierr)
#endif

  end subroutine deallocate_wannier

#else
  CONTAINS
  subroutine dummy_routine()
  end subroutine dummy_routine
#endif
endmodule m_w90_wrapper
