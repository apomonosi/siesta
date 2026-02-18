
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_occ_proj module
!! In this subroutine we compute the projection of the Kohn-Sham orbitals into
!! the states of a localized basis of LDA+U projectors,
!! \f$ \vert \phi_{m}^{I} \rangle \f$.
!!
!! Following Ref. \cite Himmetoglu:2014:RHC, Eq. (3)
!! \f{eqnarray*}{
!!       n^{I,\sigma}_{mm^\prime} & =
!!      \langle \phi_{m}^{I} \vert \hat{\rho}^{\sigma} \vert \phi_{m^\prime}^{I} \rangle
!!      \nonumber \\
!!       & = \langle \phi_{m}^{I} \vert \left( \sum_{\mu} \vert \phi_{\mu} \rangle \langle \tilde{\phi}^{\mu} \vert \right) \vert \hat{\rho}^{\sigma} \vert
!!      \left( \sum_{\nu} \vert \tilde{\phi}^{\nu} \rangle \langle \phi_{\nu} \vert \right)
!!      \vert \phi_{m^\prime}^{I} \rangle
!!      \nonumber \\
!!       & = \sum_{\mu \nu} \langle \phi_{m}^{I} \vert \phi_{\mu} \rangle
!!      \langle \tilde{\phi}^{\mu} \vert \hat{\rho}^{\sigma} \vert
!!      \tilde{\phi}^{\nu} \rangle \langle \phi_{\nu}
!!      \vert \phi_{m^\prime}^{I} \rangle
!!      \nonumber \\
!!       & = \sum_{\mu \nu} \langle \phi_{m}^{I} \vert \phi_{\mu} \rangle \rho_{\mu \nu}^{\sigma} \langle \phi_{\nu}
!!      \vert \phi_{m^\prime}^{I} \rangle,
!! \f}
!! where \f$ \vert \phi_{\mu} \rangle \f$ and \f$ \vert \phi_{\nu} \rangle \f$
!! are numerical atomic orbitals of the SIESTA basis,
!! \f$ \vert \phi^{\mu} \rangle \f$ and \f$ \vert \phi^{\nu} \rangle \f$
!! are the duals of the atomic basis, that satisfy
!! \f$ \langle \phi^{\mu} \vert \phi_{\nu} \rangle =
!!     \langle \phi_{\mu} \vert \phi^{\nu} \rangle = \delta_{\mu \nu} \f$.
!! In the previous equation, the one-particle density operator
!! \f$ \hat{\rho}^{\sigma} \f$ is defined as
!! \f{eqnarray*}{
!!   \hat{\rho}^{\sigma} = \sum_{\mathbf{k}} \sum_{i}
!!   \vert \psi_{\mathbf{k} i}^{\sigma} \rangle
!!   f_{\mathbf{k} i}^{\sigma}
!!   \langle \psi_{\mathbf{k} i}^{\sigma} \vert,
!! \f}
!! where \f$ \vert \psi_{\mathbf{k} i}^{\sigma} \rangle \f$ are the orthonormal
!! eigenstates of the Kohn-Sham one-particle Hamiltonian.
!! The LCAO density matrix is then written in terms of the dual LCAO basis,
!! \f$ \rho^{\sigma}_{\mu \nu} = \langle \phi^{\mu} \vert \hat{\rho}^{\sigma}
!! \vert \phi^{\nu} \rangle \f$.
!!
!! These occupations are stored in the pointer "occupation",
!! which has four entries,
!! that would correspond with \f$ {\rm occupation} (m, m^\prime, I, \sigma)
!! \equiv n^{I,\sigma}_{mm^\prime} \f$.
!!
!! It is important to recall that LDA+U projectors,
!! \f$ \vert \phi_{m}^{I} \rangle \f$  , should be
!! quite localized functions.  Otherwise the calculated populations
!! loose their atomic character and physical meaning.
!! In particular, they should not overlap with projectors belonging to
!! neighbour atoms, as specified in page 29 of Ref. \cite Himmetoglu:2014:RHC.
!! If this is the case, the projectors form an orthogonal set of functions:
!! between atoms, the overlap is zero. Within an atom, they are zero as
!! impossed by the inner product of the spherical harmonics.
!! Therefore, due to the orthogonality, the basis of projectors is the same
!! as the basis of their duals.
!!
!! We also compute the traces of the density matrix for every atom
!! containing correlated shells
!! \f$ {\rm tracetot} (I) \equiv n^{I} = \sum_{m} \sum_{\sigma} =
!!     n^{I \sigma}_{m}\f$,
!! where \f$ {\rm traceupup} (I) \equiv n^{I \uparrow \uparrow }_{m} =
!!    n^{I \uparrow \uparrow}_{mm}\f$, and the same for the
!!    \f$\downarrow \downarrow \f$, \f$\uparrow \downarrow \f$, and
!! \f$\downarrow \uparrow \f$.
!!
!! Finally the total number of correlated electrons and the
!! magnetic moment of every individual atom is computed.
!! Note that this magnetic moment is computed considering only the
!! the charge on the LDA+U projectors (i.e. summing only over the correlated
!! shell).
!! Therefore, it might be different from the charge computed and written in the
!! output file under the name of "spin moment:"
!! Following Ref. \cite Bousquet-10, the total number of
!! correlated electrons and the  magnetic moment
!! are computed as
!!
!! \f$ n^{I} = \sum_{m} \left(n^{I \uparrow \uparrow }_{mm} +
!!                        n^{I \downarrow \downarrow }_{mm} \right) \f$
!!
!! \f$ m^{I}_{x} = \sum_{m} \left(n^{I \uparrow \downarrow }_{mm} +
!!                                n^{I \downarrow \uparrow }_{mm} \right) \f$
!!
!! \f$ m^{I}_{y} = \sum_{m} i \left(n^{I \uparrow \downarrow }_{mm} -
!!                                  n^{I \downarrow \uparrow }_{mm} \right) \f$
!!
!! \f$ m^{I}_{z} = \sum_{m} \left(n^{I \uparrow \uparrow }_{mm} -
!!                                n^{I \downarrow \downarrow }_{mm} \right) \f$
!!
!! The computation of the occupations are parallelized,
!! and all the nodes do have access to this information.
!! The same happens with the traces and the magnetic moment
!!


module m_occ_proj

  use precision,     only : dp           ! Double precision
  complex(dp), dimension(:,:,:,:), pointer, save :: occupation => null()
                                         ! Array used to store the
                                         !   occupations of the LDA+U proj.
                                         !   In the notation of Himmetoglu et al
                                         !   Eq. (3) of Ref.
                                         !   International Journal of Quantum
                                         !   Chemistry 114, 14 (2014),
                                         !   it corresponds to
                                         !   n(m,mprime,I,sigma),
                                         !   that corresponds to the projection
                                         !   of the oneparticle density operator
                                         !   onto the states of localized
                                         !   atomic centers of LDA+U projectors,
                                         !   \vert \phi_{I,m}^{\sigma} \rangle
                                         !   First index of occupation:
                                         !      first LDA+U projector, m
                                         !   Second index of occupation:
                                         !      second LDA+U projector, mprime
                                         !   Third index of occupation:
                                         !      atomic index where the LDA+U
                                         !      projectors are centered
                                         !   Fourth index of occupation:
                                         !      spin index
  complex(dp), dimension(:,:,:,:), pointer, save :: occupation_old => null()
                                         ! Same as occupations, but in the
                                         !    previous step
  complex(dp), dimension(:), pointer, save :: traceupup => null()
                                         ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (up,up) terms of the
                                         !    density matrix
                                         !    n^{I \uparrow \uparrow} =
                                         !      \sum_{m} n^{I \up \up}_{m m}
  complex(dp), dimension(:), pointer, save :: tracedndn => null()
                                         ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (down,down) terms of the
                                         !    density matrix
                                         !    n^{I \downarrow \downarrow} =
                                         !      \sum_{m} n^{I \down \down}_{m m}
  complex(dp), dimension(:), pointer, save :: traceupdn => null()
                                         ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (up,down) terms of the
                                         !    density matrix
                                         !    n^{I \uparrow \downarrow} =
                                         !      \sum_{m} n^{I \up \down}_{m m}
  complex(dp), dimension(:), pointer, save :: tracednup => null()
                                         ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (down,up) terms of the
                                         !    density matrix
                                         !    n^{I \downarrow \uparrow} =
                                         !      \sum_{m} n^{I \down \up}_{m m}
  complex(dp), dimension(:), pointer, save :: tracetot => null()
                                         ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the trace of the
                                         !    density matrix
                                         !    n^{I} =
                                         !      \sum_{\sigma} \sum_{m}
                                         !      n^{I \sigma}_{m m}
  real(dp), dimension(:), pointer, save :: number_corr_electrons => null()
                                         ! Number of correlated electrons
                                         !    in a given atom
  real(dp), dimension(:,:), pointer, save :: magnetic_moment => null()
                                         ! Magnetic moment of a
                                         !   a given atom

CONTAINS


subroutine reset_proj_dftu()
  use alloc,         only : de_alloc
  implicit none

  call de_alloc( occupation    , 'occupation'    , 'occ_proj_dftu' )
  call de_alloc( occupation_old, 'occupation_old', 'occ_proj_dftu' )

  call de_alloc( traceupup, 'traceupup', 'occ_proj_dftu' )
  call de_alloc( tracedndn, 'tracedndn', 'occ_proj_dftu' )
  call de_alloc( traceupdn, 'traceupdn', 'occ_proj_dftu' )
  call de_alloc( tracednup, 'tracednup', 'occ_proj_dftu' )
  call de_alloc( tracetot , 'tracetot',  'occ_proj_dftu' )
  call de_alloc( number_corr_electrons, 'number_corr_electrons','occ_proj_dftu')
  call de_alloc( magnetic_moment, 'magnetic_moment', 'occ_proj_dftu' )

  nullify( occupation, occupation_old )
  nullify( traceupup, tracedndn, traceupdn, tracednup, tracetot )
  nullify( number_corr_electrons, magnetic_moment )
end subroutine reset_proj_dftu


subroutine occ_proj_dftu

  use parallel,      only : Node         ! Local processor number
  use parallel,      only : Nodes        ! Total number of processors in a
                                         !   parallel run
  use alloc,         only : re_alloc     ! Reallocation routines
  use alloc,         only : de_alloc     ! Deallocation routines
  use parallelsubs,  only : LocalToGlobalOrb
                                         ! Converts an orbital index
                                         !   in the local frame
                                         !   to the global frame
  use parallelsubs,  only : GlobalToLocalOrb
                                         ! Converts an orbital index
                                         !   in the global frame
                                         !   to the local frame
  use atm_types,     only : nspecies     ! Number of different atomic species
                                         !   in the simulation box
  use atm_types,     only : species_info ! Derived type with all the info
                                         !   about the radial functions
                                         !   (PAOs, KB projectors,
                                         !   LDA+U proj,
                                         !   VNA potentials, etc)
                                         !   for a given atomic specie
  use atm_types,     only : species      ! Actual array where the
                                         !   previous information is stored
  use atmfuncs,      only : nofis        ! Total number of basis functions for
                                         !   a given atomic species
  use atmfuncs,      only : rcut         ! Returns cutoff radius of
                                         !   Kleynman-Bylander projectors and
                                         !   atomic basis orbitals
  use atmfuncs,      only : orb_gindex   ! Returns the global index of a
                                         !   basis orbital
  use atmfuncs,      only : dftu_gindex  ! Returns the global index of a
                                         !   LDA+U projector
  use radial,        only : rad_func     ! Derived type where all the
                                         !   information related with a radial
                                         !   function is stored
  use siesta_geom,   only : isa          ! Species index of each atom
  use siesta_geom,   only : na_u         ! Number of atoms in the unit cell
  use siesta_geom,   only : na_s         ! Number of atoms in the supercell
  use siesta_geom,   only : scell        ! Lattice vector of the
                                         !   supercell in real space
                                         !   First  index: cartesian direction
                                         !   Second index: supercell vector
  use siesta_geom,   only : xa           ! Atomic positions
  use atomlist,      only : indxua       ! Index of the equivalent atom in the
                                         !   unit cell
  use atomlist,      only : no_s         ! Number of atomic orbitals in the
                                         !   supercell
  use atomlist,      only : no_l         ! Number of atomic orbitals in the
                                         !   unit cell stored in the local node
  use atomlist,      only : lasto        ! Position of last orbital
                                         !   of each atom
  use atomlist,      only : iphorb       ! Orbital index of each  orbital
                                         !   in its atom
  use m_spin,        only : spin         ! Variable containing the
                                         !   information related with the
                                         !   spin configuration
  use sys,           only : die          ! Termination routine
  use neighbour,     only : mneighb      ! Subroutine to compute the
                                         !   number of neighbours
  use neighbour,     only : maxnna       ! Maximum number of neighbours
  use neighbour,     only : iana => jan  ! Atom-index of neighbours
  use neighbour,     only : r2ki => r2ij ! Squared distances to neighbors
  use neighbour,     only : xki  => xij  ! Vector from a given atom
                                         !   to its neighbours
  use sparse_matrices, only : numh       ! Number of nonzero element of each
                                         !   row of the hamiltonian matrix
  use sparse_matrices, only : listh      ! Nonzero hamiltonian-matrix elemen
  use sparse_matrices, only : listhptr   ! Pointer to start of each row
                                         !   of the hamiltonian matrix
  use sparse_matrices, only : Dscf       ! Density matrix
  use matel_m,         only : matel      ! Returns the overlap matrix elements
                                         !   between an atomic orbital and a
                                         !   LDA+U projector
#ifdef MPI
  use m_mpi_utils, only: globalize_sum
#endif


  implicit none

  integer     :: is                      ! Counter for the loop on atomic
                                         !   species
  integer     :: io                      ! Counter for the loop on atomic
                                         !   orbitals
  integer     :: j                       ! Counter for the loop on neighbour
                                         !   atomic orbitals
  integer     :: io_local                ! Counter for the loop on atomic
                                         !   orbitals in the local node
  integer     :: ka                      ! Counter for the loop on atoms
  integer     :: ina                     ! Counter for the loop on neighbour
                                         !   atoms
  integer     :: ino                     ! Counter for the loop on neighbour
                                         !   orbitals
  integer     :: jno                     ! Counter for the second loop on
                                         !   neighbour orbitals
  integer     :: iproj                   ! Counter for the loop LDA+U projectors
  integer     :: mproj                   ! Counter for the loop LDA+U projectors
  integer     :: mprojprime              ! Counter for the loop LDA+U projectors
  integer     :: ix                      ! Counter for the loop on
                                         !   cartesian directions
  integer     :: ispin                   ! Counter for the loop on spins
  integer     :: io_global               ! Global index of the atomic orbital
  integer     :: jo                      ! Neighbour orbital
  integer     :: ioa                     ! Orbital index of the atomic orbital
                                         !   within its atom
  integer     :: ig                      ! Global index of a basis orbital
  integer     :: ind
  integer     :: kg                      ! Global index of a LDA+U projector
  integer     :: dftuidx                 ! Index of the LDA+U projector
  integer     :: kua                     ! Equivalent atom in the unit cell
  integer     :: ks                      ! Kind of atomic species for atom kua
  integer     :: ndftuproj               ! Number of LDA+U projectors for atom
                                         !  kua
  real(dp)    :: rmaxo                   ! Maximum range of the atomic orbitals
  real(dp)    :: rmaxdftu                ! Maximum range of the lda+U projectors
  real(dp)    :: rmax                    ! Total range to search the neighbour
                                         !   of a LDA+U projector
  integer     :: nneig                   ! Number of neighbour atoms
  integer     :: number_neig_orb         ! Number of neighbour orbitals of a
                                         !   given LDA+U projector
  integer     :: ia                      ! Neighbour atom
  real(dp)    :: rki                     ! Distance (in Bohr) between the two
                                         !   neighbour atoms
  real(dp)    :: rci                     ! Cutoff distance of a given atomic
                                         !   orbital
  real(dp)    :: Sik                     ! Overlap between a LDA+U projector
                                         !   and the first atomic orbital
  real(dp)    :: Sjk                     ! Overlap between a LDA+U projector
                                         !   and the second atomic orbital
  complex(dp), dimension(:), pointer :: Dij
                                         ! Value of the density matrix between
                                         !   two atomic orbitals
                                         !   If spin-orbit is activated,
                                         !   these matrix elements are
                                         !   complex (2x2) matrices
  integer, save :: max_number_dftu_proj = 0
                                         ! Maximum number of LDA+U projectors
  integer, save :: maxno = 500           ! Maximum number of basis orbitals
                                         !   overlapping a KB projector
  logical          within                ! Determines whether an atomic orbital
                                         !    and a LDA+U projector overlaps
                                         !    or not
  logical, save :: firstime = .true.     ! First time that this
                                         !   subroutine is called?
  logical, dimension(:), pointer ::  listedall
                                         ! List of all orbitals needed for
                                         !   this node
  logical, dimension(:), pointer ::  listed
                                         ! Determines whether a given orbital
                                         !   in the unit cell is neighbour of
                                         !   another orbital in the supercell
  integer, dimension(:), pointer :: iono ! The nno-eme neighbour orbital
                                         !   of the LDA+U projector is io,
                                         !   where io runs between 1 and
                                         !   the total number of orbitals
                                         !   in the supercell
  integer, dimension(:), pointer :: iano ! The nno-eme neighbour orbital
                                         !   of the LDA+U projector belongs
                                         !   to atom ia, where ia runs between
                                         !   1 and the total number of atoms
                                         !   in the supercell
  real(dp), dimension(:,:), pointer :: xno
                                         ! The relative position between
                                         !   the center of the LDA+U proj.
                                         !   and the center of the nno-eme
                                         !   neighbour orbital is xno
  real(dp), dimension(:,:), pointer :: Ski
                                         ! Overlap between a LDA+U projector
                                         !   and a neighbour atomic orbital
                                         !   First  index: index of the LDA+U pr
                                         !   Second index: nno-eme neighbour
  real(dp), dimension(:,:,:), pointer :: grSki
                                         ! Gradient of the overlap between a
                                         !   LDA+U projector
                                         !   and a neighbour atomic orbital
                                         !   First  index: cartesian direction
                                         !   Second index: index of the LDA+U pr
                                         !   Third  index: nno-eme neighbour
  real(dp), dimension(:,:), pointer :: Di
                                         ! Density matrix row of an orbital
                                         !   in the unit cell
                                         !   First  index: orbital in the
                                         !                 supercell
                                         !   Second index: spin index


  type(species_info),  pointer :: spp
  type(rad_func),      pointer :: pp

#ifdef MPI
  complex(dp), dimension(:,:), pointer :: buffer1 => null() ! Variable used
                                         ! Variable used to reduce the
                                         ! occupations and transfer it to
                                         ! all nodes
#endif


! Start time counter
  call timer( 'occ_proj_dftu', 1 )

! Nullify pointers
  nullify( listedall )
  nullify( listed    )
  nullify( Di        )
  nullify( iano      )
  nullify( iono      )
  nullify( xno       )
  nullify( Ski       )
  nullify( grSki     )

! Initialization and allocation of matrices
  if( firstime ) then

!   Find maximum number of LDA+U projectors on a given atom
    max_number_dftu_proj = 0
    do ka = 1, na_u
      is  = isa(ka)
      spp => species(is)
      max_number_dftu_proj = max(max_number_dftu_proj,spp%nprojsdftu)
    enddo

!   Allocate local array to store the occupations of the
!   LDA+U projectors
    nullify( occupation, occupation_old )
    call re_alloc( occupation, 1, max_number_dftu_proj, 1,max_number_dftu_proj,&
                   1, na_u, 1, spin%Grid, 'occupation', 'occ_proj_dftu' )

    call re_alloc( occupation_old, 1, max_number_dftu_proj, 1, &
                   max_number_dftu_proj,&
                   1, na_u, 1, spin%Grid, 'occupation_old', 'occ_proj_dftu' )
    occupation     = 0.0_dp
    occupation_old = 0.0_dp

!   Allocate array to store the traces of the occupations on the
!   LDA+U projectors
    nullify( traceupup, tracedndn, traceupdn, tracednup, tracetot )
    nullify( number_corr_electrons, magnetic_moment )
    call re_alloc( traceupup, 1 , na_u, 'traceupup', 'occ_proj_dftu' )
    call re_alloc( tracedndn, 1 , na_u, 'tracedndn', 'occ_proj_dftu' )
    call re_alloc( traceupdn, 1 , na_u, 'traceupdn', 'occ_proj_dftu' )
    call re_alloc( tracednup, 1 , na_u, 'tracednup', 'occ_proj_dftu' )
    call re_alloc( tracetot,  1 , na_u, 'tracetot',  'occ_proj_dftu' )
    call re_alloc( number_corr_electrons, 1, na_u, 'number_corr_electrons',&
                   'occ_proj_dftu' )
    call re_alloc( magnetic_moment, 1, 3, 1, na_u, 'magnetic_moment',&
                   'occ_proj_dftu' )
    firstime  = .false.
  endif

! Initialize the traces and the magnetic moments at each SCF step
  traceupup = 0.0_dp
  tracedndn = 0.0_dp
  traceupdn = 0.0_dp
  tracednup = 0.0_dp
  tracetot  = 0.0_dp
  number_corr_electrons = 0.0_dp
  magnetic_moment = 0.0_dp

! Allocate local variable for the density matrix
  nullify( Dij )
  call re_alloc( Dij, 1, spin%grid, 'Dij', 'occ_proj_dftu' )

! Make list of all orbitals needed for this node
  call re_alloc( listedall, 1, no_s,                        &
 &               name='listedall', routine='occ_proj_dftu' )
  listedall(1:no_s) = .false.

  do io_local = 1, no_l
    call LocalToGlobalOrb( io_local, Node, Nodes, io_global )
    listedall(io_global) = .true.
    do j = 1, numh(io_local)
      jo = listh(listhptr(io_local)+j)
      listedall(jo) = .true.
    enddo
  enddo

! Allocate local arrays related with the atomic orbitals that overlap
! with a LDA+U projector
  call re_alloc( iano, 1, maxno,                                  &
 &               name='iano', routine='occ_proj_dftu' )
  call re_alloc( iono, 1, maxno,                                  &
 &               name='iono', routine='occ_proj_dftu' )
  call re_alloc( xno, 1, 3, 1, maxno,                             &
 &               name='xno', routine='occ_proj_dftu' )
  call re_alloc( Ski, 1, max_number_dftu_proj, 1, maxno,          &
 &               name='Ski', routine='occ_proj_dftu' )
  call re_alloc( grSki, 1, 3, 1, max_number_dftu_proj, 1, maxno,  &
 &               name='grSki', routine='occ_proj_dftu' )

! Find maximum range of the atomic orbitals (rmaxo)
! and of the LDA+U projectors (rmaxdftu)
  rmaxo    = 0.0_dp
  rmaxdftu = 0.0_dp

  do is = 1, nspecies

!   Species orbital range
    do io = 1, nofis(is)
      rmaxo = max(rmaxo, rcut(is,io))
    enddo

!   Species DFTU range
    spp => species(is)
    do iproj = 1, spp%n_pjdftunl
      pp => spp%pjdftu(iproj)
      rmaxdftu = max(rmaxdftu, pp%cutoff)
    enddo

  enddo  ! Close loop on the number of atomic species

! Total range to search the neighbours of a LDA+U projector
  rmax = rmaxo + rmaxdftu

! Initialize arrays Di and listed only once
  call re_alloc( Di, 1, no_s, 1, spin%dm,                   &
 &               name='Di', routine='occ_proj_dftu' )
  Di = 0.0_dp

  call re_alloc( listed, 1, no_s,                           &
 &               name='listed', routine='occ_proj_dftu')
  listed = .false.

!! For debugging
!  write(6,'(a,i5)')                                         &
! &  'occ_proj_dftu: nspin = ' ,                             &
! &                  spin%dm
!  write(6,'(a,i5)')                                         &
! &  'occ_proj_dftu: maximum number of LDA+U projectors = ', &
! &                  max_number_dftu_proj
!  write(6,'(a,3f12.5)')                                     &
! &  'occ_proj_dftu: rmaxo, rmaxdftu, rmax = ' ,             &
! &                  rmaxo, rmaxdftu, rmax
!! End debugging

! Initialize occupations
  occupation = cmplx(0.0_dp,0.0_dp,kind=dp)

! Initialize neighb subroutine
  call mneighb( scell, rmax, na_s, xa, 0, 0, nneig )

! Loop on atoms with LDA+U projectors
  do ka = 1, na_s
    kua = indxua(ka)
    ks  = isa(ka)
    spp => species(ks)
    ndftuproj = spp%nprojsdftu
!   We continue only if this atom has a non-vanishing number of LDA+U projectors
    if( ndftuproj == 0 ) cycle

!   Find neighbour atoms within a range equal to the
!   maximum cutoff of the LDA+U projectors +
!   maximum cutoff of the atomic orbitals
    call mneighb( scell, rmax, na_s, xa, ka, 0, nneig )

!   Find neighbour orbitals
    number_neig_orb = 0

!   Loop on neighbour atoms
    do ina = 1, nneig
!     The neighbour atom is the atom ia...
      ia = iana(ina)
!     ... and it belongs to the species is
      is = isa(ia)
!     ... and the distance (in Bohrs) between the two atoms is rki
      rki = sqrt(r2ki(ina))

!     Loop over all the atomic orbitals of the neighbour atom
      do io = lasto(ia-1)+1, lasto(ia)

!       Only calculate if needed locally
        if (listedall(io)) then
!         Identify the index of the orbital within its atom...
          ioa = iphorb(io)
!         ... and its global index between all the radial functions ...
          ig  = orb_gindex(is,ioa)
!         ... and its cutoff distance
          rci = rcut(is,ioa)

!         Find if orbital is within range
          within = .false.
!         Let us check whether this atomic orbital overlaps with
!         any of the LDA+U projectors of the atom ka
          do mproj = 1, ndftuproj
!           Find the index of the LDA+U projector
            dftuidx = spp%pjdftu_index(mproj)
!           Point to the place where all the information on this
!           LDA+U projector is stored
            pp => spp%pjdftu(dftuidx)
!           Check whether the distance between the atom ka and the neighbour
!           atom is larger or smaller than the sum of the cutoffs of the
!           LDA+U projector plus the atomic orbital
            if ( rci + pp%cutoff > rki ) within = .true.
          enddo

!         Find the overlap between the neighbour orbitals and
!         the LDA+U projectors
          if (within) then
!         Check maxno - if too small then increase array sizes
            if ( number_neig_orb .eq. maxno ) then
              maxno = maxno + 10
              call re_alloc( iano, 1, maxno, name='iano',                      &
     &                       copy=.true., routine='occ_proj_dftu' )
              call re_alloc( iono, 1, maxno, name='iono',                      &
     &                       copy=.true., routine='occ_proj_dftu' )
              call re_alloc( xno, 1, 3, 1, maxno, name='xno',                  &
     &                       copy=.true., routine='occ_proj_dftu' )
              call re_alloc( Ski,1, max_number_dftu_proj, 1, maxno, name='Ski',&
     &                       copy=.true., routine='occ_proj_dftu' )
              call re_alloc( grSki, 1, 3, 1, max_number_dftu_proj, 1, maxno,   &
     &                       name='grSki', routine='occ_proj_dftu',            &
     &                       copy=.true. )
            endif
            number_neig_orb = number_neig_orb + 1
!           The nno-eme neighbour orbital of the LDA+U projector is io,
!           where io runs between 1 and the total number of orbitals
!           in the supercell
            iono( number_neig_orb ) = io

!           The nno-eme neighbour orbital of the LDA+U projector belongs
!           to atom ia,
!           where ia runs between 1 and the total number of atoms
!           in the supercell
            iano( number_neig_orb ) = ia

!           The relative position between the center of the LDA+U proj.
!           and the center of the nno-eme neighbour orbital
!           is xno
            do ix = 1, 3
              xno( ix, number_neig_orb ) = xki(ix,ina)
            enddo

!           Here we compute the overlap between a LDA+U projector
!           and a neighbour atomic orbital
            do mproj = 1, ndftuproj
              kg  = dftu_gindex(ks,mproj)
              call matel( 'S', kg, ig, xki(1:3,ina),        &
 &                            Ski(mproj,number_neig_orb),       &
 &                            grSki(1:3,mproj,number_neig_orb) )
            enddo

          endif  ! If on orbitals within range

        endif ! If on orbitals within the local node

      enddo ! End loop on the orbitals of the neighbour atom

    enddo ! End loop on neighbour atoms

!   Loop on neighbour orbitals
    do ino = 1, number_neig_orb

!     The nno-eme neighbour orbital of the LDA+U projector is io,
!     where io runs between 1 and the total number of orbitals
!     in the supercell
      io = iono(ino)

!     The nno-eme neighbour orbital of the LDA+U projector belongs
!     to atom ia,
!     where ia runs between 1 and the total number of atoms
!     in the supercell
      ia = iano(ino)

!     Only atoms in the unit cell
      if ( ia > na_u ) cycle

!     Only local orbitals
      call GlobalToLocalOrb( io, Node, Nodes, io_local )
      if ( io_local < 1 ) cycle

!     Scatter filter/density matrix row of orbital io
      do j = 1, numh(io_local)
        ind = listhptr(io_local) + j
        jo = listh(ind)
!       Set that the orbital io_local is a neighbour of orbital jo
        listed(jo) = .true.
        do ispin = 1, spin%dm
          Di(jo,ispin) = Di(jo,ispin) + Dscf(ind,ispin)
        enddo
      enddo

!     Find matrix elements with other neighbour orbitals
      do jno = 1, number_neig_orb

!       The nno-eme neighbour orbital of the LDA+U projector is io,
!       where io runs between 1 and the total number of orbitals
!       in the supercell
        jo = iono(jno)

!       Follow only if orbitals io and jo are linked by the density matrix
        if ( listed(jo) ) then

!         Loop on LDA+U projectors of atom ka
          do mproj = 1, ndftuproj

!           Find the overlap between projector and the first atomic orbital
            Sik = Ski(mproj,ino)

!           Second loop on LDA+U projectors of atom ka
            do mprojprime = 1, ndftuproj

!             Find the overlap between projector and the second atomic orbital
              Sjk = Ski(mprojprime,jno)

!             Transform the density matrix between the two atomic orbital
!             in a complex (2x2) matrix
!             When SO coupling is considered, the density matrix and the
!             Hamiltonian must be globally Hermitian
!             (see last sentence of Section 7 of the
!             technical SIESTA paper in JPCM 14, 2745 (2002).
!             However, in diag3k, the sign of the imaginary part of the
!             (up,down) matrix elements is changed:
!             Dnew(ind,1) = Dnew(ind,1) + real(D11,dp)
!             Dnew(ind,2) = Dnew(ind,2) + real(D22,dp)
!             Dnew(ind,3) = Dnew(ind,3) + real(D12,dp)
!             Dnew(ind,4) = Dnew(ind,4) - aimag(D12)
!             Dnew(ind,5) = Dnew(ind,5) + aimag(D11)
!             Dnew(ind,6) = Dnew(ind,6) + aimag(D22)
!             Dnew(ind,7) = Dnew(ind,7) + real(D21,dp)
!             Dnew(ind,8) = Dnew(ind,8) + aimag(D21)
!             In the subroutines to compute the corresponding DFT+U
!             matrix elements:
!             - We change locally the sign of this imaginary part, so the
!               density matrix recovers all its properties
!             - The occupations are computed with the "good" DM.
              Dij(1) =cmplx( Di(jo,1), Di(jo,5),           kind = dp )
              Dij(2) =cmplx( Di(jo,2), Di(jo,6),           kind = dp )
              Dij(3) =cmplx( Di(jo,3), -1.0_dp * Di(jo,4), kind = dp )
              Dij(4) =cmplx( Di(jo,7), Di(jo,8),           kind = dp )

              do ispin = 1, spin%grid
!               Compute the product of the density matrix times the two
!               overlaps between the LDA+U projector and the two atomic orbitals
!               as written in the documentation of the header
                occupation(mproj,mprojprime,kua,ispin) =         &
 &                occupation(mproj,mprojprime,kua,ispin) +       &
 &                Dij(ispin) * Sik * Sjk
! &                Dij * Sik * Sjk/(3.0_dp-dble(spin%grid))
              enddo

            enddo

          enddo  ! End loop on the LDA+U projectors of atom ka

        endif  ! End if orbitals io and jo are linked by the density matrix

      enddo  ! End second loop on neighbour orbitals

!     Restore Di and listed
      do j = 1, numh(io_local)
        jo = listh(listhptr(io_local)+j)
        listed(jo) = .false.
        do ispin = 1, spin%dm
          Di(jo,ispin) = 0.0_dp
        enddo
      enddo

    enddo  ! End loop on neighbour orbitals (ino)

  enddo ! End loop on the number of atoms

#ifdef MPI
! Global reduction of occupation
  call re_alloc( buffer1, 1, max_number_dftu_proj, 1, max_number_dftu_proj,  &
 &               name='buffer1', routine = 'occ_proj_dftu')
  do ka = 1, na_u
    is = isa(ka)
    spp => species(is)
    ndftuproj = spp%nprojsdftu
    if( ndftuproj == 0 ) cycle
    do ispin = 1, spin%grid
      call globalize_sum(occupation(1:ndftuproj,1:ndftuproj,ka,ispin),       &
 &                       buffer1(1:ndftuproj,1:ndftuproj))
      occupation(1:ndftuproj,1:ndftuproj,ka,ispin) =                         &
 &         buffer1(1:ndftuproj,1:ndftuproj)
    enddo
  enddo
  call de_alloc(buffer1, name='buffer1', routine = 'occ_proj_dftu')
#endif

! Compute the traces of the occupations
  do ka = 1, na_u
    ks  = isa(ka)
    spp => species(ks)
    ndftuproj = spp%nprojsdftu
    do iproj = 1, ndftuproj
      traceupup(ka)  = traceupup(ka) + occupation(iproj,iproj,ka,1)
      tracedndn(ka)  = tracedndn(ka) + occupation(iproj,iproj,ka,2)
      traceupdn(ka)  = traceupdn(ka) + occupation(iproj,iproj,ka,3)
      tracednup(ka)  = tracednup(ka) + occupation(iproj,iproj,ka,4)
      tracetot(ka)   = tracetot(ka)  + occupation(iproj,iproj,ka,1) + &
 &                                     occupation(iproj,iproj,ka,2)
    enddo
    number_corr_electrons(ka) = real(traceupup(ka) + tracedndn(ka), kind = dp)
    magnetic_moment(1,ka) = real(traceupdn(ka) + tracednup(ka), kind = dp)
    magnetic_moment(2,ka) = real((traceupdn(ka) - tracednup(ka)) * &
 &                               cmplx(0.0_dp, 1.0_dp, kind=dp) )
    magnetic_moment(3,ka) = real(traceupup(ka) - tracedndn(ka), kind = dp)
  enddo

!! For debugging
!  do ka = 1, na_u
!    is = isa(ka)
!    spp => species(is)
!    ndftuproj = spp%nprojsdftu
!    if( ndftuproj == 0 ) cycle
!    do mproj = 1, ndftuproj
!      do mprojprime = 1, ndftuproj
!        do ispin = 1, spin%grid
!          write(6,'(a,6i5,2f12.5)')                                           &
! &          'occ_proj_dftu : Node, Nodes, ka, mproj, mprojprime, spin, occ =',&
! &                           Node, Nodes, ka, mproj, mprojprime, ispin,       &
! &                           occupation(mproj,mprojprime,ka,ispin)
!        enddo
!      enddo
!    enddo
!    write(6,'(a,3i5,2f12.5)')' occ_proj_dftu: Node, Nodes, ka, traceupup = ' , &
! &    Node, Nodes, ka, traceupup(ka)
!    write(6,'(a,3i5,2f12.5)')' occ_proj_dftu: Node, Nodes, ka, tracedndn = ' , &
! &    Node, Nodes, ka, tracedndn(ka)
!    write(6,'(a,3i5,2f12.5)')' occ_proj_dftu: Node, Nodes, ka, traceupdn = ' , &
! &    Node, Nodes, ka, traceupdn(ka)
!    write(6,'(a,3i5,2f12.5)')' occ_proj_dftu: Node, Nodes, ka, tracednup = ' , &
! &    Node, Nodes, ka, tracednup(ka)
!    write(6,'(a,3i5,2f12.5)')' occ_proj_dftu: Node, Nodes, ka, tracetot  = ' , &
! &    Node, Nodes, ka, tracetot(ka)
!    write(6,'(a,3i5,f12.5)')' occ_proj_dftu: Node, Nodes, ka, number_cor= ' , &
! &    Node, Nodes, ka, number_corr_electrons(ka)
!    write(6,'(a,3i5,3f12.5)')' occ_proj_dftu: Node, Nodes, ka, magnetic  = ' , &
! &    Node, Nodes, ka, magnetic_moment(:,ka)
!  enddo
!!  call die()
!! End debugging

! Deallocate local memory
  call de_alloc( grSki, name='grSki', routine='occ_proj_dftu' )
  call de_alloc( Ski,   name='Ski', routine='occ_proj_dftu' )
  call de_alloc( xno,   name='xno', routine='occ_proj_dftu' )
  call de_alloc( iono,  name='iono', routine='occ_proj_dftu' )
  call de_alloc( iano,  name='iano', routine='occ_proj_dftu' )

  call de_alloc( listedall, name='listedall', routine='occ_proj_dftu' )
  call de_alloc( listed,    name='listed'   , routine='occ_proj_dftu' )
  call de_alloc( Dij,       name='Dij'      , routine='occ_proj_dftu' )
  call de_alloc( Di,        name='Di'       , routine='occ_proj_dftu' )

! Stop time counter
  call timer( 'occ_proj_dftu', 2 )

end subroutine occ_proj_dftu

end module m_occ_proj
