
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_pot_dftu module.
!! Here, we compute the Hamiltonian matrix elements, and the energy
!! forces, and stress contributions related with
!! the U and the J in the non-collinear version of the LSDA+U Hamiltonian.
!!
!! For the Hamiltonian matrix elements, we follow Eq. (3) and Eq.(6) 
!! of Ref. \cite Bousquet-10.
!! In the case of noncollinear magnetism, both the density and the Hamiltonian
!! matrix elements between two atomic orbitals in the SIESTA basis 
!! \f$ \mu \f$ and \f$ \nu \f$ are expressed as \f$ 2 \times 2 \f$ matrices,
!! according with the two-component spinor formulation.
!! And every entry of the matrix might be complex.
!! Therefore, within SIESTA, the dimension of the Hamiltonian and Density 
!! matrices in the case of non-collinear spin is 8 (the real and the imaginary
!! part of every of the four entries in the matrices).
!!
!! Here, we redefine locally the structure of the matrices,
!! \f{eqnarray*}{
!! \rho_{\mu \nu}=&
!! \begin{pmatrix}
!! \rho^{\uparrow\uparrow}_{\mu \nu} & \rho^{\uparrow\downarrow}_{\mu \nu} \\
!! \rho^{\downarrow\uparrow}_{\mu \nu} & \rho^{\downarrow\downarrow}_{\mu \nu}
!! \end{pmatrix}
!! \equiv \begin{pmatrix}
!! \rm{Dscf\_cmplx\_1} & \rm{Dscf\_cmplx\_3} \\
!! \rm{Dscf\_cmplx\_4} & \rm{Dscf\_cmplx\_2}
!! \end{pmatrix}
!! = \begin{pmatrix}
!! \rm{Dscf(1)} + i \rm{Dscf(5)} & \rm{Dscf(3)} - i \rm{Dscf(4)} \\
!! \rm{Dscf(7)} + i \rm{Dscf(8)} & \rm{Dscf(2)} + i \rm{Dscf(6)}
!! \end{pmatrix}
!! \f}
!!
!! The Hamiltonian matrix elements are split into a Hubbard-like term,
!! that takes the shape of
!! 
!! \f{eqnarray*}{
!! V^{\rm Hubbard, \uparrow \uparrow}_{\nu \mu}=& 
!! \sum_{I} \sum_{m m^{\prime}}
!! \sum_{m^{\prime \prime}m^{\prime \prime \prime}}
!! \left[
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} \vert m^{\prime}, m^{\prime \prime \prime} \rangle 
!! \left( 
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\downarrow \downarrow} +  
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\uparrow \uparrow} \right) -
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} \vert m^{\prime \prime \prime},m^{\prime} \rangle  
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\uparrow \uparrow} 
!! \right] 
!! \langle \phi_{\nu} \vert \phi_{m^{\prime}}^{I} \rangle
!!\langle \phi_{m}^{I} \vert \phi_{\mu} \rangle
!! \f}
!!
!! \f{eqnarray*}{
!! V^{\rm Hubbard, \downarrow \downarrow}_{\nu \mu}=& 
!! \sum_{I} \sum_{m m^{\prime}}
!! \sum_{m^{\prime \prime}m^{\prime \prime \prime}}
!! \left[
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} \vert m^{\prime}, m^{\prime \prime \prime} \rangle 
!! \left( 
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\uparrow \uparrow} +  
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\downarrow \downarrow} \right) -
!! \langle m, m^{\prime \prime} \vert V_{\rm ee} \vert m^{\prime \prime \prime},m^{\prime} \rangle  
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\downarrow \downarrow} 
!! \right] 
!! \langle \phi_{\nu} \vert \phi_{m^{\prime}}^{I} \rangle
!!\langle \phi_{m}^{I} \vert \phi_{\mu} \rangle
!! \f}
!!
!! \f{eqnarray*}{
!! V^{\rm Hubbard, \uparrow \downarrow}_{\nu \mu}=& 
!! \sum_{I} \sum_{m m^{\prime}}
!! \sum_{m^{\prime \prime}m^{\prime \prime \prime}}
!! \left( - \langle m, m^{\prime \prime} \vert V_{\rm ee} \vert m^{\prime \prime \prime}, m^{\prime} \rangle \right)
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\uparrow \downarrow} 
!! \langle \phi_{\nu} \vert \phi_{m^{\prime}}^{I} \rangle
!!\langle \phi_{m}^{I} \vert \phi_{\mu} \rangle
!! \f}
!!
!! \f{eqnarray*}{
!! V^{\rm Hubbard, \downarrow \uparrow}_{\nu \mu}=& 
!! \sum_{I} \sum_{m m^{\prime}}
!! \sum_{m^{\prime \prime}m^{\prime \prime \prime}}
!! \left( - \langle m, m^{\prime \prime} \vert V_{\rm ee} \vert m^{\prime \prime \prime}, m^{\prime} \rangle \right)
!! n_{m^{\prime \prime}m^{\prime \prime \prime}}^{I,\downarrow \uparrow} 
!! \langle \phi_{\nu} \vert \phi_{m^{\prime}}^{I} \rangle
!!\langle \phi_{m}^{I} \vert \phi_{\mu} \rangle
!! \f}
!!
!! At different points in the module we shall make use of the
!! fact that the Hamiltonian and density matrices are globally Hermitian,
!! that means 
!!
!! \f{eqnarray*}{
!!   H_{\nu\mu}^{\beta \alpha} = \left(H_{\mu\nu}^{\alpha \beta}\right)^{\ast}
!! \f}
!!
!! For instance, at the time of computing the energy
!!
!! \f{eqnarray*}{
!! E^{\rm Hubbard} = & \frac{1}{2} \mathrm{Tr} \left( \rho H \right) = \frac{1}{2} \sum_{\mu \nu} \mathrm{Tr} \left[
!! \begin{pmatrix}
!! \rho^{\uparrow\uparrow}_{\mu \nu} & \rho^{\uparrow\downarrow}_{\mu \nu} \\
!! \rho^{\downarrow\uparrow}_{\mu \nu} & \rho^{\downarrow\downarrow}_{\mu \nu}
!! \end{pmatrix}
!! \begin{pmatrix}
!! V^{\uparrow\uparrow}_{\nu \mu} & V^{\uparrow\downarrow}_{\nu \mu} \\
!! V^{\downarrow\uparrow}_{\nu \mu} & V^{\downarrow\downarrow}_{\nu \mu}
!! \end{pmatrix}
!! \right]
!! = \frac{1}{2} \sum_{\mu \nu} \left( \rho^{\uparrow\uparrow}_{\mu \nu} V^{\uparrow\uparrow}_{\nu \mu} +
!!   \rho^{\uparrow\downarrow}_{\mu \nu} V^{\downarrow\uparrow}_{\nu \mu} + 
!!   \rho^{\downarrow\uparrow}_{\mu \nu} V^{\uparrow\downarrow}_{\nu \mu} + 
!!   \rho^{\downarrow\downarrow}_{\mu \nu} V^{\downarrow\downarrow}_{\nu \mu} \right)
!! = \frac{1}{2} \sum_{\mu \nu} \left[ \left( \rho^{\uparrow\uparrow}_{\nu \mu} \right)^{\ast} V^{\uparrow\uparrow}_{\nu \mu} +
!!   \left( \rho^{\downarrow\uparrow}_{\nu \mu} \right)^{\ast}  V^{\downarrow\uparrow}_{\nu \mu} + 
!!   \left( \rho^{\uparrow\downarrow}_{\nu \mu} \right)^{\ast}  V^{\uparrow\downarrow}_{\nu \mu} + 
!!   \left( \rho^{\downarrow\downarrow}_{\nu \mu} \right)^{\ast} V^{\downarrow\downarrow}_{\nu \mu} \right]
!! \f}
!!
!! where \f$ n_{m^{\prime\prime}m^{\prime\prime\prime}}^{I,\uparrow \uparrow}\f$
!! and related terms are the occupancies computed in the module
!! m_occ_proj following Ref. \cite Himmetoglu:2014:RHC, Eq. (3)
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
!! \f$ \rho^{\sigma}_{\mu \nu} = \langle \phi^{\mu} \vert \hat{\rho}^{\sigma} \vert 
!! \phi^{\nu} \rangle \f$.
!!
!! The double counting energy is taken from Eq.(13) of Ref. \cite Bultmark-09,
!! that considers the possibility of non-collinear spins
!!
!! \f{eqnarray*}{
!!   E_{\rm dc}^{\rm FLL} = \frac{1}{4} \sum_{I}
!!     \left[ 2 U^{I} n^{I} (n^{I}-1) -2J^{I} n^{I}
!! \left( \frac{n^{I}}{2} -1 \right) - J \vec{m}^{I} \cdot \vec{m}^{I} \right],
!! \f}
!!
!! where \f$ n^{I} \f$ is the total electronic charge in the correlated shell
!! of the atom \f$ I \f$,
!! and \f$ \vec{m}^{I} \f$ the magnetization contribution.
!! In the case of non-collinear spin
!!
!! \f{eqnarray*}{
!!    n^{I} & = \sum_{m} \left( n_{mm}^{I \uparrow \uparrow} + n_{mm}^{I \downarrow \downarrow} \right),  
!!    \nonumber \\
!!    \vec{m}^{I} & = \left( 
!!    \sum_{m} \left( n_{mm}^{I \uparrow \downarrow} + n_{mm}^{I \downarrow \uparrow} \right),
!!    i \sum_{m} \left( n_{mm}^{I \uparrow \downarrow} - n_{mm}^{I \downarrow \uparrow} \right),
!!    \sum_{m} \left( n_{mm}^{I \uparrow \uparrow} - n_{mm}^{I \downarrow \downarrow} \right).
!!    \right)
!! \f}
!!
!! Taking the functional derivatives with respect the density matrix, the
!! double counting energy contributes with the following potential matrix elements
!!
!! \f{eqnarray*}{
!! V^{\rm dc, \uparrow \uparrow}_{\nu \mu}=& 
!! - \sum_{I} \sum_{m} U^{I} \left(S_{m \mu}^{I} S^{I}_{\nu m} \right) \left( n^{I} - \frac{1}{2} \right)  
!! + \sum_{I} \sum_{m} J^{I} \left(S_{m \mu}^{I} S^{I}_{\nu m} \right) 
!!  \left( \sum_{m^{\prime}} \sum_{\mu^{\prime} \nu^{\prime}} \rho_{\mu^\prime \nu^\prime}^{\uparrow \uparrow} 
!!  S^{I}_{m^\prime \mu^\prime} S^{I}_{\nu^\prime m^\prime} - \frac{1}{2}   \right)
!! \f}
!!
!! \f{eqnarray*}{
!! V^{\rm dc, \downarrow \downarrow}_{\nu \mu}=& 
!! - \sum_{I} \sum_{m} U^{I} \left(S_{m \mu}^{I} S^{I}_{\nu m} \right) \left( n^{I} - \frac{1}{2} \right)  
!! + \sum_{I} \sum_{m} J^{I} \left(S_{m \mu}^{I} S^{I}_{\nu m} \right) 
!!  \left( \sum_{m^{\prime}} \sum_{\mu^{\prime} \nu^{\prime}} \rho_{\mu^\prime \nu^\prime}^{\downarrow \downarrow} 
!!  S^{I}_{m^\prime \mu^\prime} S^{I}_{\nu^\prime m^\prime} - \frac{1}{2}   \right)
!! \f}
!!
!! \f{eqnarray*}{
!! V^{\rm dc, \uparrow \downarrow}_{\nu \mu}=& 
!! + \sum_{I} \sum_{m} J^{I} \left(S_{m \mu}^{I} S^{I}_{\nu m} \right) 
!!  \left( \sum_{m^{\prime}} \sum_{\mu^{\prime} \nu^{\prime}} \rho_{\mu^\prime \nu^\prime}^{\uparrow \downarrow} 
!!  S^{I}_{m^\prime \mu^\prime} S^{I}_{\nu^\prime m^\prime}  \right)
!! \f}
!!
!! \f{eqnarray*}{
!! V^{\rm dc, \downarrow \uparrow}_{\nu \mu}=& 
!! + \sum_{I} \sum_{m} J^{I} \left(S_{m \mu}^{I} S^{I}_{\nu m} \right) 
!!  \left( \sum_{m^{\prime}} \sum_{\mu^{\prime} \nu^{\prime}} \rho_{\mu^\prime \nu^\prime}^{\downarrow \uparrow} 
!!  S^{I}_{m^\prime \mu^\prime} S^{I}_{\nu^\prime m^\prime}  \right)
!! \f}
!!
!! 

module m_pot_dftu

  use precision,     only : dp           ! Double precision

CONTAINS

subroutine dftu_so_hamil_2( H_dftu_so, fal, stressl )
  
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
  use atm_types,     only : dftu_so_integrals_type
                                         ! Derived type for the
                                         !   definition of the on-site
                                         !   four-center-integrals
                                         !   required for LDA+U+Spin orbit
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
  use sparse_matrices, only : maxnh      ! Maximum number of orbitals
                                         !   interacting
                                         !   NOTE: In parallel runs,
                                         !   maxnh changes from node to node
  use matel_m,         only : matel      ! Returns the overlap matrix elements
                                         !   between an atomic orbital and a
                                         !   LDA+U projector 
  use m_energies,      only : E_dftu_so  ! Contribution to the energy coming
                                         !   from non-collinear LDA+U
  use m_energies,      only : E_correc_dc   
                                         ! Correction energy required for
                                         !   the LDA+U+SO calculations
  use m_occ_proj,      only : occupation ! Projection of the Kohn-Sham orbitals
                                         !   into the states of a 
                                         !   localized basis of LDA+U projectors
  use m_occ_proj,      only : traceupup  ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (up,up) terms of the
                                         !    density matrix
                                         !    n^{I \uparrow \uparrow} =
                                         !      \sum_{m} n^{I \up \up}_{m m}
  use m_occ_proj,      only : tracedndn  ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (down,down) terms of the
                                         !    density matrix
                                         !    n^{I \downarrow \downarrow} =
                                         !      \sum_{m} n^{I \down \down}_{m m}
  use m_occ_proj,      only : traceupdn  ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (up,down) terms of the
                                         !    density matrix
                                         !    n^{I \uparrow \downarrow} =
                                         !      \sum_{m} n^{I \up \down}_{m m}
  use m_occ_proj,      only : tracednup  ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the (down,up) terms of the
                                         !    density matrix
                                         !    n^{I \downarrow \uparrow} =
                                         !      \sum_{m} n^{I \down \up}_{m m}
  use m_occ_proj,      only : tracetot   ! Sum over all the magnetic quantum
                                         !    number of the correlated shell
                                         !    of the trace of the
                                         !    density matrix
                                         !    n^{I} =
                                         !      \sum_{\sigma} \sum_{m}
                                         !      n^{I \sigma}_{m m}
  use m_occ_proj,      only : number_corr_electrons   
                                         ! Number of correlated electrons 
                                         !    on a given atom
  use units, only: eV, Ang

#ifdef MPI
  use m_mpi_utils, only: globalize_sum
#endif


  implicit none

  complex(dp), intent(inout) :: H_dftu_so(maxnh,spin%Grid)
                                         ! Hamiltonian matrix elements between
                                         !    atomic orbitals in the basis 
                                         !    set of SIESTA
                                         !    For a given pair of orbitals,
                                         !    mu and nu, these matrix elements
                                         !    might be:
                                         !    Within Spin-Orbit: A matrix (2x2)
                                         !       that means, four entries
                                         !       that are complex numbers
                                         !       In this case, spin%Grid = 4
                                         !    Without Spin-Orbit, but spin pol
                                         !       the off-diagonal terms of the 
                                         !       previous matrices are zero
                                         !       and we keep only the diagonal
                                         !       part. 
                                         !       Only the component up and down
                                         !       will be required, and they are
                                         !       real
                                         !       In this case, spin%Grid = 2
                                         !    Without Spin-Orbit, and no-spin 
                                         !       the diagonal terms are equal.
                                         !       Only one real number required
                                         !       In this case, spin%Grid = 1
  real(dp), intent(inout) :: fal(3,na_u) ! Local copy of atomic forces
  real(dp), intent(inout) :: stressl(3,3)! Local copy of stress tensor
  integer     :: is                      ! Counter for the loop on atomic 
                                         !   species
  integer     :: io                      ! Counter for the loop on atomic
                                         !   orbitals
  integer     :: jneig                   ! Counter for the loop on neighbour
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
  integer     :: l_correlated            ! Angular momentum of the correlated
                                         !    shell
  integer     :: mproj                   ! Counter for the loop LDA+U projectors
  integer     :: mprojprime              ! Counter for the loop LDA+U projectors
  integer     :: mproj2prime             ! Counter for the loop LDA+U projectors
  integer     :: mproj3prime             ! Counter for the loop LDA+U projectors
  integer     :: m                       ! Counter for the loop LDA+U projectors
  integer     :: mprime                  ! Counter for the loop LDA+U projectors
  integer     :: m2prime                 ! Counter for the loop LDA+U projectors
  integer     :: m3prime                 ! Counter for the loop LDA+U projectors
  integer     :: ix                      ! Counter for the loop on 
                                         !   cartesian directions
  integer     :: jx                      ! Counter for the loop on 
                                         !   cartesian directions
  integer     :: idftuso                 ! Counter for loop on 4 center integral
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
  integer     :: ja                      ! Neighbour atom
  integer     :: jua                     ! Equivalent of the neighbour atom in the unit cell
  real(dp)    :: rki                     ! Distance (in Bohr) between the two
                                         !   neighbour atoms
  real(dp)    :: rci                     ! Cutoff distance of a given atomic 
                                         !   orbital
  real(dp)    :: Sik                     ! Overlap between a LDA+U projector
                                         !   and the first atomic orbital
  real(dp)    :: Sjk                     ! Overlap between a LDA+U projector
                                         !   and the second atomic orbital
  real(dp)    :: U                       ! Value of U
  real(dp)    :: J                       ! Value of J
  real(dp)    :: direct_int              ! Direct four center integrals
  real(dp)    :: fock_int                ! Fock (exchange) four center integrals
  real(dp)    :: Cijk                    ! Auxiliary variable to compute the 
                                         !   forces
  real(dp)    :: fik                     ! Auxiliary variable to compute the 
                                         !   forces
  real(dp)    :: volume                  ! Unit cell volume
  complex(dp) ::  Dscf_cmplx_1           ! (up,up) matrix element of the density
                                         !   matrix between two orbital mu,nu
  complex(dp) ::  Dscf_cmplx_2           ! (down,down) matrix element of the 
                                         !   density matrix between two orbital
                                         !   mu,nu
  complex(dp) ::  Dscf_cmplx_3           ! (up,down) matrix element of the 
                                         !   density matrix between two orbital
                                         !   mu,nu
  complex(dp) ::  Dscf_cmplx_4           ! (down,up) matrix element of the 
                                         !   density matrix between two orbital
                                         !   mu,nu
  complex(dp) ::  H_dftu_so_Hubbard(maxnh,4)
                                         ! Hubbard contribution to the 
                                         !   Hamiltonian matrix elements
  complex(dp) ::  H_dftu_so_dc(maxnh,4)
                                         ! Double counting contribution to the 
                                         !   Hamiltonian matrix elements
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
                                         !   First  index: neighbour orbital 
                                         !                 in the supercell
                                         !   Second index: spin index
  complex(dp), dimension(:,:), pointer :: Vi_Hubbard
                                         ! Hubbard contribution to the 
                                         !   Hamiltonian matrix row of an 
                                         !   orbital in the unit cell
                                         !   First  index: neighbour orbital  
                                         !                 in the supercell
                                         !   Second index: spin index
  complex(dp), dimension(:,:), pointer :: Vi_dc
                                         ! Double counting contribution to the 
                                         !   Hamiltonian matrix row of an 
                                         !   orbital in the unit cell
                                         !   First  index: neighbour orbital  
                                         !                 in the supercell
                                         !   Second index: spin index
  real(dp) :: volcel                     ! Function to compute the volume of 
                                         !   unit cell
  external :: timer, volcel

  type(species_info),  pointer :: spp
  type(dftu_so_integrals_type), pointer :: dftuintegrals
  type(rad_func),      pointer :: pp

#ifdef MPI
  real(dp), dimension(:,:), pointer :: buffer1 => null() ! Variable used 
                                         ! Variable used to reduce the 
                                         ! occupations and transfer it to 
                                         ! all nodes
#endif


! Start time counter
  call timer( 'dftu_so_hamil_2', 1 )

! Nullify pointers
  nullify( listedall  )
  nullify( listed     )
  nullify( Di         )
  nullify( Vi_Hubbard )
  nullify( Vi_dc      )
  nullify( iano       )
  nullify( iono       )
  nullify( xno        )
  nullify( Ski        )
  nullify( grSki      )

! First of all, we initialize the LDA+U Hamiltonian and the
! correction energy from the double counting term at every SCF step
  H_dftu_so         = cmplx( 0.0_dp, 0.0_dp, kind=dp )
  H_dftu_so_Hubbard = cmplx( 0.0_dp, 0.0_dp, kind=dp )
  H_dftu_so_dc      = cmplx( 0.0_dp, 0.0_dp, kind=dp )
  E_correc_dc = 0.0_dp

! Find unit cell volume
  volume = volcel( scell ) * na_u / na_s

! Initialization and allocation of matrices
  if( firstime ) then

!   Find maximum number of LDA+U projectors on a given atom
    max_number_dftu_proj = 0
    do ka = 1, na_u
      is  = isa(ka)
      spp => species(is)
      max_number_dftu_proj = max(max_number_dftu_proj,spp%nprojsdftu)
    enddo

    firstime  = .false.
  endif

! Make list of all orbitals needed for this node
  call re_alloc( listedall, 1, no_s,                        &
 &               name='listedall', routine='dftu_so_hamil_2' )
  listedall(1:no_s) = .false.

  do io_local = 1, no_l
    call LocalToGlobalOrb( io_local, Node, Nodes, io_global )
    listedall(io_global) = .true.
    do jneig = 1, numh(io_local)
      jo = listh(listhptr(io_local)+jneig)
      listedall(jo) = .true.
    enddo
  enddo

! Allocate local arrays related with the atomic orbitals that overlap
! with a LDA+U projector
  call re_alloc( iano, 1, maxno,                                  &
 &               name='iano', routine='dftu_so_hamil_2' )
  call re_alloc( iono, 1, maxno,                                  &
 &               name='iono', routine='dftu_so_hamil_2' )
  call re_alloc( xno, 1, 3, 1, maxno,                             &
 &               name='xno', routine='dftu_so_hamil_2' )
  call re_alloc( Ski, 1, max_number_dftu_proj, 1, maxno,          &
 &               name='Ski', routine='dftu_so_hamil_2' )
  call re_alloc( grSki, 1, 3, 1, max_number_dftu_proj, 1, maxno,  &
 &               name='grSki', routine='dftu_so_hamil_2' )

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
 &               name='Di', routine='dftu_so_hamil_2' )
  Di = 0.0_dp

! Initialize arrays Vi_Hubbard, Vi_dc, and listed only once
  call re_alloc( Vi_Hubbard, 1, no_s, 1, spin%grid,         &
 &               name='Vi_Hubbard', routine='dftu_so_hamil_2' )
  Vi_Hubbard = cmplx(0.0_dp, 0.0_dp, kind = dp)

  call re_alloc( Vi_dc, 1, no_s, 1, spin%grid,              &
 &               name='Vi_dc', routine='dftu_so_hamil_2' )
  Vi_dc = cmplx(0.0_dp, 0.0_dp, kind = dp)

  call re_alloc( listed, 1, no_s,                           &
 &               name='listed', routine='dftu_so_hamil_2')
  listed = .false.

!! For debugging
!  write(6,'(a,i5)')                                           &
! &  'dftu_so_hamil_2: nspin       = ' ,                       &
! &                    spin%dm
!  write(6,'(a,i5)')                                           &
! &  'dftu_so_hamil_2: nspin DM    = ' ,                       &
! &                    spin%DM  
!  write(6,'(a,i5)')                                           &
! &  'dftu_so_hamil_2: nspin H     = ' ,                       &
! &                    spin%H   
!  write(6,'(a,i5)')                                           &
! &  'dftu_so_hamil_2: nspin Grid  = ' ,                       &
! &                    spin%grid   
!  write(6,'(a,i5)')                                           &
! &  'dftu_so_hamil_2: nspin Spinor= ' ,                       &
! &                    spin%spinor 
!  write(6,'(a,i5)')                                           &
! &  'dftu_so_hamil_2: maximum number of LDA+U projectors = ', &
! &                    max_number_dftu_proj
!  write(6,'(a,2i5,3f12.5)')                                   &
! &  'dftu_so_hamil_2: Node, Nodes, rmaxo, rmaxdftu, rmax = ', &
! &                    Node, Nodes, rmaxo, rmaxdftu, rmax 
!  call die()
!! End debugging

! Initialize neighb subroutine
  call mneighb( scell, rmax, na_s, xa, 0, 0, nneig )

! Loop on atoms with LDA+U projectors
  do ka = 1, na_s
    kua = indxua(ka)
    ks  = isa(ka)
    spp => species(ks)
    ndftuproj = spp%nprojsdftu

!!   For debugging
!    write(6,'(a,5i5)')                                                        &
! &    'dftu_so_hamil_2: Node, Nodes, atom, species, number LDA+U projectors: ',&
! &                      Node, Nodes, ka, ks, ndftuproj
!    do iproj = 1, ndftuproj
!      write(6,'(a,8i5,2f12.5)')                                               &
! &      'dftu_so_hamil_2: Node, Nodes, iproj, index, n, l, m, gindex: ',      &
! &      Node, Nodes, iproj, spp%pjdftu_index(iproj), spp%pjdftu_n(iproj),     &
! &      spp%pjdftu_l(iproj), spp%pjdftu_m(iproj), spp%pjdftu_gindex(iproj),   &
! &      spp%pjdftunl_U(spp%pjdftu_index(iproj)),                              &
! &      spp%pjdftunl_J(spp%pjdftu_index(iproj))
!      l_correlated = spp%pjdftunl_l(spp%pjdftu_index(iproj))
!      dftuintegrals => spp%dftu_so_integrals(spp%pjdftu_index(iproj))
!      do idftuso = 0, 2 * l_correlated
!         write(6,'(a,3i5,f12.5)')                                             &
! &             'dftu_so_hamil_2: Node, Nodes, index, Slater = ',              &
! &             Node, Nodes, idftuso, dftuintegrals%Slater_F(idftuso)
!      enddo
!
!      do m = 1, 2*l_correlated+1
!        do mprime = 1, 2*l_correlated+1
!          do m2prime = 1, 2*l_correlated+1
!            do m3prime = 1, 2*l_correlated+1
!            write(6,'(a,6i5,f12.5)')                                          &
! &  'dftu_so_hamil_2: Node, Nodes, m, mprime, m2prime, m3prime, vee_integral_real = ',   &
! &            Node, Nodes, m, mprime, m2prime, m3prime,                       &
! &            dftuintegrals%vee_4center_integrals(m, mprime, m2prime, m3prime)/eV
!            enddo
!          enddo
!        enddo
!      enddo
!
!    enddo
!    call die()
!!   End debugging

!   We continue only if this atom has a non-vanishing number of LDA+U projectors
    if( ndftuproj == 0 ) cycle

!   Compute the correction to the energy coming from the double-counting term
!   Loop over all the correlated shells of the atom ka 
!   we are not counting here the "m copies"
    do iproj = 1, spp%n_pjdftunl

!     Identify the U and the J for this shell
      U = spp%pjdftunl_U(iproj)
      J = spp%pjdftunl_J(iproj)

      if (Node .eq. 0 ) then
        if( ka .le. na_u ) then
          E_correc_dc = E_correc_dc +                                    &
 &           0.25_dp * (U-J) * number_corr_electrons(ka)
        endif
      else
        E_correc_dc = 0.0_dp
      endif
    enddo

!!   For debugging
!    write(6,'(a,2i5,f12.5)')                              &
! &    'dftu_so_hamil_2: Node, Nodes, E_correc_dc = ',     &
! &                      Node, Nodes, E_correc_dc/eV
!    call die()
!!   End debugging

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
     &                       copy=.true., routine='dftu_so_hamil_2' )
              call re_alloc( iono, 1, maxno, name='iono',                      &
     &                       copy=.true., routine='dftu_so_hamil_2' )
              call re_alloc( xno, 1, 3, 1, maxno, name='xno',                  &
     &                       copy=.true., routine='dftu_so_hamil_2' )
              call re_alloc( Ski,1, max_number_dftu_proj, 1, maxno, name='Ski',&
     &                       copy=.true., routine='dftu_so_hamil_2' )
              call re_alloc( grSki, 1, 3, 1, max_number_dftu_proj, 1, maxno,   &
     &                       name='grSki', routine='dftu_so_hamil_2',          &
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
      do jneig = 1, numh(io_local)
        ind = listhptr(io_local) + jneig
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

!       The nno-eme neighbour orbital of the LDA+U projector belongs
!       to atom ja,
!       where ia runs between 1 and the total number of atoms
!       in the supercell
        ja  = iano(jno)
        jua = indxua(ja)

!       Follow only if orbitals io and jo are linked by the density matrix
        if ( listed(jo) ) then

!         Store in complex format the density matrix between
!         orbitals mu and nu
!         When SO coupling is considered, the density matrix and the Hamiltonian
!         must be globally Hermitian (see last sentence of Section 7 of the
!         technical SIESTA paper in JPCM 14, 2745 (2002).
!         However, in diag3k, the sign of the imaginary part of the (up,down)
!         matrix elements is changed:
!         Dnew(ind,1) = Dnew(ind,1) + real(D11,dp)
!         Dnew(ind,2) = Dnew(ind,2) + real(D22,dp)
!         Dnew(ind,3) = Dnew(ind,3) + real(D12,dp)
!         Dnew(ind,4) = Dnew(ind,4) - aimag(D12)
!         Dnew(ind,5) = Dnew(ind,5) + aimag(D11)
!         Dnew(ind,6) = Dnew(ind,6) + aimag(D22)
!         Dnew(ind,7) = Dnew(ind,7) + real(D21,dp)
!         Dnew(ind,8) = Dnew(ind,8) + aimag(D21)
!         In the subroutines to compute the corresponding DFT+U matrix elements:
!         - We change locally the sign of this imaginary part, so the
!           density matrix recovers all its properties
!         - The Potential matrix elements are computed with the "good" DM.
 
          Dscf_cmplx_1 = cmplx(Di(jo,1),Di(jo,5),  dp)
          Dscf_cmplx_2 = cmplx(Di(jo,2),Di(jo,6),  dp)
          Dscf_cmplx_3 = cmplx(Di(jo,3),-Di(jo,4), dp)
          Dscf_cmplx_4 = cmplx(Di(jo,7),Di(jo,8),  dp)

!!         For debugging
!          write(6,'(a,4i5,2f12.5)')                                        &
! &          'dftu_so_hamil_2: Node, Nodes, io, jo, Dscf_cmplx 1 = ',       &
! &            Node, Nodes, io, jo, Dscf_cmplx_1
!          write(6,'(a,4i5,2f12.5)')                                        &
! &          'dftu_so_hamil_2: Node, Nodes, io, jo, Dscf_cmplx 2 = ',       &
! &            Node, Nodes, io, jo, Dscf_cmplx_2
!          write(6,'(a,4i5,2f12.5)')                                        &
! &          'dftu_so_hamil_2: Node, Nodes, io, jo, Dscf_cmplx 3 = ',       &
! &            Node, Nodes, io, jo, Dscf_cmplx_3
!          write(6,'(a,4i5,2f12.5)')                                        &
! &          'dftu_so_hamil_2: Node, Nodes, io, jo, Dscf_cmplx 4 = ',       &
! &            Node, Nodes, io, jo, Dscf_cmplx_4
!!         End debugging

!         Loop on LDA+U projectors of atom ka
          do mproj = 1, ndftuproj      

!           Point to the memory array where the four center integrals
!           are stored
            dftuintegrals => spp%dftu_so_integrals(spp%pjdftu_index(mproj))

!           Set the values for U and J for this correlated shell
            U = spp%pjdftunl_U(spp%pjdftu_index(mproj))
            J = spp%pjdftunl_J(spp%pjdftu_index(mproj))

!           Find the overlap between projector and the first atomic orbital
            Sik = Ski(mproj,ino)

!           Compute the contribution to the potential that comes from the
!           double counting term
!           For the (up,up) term
            Vi_dc(jo,1) = Vi_dc(jo,1) -              &
 &              (U * (tracetot(kua)  - 0.5_dp)  -    & 
 &               J * (traceupup(kua) - 0.5_dp) ) * Sik * Ski(mproj,jno) 
!           Compute the contributions to the forces and stress from the 
!           (up,up) of the double counting energy           
            Cijk = -2.0_dp * conjg(Dscf_cmplx_1) *   &
 &              (U * (tracetot(kua)  - 0.5_dp)  -    & 
 &               J * (traceupup(kua) - 0.5_dp) ) * Ski(mproj,jno) 
            do ix = 1, 3
              fik = Cijk * grSki(ix,mproj,ino)
              fal(ix,ia)  = fal(ix,ia) - fik
              fal(ix,kua) = fal(ix,kua) + fik
              do jx = 1, 3
                stressl(jx,ix) = stressl(jx,ix) +    &
 &                              xno(jx,ino)*fik/volume
              enddo
            enddo

!           For the (down,down) term
            Vi_dc(jo,2) = Vi_dc(jo,2) -              &
 &              (U * (tracetot(kua)  - 0.5_dp)  -    & 
 &               J * (tracedndn(kua) - 0.5_dp) ) * Sik * Ski(mproj,jno) 
!           Compute the contributions to the forces and stress from the 
!           (down,down) of the double counting energy           
            Cijk = -2.0_dp * conjg(Dscf_cmplx_2) *   &
 &              (U * (tracetot(kua)  - 0.5_dp)  -    & 
 &               J * (tracedndn(kua) - 0.5_dp) ) * Ski(mproj,jno) 
            do ix = 1, 3
              fik = Cijk * grSki(ix,mproj,ino)
              fal(ix,ia)  = fal(ix,ia) - fik
              fal(ix,kua) = fal(ix,kua) + fik
              do jx = 1, 3
                stressl(jx,ix) = stressl(jx,ix) +    &
 &                              xno(jx,ino)*fik/volume
              enddo
            enddo

!           For the (up,down) term
            Vi_dc(jo,3) = Vi_dc(jo,3) +                             &
 &            J * traceupdn(kua) * Sik * Ski(mproj,jno) 
!           Compute the contributions to the forces and stress from the 
!           (up,down) of the double counting energy           
            Cijk = 2.0_dp * conjg(Dscf_cmplx_3) *    &
 &            J * traceupdn(kua) * Ski(mproj,jno) 
            do ix = 1, 3
              fik = Cijk * grSki(ix,mproj,ino)
              fal(ix,ia)  = fal(ix,ia) - fik
              fal(ix,kua) = fal(ix,kua) + fik
              do jx = 1, 3
                stressl(jx,ix) = stressl(jx,ix) +    &
 &                              xno(jx,ino)*fik/volume
              enddo
            enddo

!           For the (down,up) term
            Vi_dc(jo,4) = Vi_dc(jo,4) +                             &
 &            J * tracednup(kua) * Sik * Ski(mproj,jno) 
!           Compute the contributions to the forces and stress from the 
!           (up,down) of the double counting energy           
            Cijk = 2.0_dp * conjg(Dscf_cmplx_4) *    &
 &            J * tracednup(kua) * Ski(mproj,jno) 
            do ix = 1, 3
              fik = Cijk * grSki(ix,mproj,ino)
              fal(ix,ia)  = fal(ix,ia) - fik
              fal(ix,kua) = fal(ix,kua) + fik
              do jx = 1, 3
                stressl(jx,ix) = stressl(jx,ix) +      &
 &                              xno(jx,ino)*fik/volume
              enddo
            enddo

!           Second loop on LDA+U projectors of atom ka
            do mprojprime = 1, ndftuproj

!             Find the overlap between projector and the second atomic orbital
              Sjk = Ski(mprojprime,jno)

!               Third loop on LDA+U projectors of atom ka
                do mproj2prime = 1, ndftuproj

!                 Fourth loop on LDA+U projectors of atom ka
                  do mproj3prime = 1, ndftuproj

                    direct_int =                                               &
 &                    dftuintegrals%vee_4center_integrals( mproj,              &
 &                                                         mprojprime,         &
 &                                                         mproj2prime,        & 
 &                                                         mproj3prime  )
                    fock_int   =                                               &
 &                    dftuintegrals%vee_4center_integrals( mproj,              &
 &                                                         mproj3prime,        &
 &                                                         mproj2prime,        &
 &                                                         mprojprime)

!                   Define the component (up,up) of the Hubbard potential
                    Vi_Hubbard(jo,1) = Vi_Hubbard(jo,1) +                      &
 &                    (direct_int*(occupation(mproj2prime,mproj3prime,kua,2) + &
 &                                 occupation(mproj2prime,mproj3prime,kua,1))- &
 &                     fock_int  * occupation(mproj2prime,mproj3prime,kua,1))* &
 &                     Sik * Sjk
!                   Compute the contributions to the forces and stress from the 
!                   (up,up) of the double counting energy           
                    Cijk = 2.0_dp * conjg(Dscf_cmplx_1) *                      &
 &                    (direct_int*(occupation(mproj2prime,mproj3prime,kua,2) + &
 &                                 occupation(mproj2prime,mproj3prime,kua,1))- &
 &                     fock_int  * occupation(mproj2prime,mproj3prime,kua,1))* &
 &                     Sjk
                    do ix = 1, 3
                      fik = Cijk * grSki(ix,mproj,ino)
                      fal(ix,ia)  = fal(ix,ia) - fik
                      fal(ix,kua) = fal(ix,kua) + fik
                      do jx = 1, 3
                        stressl(jx,ix) = stressl(jx,ix) +                        &
 &                                      xno(jx,ino)*fik/volume
                      enddo
                    enddo

!                   Define the component (down,down) of the Hubbard potential
                    Vi_Hubbard(jo,2) = Vi_Hubbard(jo,2) +                      &
 &                    (direct_int*(occupation(mproj2prime,mproj3prime,kua,2) + &
 &                                 occupation(mproj2prime,mproj3prime,kua,1))- &
 &                     fock_int  * occupation(mproj2prime,mproj3prime,kua,2))* &
 &                     Sik * Sjk
!                   Compute the contributions to the forces and stress from the 
!                   (down,down) of the double counting energy           
                    Cijk = 2.0_dp * conjg(Dscf_cmplx_2) *                      &
 &                    (direct_int*(occupation(mproj2prime,mproj3prime,kua,2) + &
 &                                 occupation(mproj2prime,mproj3prime,kua,1))- &
 &                     fock_int  * occupation(mproj2prime,mproj3prime,kua,2))* &
 &                     Sjk
                    do ix = 1, 3
                      fik = Cijk * grSki(ix,mproj,ino)
                      fal(ix,ia)  = fal(ix,ia) - fik
                      fal(ix,kua) = fal(ix,kua) + fik
                      do jx = 1, 3
                        stressl(jx,ix) = stressl(jx,ix) +                      &
 &                                      xno(jx,ino)*fik/volume
                      enddo
                    enddo

!                   Define the component (up,down) of the Hubbard potential
                    Vi_Hubbard(jo,3) = Vi_Hubbard(jo,3) -                      &
 &                     fock_int  * occupation(mproj2prime,mproj3prime,kua,3) * &
 &                     Sik * Sjk
!                   Compute the contributions to the forces and stress from the 
!                   (up,down) of the double counting energy           
                    Cijk = 2.0_dp * conjg(Dscf_cmplx_3) *                       &
 &                     (-fock_int * occupation(mproj2prime,mproj3prime,kua,3))* &
 &                     Sjk
                    do ix = 1, 3
                      fik = Cijk * grSki(ix,mproj,ino)
                      fal(ix,ia)  = fal(ix,ia) - fik
                      fal(ix,kua) = fal(ix,kua) + fik
                      do jx = 1, 3
                        stressl(jx,ix) = stressl(jx,ix) +                      &
 &                                      xno(jx,ino)*fik/volume
                      enddo
                    enddo

!                   Define the component (down,up) of the Hubbard potential
                    Vi_Hubbard(jo,4) = Vi_Hubbard(jo,4) -                      &
 &                     fock_int  * occupation(mproj2prime,mproj3prime,kua,4) * &
 &                     Sik * Sjk
!                   Compute the contributions to the forces and stress from the 
!                   (down,up) of the double counting energy           
                    Cijk = 2.0_dp * conjg(Dscf_cmplx_4) *                       &
 &                     (-fock_int * occupation(mproj2prime,mproj3prime,kua,4))* &
 &                     Sjk
                    do ix = 1, 3
                      fik = Cijk * grSki(ix,mproj,ino)
                      fal(ix,ia)  = fal(ix,ia) - fik
                      fal(ix,kua) = fal(ix,kua) + fik
                      do jx = 1, 3
                        stressl(jx,ix) = stressl(jx,ix) +                      &
 &                                      xno(jx,ino)*fik/volume
                      enddo
                    enddo

                  enddo ! End loop on mproj3prime

                enddo   ! End loop on mproj2prime

            enddo       ! End loop on mprojprime

          enddo         ! End loop on mproj (the LDA+U projectors of atom ka)

        endif  ! End if orbitals io and jo are linked by the density matrix

      enddo  ! End second loop on neighbour orbitals

!     Restore Di and listed
      do jneig = 1, numh(io_local)
        ind = listhptr(io_local) + jneig
        jo = listh(ind)
        listed(jo) = .false.
        do ispin = 1, spin%dm
          Di(jo,ispin) = 0.0_dp
        enddo

!       In the equations before we have computed V_{\nu \mu},
!       and we want to store them in the matrix elements of H_{\mu \nu}
!       The transformation explained in the doxygen documentation
!       related with the globally Hermitian behavior of the Hamiltonian
!       and density matrix has to be performed
!       See the end of Sec. 7 of the SIESTA paper
        H_dftu_so_dc(ind,1)      = H_dftu_so_dc(ind,1)      +    &
  &                                conjg(Vi_dc(jo,1))      
        H_dftu_so_Hubbard(ind,1) = H_dftu_so_Hubbard(ind,1) +    &
  &                                conjg(Vi_Hubbard(jo,1))       
        H_dftu_so(ind,1)         = H_dftu_so(ind,1)         +    &
  &                                conjg(Vi_dc(jo,1))       +    &
  &                                conjg(Vi_Hubbard(jo,1))       

        H_dftu_so_dc(ind,2)      = H_dftu_so_dc(ind,2)      +    &
  &                                conjg(Vi_dc(jo,2))      
        H_dftu_so_Hubbard(ind,2) = H_dftu_so_Hubbard(ind,2) +    &
  &                                conjg(Vi_Hubbard(jo,2))       
        H_dftu_so(ind,2)         = H_dftu_so(ind,2)         +    &
  &                                conjg(Vi_dc(jo,2))       +    &
  &                                conjg(Vi_Hubbard(jo,2))       

        H_dftu_so_dc(ind,3)      = H_dftu_so_dc(ind,3)      +    &
  &                                conjg(Vi_dc(jo,4))      
        H_dftu_so_Hubbard(ind,3) = H_dftu_so_Hubbard(ind,3) +    &
  &                                conjg(Vi_Hubbard(jo,4))       
        H_dftu_so(ind,3)         = H_dftu_so(ind,3)         +    &
  &                                conjg(Vi_dc(jo,4))       +    &
  &                                conjg(Vi_Hubbard(jo,4))       

        H_dftu_so_dc(ind,4)      = H_dftu_so_dc(ind,4)      +    &
  &                                conjg(Vi_dc(jo,3))      
        H_dftu_so_Hubbard(ind,4) = H_dftu_so_Hubbard(ind,4) +    &
  &                                conjg(Vi_Hubbard(jo,3))       
        H_dftu_so(ind,4)         = H_dftu_so(ind,4)         +    &
  &                                conjg(Vi_dc(jo,3))       +    &
  &                                conjg(Vi_Hubbard(jo,3))       

        Vi_dc(jo,:)      = cmplx( 0.0_dp, 0.0_dp, kind=dp )
        Vi_Hubbard(jo,:) = cmplx( 0.0_dp, 0.0_dp, kind=dp )
      enddo

    enddo  ! End loop on neighbour orbitals (ino)

  enddo ! End loop on the number of atoms

  E_dftu_so   = 0.0_dp
  do ind = 1, maxnh
!   Note the change in the imaginary part of the (up,down) component,
!   as discussed above
    Dscf_cmplx_1 = cmplx(Dscf(ind,1),Dscf(ind,5), dp)
    Dscf_cmplx_2 = cmplx(Dscf(ind,2),Dscf(ind,6), dp)
    Dscf_cmplx_3 = cmplx(Dscf(ind,3),-Dscf(ind,4), dp)
    Dscf_cmplx_4 = cmplx(Dscf(ind,7),Dscf(ind,8), dp)
!   Compute the energy according to the Equations developed in
!   the doxygen documentation
!   Here, we are computing the trace of the the potential times the 
!   density matrix, considering that both of them are (2x2) matrices
!   and the global hermitian of the density matrix
    E_dftu_so = E_dftu_so +                                                   &
 &    0.5_dp * ( real( H_dftu_so_Hubbard(ind,1)*conjg(Dscf_cmplx_1), dp)   +  &
 &               real( H_dftu_so_Hubbard(ind,2)*conjg(Dscf_cmplx_2), dp)   +  &
 &               real( H_dftu_so_Hubbard(ind,4)*conjg(Dscf_cmplx_4), dp)   +  &
 &               real( H_dftu_so_Hubbard(ind,3)*conjg(Dscf_cmplx_3), dp) ) +  &
 &    0.5_dp * ( real( H_dftu_so_dc(ind,1)*conjg(Dscf_cmplx_1), dp)        +  &
 &               real( H_dftu_so_dc(ind,2)*conjg(Dscf_cmplx_2), dp)        +  &
 &               real( H_dftu_so_dc(ind,4)*conjg(Dscf_cmplx_4), dp)        +  &
 &               real( H_dftu_so_dc(ind,3)*conjg(Dscf_cmplx_3), dp) )
  enddo 

! For debugging
!  do io_local = 1, no_l
!    do jneig = 1, numh(io_local)
!      ind = listhptr(io_local) + jneig
!      jo = listh(ind)
!      do ispin = 1, spin%grid
!        if( abs(real(H_dftu_so(ind,ispin))) .gt. 1.d-4 .or.               &
! &          abs(aimag(H_dftu_so(ind,ispin))) .gt. 1.d-4 )                 &
! &      write(6,'(a,6i7,2f12.5)')                                         &
! &      'dftu_so_hamil_2: Node, Nodes, io, ind, jo, ispin, H_dftu_so = ', &
! &       Node, Nodes, io_local, ind, jo, ispin, H_dftu_so(ind,ispin)
!      enddo 
!    enddo 
!  enddo 
!
!   do ia = 1, 3
!     write(6,'(a,i5,3f12.5)')            &
! &     'dftu_so_hamil_2: ia, stress = ',     &
! &      ia, stressl(:,ia)
!   enddo
!
!   write(6,'(a,2i5,f12.5)')                           &
! &   'dftu_so_hamil_2: Node, Nodes, E_dftu_so    = ', &
! &                     Node, Nodes, E_dftu_so / eV
!
!   write(6,'(a,2i5,f12.5)')                           &
! &   'dftu_so_hamil_2: Node, Nodes, E_correc_dc  = ', &
! &                     Node, Nodes, E_correc_dc / eV
!   
!   do ia = 1, na_u
!     write(6,'(a,3i5,3f12.5)')                          &
! &     'dftu_so_hamil_2: Node, Nodes, ia, fa = ',      &
! &     Node, Nodes, ia, fa(:,ia)
!   enddo 
!  call die()
!! End debugging

! Deallocate local memory
  call de_alloc( grSki, name='grSki' )
  call de_alloc( Ski,   name='Ski' )
  call de_alloc( xno,   name='xno' )
  call de_alloc( iono,  name='iono' )
  call de_alloc( iano,  name='iano' )

  call de_alloc( listedall,  name='listedall'  )
  call de_alloc( listed,     name='listed'     )
  call de_alloc( Di,         name='Di'         )
  call de_alloc( Vi_Hubbard, name='Vi_Hubbard' )
  call de_alloc( Vi_dc,      name='Vi_dc'      )


! Stop time counter
  call timer( 'dftu_so_hamil_2', 2 )

end subroutine dftu_so_hamil_2

end module m_pot_dftu

