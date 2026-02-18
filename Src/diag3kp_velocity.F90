! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

#ifdef MPI
subroutine diag3kp_velocity( spin, no_l, no_u, no_s, nnz, &
    ncol, ptr, col, H, S, getD, &
    qtot, temp, e1, e2, xij, &
    nk, kpoint, wk, eo, qo, DM, EDM, ef, &
    Entropy, Hk, Sk, psi, aux, &
    occtol, iscf, neigwanted)
  
  ! *********************************************************************
  ! Subroutine to calculate the eigenvalues and eigenvectors, density
  ! and energy-density matrices, and occupation weights of each 
  ! eigenvector, for given Hamiltonian and Overlap matrices with spin-orbit
  ! K-sampling version.
  ! This routine also calculates the velocities for the bands
  ! and shifts the eigenvalues according to the 
  ! Created by Nick Papior, 2020
  ! Uses parallelisation over K points instead of parallelisation 
  ! within them.
  ! **************************** INPUT **********************************
  ! type(t_spin) spin           : Spin type
  ! integer no_l                : Number of basis orbitals in unit cell
  !                               local to this processor
  ! integer no_u                : Number of basis orbitals in unit cell
  ! integer no_s                : Number of basis orbitals in supercell
  ! integer nnz                 : Maximum number of orbitals interacting
  ! integer ncol(no_l)          : Number of nonzero elements of each row 
  !                               of hamiltonian/density matrix locally
  ! integer ptr(no_l)           : Pointer to each row (-1) of the
  !                               hamiltonian matrix locally
  ! integer col(nnz)            : Nonzero hamiltonian-matrix element  
  !                               column indexes for each matrix row
  ! real*8  H(nnz,spin%H)       : Hamiltonian in sparse form
  ! real*8  S(nnz)              : Overlap in sparse form
  ! logical getD                : Find occupations and density matrices?
  ! real*8  qtot                : Number of electrons in unit cell
  ! real*8  temp                : Electronic temperature 
  ! real*8  e1, e2              : Energy range for density-matrix states
  !                               (to find local density of states)
  !                               Not used if e1 > e2
  ! real*8  xij(3,nnz)          : Vectors between orbital centers (sparse)
  !                               (not used if only gamma point)
  ! integer nk                  : Number of k points
  ! real*8  kpoint(3,nk)        : k point vectors
  ! real*8  wk(nk)              : k point weights (must sum one)
  ! real*8  occtol              : Occupancy threshold for DM build
  ! integer neigwanted          : maximum number of eigenvalues wanted
  ! *************************** OUTPUT **********************************
  ! real*8 eo(maxo*spin%spinor,nk)   : Eigenvalues
  ! ******************** OUTPUT (only if getD=.true.) *******************
  ! real*8 qo(maxo*spin%spinor,nk)   : Occupations of eigenstates
  ! real*8 Dnew(maxnd,spin%DM)    : Output Density Matrix
  ! real*8 Enew(maxnd,spin%EDM)    : Output Energy-Density Matrix
  ! real*8 ef                   : Fermi energy
  ! real*8 Entropy              : Electronic entropy
  ! *************************** AUXILIARY *******************************
  ! complex*16 Hk(2,no_u,2,no_u) : Auxiliary space for the hamiltonian matrix
  ! complex*16 Sk(2,no_u,2,no_u) : Auxiliary space for the overlap matrix
  ! complex*16 psi(2,no_u,no_u*2) : Auxiliary space for the eigenvectors
  ! real*8 aux(2,no_u*6)          : Extra auxiliary space
  ! *************************** UNITS ***********************************
  ! xij and kpoint must be in reciprocal coordinates of each other.
  ! temp and H must be in the same energy units.
  ! eo, Enew and ef returned in the units of H.
  ! *************************** PARALLEL ********************************
  ! The auxiliary arrays are now no longer symmetric and so the order
  ! of referencing has been changed in several places to reflect this.
  ! Note : It is assumed in a couple of places that the sparsity of
  ! H/S and Dscf are the same (which is the case in the code at present)
  ! and therefore the only one set of pointer arrays are globalised.
  ! *********************************************************************
  !
  !  Modules
  !
  use precision
  use sys
  use units, only: eV, Ang
  use parallel,     only : Node, Nodes
  use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
  use mpi_siesta
  use m_fermid,     only : fermid, stepf
  use alloc,        only : re_alloc, de_alloc
  use t_spin, only: tSpin
  use velocity_shift_m, only: velocity_shift, calc_velocity_current
  use velocity_shift_m, only: velocity_results, velocity_results_print
  use siesta_geom, only: ucell, cell_periodic
  use intrinsic_missing, only: MODP

  implicit none

  integer :: MPIerror

  type(tSpin), intent(in) :: spin

  ! Pass system size information
  ! Local, unit-cell, supercell
  integer, intent(in) :: no_l, no_u, no_s

  ! K-point information
  integer, intent(in) :: nk
  real(dp), intent(in) :: kpoint(3,nk), wk(nk)

  ! Now pass sparse patterns
  integer, intent(in) :: nnz
  integer, intent(in) :: ncol(no_l), ptr(no_l), col(nnz)

  ! Matrices
  real(dp), intent(in) :: H(nnz,spin%H), S(nnz), xij(3,nnz)
  real(dp), intent(inout) :: DM(nnz,spin%DM), EDM(nnz,spin%EDM)

  ! Sought charges per spin
  real(dp), intent(in) :: Qtot
  ! Fermi-levels
  real(dp), intent(inout) :: Ef
  ! Energy range of DM
  real(dp), intent(in) :: E1, E2
  
  ! Eigenvalues and occupations (charges)
  real(dp), intent(inout) :: qo(no_u*spin%spinor,nk), eo(no_u*spin%spinor,nk)
  real(dp), intent(in) :: Occtol, Temp
  ! Calculated quantities
  real(dp), intent(inout) :: Entropy

  ! Control whether DM/EDM is calculated
  logical, intent(in) :: getD

  ! Current SCF step
  integer, intent(in) :: iSCF

  ! Number of calculated eigenstates
  integer, intent(in) :: neigwanted

  ! Auxiliary arrays
  complex(dp), intent(inout), target :: Hk(2,no_u,2,no_u), Sk(2,no_u,2,no_u)
  complex(dp), intent(inout), target :: psi(2,no_u,no_u*2)
  real(dp), intent(inout), target :: aux(2,no_u,6)

  ! Internal variables
  integer :: BNode, ie, ierror, ik, is, lo
  integer :: io, iio, jo, ind, neigneeded
  integer :: no_u2, neigwanted2
  real(dp) :: kxij, t
  complex(dp) :: kph, D11, D22, D12, D21, cp

  !< Current calculation
  real(dp) :: BB_res(3,5) ! I[[x, y, z]] + qd[dq, 0]
#ifdef MPI
  real(dp) :: mpiBB_res(3,5) ! for reduction of I + qd
#endif

  ! Arrays for figuring out the degenerate states
  real(dp), parameter :: deg_EPS = 7.3498067e-06_dp ! 1e-4 eV
  integer :: ndeg
  logical :: deg_setup
  real(dp), pointer :: vdeg(:,:) => null(), Identity(:,:) => null()

  ! Globalized matrices
  integer :: g_nnz
  integer, pointer :: g_ncol(:), g_ptr(:), g_col(:)
  real(dp), pointer :: g_H(:,:), g_S(:), g_xij(:,:)
  real(dp), pointer :: g_DM(:,:), g_EDM(:,:)
  real(dp), pointer :: v(:,:)

  ! Not allocated, only pointers
  complex(dp), pointer :: Dk(:,:,:,:), Ek(:,:,:,:)
  complex(dp), pointer :: dH(:,:,:,:), dS(:,:,:,:)
  real(dp), pointer :: eig_aux(:), SJ(:,:)

#ifdef DEBUG
  call write_debug( '    PRE diag3kp_velocity' )
#endif

  no_u2 = no_u * 2
  neigwanted2 = neigwanted * 2

  call diag3kp_velocity_1d(no_u2, aux(1,1,3), eig_aux)
  call diag3kp_velocity_2d(3,no_u2, aux(1,1,4), SJ)

  ! Globalise sparsity pattern
  call MPI_AllReduce(nnz, g_nnz, 1, MPI_Integer, MPI_Sum, &
      MPI_Comm_World, MPIerror)

  ! Nullify arrays
  nullify(g_ncol, g_ptr, g_col, g_H, g_S, g_xij, g_DM, g_EDM)
  nullify(v)
  
  ! Allocate local memory for global list arrays
  call re_alloc( g_ncol, 1, no_u, name='g_ncol', routine= 'diag3kp_velocity' )
  call re_alloc( g_ptr, 1, no_u, name='g_ptr', routine= 'diag3kp_velocity' )
  call re_alloc( g_col, 1, g_nnz, name='g_col', routine= 'diag3kp_velocity' )
  call re_alloc( g_H, 1, spin%H, 1, g_nnz, name='g_H', routine= 'diag3kp_velocity' )
  call re_alloc( g_S, 1, g_nnz, name='g_S', routine= 'diag3kp_velocity' )
  call re_alloc( g_xij, 1, 3, 1, g_nnz, name='g_xij', routine= 'diag3kp_velocity' )
  call re_alloc( v, 1, 3, 1, no_u2, name='v', routine= 'diag3kp_velocity' )

  ! Create pointers for d* arrays
  ! Note that these arrays may then not be inter-used
  dH => Hk
  dS => Sk
  Dk => Hk
  Ek => Sk

  ! Globalize arrays
  g_ptr(1) = 0
  io = 1
  do while ( io <= no_u )
    
    call WhichNodeOrb(io,Nodes,BNode)

        ! Number of consecutive elements
    ! Generally this should be the block size...
    ik = count_consecutive(io, no_u)
    
    if ( Node == BNode ) then

      ! Do copy on hosting node
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      ind = 0
      do lo = 1, ik
        ! Global node
        jo = io-1+lo
        g_ncol(jo) = ncol(iio-1+lo)
        ind = ind + g_ncol(jo)
        if ( jo < no_u ) g_ptr(jo+1) = g_ptr(jo) + g_ncol(jo)
      end do
      ! Copy data
      g_col(g_ptr(io)+1:g_ptr(io)+ind) = col(ptr(iio)+1:ptr(iio)+ind)
      g_H(:,g_ptr(io)+1:g_ptr(io)+ind) = transpose(H(ptr(iio)+1:ptr(iio)+ind,:))
      g_S(g_ptr(io)+1:g_ptr(io)+ind) = S(ptr(iio)+1:ptr(iio)+ind)
      g_xij(:,g_ptr(io)+1:g_ptr(io)+ind) = xij(:,ptr(iio)+1:ptr(iio)+ind)

    end if
    
    call MPI_Bcast(g_ncol(io),ik,MPI_Integer, BNode, &
        MPI_Comm_World,MPIerror)

    if ( Node /= BNode ) then
      ind = 0
      do jo = io, io + ik - 1
        ind = ind + g_ncol(jo)
        if ( jo < no_u ) g_ptr(jo+1) = g_ptr(jo) + g_ncol(jo)
      end do
    end if

    call MPI_Bcast(g_col(g_ptr(io)+1),ind,MPI_Integer, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_H(1,g_ptr(io)+1),ind*spin%H,MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_S(g_ptr(io)+1),ind,MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_xij(1,g_ptr(io)+1),3*ind,MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)

    io = io + ik
    
  end do

  ! Perform diagonalization loop
  do ik = 1 + Node, nk, Nodes
    
    call setup_k(kpoint(:,ik))

    ! Since we want to calculate the velocities as well we do need the eigenstates
    call cdiag(Hk,Sk,no_u2,no_u2,no_u2,eo(1,ik),psi, &
        neigwanted2,iscf,ierror, -1)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    else if ( ierror < 0 ) then
      call setup_k(kpoint(:,ik))
      call cdiag(Hk,Sk,no_u2,no_u2,no_u2,eo(1,ik),psi, &
          neigwanted2,iscf,ierror, -1)
    end if

    ! Flag we need to setup the sum(dH * dir)
    deg_setup = .true.

    ! Figure out the degenerate states
    ndeg = 1
    do ie = 2, neigwanted2
        
      ! Eigenvalues are sorted, so this will work fine
      if ( eo(ie,ik) - eo(ie-1,ik) < deg_EPS ) then

        ndeg = ndeg + 1
          
      else if ( ndeg > 1 ) then

        ! Decouple ndeg from ie - 1 and ndeg back
        call degenerate_decouple(neigwanted2, eo(1,ik), ndeg, ie - 1)

        ! To debug, one may calculate the velocities
        ! before and after the degenerate decoupling
        ndeg = 1
          
      end if
        
    end do

    ! Decouple the last elements
    if ( ndeg > 1 ) &
        call degenerate_decouple(neigwanted2, eo(1,ik), ndeg, neigwanted2)

    ! After having decoupled the degenerate states we can calculate
    ! the velocities
    call calculate_velocity(neigwanted2, eo(:,ik))

    ! Shift eigenvalues according to the velocity projection onto the field
    call velocity_shift(+1, neigwanted2, eo(:,ik), v)

  end do

  ! Globalise eigenvalues
  do ik = 1,nk
    BNode = mod(ik-1, Nodes)
    call MPI_Bcast(eo(1,ik), neigwanted2, MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
  end do

  ! Check if we are done ................................................
  if ( .not. getD ) goto 999

  ! Find new Fermi energy and occupation weights ........................
  call fermid(2, 1, nk, wk, no_u2, &
      neigwanted2, eo, Temp, qtot, qo, ef, Entropy )

  ! Initialize current
  BB_res(:,:) = 0._dp

  ! Allocate globalized DM and EDM
  call re_alloc( g_DM, 1, g_nnz, 1, spin%DM, name='g_DM', routine= 'diag3kp_velocity' )
  call re_alloc( g_EDM, 1, g_nnz, 1, spin%EDM, name='g_EDM', routine= 'diag3kp_velocity' )

  ! Find weights for local density of states ............................
  if ( e1 < e2 ) then
    
    t = max( temp, 1.d-6 )
!$OMP parallel do default(shared), private(ik,io), firstprivate(t)
    do ik = 1,nk
      do io = 1, neigwanted2
        qo(io,ik) = wk(ik) * &
            ( stepf((eo(io,ik)-e2)/t) - stepf((eo(io,ik)-e1)/t) )
      end do
    end do
!$OMP end parallel do

  end if

  ! Initialize to 0
  g_DM(:,:) = 0._dp
  g_EDM(:,:) = 0._dp

  do ik = 1 + Node, nk, Nodes

    ! In this case we shouldn't determine number of states
    ! depending on occupations. To correctly calculate the current
    ! we really need the unoccupied states.
    call setup_k(kpoint(:,ik))
    call cdiag(Hk,Sk,no_u2,no_u2,no_u2,eig_aux,psi,neigwanted2,iscf,ierror, -1)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    else if ( ierror < 0 ) then
      call setup_k(kpoint(:,ik))
      call cdiag(Hk,Sk,no_u2,no_u2,no_u2,eig_aux,psi,neigwanted2,iscf,ierror, -1)
    end if

    ! Flag we need to setup the sum(dH * dir)
    deg_setup = .true.

    ! Before expanding eigenvectors we need to decouple the degenerate
    ! states so that we don't populate the degeneracies wrongly
    ! This assumes that the decoupling is stable.
    ndeg = 1
    do ie = 2, neigwanted2
        
      ! Eigenvalues are sorted, so this will work fine
      if ( eig_aux(ie) - eig_aux(ie-1) < deg_EPS ) then

        ndeg = ndeg + 1
          
      else if ( ndeg > 1 ) then

        ! Decouple ndeg from ie - 1 and ndeg back
        call degenerate_decouple(neigwanted2, eig_aux, ndeg, ie - 1)

        ! To debug, one may calculate the velocities
        ! before and after the degenerate decoupling
        ndeg = 1
          
      end if
        
    end do

    ! Decouple the last elements
    if ( ndeg > 1 ) &
        call degenerate_decouple(neigwanted2, eig_aux, ndeg, neigwanted2)

    ! If we want to calculate the current and print it, we can do so here
    if ( calc_velocity_current ) then
      call calculate_velocity(neigwanted2, eig_aux)
      call calculate_spin_moment(neigwanted2, SJ)
      call velocity_results(neigwanted2, eig_aux, qo(:,ik), &
          v, SJ, Ef, wk(ik), Temp, BB_res)
    end if

    ! Find maximum eigenvector that is required for this k point and spin
    ! Note that since eo has averaged out the degeneracy eigenvalues this below
    ! block will also group *all* degenerate eigenvalues!
    neigneeded = 1
    do ie = neigwanted2, 1, -1
      if ( abs(qo(ie,ik)) > occtol ) then
        neigneeded = ie
        exit
      end if
    end do

    ! Expand the eigenvectors to the density matrix
    Dk = cmplx(0._dp, 0._dp, dp)
    Ek = cmplx(0._dp, 0._dp, dp)
      
!$OMP parallel default(shared), &
!$OMP&private(ie,io,jo,ind), &
!$OMP&private(kxij,kph,D11,D22,D12,D21,cp)

    ! Add contribution to density matrices of unit-cell orbitals

    ! Global operation to form new density matrix
    do ie = 1, neigneeded
        
!$OMP do
      do io = 1, no_u
        D11 = qo(ie,ik) * psi(1,io,ie)
        D22 = qo(ie,ik) * psi(2,io,ie)
        D12 = D11 * eig_aux(ie)
        D21 = D22 * eig_aux(ie)
        do jo = 1, no_u

          cp = conjg(psi(1,jo,ie))
          Dk(1,jo,1,io) = Dk(1,jo,1,io) + D11 * cp
          Ek(1,jo,1,io) = Ek(1,jo,1,io) + D12 * cp          
          Dk(2,jo,1,io) = Dk(2,jo,1,io) + D22 * cp
          Ek(2,jo,1,io) = Ek(2,jo,1,io) + D21 * cp

          cp = conjg(psi(2,jo,ie))
          Dk(1,jo,2,io) = Dk(1,jo,2,io) + D11 * cp
          Ek(1,jo,2,io) = Ek(1,jo,2,io) + D12 * cp
          Dk(2,jo,2,io) = Dk(2,jo,2,io) + D22 * cp
          Ek(2,jo,2,io) = Ek(2,jo,2,io) + D21 * cp

        end do
      end do
!$OMP end do

    end do
      
!$OMP do
    do io = 1, no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = kpoint(1,ik) * g_xij(1,ind) + kpoint(2,ik) * g_xij(2,ind) + kpoint(3,ik) * g_xij(3,ind)
        kph = exp(cmplx(0._dp, - kxij, dp))

        D11 = Dk(1,jo,1,io) * kph
        D22 = Dk(2,jo,2,io) * kph
        D12 = Dk(1,jo,2,io) * kph
        D21 = Dk(2,jo,1,io) * kph
          
        g_DM(ind,1) = g_DM(ind,1) + real(D11, dp)
        g_DM(ind,2) = g_DM(ind,2) + real(D22, dp)
        g_DM(ind,3) = g_DM(ind,3) + real(D12, dp)
        g_DM(ind,4) = g_DM(ind,4) - aimag(D12)
        g_DM(ind,5) = g_DM(ind,5) + aimag(D11)
        g_DM(ind,6) = g_DM(ind,6) + aimag(D22)
        g_DM(ind,7) = g_DM(ind,7) + real(D21, dp)
        g_DM(ind,8) = g_DM(ind,8) + aimag(D21)

        D11 = Ek(1,jo,1,io) * kph
        D22 = Ek(2,jo,2,io) * kph
        D12 = Ek(1,jo,2,io) * kph
        D21 = Ek(2,jo,1,io) * kph

        g_EDM(ind,1) = g_EDM(ind,1) + real(D11, dp)
        g_EDM(ind,2) = g_EDM(ind,2) + real(D22, dp)
        g_EDM(ind,3) = g_EDM(ind,3) + real(D12, dp)
        g_EDM(ind,4) = g_EDM(ind,4) - aimag(D12)
          
      end do
    end do
!$OMP end do
    
!$OMP end parallel

  end do

  ! Now do bcast of the respective arrays

  do io = 1, no_u
    
    call WhichNodeOrb(io,Nodes,BNode)

    if ( BNode == Node ) then
      call GlobalToLocalOrb(io,BNode,Nodes,iio)
      
      do is = 1, spin%DM
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            DM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      do is = 1, spin%EDM
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            EDM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do

    else

      ! This is because this is *NOT* the receiving node
      ! and hence the recv_buff is not important
      do is = 1, spin%DM
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            eig_aux(1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      do is = 1, spin%EDM
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            eig_aux(1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      
    end if
  end do

  call de_alloc( g_DM, name='g_DM', routine= 'diag3kp_velocity' )
  call de_alloc( g_EDM, name='g_EDM', routine= 'diag3kp_velocity' )

  if ( calc_velocity_current ) then
    call MPI_Reduce(BB_res(1,1), mpiBB_res(1,1), 14, MPI_Double_Precision, &
        MPI_sum,0,MPI_Comm_World,MPIerror)
    call velocity_results_print(ucell, cell_periodic, mpiBB_res)
  end if

  ! Exit point 
999 continue

  ! Clean-up the arrays required for degenerate decoupling
  call de_alloc(vdeg, name='vdeg', routine= 'diag3kp_velocity' )
  call de_alloc(Identity, name='identity', routine= 'diag3kp_velocity' )

  ! Clean-up the velocities
  call de_alloc( v, name='v', routine= 'diag3kp_velocity' )

  ! Clean up memory
  call de_alloc( g_ncol, name='g_ncol', routine= 'diag3kp_velocity' )
  call de_alloc( g_ptr, name='g_ptr', routine= 'diag3kp_velocity' )
  call de_alloc( g_col, name='g_col', routine= 'diag3kp_velocity' )
  call de_alloc( g_H, name='g_H', routine= 'diag3kp_velocity' )
  call de_alloc( g_S, name='g_S', routine= 'diag3kp_velocity' )
  call de_alloc( g_xij, name='g_xij', routine= 'diag3kp_velocity' )

#ifdef DEBUG
  call write_debug( '    POS diag3kp_velocity' )
#endif

contains
  
  subroutine setup_k(k)
    real(dp), intent(in) :: k(3)

    integer :: io, ind, jo

    Hk = cmplx(0._dp, 0._dp, dp)
    Sk = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,kph)
    do io = 1,no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)

        ! Calculate the complex phase
        kph = exp(cmplx(0._dp, - kxij, dp))

        Sk(1,jo,1,io) = Sk(1,jo,1,io) + g_S(ind) * kph
        Sk(2,jo,2,io) = Sk(2,jo,2,io) + g_S(ind) * kph
        Hk(1,jo,1,io) = Hk(1,jo,1,io) + cmplx(g_H(1,ind), g_H(5,ind), dp) * kph
        Hk(2,jo,2,io) = Hk(2,jo,2,io) + cmplx(g_H(2,ind), g_H(6,ind), dp) * kph
        Hk(1,jo,2,io) = Hk(1,jo,2,io) + cmplx(g_H(3,ind), -g_H(4,ind), dp) * kph
        Hk(2,jo,1,io) = Hk(2,jo,1,io) + cmplx(g_H(7,ind), g_H(8,ind), dp) * kph

      end do
    end do
!$OMP end parallel do

  end subroutine setup_k

  subroutine setup_Sk(k)
    real(dp), intent(in) :: k(3)

    integer :: io, ind, jo

    Sk = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,kph)
    do io = 1,no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)

        ! Calculate the complex phase
        kph = exp(cmplx(0._dp, - kxij, dp))

        Sk(1,jo,1,io) = Sk(1,jo,1,io) + g_S(ind) * kph
        Sk(2,jo,2,io) = Sk(2,jo,2,io) + g_S(ind) * kph

      end do
    end do
!$OMP end parallel do

  end subroutine setup_Sk

  subroutine calculate_spin_moment(neig, SJ)
    integer, intent(in) :: neig
    ! Spin texture S[x, y, z]
    real(dp), intent(inout) :: SJ(3,neig)
    complex(dp) :: Sbox(2,2)

    integer :: ie

    call setup_Sk(kpoint(:,ik))

    do ie = 1, neig

      ! S(k) | psi_i >
      call zgemv('N', no_u2, no_u2, cmplx(1._dp, 0._dp, dp), &
          Sk(1,1,1,1), no_u2, psi(1,1,ie), 1, &
          cmplx(0._dp, 0._dp, dp), aux(1,1,1), 1)

      ! Now calculate spin-box
      ! Since the arrays are transposed (they should be < psi | aux
      ! we need to back-transpose in the following calculation
      call zgemm('N', 'C', 2, 2, no_u, cmplx(1._dp, 0._dp, dp), &
          aux(1,1,1), 2, psi(1,1,ie), 2, &
          cmplx(0._dp, 0._dp, dp), Sbox, 2)

      SJ(1,ie) = real(Sbox(1,2) + Sbox(2,1), dp)
      SJ(2,ie) = aimag(Sbox(1,2) - Sbox(2,1))
      SJ(3,ie) = real(Sbox(1,1) - Sbox(2,2), dp)

    end do

  end subroutine calculate_spin_moment
  
  subroutine setup_dk(ix, k)
    integer, intent(in) :: ix
    real(dp), intent(in) :: k(3)

    integer :: io, ind, jo

    dS = cmplx(0._dp, 0._dp, dp)
    dH = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,kph)
    do io = 1, no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)

        ! Calculate the complex phase
        kph = cmplx(0._dp, - g_xij(ix,ind), dp) * exp(cmplx(0._dp, - kxij, dp))

        dS(1,jo,1,io) = dS(1,jo,1,io) + g_S(ind) * kph
        dS(2,jo,2,io) = dS(2,jo,2,io) + g_S(ind) * kph
        dH(1,jo,1,io) = dH(1,jo,1,io) + cmplx(g_H(1,ind), g_H(5,ind), dp) * kph
        dH(2,jo,2,io) = dH(2,jo,2,io) + cmplx(g_H(2,ind), g_H(6,ind), dp) * kph
        dH(1,jo,2,io) = dH(1,jo,2,io) + cmplx(g_H(3,ind), -g_H(4,ind), dp) * kph
        dH(2,jo,1,io) = dH(2,jo,1,io) + cmplx(g_H(7,ind), g_H(8,ind), dp) * kph
        
      end do
    end do
!$OMP end parallel do

  end subroutine setup_dk

  subroutine setup_dk_direction(k)
    use velocity_shift_m, only: dir => velocity_dir
    real(dp), intent(in) :: k(3)

    real(dp) :: fac
    integer :: io, ind, jo

    dS = cmplx(0._dp, 0._dp, dp)
    dH = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,kph,fac)
    do io = 1, no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)

        fac = g_xij(1,ind)*dir(1) + g_xij(2,ind)*dir(2) + g_xij(3,ind)*dir(3)
        ! Calculate the complex phase
        kph = cmplx(0._dp, -fac, dp) * exp(cmplx(0._dp, -kxij, dp))

        dS(1,jo,1,io) = dS(1,jo,1,io) + g_S(ind) * kph
        dS(2,jo,2,io) = dS(2,jo,2,io) + g_S(ind) * kph
        dH(1,jo,1,io) = dH(1,jo,1,io) + cmplx(g_H(1,ind), g_H(5,ind), dp) * kph
        dH(2,jo,2,io) = dH(2,jo,2,io) + cmplx(g_H(2,ind), g_H(6,ind), dp) * kph
        dH(1,jo,2,io) = dH(1,jo,2,io) + cmplx(g_H(3,ind), -g_H(4,ind), dp) * kph
        dH(2,jo,1,io) = dH(2,jo,1,io) + cmplx(g_H(7,ind), g_H(8,ind), dp) * kph

      end do
    end do
!$OMP end parallel do

  end subroutine setup_dk_direction

  ! This subroutine decouples consecutive degenerate states
  !
  ! NOTE this routine may not use the global variable ie, is or ik
  subroutine degenerate_decouple(neig, eig, ndeg, io_end)
    integer, intent(in) :: neig
    real(dp), intent(inout) :: eig(neig)
    integer, intent(in) :: ndeg, io_end
    integer :: io_start, ix, io, jo
    real(dp) :: e_avg

    ! Setup the sum(dH(1:3) * dir(:))
    if ( deg_setup ) then
      call setup_dk_direction(kpoint(:,ik))
      deg_setup = .false.
    end if

    call re_alloc(vdeg, 1, 2, 1, ndeg**2, name='vdeg', routine= 'diag3kp_velocity', &
        copy=.false., shrink=.false.)
    call re_alloc(Identity, 1, 2, 1, ndeg**2, name='identity', routine= 'diag3kp_velocity', &
        copy=.false., shrink=.false.)

    ! Figure out the start and end of the orbitals
    io_start = io_end - ndeg + 1

    ! Calculate the average eigenvalue
    e_avg = 0._dp
    do io = io_start, io_end
      e_avg = e_avg + eig(io)
    end do
    
    ! Average the eigenvalues
    e_avg = e_avg / ndeg

    ! Copy them back to the energy so that they have a unified energy
    ! A decoupling forces them to be aligned since we "mix" the eigenvalues
    do io = io_start, io_end
      eig(io) = e_avg
    end do

    ! Calculate < psi_i | dH - e dS | psi_j >
    Identity(:,1:ndeg**2) = 0._dp
    do io = 1, ndeg
      jo = io_start + io - 1

      ! dH | psi_i >
      call zgemv('N', no_u2, no_u2, cmplx(1._dp, 0._dp, dp), &
          dH(1,1,1,1), no_u2, psi(1,1,jo), 1, &
          cmplx(0._dp, 0._dp, dp), aux(1,1,1), 1)
      ! dH - e dS | psi_i >
      call zgemv('N', no_u2, no_u2, cmplx(-e_avg, 0._dp, dp), &
          dS(1,1,1,1), no_u2, psi(1,1,jo), 1, &
          cmplx(1._dp, 0._dp, dp), aux(1,1,1), 1)

      ! Calculte < psi_: | dH - e dS | psi_i >
      call zgemv('C', no_u2, ndeg, cmplx(1._dp, 0._dp, dp), &
          psi(1,1,io_start), no_u2, aux(1,1,1), 1, &
          cmplx(0._dp, 0._dp, dp), vdeg(1,(io-1)*ndeg+1), 1)

      ! Initialize identity matrix
      Identity(1,(io-1)*ndeg+io) = 1._dp

    end do

    ! Now we have the < psi_j | dH - e dS | psi_i > matrix, we need to diagonalize to find
    ! the linear combinations of the new states
    call cdiag(vdeg, Identity, ndeg, ndeg, ndeg, aux, Hk, ndeg, 1, ierror, -1)
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation in degenerate sub-space')
    end if
      
    ! Re-create the new states
    ! Now Hk contains the linear combination of the states that
    ! should decouple the degenerate subspace
    call zgemm('N', 'N', no_u2, ndeg, ndeg, cmplx(1._dp, 0._dp, dp), &
        psi(1,1,io_start), no_u2, Hk(1,1,1,1), ndeg, &
        cmplx(0._dp, 0._dp, dp), Sk(1,1,1,1), no_u2)

    ! Copy back the new states
    call zcopy(ndeg*no_u2, Sk(1,1,1,1), 1, psi(1,1,io_start), 1)
      
  end subroutine degenerate_decouple

  ! This routine calculates the velocities of all the states
  ! It uses the variable ie, so don't use this outside
  subroutine calculate_velocity(neig, eig)
    integer, intent(in) :: neig
    real(dp), intent(in) :: eig(neig)
    integer :: ix, ie
    
    complex(dp), external :: zdotc

    do ix = 1, 3
      
      ! Calculate velocities using the analytic d H(k) / d k_alpha where
      ! alpha is the Cartesian direction
      call setup_dk(ix, kpoint(:,ik))
      
      do ie = 1, neig

        ! dH | psi >
        call zgemv('N', no_u2, no_u2, cmplx(1._dp, 0._dp, dp), &
            dH(1,1,1,1), no_u2, psi(1,1,ie), 1, &
            cmplx(0._dp, 0._dp, dp), aux(1,1,1), 1)
        ! dH - e dS | psi >
        call zgemv('N', no_u2, no_u2, cmplx(-eig(ie), 0._dp, dp), &
            dS(1,1,1,1), no_u2, psi(1,1,ie), 1, &
            cmplx(1._dp, 0._dp, dp), aux(1,1,1), 1)

        ! Now calculate the velocity along ix
        ! < psi | H - e dS | psi >
        v(ix,ie) = real(zdotc(no_u2,psi(1,1,ie),1,aux(1,1,1),1), dp)

      end do
      
    end do

!!$    do ie = 1, neigwanted2
!!$      print '(i5,3(tr1,f10.4))', ik, v(:,ie) / Ang / hbar_Rys * 1.e-12_dp
!!$    end do

  end subroutine calculate_velocity

  subroutine diag3kp_velocity_1d(n,from,to)
    integer, intent(in) :: n
    real(dp), intent(in), target :: from(n)
    real(dp), pointer :: to(:)
    to => from(:)
  end subroutine diag3kp_velocity_1d

  subroutine diag3kp_velocity_2d(n1,n2,from,to)
    integer, intent(in) :: n1,n2
    real(dp), intent(in), target :: from(n1,n2)
    real(dp), pointer :: to(:,:)
    to => from(:,:)
  end subroutine diag3kp_velocity_2d

  ! Returns a consecutive number of contributions
  ! starting from the specified index
  function count_consecutive(no_u,io) result(n)
    integer, intent(in) :: no_u, io
    integer :: n
    ! Local variables
    integer :: io_node, i, i_node

    n = 1
    call WhichNodeOrb(io,Nodes,io_node)
    do i = io + 1 , no_u
      ! if the idx is not present, just return
      call WhichNodeOrb(i,Nodes,i_node)
      if ( io_node /= i_node ) return
      n = n + 1
    end do

  end function count_consecutive

end subroutine diag3kp_velocity


#else

subroutine diag3kp_velocity_dummy()
end subroutine diag3kp_velocity_dummy

#endif
