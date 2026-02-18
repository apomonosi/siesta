! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

#ifdef MPI
subroutine diagkp_velocity( spin, no_l, no_u, no_s, nnz, &
    ncol, ptr, col, H, S, getD, &
    fixspin, qtot, qs, temp, e1, e2, xij, &
    nk, kpoint, wk, eo, qo, DM, EDM, ef, efs, &
    Entropy, Hk, Sk, psi, aux, &
    occtol, iscf, neigwanted)
  
  ! *********************************************************************
  ! Subroutine to calculate the eigenvalues and eigenvectors, density
  ! and energy-density matrices, and occupation weights of each 
  ! eigenvector, for given Hamiltonian and Overlap matrices (including
  ! spin polarization). K-sampling version.
  ! This routine also calculates the velocities for the bands
  ! and shifts the eigenvalues according to the 
  ! Created from diagkp, written by Nick Papior, 2018
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
  ! logical fixspin             : Spin-dependent fermi-levels
  ! real*8  qtot                : Number of electrons in unit cell
  ! real*8  qs                  : Number of electrons in unit cell per spin
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
  ! real*8 eo(maxo,spin%spinor,nk)   : Eigenvalues
  ! ******************** OUTPUT (only if getD=.true.) *******************
  ! real*8 qo(maxo,spin%spinor,nk)   : Occupations of eigenstates
  ! real*8 Dnew(maxnd,spin%DM)    : Output Density Matrix
  ! real*8 Enew(maxnd,spin%EDM)    : Output Energy-Density Matrix
  ! real*8 ef                   : Fermi energy
  ! real*8 Entropy              : Electronic entropy
  ! *************************** AUXILIARY *******************************
  ! real*8 Hk(2,no_u,no_l) : Auxiliary space for the hamiltonian matrix
  ! real*8 Sk(2,no_u,no_l) : Auxiliary space for the overlap matrix
  ! real*8 psi(2,no_u,no_l)  : Auxiliary space for the eigenvectors
  ! real*8 aux(2,no_u)       : Extra auxiliary space
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
  use parallel,     only : Node, Nodes
  use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
  use mpi_siesta
  use m_fermid,     only : fermid, fermispin, stepf
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
  real(dp), intent(in) :: Qs(spin%spinor), Qtot
  ! Fermi-levels
  real(dp), intent(inout) :: Ef, Efs(spin%spinor)
  ! Energy range of DM
  real(dp), intent(in) :: E1, E2
  
  ! Eigenvalues and occupations (charges)
  real(dp), intent(inout) :: qo(no_u,spin%spinor,nk), eo(no_u,spin%spinor,nk)
  real(dp), intent(in) :: Occtol, Temp
  ! Calculated quantities
  real(dp), intent(inout) :: Entropy

  ! Control whether DM/EDM is calculated and for fixed spin
  logical, intent(in) :: getD, fixspin

  ! Current SCF step
  integer, intent(in) :: iSCF

  ! Number of calculated eigenstates
  integer, intent(in) :: neigwanted

  ! Auxiliary arrays
  real(dp), intent(inout), target :: Hk(2,no_u,no_u), Sk(2,no_u,no_u)
  real(dp), intent(inout), target :: psi(2,no_u,no_u), aux(2,no_u*2)

  ! Internal variables
  integer :: BNode, ie, ierror, ik, is, lo
  integer :: io, iio, jo, ind, neigneeded
  real(dp) :: ee, qe, kxij, ckxij, skxij, t
  real(dp) :: pipj1, pipj2

  !< Current calculation
  real(dp) :: BB_res(5,spin%spinor) ! I[[x, y, z]] + qd[dq, 0]
#ifdef MPI
  real(dp) :: mpiBB_res(5,spin%spinor) ! for reduction of I + qd
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
  real(dp), pointer :: Dk(:,:,:), Ek(:,:,:)
  real(dp), pointer :: dH(:,:,:), dS(:,:,:)
  real(dp), pointer :: eig_aux(:)

#ifdef DEBUG
  call write_debug( '    PRE diagkp_velocity' )
#endif

  ! Globalise sparsity pattern
  call MPI_AllReduce(nnz, g_nnz, 1, MPI_Integer, MPI_Sum, &
      MPI_Comm_World, MPIerror)

  call diagkp_velocity_2d1d(no_u, aux(1,no_u+1), eig_aux)

  ! Nullify arrays
  nullify(g_ncol, g_ptr, g_col, g_H, g_S, g_xij, g_DM, g_EDM)
  nullify(v)
  
  ! Allocate local memory for global list arrays
  call re_alloc( g_ncol, 1, no_u, name='g_ncol', routine= 'diagkp_velocity' )
  call re_alloc( g_ptr, 1, no_u, name='g_ptr', routine= 'diagkp_velocity' )
  call re_alloc( g_col, 1, g_nnz, name='g_col', routine= 'diagkp_velocity' )
  call re_alloc( g_H, 1, g_nnz, 1, spin%spinor, name='g_H', routine= 'diagkp_velocity' )
  call re_alloc( g_S, 1, g_nnz, name='g_S', routine= 'diagkp_velocity' )
  call re_alloc( g_xij, 1, 3, 1, g_nnz, name='g_xij', routine= 'diagkp_velocity' )
  call re_alloc( v, 1, 3, 1, no_u, name='v', routine= 'diagkp_velocity' )

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

    ! Determine the hosting node
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
      g_H(g_ptr(io)+1:g_ptr(io)+ind,:) = H(ptr(iio)+1:ptr(iio)+ind,:)
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
    do is = 1, spin%H
      call MPI_Bcast(g_H(g_ptr(io)+1,is),ind,MPI_Double_Precision, &
          BNode,MPI_Comm_World,MPIerror)
    end do
    call MPI_Bcast(g_S(g_ptr(io)+1),ind,MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_xij(1,g_ptr(io)+1),3*ind,MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)

    io = io + ik
    
  end do

  ! Perform diagonalization loop
  do ik = 1 + Node, nk, Nodes
    do is = 1, spin%spinor
      
      call setup_k(kpoint(:,ik))

      ! Since we want to calculate the velocities as well we do need the eigenstates
      call cdiag(Hk,Sk,no_u,no_u,no_u,eo(1,is,ik),psi,neigwanted,iscf,ierror, -1)

      ! Check error flag and take appropriate action
      if ( ierror > 0 ) then
        call die('Terminating due to failed diagonalisation')
      else if ( ierror < 0 ) then
        call setup_k(kpoint(:,ik))
        call cdiag(Hk,Sk,no_u,no_u,no_u,eo(1,is,ik),psi,neigwanted,iscf,ierror, -1)
      end if

      ! Flag we need to setup the sum(dH * dir)
      deg_setup = .true.

      ! Figure out the degenerate states
      ndeg = 1
      do ie = 2, neigwanted
        
        ! Eigenvalues are sorted, so this will work fine
        if ( eo(ie,is,ik) - eo(ie-1,is,ik) < deg_EPS ) then

          ndeg = ndeg + 1
          
        else if ( ndeg > 1 ) then

          ! Decouple ndeg from ie - 1 and ndeg back
          call degenerate_decouple(neigwanted, eo(1,is,ik), ndeg, ie - 1)

          ! To debug, one may calculate the velocities
          ! before and after the degenerate decoupling
          ndeg = 1
          
        end if
        
      end do

      ! Decouple the last elements
      if ( ndeg > 1 ) &
          call degenerate_decouple(neigwanted, eo(1,is,ik), ndeg, neigwanted)

      ! After having decoupled the degenerate states we can calculate
      ! the velocities
      call calculate_velocity(neigwanted, eo(:,is,ik))

      ! Shift eigenvalues according to the velocity projection onto the field
      call velocity_shift(+1, neigwanted, eo(:,is,ik), v)

    end do
    
  end do

  ! Globalise eigenvalues
  if ( neigwanted /= no_u ) then
    do ik = 1,nk
      do is = 1,spin%spinor
        BNode = mod(ik-1, Nodes)
        call MPI_Bcast(eo(1,is,ik), neigwanted, MPI_Double_Precision, &
            BNode,MPI_Comm_World,MPIerror)
      end do
    end do
  else
    io = no_u * spin%spinor
    do ik = 1,nk
      BNode = mod(ik-1, Nodes)
      call MPI_Bcast(eo(1,1,ik), io, MPI_Double_Precision, &
          BNode,MPI_Comm_World,MPIerror)
    end do
  end if

  ! Check if we are done ................................................
  if ( .not. getD ) goto 999

  ! Find new Fermi energy and occupation weights ........................
  if ( fixspin ) then
    call fermispin( spin%spinor, spin%spinor, nk, wk, no_u, &
        neigwanted, eo, Temp, qs, qo, efs, Entropy )
  else
    call fermid( spin%spinor, spin%spinor, nk, wk, no_u, &
        neigwanted, eo, Temp, qtot, qo, ef, Entropy )
    ! This just means we can use efs independently on the used method
    efs(:) = ef
  end if

  ! Initialize current
  BB_res(:,:) = 0._dp

  ! Allocate globalized DM and EDM
  call re_alloc( g_DM, 1, g_nnz, 1, spin%DM, name='g_DM', routine= 'diagkp_velocity' )
  call re_alloc( g_EDM, 1, g_nnz, 1, spin%EDM, name='g_EDM', routine= 'diagkp_velocity' )

  ! Find weights for local density of states ............................
  if ( e1 < e2 ) then
    
    t = max( temp, 1.d-6 )
!$OMP parallel do default(shared), private(ik,is,io), firstprivate(t)
    do ik = 1,nk
      do is = 1,spin%spinor
        do io = 1,neigwanted
          qo(io,is,ik) = wk(ik) * &
              ( stepf((eo(io,is,ik)-e2)/t) - &
              stepf((eo(io,is,ik)-e1)/t)) * 2.0d0 / spin%spinor
        end do
      end do
    end do
!$OMP end parallel do

  end if

  ! Initialize to 0
  g_DM(:,:) = 0._dp
  g_EDM(:,:) = 0._dp

  do ik = 1 + Node, nk, Nodes
    do is = 1, spin%spinor

      ! In this case we shouldn't determine number of states
      ! depending on occupations. To correctly calculate the current
      ! we really need the unoccupied states.
      call setup_k(kpoint(:,ik))
      call cdiag(Hk,Sk,no_u,no_u,no_u,eig_aux,psi,neigwanted,iscf,ierror, -1)

      ! Check error flag and take appropriate action
      if ( ierror > 0 ) then
        call die('Terminating due to failed diagonalisation')
      else if ( ierror < 0 ) then
        call setup_k(kpoint(:,ik))
        call cdiag(Hk,Sk,no_u,no_u,no_u,eig_aux,psi,neigwanted,iscf,ierror, -1)
      end if

      ! Flag we need to setup the sum(dH * dir)
      deg_setup = .true.

      ! Before expanding eigenvectors we need to decouple the degenerate
      ! states so that we don't populate the degeneracies wrongly
      ! This assumes that the decoupling is stable.
      ndeg = 1
      do ie = 2, neigwanted
        
        ! Eigenvalues are sorted, so this will work fine
        if ( eig_aux(ie) - eig_aux(ie-1) < deg_EPS ) then

          ndeg = ndeg + 1
          
        else if ( ndeg > 1 ) then

          ! Decouple ndeg from ie - 1 and ndeg back
          call degenerate_decouple(neigwanted, eig_aux, ndeg, ie - 1)

          ! To debug, one may calculate the velocities
          ! before and after the degenerate decoupling
          ndeg = 1
          
        end if
        
      end do

      ! Decouple the last elements
      if ( ndeg > 1 ) &
          call degenerate_decouple(neigwanted, eig_aux, ndeg, neigwanted)

      ! If we want to calculate the current and print it, we can do so here
      if ( calc_velocity_current ) then
        call calculate_velocity(neigwanted, eig_aux)
        call velocity_results(neigwanted, eig_aux, qo(:,is,ik), &
            v, Efs(is), wk(ik), Temp, BB_res(:,is))
      end if

      ! Find maximum eigenvector that is required for this k point and spin
      ! Note that since eo has averaged out the degeneracy eigenvalues this below
      ! block will also group *all* degenerate eigenvalues!
      neigneeded = 1
      do ie = neigwanted, 1, -1
        if ( abs(qo(ie,is,ik)) > occtol ) then
          neigneeded = ie
          exit
        end if
      end do

      ! Expand the eigenvectors to the density matrix

      Dk = cmplx(0._dp, 0._dp, dp)
      Ek = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared), &
!$OMP&private(ie,qe,ee,io,jo,pipj1,pipj2,ind), &
!$OMP&private(kxij,ckxij,skxij)

      ! Add contribution to density matrices of unit-cell orbitals
      
      ! Global operation to form new density matrix
      do ie = 1, neigneeded
        
        qe = qo(ie,is,ik)
        ee = qe * eig_aux(ie)
        
!$OMP do
        do io = 1, no_u
          do jo = 1, no_u
            pipj1 = psi(1,io,ie) * psi(1,jo,ie) + psi(2,io,ie) * psi(2,jo,ie)
            pipj2 = psi(1,io,ie) * psi(2,jo,ie) - psi(2,io,ie) * psi(1,jo,ie)
            Dk(1,jo,io) = Dk(1,jo,io) + qe * pipj1
            Dk(2,jo,io) = Dk(2,jo,io) + qe * pipj2
            Ek(1,jo,io) = Ek(1,jo,io) + ee * pipj1
            Ek(2,jo,io) = Ek(2,jo,io) + ee * pipj2
          end do
        end do
!$OMP end do

      end do
      
!$OMP do
      do io = 1, no_u
        do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
          jo = modp(g_col(ind), no_u)
          kxij = kpoint(1,ik) * g_xij(1,ind) + kpoint(2,ik) * g_xij(2,ind) + kpoint(3,ik) * g_xij(3,ind)
          ckxij = cos(kxij)
          skxij = sin(kxij)
          
          g_DM(ind,is) = g_DM(ind,is) + Dk(1,jo,io)*ckxij - Dk(2,jo,io)*skxij
          g_EDM(ind,is) = g_EDM(ind,is) + Ek(1,jo,io)*ckxij - Ek(2,jo,io)*skxij
          
        end do
      end do
!$OMP end do
    
!$OMP end parallel

    end do

  end do

  ! Now do bcast of the respective arrays

  do io = 1, no_u
    
    call WhichNodeOrb(io,Nodes,BNode)

    if ( BNode == Node ) then
      call GlobalToLocalOrb(io,BNode,Nodes,iio)
      
      do is = 1, spin%spinor
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            DM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            EDM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do

    else

      ! This is because this is *NOT* the receiving node
      ! and hence the recv_buff is not important
      do is = 1, spin%spinor
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            aux(1,1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            aux(1,1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      
    end if
  end do

  call de_alloc( g_DM, name='g_DM', routine= 'diagkp_velocity' )
  call de_alloc( g_EDM, name='g_EDM', routine= 'diagkp_velocity' )

  if ( calc_velocity_current ) then
    call MPI_Reduce(BB_res(1,1), mpiBB_res(1,1), 5*spin%spinor, MPI_Double_Precision, &
        MPI_sum,0,MPI_Comm_World,MPIerror)
    call velocity_results_print(spin, ucell, cell_periodic, mpiBB_res)
  end if

  ! Exit point 
999 continue

  ! Clean-up the arrays required for degenerate decoupling
  call de_alloc(vdeg, name='vdeg', routine= 'diagkp_velocity' )
  call de_alloc(Identity, name='identity', routine= 'diagkp_velocity' )

  ! Clean-up the velocities
  call de_alloc( v, name='v', routine= 'diagkp_velocity' )

  ! Clean up memory
  call de_alloc( g_ncol, name='g_ncol', routine= 'diagkp_velocity' )
  call de_alloc( g_ptr, name='g_ptr', routine= 'diagkp_velocity' )
  call de_alloc( g_col, name='g_col', routine= 'diagkp_velocity' )
  call de_alloc( g_H, name='g_H', routine= 'diagkp_velocity' )
  call de_alloc( g_S, name='g_S', routine= 'diagkp_velocity' )
  call de_alloc( g_xij, name='g_xij', routine= 'diagkp_velocity' )

#ifdef DEBUG
  call write_debug( '    POS diagkp_velocity' )
#endif

contains
  
  subroutine setup_k(k)
    real(dp), intent(in) :: k(3)

    integer :: io, ind, jo

    Hk = cmplx(0._dp, 0._dp, dp)
    Sk = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,ckxij,skxij)
    do io = 1,no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)
        
        ckxij = cos(kxij)
        skxij = sin(kxij)
        
        ! Note : sign of complex part changed to match change in order of io/jo
        Sk(1,jo,io) = Sk(1,jo,io) + g_S(ind) * ckxij
        Sk(2,jo,io) = Sk(2,jo,io) - g_S(ind) * skxij
        Hk(1,jo,io) = Hk(1,jo,io) + g_H(ind,is) * ckxij
        Hk(2,jo,io) = Hk(2,jo,io) - g_H(ind,is) * skxij
        
      end do
    end do
!$OMP end parallel do

  end subroutine setup_k
  
  subroutine setup_dk(ix, k)
    integer, intent(in) :: ix
    real(dp), intent(in) :: k(3)

    integer :: io, ind, jo

    dS = cmplx(0._dp, 0._dp, dp)
    dH = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,ckxij,skxij)
    do io = 1, no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)
        
        ckxij = cos(kxij) * g_xij(ix,ind)
        skxij = sin(kxij) * g_xij(ix,ind)
        
        ! Note : sign of complex part changed to match change in order of io/jo
        ! Since the phases are
        !    exp(-i k xij)
        ! we get:
        !  - i xij exp(-i k xij)
        dS(1,jo,io) = dS(1,jo,io) - g_S(ind) * skxij
        dS(2,jo,io) = dS(2,jo,io) - g_S(ind) * ckxij
        dH(1,jo,io) = dH(1,jo,io) - g_H(ind,is) * skxij
        dH(2,jo,io) = dH(2,jo,io) - g_H(ind,is) * ckxij
        
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

!$OMP parallel do default(shared), private(io,jo,ind,kxij,ckxij,skxij,fac)
    do io = 1, no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)

        fac = g_xij(1,ind)*dir(1) + g_xij(2,ind)*dir(2) + g_xij(3,ind)*dir(3)
        ckxij = cos(kxij) * fac
        skxij = sin(kxij) * fac

        ! Note : sign of complex part changed to match change in order of io/jo
        ! Since the phases are
        !    exp(-i k xij)
        ! we get:
        !  - i xij exp(-i k xij)
        dS(1,jo,io) = dS(1,jo,io) - g_S(ind) * skxij
        dS(2,jo,io) = dS(2,jo,io) - g_S(ind) * ckxij
        dH(1,jo,io) = dH(1,jo,io) - g_H(ind,is) * skxij
        dH(2,jo,io) = dH(2,jo,io) - g_H(ind,is) * ckxij

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

    call re_alloc(vdeg, 1, 2, 1, ndeg**2, name='vdeg', routine= 'diagkp_velocity', &
        copy=.false., shrink=.false.)
    call re_alloc(Identity, 1, 2, 1, ndeg**2, name='identity', routine= 'diagkp_velocity', &
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
      call zgemv('N', no_u, no_u, cmplx(1._dp, 0._dp, dp), &
          dH(1,1,1), no_u, psi(1,1,jo), 1, &
          cmplx(0._dp, 0._dp, dp), aux(1,1), 1)
      ! dH - e dS | psi_i >
      call zgemv('N', no_u, no_u, cmplx(-e_avg, 0._dp, dp), &
          dS(1,1,1), no_u, psi(1,1,jo), 1, &
          cmplx(1._dp, 0._dp, dp), aux(1,1), 1)

      ! Calculte < psi_: | dH - e dS | psi_i >
      call zgemv('C', no_u, ndeg, cmplx(1._dp, 0._dp, dp), &
          psi(1,1,io_start), no_u, aux(1,1), 1, &
          cmplx(0._dp, 0._dp, dp), vdeg(1,(io-1)*ndeg+1), 1)

      ! Initialize identity matrix
      Identity(1,(io-1)*ndeg+io) = 1._dp

    end do

    ! Now we have the < psi_j | dH - e dS | psi_i > matrix, we need to diagonalize to find
    ! the linear combinations of the new states
    call cdiag(vdeg(1,1), Identity(1,1), ndeg, ndeg, ndeg, &
        aux, Hk(1,1,1), ndeg, 1, ierror, -1)
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation in degenerate sub-space')
    end if
      
    ! Re-create the new states
    ! Now Hk contains the linear combination of the states that
    ! should decouple the degenerate subspace
    call zgemm('N', 'N', no_u, ndeg, ndeg, cmplx(1._dp, 0._dp, dp), &
        psi(1,1,io_start), no_u, Hk(1,1,1), ndeg, &
        cmplx(0._dp, 0._dp, dp), Sk(1,1,1), no_u)

    ! Copy back the new states
    call zcopy(ndeg*no_u, Sk(1,1,1), 1, psi(1,1,io_start), 1)
      
  end subroutine degenerate_decouple

  ! This routine calculates the velocities of all the states
  ! It uses the variable ie, so don't use this outside
  subroutine calculate_velocity(neig, eig)
    integer, intent(in) :: neig
    real(dp), intent(in) :: eig(neig)
    integer :: ix, ie, io

    real(dp), external :: ddot

    do ix = 1, 3
      
      ! Calculate velocities using the analytic d H(k) / d k_alpha where
      ! alpha is the Cartesian direction
      call setup_dk(ix, kpoint(:,ik))
      
      do ie = 1, neig

        ! dH | psi >
        call zgemv('N', no_u, no_u, cmplx(1._dp, 0._dp, dp), &
            dH(1,1,1), no_u, psi(1,1,ie), 1, &
            cmplx(0._dp, 0._dp, dp), aux(1,1), 1)
        ! dH - e dS | psi >
        call zgemv('N', no_u, no_u, cmplx(-eig(ie), 0._dp, dp), &
            dS(1,1,1), no_u, psi(1,1,ie), 1, &
            cmplx(1._dp, 0._dp, dp), aux(1,1), 1)

        ! Now calculate the velocity along ix
        ! < psi | H - e dS | psi >
        v(ix,ie) = ddot(no_u*2,psi(1,1,ie),1,aux(1,1),1)

      end do
      
    end do

!!$    do ie = 1, neigwanted
!!$      print '(i5,3(tr1,f10.4))', ik, v(:,ie) / Ang / hbar_Rys * 1.e-12_dp
!!$    end do

  end subroutine calculate_velocity

  subroutine diagkp_velocity_2d1d(n,from,to)
    integer, intent(in) :: n
    real(dp), intent(in), target :: from(n)
    real(dp), pointer :: to(:)
    to => from(:)
  end subroutine diagkp_velocity_2d1d

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

end subroutine diagkp_velocity


#else

subroutine diagkp_velocity_dummy()
end subroutine diagkp_velocity_dummy

#endif
