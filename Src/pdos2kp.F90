! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
#ifdef MPI
subroutine pdos2kp( nuo, no, maxuo, maxnh, &
    nuotot, numh, listhptr, listh, H, S, &
    E1, E2, nhist, sigma, &
    xij, indxuo, nk, kpoint, wk, &
    Haux, Saux, psi, eo, dtot, dpr )

  ! **********************************************************************
  ! Find the density of states projected onto the atomic orbitals
  !     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
  ! where n run over all the bands between two given energies
  ! Written by J. Junquera and E. Artacho. Nov' 99
  ! Spin-orbit coupling version by J. Ferrer, October 2007
  ! Huge performance increase by N. Papior, 2022
  ! ****  INPUT  *********************************************************
  ! INTEGER NUO               : Number of atomic orbitals in the unit cell
  ! INTEGER NO                : Number of atomic orbitals in the supercell
  ! INTEGER MAXUO             : Maximum number of atomic orbitals in the unit cell
  ! INTEGER MAXNH             : Maximum number of orbitals interacting
  !                             with any orbital
  ! INTEGER NUOTOT            : Total number of orbitals per unit cell
  ! INTEGER NUMH(NUO)         : Number of nonzero elements of each row
  !                             of hamiltonian matrix
  ! INTEGER LISTHPTR(NUO)     : Pointer to each row (-1) of the
  !                             hamiltonian matrix
  ! INTEGER LISTH(MAXNH)      : Nonzero hamiltonian-matrix element
  !                             column indexes for each matrix row
  ! REAL*8  H(MAXNH,4)    : Hamiltonian in sparse format
  ! REAL*8  S(MAXNH)          : Overlap in sparse format
  ! REAL*8  E1, E2            : Energy range for density-matrix states
  !                             (to find local density of states)
  !                             Not used if e1 > e2
  ! INTEGER NHIST             : Number of the subdivisions of the histogram
  ! REAL*8  SIGMA             : Width of the gaussian to expand the eigenvectors
  ! REAL*8  XIJ(3,MAXNH)      : Vectors between orbital centers (sparse)
  !                             (not used if only gamma point)
  ! INTEGER INDXUO(NO)        : Index of equivalent orbital in unit cell
  ! INTEGER NK                : Number of k points
  ! REAL*8  KPOINT(3,NK)      : k point vectors
  ! REAL*8  WK(NK)            : Weights for k points
  ! ****  AUXILIARY  *****************************************************
  ! complex*16  HAUX(2,NUOTOT,2,NUOTOT) : Auxiliary space for the hamiltonian matrix
  ! complex*16  SAUX(2,NUOTOT,2,NUOTOT) : Auxiliary space for the overlap matrix
  ! complex*16  PSI(2,NUOTOT,2,NUOTOT)  : Auxiliary space for the eigenvectors
  ! real*8      eo(nuotot*2)         : Auxiliary space for the eigenvalues
  ! ****  OUTPUT  ********************************************************
  ! REAL*8  DTOT(4,NHIST)      : Total density of states
  ! REAL*8  DPR(4,nuotot,NHIST): Projected density of states
  ! **********************************************************************

  use precision
  use parallel,     only : Node, Nodes
  use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
  use units,        only : pi
  use alloc,        only : re_alloc, de_alloc
#ifdef MPI
  use mpi_siesta
#endif
  use sys, only : die
  use intrinsic_missing, only: modp  

  implicit none

  integer :: nuo, no, maxuo, maxnh, nk, nhist, nuotot
  integer :: numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)

  real(dp) :: H(maxnh,4), S(maxnh), E1, E2, sigma, &
      xij(3,maxnh), kpoint(3,nk), &
      dtot(4,nhist), dpr(4,nuotot,nhist), wk(nk)
  complex(dp) :: Haux(2,nuotot,2,nuotot), Saux(2,nuotot,2,nuotot)
  complex(dp) :: psi(2,nuotot,2,nuotot)
  real(dp) :: eo(nuotot*2)

  ! Internal variables ---------------------------------------------------
  integer :: ik, ispin, iuo, io, iio, juo, j, jo, ihist, iband, ind, ierror
  integer :: iEmin, iEmax
  integer :: nuotot2, BNode

  ! Globalized matrices
  integer :: maxnhg
  integer, pointer :: numhg(:), listhptrg(:), listhg(:)
  real(dp), pointer :: Hg(:,:), Sg(:), xijg(:,:)

  real(dp) :: delta, ener, diff, gauss, norm, wksum
  real(dp) :: limit, inv_sigma2
  real(dp) :: D1, D2

  integer ::  MPIerror
  real(dp), dimension(:,:,:), pointer :: aux_red => null()
  
  external :: cdiag

  ! START -----------------------------------------------------------------

  ! Initialize some variables
  delta = (e2 - e1)/(nhist-1)
  inv_sigma2 = 1._dp / sigma**2
  ! Limit is exp(-20) ~ 2e-9
  limit = sqrt(20._dp) * sigma
  nuotot2 = nuotot * 2

  ! Globalise sparsity pattern
  call MPI_AllReduce(maxnh, maxnhg, 1, MPI_Integer, MPI_Sum, &
      MPI_Comm_World, MPIerror)

  ! Nullify arrays
  nullify(numhg, listhptrg, listhg, Hg, Sg, xijg)
  
  ! Allocate local memory for global list arrays
  call re_alloc( numhg, 1, nuotot, name='numhg', routine= 'pdos2kp' )
  call re_alloc( listhptrg, 1, nuotot, name='listhptrg', routine= 'pdos2kp' )
  call re_alloc( listhg, 1, maxnhg, name='listhg', routine= 'pdos2kp' )
  call re_alloc( Hg, 1, 4, 1, maxnhg, name='Hg', routine= 'pdos2kp' )
  call re_alloc( Sg, 1, maxnhg, name='Sg', routine= 'pdos2kp' )
  call re_alloc( xijg, 1, 3, 1, maxnhg, name='xijg', routine= 'pdos2kp' )

  ! Globalize arrays
  listhptrg(1) = 0
  do io = 1, nuotot
    
    call WhichNodeOrb(io,Nodes,BNode)
    
    if ( Node == BNode ) then
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      numhg(io) = numh(iio)
      do jo = 1, numh(iio)
        listhg(listhptrg(io)+jo) = listh(listhptr(iio)+jo)
        Hg(:,listhptrg(io)+jo) = H(listhptr(iio)+jo,:)
        Sg(listhptrg(io)+jo) = S(listhptr(iio)+jo)
        xijg(:,listhptrg(io)+jo) = xij(:,listhptr(iio)+jo)
      end do
    end if
    
    call MPI_Bcast(numhg(io),1,MPI_Integer, BNode, &
        MPI_Comm_World,MPIerror)
    call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_Integer, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(Hg(1,listhptrg(io)+1),4*numhg(io),MPI_Double_Precision, &
          BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(Sg(listhptrg(io)+1),numhg(io),MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(xijg(1,listhptrg(io)+1),3*numhg(io),MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)

    ! Update list-pointer
    if ( io < nuotot ) listhptrg(io+1) = listhptrg(io) + numhg(io)
    
  end do
  
  ! Find eigenvalues at every k point
  do ik = 1 + Node, nk, Nodes

    call setup_k(kpoint(:,ik))

    ! Diagonalize for each k point. Note duplication of problem size
    call cdiag( Haux, Saux, nuotot2, nuotot2, nuotot2, &
        eo, psi, nuotot2, 1, ierror, -1 ) ! dummy BlockSize

    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')

    elseif ( ierror < 0 ) then
      
      ! Repeat diagonalisation with increased memory to handle clustering
      call setup_k(kpoint(:,ik))

      call cdiag( Haux, Saux, nuotot2, nuotot2, nuotot2, &
          eo, psi, nuotot2, 1, ierror, -1 ) ! dummy BlockSize
    endif

    ! Figure out the minimum and maximum eigenstates that will contribute
    ! This ensures we calculate fewer columns of the psi basis
    ! Note, eo *MUST* be sorted. This is ensured by lapack/scalapack.
    iEmin = nuotot2
    do jo = 1, nuotot2
      if ( (E1 - limit) < EO(jo) .and. EO(jo) < (E2 + limit) ) then
        iEmin = jo
        exit
      end if
    end do

    iEmax = 1
    do jo = nuotot2, 1, -1
      if ( (E1 - limit) < EO(jo) .and. EO(jo) < (E2 + limit) ) then
        iEmax = jo
        exit
      end if
    end do
    ! Ensure that they are in full sections (makes the below easier)
    iEmin = iEmin - mod(iEmin-1,2)
    iEmax = iEmax + mod(iEmax,2)

    ! correct wrong cases, should probably never be found?
    if ( iEmin > iEmax ) then
      iEmin = 1
      iEmax = 0
    end if

    ! Recalculate again the overlap matrix in k-space
    call setup_Sk(kpoint(:,ik))

    ! Convert iEmin/iEmax to local indices (not factor 2)
    iEmin = (iEmin + 1)/2
    iEmax = iEmax / 2
    ! Total number of elements calculated (jo is allowed to be 0)
    jo = (iEmax - iEmin + 1) * 2

    ! Now perform the matrix-multiplications
    ! This is: S | psi >

    call zgemm('N','N',nuotot2, jo, nuotot2, cmplx(1._dp, 0._dp, dp), &
        Saux(1,1,1,1),nuotot2, psi(1,1,1,iEmin),nuotot2, cmplx(0._dp, 0._dp, dp), &
        Haux(1,1,1,iEmin), nuotot2)

!!$OMP parallel default(none) shared(Haux,psi,dtot,dpr,iEmin,iEmax,inv_sigma2,D1,D2) &
!!$OMP& shared(nhist,Node,Nodes,eo,limit,wk,nuotot,nuo,ik,delta,e1) &
!!$OMP& private(jo,iuo,ihist,ener,iband,diff,gauss,j)

    ! Ensure we multiply with the local nodes complex conjugate
    ! This is the final step of < psi | S | psi >
    ! but doing it element wise, rather than a dot-product
    ! Now we will do some magic.
    ! A complex array is equivalent to two reals next to each other.
    ! So to utilize daxpy below, we just ensure the following:
    !   UU = real(UU)
    !   DD = aimag(UU)
    !   X  = real(UD)
    !   Y  = aimag(UD)
    ! This enables us efficient usage of BLAS while obfuscating things
    ! a bit.
    ! Additionally consider that we want the correct sign to use daxpy.
    ! So basically we need to calculate D21 + conjg(D12) which
    ! then has the correct signs.

!!$OMP do schedule(static)
    do jo = iEmin, iEmax
      do io = 1, nuotot
        D1 = real(conjg(psi(1,io,1,jo)) * Haux(1,io,1,jo), dp)
        D2 = real(conjg(psi(2,io,1,jo)) * Haux(2,io,1,jo), dp)
        Haux(2,io,1,jo) = conjg(psi(1,io,1,jo)) * Haux(2,io,1,jo) &
            + psi(2,io,1,jo) * conjg(Haux(1,io,1,jo))
        Haux(1,io,1,jo) = cmplx(D1, D2, dp)
      end do
      do io = 1, nuotot
        D1 = real(conjg(psi(1,io,2,jo)) * Haux(1,io,2,jo), dp)
        D2 = real(conjg(psi(2,io,2,jo)) * Haux(2,io,2,jo), dp)
        Haux(2,io,2,jo) = conjg(psi(1,io,2,jo)) * Haux(2,io,2,jo) &
            + psi(2,io,2,jo) * conjg(Haux(1,io,2,jo))
        Haux(1,io,2,jo) = cmplx(D1, D2, dp)
      end do
    end do
!!$OMP end do

!!$OMP do schedule(static,16)
    do ihist = 1, nhist
      ener = E1 + (ihist - 1) * delta
      do iband = iEmin, iEmax
        ! the energy comes from the global array
        diff = abs(ener - eo(iband*2-1))
          
        if ( diff < limit ) then
          gauss = exp(-diff**2*inv_sigma2) * wk(ik)
          ! See discussion about daxpy + OMP usage in pdosg.F90
          call daxpy(nuotot*4,gauss,Haux(1,1,1,iband),1,dpr(1,1,ihist),1)
        end if
          
        diff = abs(ener - eo(iband*2))
        if ( diff < limit ) then
          gauss = exp(-diff**2*inv_sigma2) * wk(ik)
          ! See discussion about daxpy + OMP usage in pdosg.F90
          call daxpy(nuotot*4,gauss,Haux(1,1,2,iband),1,dpr(1,1,ihist),1)
        end if
          
      end do
    end do
!!$OMP end do nowait

!!$OMP end parallel

  end do ! nk

  ! Free local memory from computation of dpr
  call de_alloc( xijg, name='xijg', routine= 'pdos2kp' )
  call de_alloc( Hg, name='Hg', routine= 'pdos2kp' )
  call de_alloc( Sg, name='Sg', routine= 'pdos2kp' )
  call de_alloc( listhg, name='listhg', routine= 'pdos2kp' )
  call de_alloc( listhptrg, name='listhptrg', routine= 'pdos2kp' )
  call de_alloc( numhg, name='numhg', routine= 'pdos2kp' )

  ! Allocate workspace array for global reduction
  call re_alloc( aux_red, 1, 4, 1, nuotot, 1, nhist, &
      name='aux_red', routine='pdos2kp')
  
  ! Global reduction of dpr matrix
  call MPI_Reduce(dpr(1,1,1),aux_red(1,1,1),4*nuotot*nhist, &
      MPI_double_precision,MPI_sum,0,MPI_Comm_World,ierror)
  dpr(:,:,:) = aux_red(:,:,:)

  call de_alloc(aux_red, name='aux_red', routine='pdos2kp')
  
  wksum = 0.0d0
  do ik = 1, nk
    wksum = wksum + wk(ik)
  end do
  norm = 1._dp / (sigma * sqrt(pi) * wksum)
  call dscal(4*nuotot*nhist,norm,dpr(1,1,1),1)

  do ihist = 1, nhist
    dtot(1,ihist) = sum(dpr(1,:,ihist))
    dtot(2,ihist) = sum(dpr(2,:,ihist))
    dtot(3,ihist) = sum(dpr(3,:,ihist))
    dtot(4,ihist) = sum(dpr(4,:,ihist))
  end do

contains

  subroutine setup_k(k)
    real(dp), intent(in) :: k(3)
    real(dp) :: kxij
    complex(dp) :: kphs

    Saux(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
    Haux(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do io = 1,nuotot
      do ind = listhptrg(io) + 1, listhptrg(io) + numhg(io)
        jo = modp(listhg(ind), nuotot)
        kxij = k(1) * xijg(1,ind) + k(2) * xijg(2,ind) + k(3) * xijg(3,ind)

        ! Calculate the complex phase
        kphs = exp(cmplx(0._dp, -kxij, dp))
        
        Saux(1,jo,1,io) = Saux(1,jo,1,io) + Sg(ind) * kphs
        Saux(2,jo,2,io) = Saux(2,jo,2,io) + Sg(ind) * kphs
        Haux(1,jo,1,io) = Haux(1,jo,1,io) + cmplx(Hg(1,ind), 0._dp, dp) * kphs
        Haux(2,jo,2,io) = Haux(2,jo,2,io) + cmplx(Hg(2,ind), 0._dp, dp) * kphs
        Haux(1,jo,2,io) = Haux(1,jo,2,io) + cmplx(Hg(3,ind), -Hg(4,ind), dp) * kphs
        Haux(2,jo,1,io) = Haux(2,jo,1,io) + cmplx(Hg(3,ind), Hg(4,ind), dp) * kphs

      end do
    end do

  end subroutine setup_k


  subroutine setup_Sk(k)
    real(dp), intent(in) :: k(3)
    real(dp) :: kxij
    complex(dp) :: kphs

    Saux(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do io = 1,nuotot
      do ind = listhptrg(io) + 1, listhptrg(io) + numhg(io)
        jo = modp(listhg(ind), nuotot)
        kxij = k(1) * xijg(1,ind) + k(2) * xijg(2,ind) + k(3) * xijg(3,ind)

        ! Calculate the complex phase
        kphs = exp(cmplx(0._dp, -kxij, dp))
        
        Saux(1,jo,1,io) = Saux(1,jo,1,io) + Sg(ind) * kphs
        Saux(2,jo,2,io) = Saux(2,jo,2,io) + Sg(ind) * kphs

      end do
    end do
    
  end subroutine setup_Sk
  
end subroutine pdos2kp
#endif
