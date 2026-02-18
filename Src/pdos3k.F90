! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine pdos3k( nuo, no, maxuo, maxnh, &
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
  ! REAL*8  H(MAXNH,8)    : Hamiltonian in sparse format
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
  ! complex*16  HAUX(2,NUOTOT,2,NUO) : Auxiliary space for the hamiltonian matrix
  ! complex*16  SAUX(2,NUOTOT,2,NUO) : Auxiliary space for the overlap matrix
  ! complex*16  PSI(2,NUOTOT,2,NUO)  : Auxiliary space for the eigenvectors
  ! real*8      eo(nuotot*2)         : Auxiliary space for the eigenvalues
  ! ****  OUTPUT  ********************************************************
  ! REAL*8  DTOT(4,NHIST)      : Total density of states
  ! REAL*8  DPR(4,nuotot,NHIST): Projected density of states
  ! **********************************************************************

  use precision
  use parallel,     only : Node, Nodes, BlockSize
  use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
  use units,        only : pi
  use alloc,        only : re_alloc, de_alloc
#ifdef MPI
  use mpi_siesta
#endif
  use sys, only : die
  use m_diag, only: diag_get_1d_context, diag_descinit
  use m_diag_option, only: Serial

  implicit none

  integer :: nuo, no, maxuo, maxnh, nk, nhist, nuotot
  integer :: numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)

  real(dp) :: H(maxnh,8), S(maxnh), E1, E2, sigma, &
      xij(3,maxnh), kpoint(3,nk), &
      dtot(4,nhist), dpr(4,nuotot,nhist), wk(nk)
  complex(dp) :: Haux(2,nuotot,2,nuo), Saux(2,nuotot,2,nuo)
  complex(dp) :: psi(2,nuotot,2,nuo)
  real(dp) :: eo(nuotot*2)

  ! Internal variables ---------------------------------------------------
  integer :: ik, ispin, iuo, io, juo, j, jo, ihist, iband, ind, ierror
  integer :: iEmin, iEmax
  integer :: nuo2, nuotot2, BlockSize2

  real(dp) :: delta, ener, diff, gauss, norm, wksum
  real(dp) :: limit, inv_sigma2
  real(dp) :: D1, D2

#ifdef MPI
  ! All of our matrices are described in the same manner
  ! So we only need one descriptor
  integer :: ctxt, desc(9)
  real(dp), dimension(:,:,:), pointer :: aux_red => null()
#endif

  external :: cdiag

  ! START -----------------------------------------------------------------

  ! Initialize some variables
  delta = (e2 - e1)/(nhist-1)
  inv_sigma2 = 1._dp / sigma**2
  ! Limit is exp(-20) ~ 2e-9
  limit = sqrt(20._dp) * sigma
  nuo2 = nuo * 2
  nuotot2 = nuotot * 2
  BlockSize2 = BlockSize * 2

  ! Find eigenvalues at every k point
  do ik = 1,nk

    call setup_k(kpoint(:,ik))

    ! Diagonalize for each k point. Note duplication of problem size
    call cdiag( Haux, Saux, nuotot2, nuo2, nuotot2, &
        eo, psi, nuotot2, 1, ierror, BlockSize2 )

    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')

    elseif ( ierror < 0 ) then
      
      ! Repeat diagonalisation with increased memory to handle clustering
      call setup_k(kpoint(:,ik))

      call cdiag( Haux, Saux, nuotot2, nuo2, nuotot2, &
          eo, psi, nuotot2, 1, ierror, BlockSize2 )
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

    if ( Serial ) then
      call zgemm('N','N',nuotot2, jo, nuotot2, cmplx(1._dp, 0._dp, dp), &
          Saux(1,1,1,1),nuotot2, psi(1,1,1,iEmin),nuotot2, cmplx(0._dp, 0._dp, dp), &
          Haux(1,1,1,iEmin), nuotot2)
    end if

#ifdef MPI

    if ( .not. Serial ) then
      ! We need to define the contexts and matrix descriptors
      ctxt = diag_get_1d_context()
      call diag_descinit(nuotot2, nuotot2, BlockSize2, desc, ctxt)

      call pzgemm('N', 'N', nuotot2, jo, nuotot2, cmplx(1._dp, 0._dp, dp), &
          Saux(1,1,1,1), 1, 1, desc, psi(1,1,1,1), 1, iEmin*2-1, desc, &
          cmplx(0._dp, 0._dp, dp), Haux(1,1,1,1), 1, iEmin*2-1, desc)
    end if

    ! reset
    iuo = iEmin
    juo = iEmax
    iEmin = 1
    iEmax = 0
    do jo = 1, nuo
      call LocalToGlobalOrb(jo, Node, Nodes, j)
      if ( iuo <= j ) then
        ! lowest point where we have this orbital
        iEmin = jo
        exit
      end if
    end do
    do jo = nuo, 1, -1
      call LocalToGlobalOrb(jo, Node, Nodes, j)
      if ( j <= juo ) then
        iEmax = jo
        exit
      end if
    end do

    ! Now iEmin, iEmax are local indices

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
        call LocalToGlobalOrb(iband, Node, Nodes, j)
        diff = abs(ener - eo(j*2-1))

        ! TODO, this *could* be merged into a single zgemm call with gauss(2), but...
        ! In fact, all of iEmin ... iEmax could be merged to do everything *once*
        ! But a bit more complicated.
        
        if ( diff < limit ) then
          gauss = exp(-diff**2*inv_sigma2) * wk(ik)
          ! See discussion about daxpy + OMP usage in pdosg.F90
          call daxpy(nuotot*4,gauss,Haux(1,1,1,iband),1,dpr(1,1,ihist),1)
        end if
          
        diff = abs(ener - eo(j*2))
        if ( diff < limit ) then
          gauss = exp(-diff**2*inv_sigma2) * wk(ik)
          ! See discussion about daxpy + OMP usage in pdosg.F90
          call daxpy(nuotot*4,gauss,Haux(1,1,2,iband),1,dpr(1,1,ihist),1)
        end if
          
      end do
    end do
!!$OMP end do nowait

!!$OMP end parallel

#else

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

#endif

  end do ! nk

#ifdef MPI
  ! Allocate workspace array for global reduction
  call re_alloc( aux_red, 1, 4, 1, nuotot, 1, nhist, &
      name='aux_red', routine='pdos3k')
  
  ! Global reduction of dpr matrix
  call MPI_Reduce(dpr(1,1,1),aux_red(1,1,1),4*nuotot*nhist, &
      MPI_double_precision,MPI_sum,0,MPI_Comm_World,ierror)
  dpr(:,:,:) = aux_red(:,:,:)

  call de_alloc(aux_red, name='aux_red', routine='pdos3k')

#endif

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

    do iuo = 1,nuo
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        kxij = k(1) * xij(1,ind) + k(2) * xij(2,ind) + k(3) * xij(3,ind)
        kphs = exp(cmplx(0.0_dp, -kxij, dp))

        Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind) * kphs
        Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind) * kphs
        Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + cmplx(H(ind,1), H(ind,5),dp) * kphs
        Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + cmplx(H(ind,2), H(ind,6),dp) * kphs
        Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo) + cmplx(H(ind,3), - H(ind,4),dp) * kphs
        Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo) + cmplx(H(ind,7), + H(ind,8),dp) * kphs

      end do
    end do

  end subroutine setup_k


  subroutine setup_Sk(k)
    real(dp), intent(in) :: k(3)
    real(dp) :: kxij
    complex(dp) :: kphs

    Saux(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do iuo = 1,nuo
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        kxij = k(1) * xij(1,ind) + k(2) * xij(2,ind) + k(3) * xij(3,ind)
        kphs = exp(cmplx(0.0_dp, -kxij, dp))

        Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind) * kphs
        Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind) * kphs

      end do
    end do

  end subroutine setup_Sk
  
end subroutine pdos3k
