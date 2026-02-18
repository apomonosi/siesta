! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
#ifdef MPI
subroutine pdoskp(nspin, nuo, no, maxnh, &
    nuotot, numh, listhptr, listh, H, S, &
    E1, E2, nhist, sigma, &
    xij, indxuo, nk, kpoint, wk, &
    Haux, Saux, psi, eo, dtot, dpr)

  ! **********************************************************************
  ! Find the density of states projected onto the atomic orbitals
  !     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
  ! where n run over all the bands between two given energies
  ! Written by J. Junquera and E. Artacho. Nov' 99
  ! Modified version for parallel execution over K points by J.D. Gale
  ! March 2005
  ! ****  INPUT  *********************************************************
  ! integer nspin             : Number of spin components (1 or 2)
  ! integer nuo               : Number of atomic orbitals in the unit cell
  ! integer NO                : Number of atomic orbitals in the supercell
  ! integer maxnh             : Maximum number of orbitals interacting
  !                             with any orbital
  ! integer nuotot            : Total number of orbitals per unit cell
  ! integer numh(nuo)         : Number of nonzero elements of each row
  !                             of hamiltonian matrix
  ! integer listhptr(nuo)     : Pointer to each row (-1) of the
  !                             hamiltonian matrix
  ! integer listh(maxnh)      : Nonzero hamiltonian-matrix element
  !                             column indexes for each matrix row
  ! real*8  H(maxnh,nspin)    : Hamiltonian in sparse format
  ! real*8  S(maxnh)          : Overlap in sparse format
  ! real*8  E1, E2            : Energy range for density-matrix states
  !                             (to find local density of states)
  !                             Not used if e1 > e2
  ! integer nhist             : Number of the subdivisions of the histogram
  ! real*8  sigma             : Width of the gaussian to expand the eigenvectors
  ! real*8  xij(3,maxnh)      : Vectors between orbital centers (sparse)
  !                             (not used if only gamma point)
  ! integer indxuo(no)        : Index of equivalent orbital in unit cell
  ! integer nk                : Number of k points
  ! real*8  kpoint(3,nk)      : k point vectors
  ! real*8  wk(nk)            : Weights for k points
  ! ****  AUXILIARY  *****************************************************
  ! real*8  Haux(2,nuotot,nuo) : Auxiliary space for the hamiltonian matrix
  ! real*8  Saux(2,nuotot,nuo) : Auxiliary space for the overlap matrix
  ! real*8  psi(2,nuotot,nuo)  : Auxiliary space for the eigenvectors
  ! real*8  eo(nuotot)         : Auxiliary space for the eigenvalues
  ! ****  OUTPUT  ********************************************************
  ! real*8  dtot(nhist,2)   : Total density of states
  ! real*8  dpr(nuotot,nhist,2): Proyected density of states
  ! **********************************************************************

  use precision
  use parallel,     only : Node, Nodes
  use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
  use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
  use units,        only : pi
  use alloc,        only : re_alloc, de_alloc
#ifdef MPI
  use mpi_siesta
#endif
  use sys,          only : die

  implicit none

  integer :: nspin, nuo, no, maxspn, maxnh, NK, nhist, nuotot

  integer :: numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)

  real(dp) :: H(maxnh,nspin), S(maxnh), E1, E2, sigma, &
      xij(3,maxnh), kpoint(3,nk), &
      Haux(2,nuotot,nuotot), Saux(2,nuotot,nuotot), &
      psi(2,nuotot,nuotot), &
      dtot(nhist,nspin), dpr(nuotot,nhist,nspin), wk(nk)
  real(dp) :: eo(nuotot)

  ! Internal variables ---------------------------------------------------
  integer :: ispin, ik, iio, io, iuo, juo, j, jo, ihist, iband, ind, ierror
  integer :: maxnhg, nuog, BNode
  integer :: iEmin, iEmax

  integer, dimension(:), pointer :: numhg, listhptrg, listhg

  real(dp) :: delta, ener, diff
  real(dp) :: gauss, norm, wksum
  real(dp) :: limit, inv_sigma2

  real(dp), dimension(:), pointer :: Sg
  real(dp), dimension(:,:), pointer :: Hg, xijg

#ifdef MPI
  integer ::  MPIerror
  real(dp), dimension(:,:,:), pointer :: aux_red => null()
#endif

  external :: cdiag

#ifdef DEBUG
  call write_debug( '    PRE pdoskp' )
#endif

  ! Initialize some variables
  delta = (E2 - E1)/(nhist-1)
  inv_sigma2 = 1._dp / sigma**2
  ! Limit is exp(-20) ~ 2e-9
  limit = sqrt(20._dp) * sigma

  ! Globalise list arrays - assumes listh and listd are the same

  ! Allocate local memory for global list arrays
  nullify( numhg )
  call re_alloc( numhg, 1, nuotot, name='numhg', routine='pdoskp')
  nullify( listhptrg )
  call re_alloc( listhptrg, 1, nuotot, name='listhptrg', routine='pdoskp')

  ! Globalise numh
  do io = 1,nuotot
    call WhichNodeOrb(io,Nodes,BNode)
    if (Node.eq.BNode) then
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      numhg(io) = numh(iio)
    endif
#ifdef MPI
    call MPI_Bcast(numhg(io),1,MPI_integer,BNode, MPI_Comm_World,MPIerror)
#endif
  enddo

  ! Build global listhptr
  listhptrg(1) = 0
  do io = 2,nuotot
    listhptrg(io) = listhptrg(io-1) + numhg(io-1)
  end do

  ! Globalse listh
  maxnhg = listhptrg(nuotot) + numhg(nuotot)
  nullify( listhg )
  call re_alloc( listhg, 1, maxnhg, name='listhg', routine='pdoskp')
  do io = 1,nuotot
    call WhichNodeOrb(io,Nodes,BNode)
    if (Node.eq.BNode) then
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      do jo = 1,numhg(io)
        listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
            listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
      enddo
    endif
#ifdef MPI
    call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer, &
        BNode,MPI_Comm_World,MPIerror)
#endif
  enddo

  ! Create new distribution of H and S
  nuog = nuotot

  nullify( Sg )
  call re_alloc( Sg, 1, maxnhg, name='Sg', routine='pdoskp')
  nullify( Hg )
  call re_alloc( Hg, 1, maxnhg, 1, nspin, name='Hg', routine='pdoskp')
  nullify( xijg )
  call re_alloc( xijg, 1, 3, 1, maxnhg, name='xijg', routine='pdoskp')

  do io = 1,nuotot
    call WhichNodeOrb(io,Nodes,BNode)
    if (Node.eq.BNode) then
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      do ispin = 1,nspin
        do jo = 1,numh(iio)
          Hg(listhptrg(io)+jo,ispin) = H(listhptr(iio)+jo,ispin)
        enddo
      enddo
      do jo = 1,numh(iio)
        Sg(listhptrg(io)+jo) = S(listhptr(iio)+jo)
      end do
      do jo = 1,numh(iio)
        xijg(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
      enddo
    endif
#ifdef MPI
    do ispin = 1,nspin
      call MPI_Bcast(Hg(listhptrg(io)+1,ispin),numhg(io), &
          MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
    end do
    call MPI_Bcast(Sg(listhptrg(io)+1),numhg(io), &
        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(xijg(1,listhptrg(io)+1),3*numhg(io), &
        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
#endif
  end do

  ! Solve eigenvalue problem for each k-point
  do ispin = 1,nspin

    do ik = 1+Node,nk,Nodes

      call setup_k(kpoint(:,ik))

      ! Diagonalize for each k point
      call cdiag( Haux, Saux, nuotot, nuog, nuotot, &
          eo, psi, nuotot, 1, ierror, -1) !dummy blocksize

      ! Check error flag and take appropriate action
      if ( ierror > 0 ) then
        call die('Terminating due to failed diagonalisation')
      elseif ( ierror < 0 ) then

        ! Repeat diagonalisation with increased memory to handle clustering
        call setup_k(kpoint(:,ik))
        call cdiag( Haux, Saux, nuotot, nuog, nuotot, &
            eo, psi, nuotot, 1, ierror, -1) !dummy blocksize
      end if

      ! Figure out the minimum and maximum eigenstates that will contribute
      ! This ensures we calculate fewer columns of the psi basis
      ! Note, eo *MUST* be sorted. This is ensured by lapack/scalapack.
      iEmin = nuotot
      do jo = 1, nuotot
        if ( (E1 - limit) < EO(jo) .and. EO(jo) < (E2 + limit) ) then
          iEmin = jo
          exit
        end if
      end do

      iEmax = 1
      do jo = nuotot, 1, -1
        if ( (E1 - limit) < EO(jo) .and. EO(jo) < (E2 + limit) ) then
          iEmax = jo
          exit
        end if
      end do

      ! correct wrong cases, should probably never be found?
      if ( iEmin > iEmax ) then
        iEmin = 1
        iEmax = 0
      end if

      ! Recalculate again the overlap matrix in k-space
      call setup_Sk(kpoint(:,ik))

      ! Total number of elements calculated (jo is allowed to be 0)
      jo = iEmax - iEmin + 1

      ! Now perform the matrix-multiplications
      ! This is: S | psi >
      call zgemm('N','N',nuotot, jo, nuotot, cmplx(1._dp, 0._dp, dp), &
          Saux(1,1,1),nuotot, psi(1,1,iEmin),nuotot,cmplx(0._dp, 0._dp, dp), &
          Haux(1,1,iEmin), nuotot)


      !!$OMP parallel default(none) shared(Haux,psi,dtot,dpr,iEmin,iEmax,inv_sigma2) &
!!$OMP& shared(nhist,Node,Nodes,eo,limit,wk,nuotot,nuo,ispin,ik,delta,e1) &
!!$OMP& private(jo,iuo,ihist,ener,iband,diff,gauss,j)

      ! Ensure we multiply with the local nodes complex conjugate
      ! This is the final step of < psi | S | psi >
      ! but doing it element wise, rather than a dot-product
!!$OMP do schedule(static)
      do jo = iEmin, iEmax
        do iuo = 1, nuotot
          Haux(1,iuo,jo) = psi(1,iuo,jo) * Haux(1,iuo,jo) + psi(2,iuo,jo) * Haux(2,iuo,jo)
        end do
      end do
!!$OMP end do

!!$OMP do schedule(static,16)
      do ihist = 1, nhist
        ener = E1 + (ihist - 1) * delta
        do iband = iEmin, iEmax
          diff = abs(ener - eo(iband))
          
          if ( diff < limit ) then
            gauss = exp(-diff**2*inv_sigma2) * wk(ik)
            dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
            ! See discussion about daxpy + OMP usage in pdosg.F90
            call daxpy(nuotot,gauss,Haux(1,1,iband),2,dpr(1,ihist,ispin),1)
          end if
          
        end do
      end do
!!$OMP end do nowait

!!$OMP end parallel

    end do

  end do

  ! Free local memory from computation of dpr
  call de_alloc( xijg, name='xijg' )
  call de_alloc( Hg, name='Hg' )
  call de_alloc( Sg, name='Sg' )
  call de_alloc( listhg, name='listhg' )
  call de_alloc( listhptrg, name='listhptrg' )
  call de_alloc( numhg, name='numhg' )

#ifdef MPI
  ! Allocate workspace array for global reduction
  call re_alloc( aux_red, 1, nuotot, 1, nhist, 1, nspin, &
      name='aux_red_dpr', routine='pdosk' )

  ! Global reduction of dpr matrix
  call MPI_Reduce(dpr(1,1,1),aux_red(1,1,1),nuotot*nhist*nspin, &
      MPI_double_precision,MPI_sum,0,MPI_Comm_World,ierror)
  dpr(:,:,:) = aux_red(:,:,:)

  call de_alloc(aux_red, name='aux_red_dpr', routine='pdosk' )

  call re_alloc( aux_red, 1, nhist, 1, nspin, 1, 1, &
      name='aux_red_dtot', routine='pdosk' )

  ! Global reduction of dtot matrix
  call MPI_Reduce(dtot(1,1),aux_red(1,1,1),nhist*nspin, &
      MPI_double_precision,MPI_sum,0,MPI_Comm_World,ierror)
  dtot(:,:) = aux_red(:,:,1)

  call de_alloc(aux_red, name='aux_red_dtot', routine='pdosk' )

#endif

  wksum = 0.0d0
  do ik = 1,nk
    wksum = wksum + wk(ik)
  enddo

  norm = 1._dp / (sigma * sqrt(pi) * wksum)

  do ispin = 1,nspin
    do ihist = 1,nhist
      dtot(ihist,ispin) = dtot(ihist,ispin) * norm
      do iuo = 1,nuotot
        dpr(iuo,ihist,ispin) = dpr(iuo,ihist,ispin) * norm
      enddo
    enddo
  enddo

#ifdef DEBUG
  call write_debug( '    POS pdoskp' )
#endif

contains

  subroutine setup_k(k)
    real(dp), intent(in) :: k(3)
    real(dp) :: kxij, Ckxij, Skxij

    ! Initialize auxiliary variables
    do iuo = 1,nuog
      do juo = 1,nuotot
        Saux(1,juo,iuo) = 0.0d0
        Saux(2,juo,iuo) = 0.0d0
        Haux(1,juo,iuo) = 0.0d0
        Haux(2,juo,iuo) = 0.0d0
      enddo
    enddo

    do iuo = 1, nuog
      do j = 1, numhg(iuo)
        ind = listhptrg(iuo) + j
        jo  = listhg(ind)
        juo = indxuo(jo)
        ! Calculate the phases k*r_ij
        kxij = k(1) * xijg(1,ind) + k(2) * xijg(2,ind) + k(3) * xijg(3,ind)
        Ckxij = cos(kxij)
        Skxij = sin(kxij)
        ! Calculate the Hamiltonian and the overlap in k space
        ! H(k) = Sum(R) exp(i*k*R) * H(R)
        Saux(1,juo,iuo) = Saux(1,juo,iuo) + Sg(ind) * Ckxij
        Saux(2,juo,iuo) = Saux(2,juo,iuo) - Sg(ind) * Skxij
        Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hg(ind,ispin) * Ckxij
        Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hg(ind,ispin) * Skxij
      end do
    end do

  end subroutine setup_k

  subroutine setup_Sk(k)
    real(dp), intent(in) :: k(3)
    real(dp) :: kxij, Ckxij, Skxij

    ! Initialize auxiliary variables 
    do iuo = 1,nuog
      do juo = 1,nuotot
        Saux(1,juo,iuo) = 0.0d0
        Saux(2,juo,iuo) = 0.0d0
      end do
    end do

    do iuo = 1, nuog
      do j = 1, numhg(iuo)
        ind = listhptrg(iuo) + j
        jo  = listhg(ind)
        juo = indxuo(jo)
        ! Calculate the phases k*r_ij
        kxij = k(1) * xijg(1,ind) + k(2) * xijg(2,ind) + k(3) * xijg(3,ind)
        Ckxij = cos(kxij)
        Skxij = sin(kxij)
        ! Calculate the Hamiltonian and the overlap in k space
        ! H(k) = Sum(R) exp(i*k*R) * H(R)
        Saux(1,juo,iuo) = Saux(1,juo,iuo) + Sg(ind) * Ckxij
        Saux(2,juo,iuo) = Saux(2,juo,iuo) - Sg(ind) * Skxij
      enddo
    enddo

  end subroutine setup_Sk

end subroutine pdoskp
#endif
