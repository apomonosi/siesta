! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine pdosg( nspin, nuo, no, maxnh, &
    nuotot, numh, listhptr, listh, H, S, &
    E1, E2, nhist, sigma, indxuo, &
    haux, saux, psi, eo, dtot, dpr)
  ! **********************************************************************
  ! Find the density of states projected onto the atomic orbitals
  !     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
  ! where n run over all the bands between two given energies
  ! Written by J. Junquera and E. Artacho. Nov' 99
  ! Gamma point version adapted from PDOSK by Julian Gale. Feb' 03
  ! Huge performance increase by N. Papior, 2022
  ! ****  INPUT  *********************************************************
  ! integer nspin             : Number of spin components (1 or 2)
  ! integer nuo               : Number of atomic orbitals in the unit cell
  ! integer no                : Number of atomic orbitals in the supercell
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
  ! integer indxuo(no)        : Index of equivalent orbital in unit cell
  ! ****  AUXILIARY  *****************************************************
  ! real*8  haux(nuo,nuo)     : Auxiliary space for the hamiltonian matrix
  ! real*8  saux(nuo,nuo)     : Auxiliary space for the overlap matrix
  ! real*8  psi(nuo,nuo)      : Auxiliary space for the eigenvectors
  ! real*8  eo(nuotot)        : Auxiliary space for the eigenvalues
  ! ****  OUTPUT  ********************************************************
  ! real*8  dtot(nhist,2)   : Total density of states
  ! real*8  dpr(nuo,nhist,2): Projected density of states
  ! **********************************************************************

  use precision
  use parallel,     only : Node, Nodes, BlockSize
  use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
  use units,        only : pi
  use alloc,        only : re_alloc, de_alloc
#ifdef MPI
  use mpi_siesta
#endif
  use sys,          only : die
  use m_diag, only: diag_get_1d_context, diag_descinit
  use m_diag_option, only: Serial

  implicit none

  integer :: nspin, nuo, no, maxnh, nhist, nuotot

  integer :: numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)

  real(dp) :: H(maxnh,nspin), S(maxnh), E1, E2, sigma, &
      haux(nuotot,nuo), saux(nuotot,nuo), psi(nuotot,nuo), &
      dtot(nhist,nspin), dpr(nuotot,nhist,nspin)
  real(dp) :: eo(nuotot)

  ! Internal variables ---------------------------------------------------
  integer :: ispin, iuo, juo, j, jo, ihist, iband, ind, ierror
  integer :: iEmin, iEmax

  real(dp) :: delta, ener, diff, gauss, norm
  real(dp) :: limit, inv_sigma2

#ifdef MPI
  ! All of our matrices are described in the same manner
  ! So we only need one descriptor
  integer :: ctxt, desc(9)
  
  real(dp), dimension(:,:,:), pointer :: aux_red => null()
#endif

  external :: rdiag

  ! Initialize some variables
  delta = (E2 - E1)/(nhist-1)

  inv_sigma2 = 1._dp / sigma**2
  ! Limit is exp(-20) ~ 2e-9
  limit = sqrt(20._dp) * sigma

  ! Solve eigenvalue problem for each k-point
  do ispin = 1, nspin

    call setup_Gamma()

    ! Diagonalize at the Gamma point
    call rdiag( haux, saux, nuotot, nuo, nuotot, eo, &
        psi, nuotot, 1, ierror, BlockSize)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    elseif ( ierror < 0 ) then
      call setup_Gamma()
      call rdiag( haux, saux, nuotot, nuo, nuotot, eo, &
          psi, nuotot, 1, ierror, BlockSize)
    endif

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

    ! Setup overlap matrix
    call setup_S()

    ! Total number of elements calculated (jo is allowed to be 0)
    jo = iEmax - iEmin + 1

    ! Now perform the matrix-multiplications
    ! This is: S | psi >

    if ( Serial ) then
      call dgemm('N','N',nuotot, jo, nuotot, 1._dp, &
          Saux(1,1), nuotot, psi(1,iEmin),nuotot, 0._dp, &
          Haux(1,iEmin), nuotot)
    end if

#ifdef MPI

    if ( .not. Serial ) then
      ! We need to define the contexts and matrix descriptors
      ctxt = diag_get_1d_context()
      write(6,*) 'pdosg.F90: ctxt=', ctxt
      call diag_descinit(nuotot, nuotot, BlockSize, desc, ctxt)

      call pdgemm('N', 'N', nuotot, jo, nuotot, 1._dp, &
          Saux(1,1), 1, 1, desc, psi(1,1), 1, iEmin, desc, &
          0._dp, Haux(1,1), 1, iEmin, desc)
    end if

    ! Convert iEmin/iEmax to local indices
    iuo = iEmin
    juo = iEmax
    ! reset
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

!!$OMP parallel default(none) shared(Haux,psi,dtot,dpr,iEmin,iEmax,inv_sigma2) &
!!$OMP& shared(nhist,Node,Nodes,eo,limit,nuotot,nuo,ispin,delta,e1) &
!!$OMP& private(jo,iuo,ihist,ener,iband,diff,gauss,j)

    ! Ensure we multiply with the local nodes complex conjugate
    ! This is the final step of < psi | S | psi >
    ! but doing it element wise, rather than a dot-product
!!$OMP do schedule(static)
    do jo = iEmin, iEmax
      do iuo = 1, nuotot
        Haux(iuo,jo) = psi(iuo,jo) * Haux(iuo,jo)
      end do
    end do
!!$OMP end do

!!$OMP do schedule(static,16)
    do ihist = 1, nhist
      ener = E1 + (ihist - 1) * delta
      do iband = iEmin, iEmax
        ! the energy comes from the global array
        call LocalToGlobalOrb(iband, Node, Nodes, j)
        diff = abs(ener - eo(j))
          
        if ( diff < limit ) then
          gauss = exp(-diff**2*inv_sigma2)
          dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
          ! When OMP is inserted, this should still work fine since
          ! BLAS does not fork threads that are already forked.
          ! However, this should be carefully checked, and possibly be able
          ! to switch off the behaviour. I do know that ELPA does the same,
          ! so perhaps it is fine in generic situations after all.
          call daxpy(nuotot,gauss,Haux(1,iband),1,dpr(1,ihist,ispin),1)
        end if
          
      end do
    end do
!!$OMP end do nowait

!!$OMP end parallel

#else

!!$OMP parallel default(none) shared(Haux,psi,dtot,dpr,iEmin,iEmax,inv_sigma2) &
!!$OMP& shared(nhist,Node,Nodes,eo,limit,nuotot,nuo,ispin,delta,e1) &
!!$OMP& private(jo,iuo,ihist,ener,iband,diff,gauss,j)

    ! Ensure we multiply with the local nodes complex conjugate
    ! This is the final step of < psi | S | psi >
    ! but doing it element wise, rather than a dot-product
!!$OMP do schedule(static)
    do jo = iEmin, iEmax
      do iuo = 1, nuotot
        Haux(iuo,jo) = psi(iuo,jo) * Haux(iuo,jo)
      end do
    end do
!!$OMP end do

!!$OMP do schedule(static,16)
    do ihist = 1, nhist
      ener = E1 + (ihist - 1) * delta
      do iband = iEmin, iEmax
        diff = abs(ener - eo(iband))
          
        if ( diff < limit ) then
          gauss = exp(-diff**2*inv_sigma2)
          dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
          ! When OMP is inserted, this should still work fine since
          ! BLAS does not fork threads that are already forked.
          ! However, this should be carefully checked, and possibly be able
          ! to switch off the behaviour. I do know that ELPA does the same,
          ! so perhaps it is fine in generic situations after all.
          call daxpy(nuotot,gauss,Haux(1,iband),1,dpr(1,ihist,ispin),1)
        end if
          
      end do
    end do
!!$OMP end do nowait

!!$OMP end parallel
#endif

  enddo

#ifdef MPI
  ! Allocate workspace array for global reduction
  call re_alloc( aux_red, 1, nuotot, 1, nhist, 1, nspin, &
      name='aux_red_dpr', routine='pdosg' )

  ! Global reduction of dpr matrix
  call MPI_Reduce(dpr(1,1,1),aux_red(1,1,1),nuotot*nhist*nspin, &
      MPI_double_precision,MPI_sum,0,MPI_Comm_World,ierror)
  dpr(:,:,:) = aux_red(:,:,:)

  call de_alloc(aux_red, name='aux_red_dpr', routine='pdosg' )

  call re_alloc( aux_red, 1, nhist, 1, nspin, 1, 1, &
      name='aux_red_dtot', routine='pdosg' )

  ! Global reduction of dtot matrix
  call MPI_Reduce(dtot(1,1),aux_red(1,1,1),nhist*nspin, &
      MPI_double_precision,MPI_sum,0,MPI_Comm_World,ierror)
  dtot(:,:) = aux_red(:,:,1)

  call de_alloc(aux_red, name='aux_red_dtot', routine='pdosg' )
#endif

  norm = 1._dp /(sigma * sqrt(pi))

  do ispin = 1, nspin
    do ihist = 1, nhist
      dtot(ihist,ispin) = dtot(ihist,ispin) * norm
    enddo
  enddo

  do ispin = 1, nspin
    do ihist = 1, nhist
      do iuo = 1, nuotot
        dpr(iuo,ihist,ispin) = dpr(iuo,ihist,ispin) * norm
      enddo
    enddo
  enddo

contains

  subroutine setup_Gamma()

    ! Initialize auxiliary variables
    do iuo = 1,nuo
      do juo = 1,nuotot
        Saux(juo,iuo) = 0.0d0
        Haux(juo,iuo) = 0.0d0
      enddo
    enddo

    do iuo = 1, nuo
      do j = 1, numh(iuo)
        ind = listhptr(iuo) + j
        jo  = listh(ind)
        juo = indxuo(jo)
        ! Calculate the Hamiltonian and the overlap
        Saux(juo,iuo) = Saux(juo,iuo) + S(ind)
        Haux(juo,iuo) = Haux(juo,iuo) + H(ind,ispin)
      enddo
    enddo

  end subroutine setup_Gamma

  subroutine setup_S()

    ! Initialize auxiliary variables 
    do iuo = 1,nuo
      do juo = 1,nuotot
        Saux(juo,iuo) = 0.0d0
      end do
    end do

    do iuo = 1, nuo
      do j = 1, numh(iuo)
        ind = listhptr(iuo) + j
        jo  = listh(ind)
        juo = indxuo(jo)
        Saux(juo,iuo) = Saux(juo,iuo) + S(ind)
      enddo
    enddo

  end subroutine setup_S

end subroutine pdosg
