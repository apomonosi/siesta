! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module io_hsx_m

  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData1D
  use class_dSpData2D

  implicit none

  public :: write_hsx
  public :: write_hs_formatted

  private
  
contains

  subroutine write_hsx(H, S, Ef, qtot, temp, prec)
! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  Ef                  : Fermi level
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! integer prec                : The precision that stores the data (sp|dp)

    use precision, only: dp
    use io_sparse_m, only: io_write, io_write_r
    use parallel, only : Node, Nodes
    use atm_types, only : nspecies
    use atomlist, only : iphorb, iaorb, lasto
    use siesta_geom, only: na_u, ucell, xa, isa, nsc, isc_off
    use atmfuncs, only : nofis, labelfis, zvalfis
    use atmfuncs, only : cnfigfio, lofio, zetafio
    use fdf
    use files, only : slabel
    use sys, only : die
#ifdef MPI
    use mpi_siesta
#endif

    type(dSpData2D), intent(inout) :: H
    type(dSpData1D), intent(inout) :: S
    real(dp) :: Ef, qtot, temp
    integer, intent(in), optional :: prec

    external :: io_assign, io_close

    ! Internal variables and arrays
    integer :: no_u
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    integer :: iu, is, io, nspin
    logical :: is_dp
    integer, allocatable, target :: gncol(:)

    call timer("writeHSX",1)

    dit => dist(H)
    sp => spar(H)
    ! total number of rows
    no_u = nrows_g(sp)
    nspin = size(H, 2)

    ! Check which precision
    is_dp = .true.
    if ( present(prec) ) is_dp = prec == dp

    if ( Node == 0 ) then

      ! Open file
      call io_assign( iu )
      open( iu, file=trim(slabel)//'.HSX', form='unformatted', status='unknown' )

      ! Write version specification (to easily distinguish between different versions)
      write(iu) 1
      ! And what precision
      write(iu) is_dp

      ! Write overall data
      write(iu) na_u, no_u, nspin, nspecies, nsc
      write(iu) ucell, Ef, qtot, temp

      write(iu) isc_off, xa(:,1:na_u), isa(1:na_u), lasto(1:na_u)

      ! Write other useful info
      write(iu) (labelfis(is),zvalfis(is),nofis(is),is=1,nspecies)
      do is = 1, nspecies
        write(iu) (cnfigfio(is,io), lofio(is,io), zetafio(is,io),io=1,nofis(is))
      end do

    end if

    allocate(gncol(no_u))
    ! Here we signal that the gncol should be globalized
    ! to the 1st node
    gncol(1) = -1
    
    ! Write sparsity pattern...
    call io_write(iu, sp, dit=dit, gncol=gncol)
    
    ! Write H and overlap
    if ( is_dp ) then
      call io_write(iu, H, gncol=gncol)
      call io_write(iu, S, gncol=gncol)
    else
      call io_write_r(iu, H, gncol=gncol)
      call io_write_r(iu, S, gncol=gncol)
    end if

    deallocate(gncol)

    if ( Node == 0 ) then
      call io_close( iu )
    end if

    call timer("writeHSX",2)

    end subroutine write_hsx

!-----------------------------------------------------------------
    subroutine write_hs_formatted(no_u, nspin, &
        maxnh, numh, listhptr, listh, H, S)
! *********************************************************************
! Saves the hamiltonian and overlap matrices in formatted form
! ONLY for gamma case
! ******************** INPUT
! integer no_u                : Number of basis orbitals per unit cell
! integer nspin               : Spin polarization (1 or 2)
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form

      use precision, only: dp, sp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb, &
          GlobalToLocalOrb, GetNodeOrbs
      use atm_types,    only : nspecies
#ifdef MPI
      use mpi_siesta
#endif

      integer           maxnh, no_u, nspin
      integer           listh(maxnh), numh(*), listhptr(*)
      real(dp)          H(maxnh,nspin), S(maxnh)
      external          io_assign, io_close


      integer    im, is, iu, ius, ju, k, mnh, ns, ia, io
      integer    ih,hl,nuo,maxnhtot,maxhg
      integer, dimension(:), allocatable :: numhg, hg_ptr
#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
#endif


      call timer("write_HS_fmt",1)

      ! Find total numbers over all Nodes
#ifdef MPI
      call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum, &
          MPI_Comm_World,MPIerror)
#else
      maxnhtot = maxnh
#endif

      if (Node.eq.0) then
        ! Open file
        call io_assign( iu )
        call io_assign( ius )
        open( iu, file="H.matrix", form='formatted', status='unknown', &
            position="rewind")      
        open( ius,file="S.matrix", form='formatted', status='unknown', &
            position="rewind")      

        ! Write overall data
        write(iu,*) no_u, no_u, maxnhtot
        write(ius,*) no_u, no_u, maxnhtot

        ! Allocate local array for global numh
        allocate(numhg(no_u))
        allocate(hg_ptr(no_u+1))
      endif

      ! Create globalised numh
      do ih = 1,no_u
#ifdef MPI
        call WhichNodeOrb(ih,Nodes,BNode)
        if (BNode.eq.0.and.Node.eq.BNode) then
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
          hl = ih
#endif
          numhg(ih) = numh(hl)
#ifdef MPI
        elseif (Node.eq.BNode) then
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
          call MPI_ISend(numh(hl),1,MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
        elseif (Node.eq.0) then
          call MPI_IRecv(numhg(ih),1,MPI_integer, &
              BNode,1,MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
        endif
        if (BNode.ne.0) then
          call MPI_Barrier(MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

      if (Node.eq.0) then
        ! Write row pointers
        maxhg = 0
        hg_ptr(1) = 1
        do ih = 1,no_u
          maxhg = max(maxhg,numhg(ih))
          hg_ptr(ih+1) = hg_ptr(ih) + numhg(ih)
        enddo
        write(iu,*) (hg_ptr(ih),ih=1,no_u+1)
        write(ius,*) (hg_ptr(ih),ih=1,no_u+1)
        deallocate(hg_ptr)
#ifdef MPI
        allocate(buffer(maxhg))
        call memory('A','D',maxhg,'iohs')
        allocate(ibuffer(maxhg))
        call memory('A','I',maxhg,'iohs')
#endif
      endif

      ! Write listh
      do ih = 1,no_u
#ifdef MPI
        call WhichNodeOrb(ih,Nodes,BNode)
        if (BNode.eq.0.and.Node.eq.BNode) then
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
          hl = ih
#endif
          write(iu,*) (listh(listhptr(hl)+im),im = 1,numh(hl))
          write(ius,*) (listh(listhptr(hl)+im),im = 1,numh(hl))
#ifdef MPI
        elseif (Node.eq.0) then
          call MPI_IRecv(ibuffer,numhg(ih),MPI_integer,BNode,1, &
              MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
        elseif (Node.eq.BNode) then
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
          call MPI_ISend(listh(listhptr(hl)+1),numh(hl),MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
        endif
        if (BNode.ne.0) then
          call MPI_Barrier(MPI_Comm_World,MPIerror)
          if (Node.eq.0) then
            write(iu,*) (ibuffer(im),im = 1,numhg(ih))
            write(ius,*) (ibuffer(im),im = 1,numhg(ih))
          endif
        endif
#endif
      enddo

#ifdef MPI
      if (Node.eq.0) then
        call memory('D','I',size(ibuffer),'iohs')
        deallocate(ibuffer)
      endif
#endif

      ! Write Hamiltonian
      do is=1,nspin
        do ih=1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            write(iu,*) (real(H(listhptr(hl)+im,is),kind=sp), &
                im=1,numh(hl))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(H(listhptr(hl)+1,is),numh(hl), &
                MPI_double_precision,0,1,MPI_Comm_World, &
                Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
              write(iu,*) (real(buffer(im),kind=sp),im=1,numhg(ih))
            endif
          endif
#endif
        enddo
      enddo

      if (node == 0) then
        call io_close(iu)
      endif

      ! Write Overlap matrix
      do ih = 1,no_u
#ifdef MPI
        call WhichNodeOrb(ih,Nodes,BNode)
        if (BNode.eq.0.and.Node.eq.BNode) then
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
          hl = ih
#endif
          write(ius,*) (real(S(listhptr(hl)+im),kind=sp), im = 1,numh(hl))
#ifdef MPI
        elseif (Node.eq.0) then
          call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
              BNode,1,MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
        elseif (Node.eq.BNode) then
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
          call MPI_ISend(S(listhptr(hl)+1),numh(hl), &
              MPI_double_precision,0,1,MPI_Comm_World, &
              Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
        endif
        if (BNode.ne.0) then
          call MPI_Barrier(MPI_Comm_World,MPIerror)
          if (Node.eq.0) then
            write(ius,*) (real(buffer(im),kind=sp),im=1,numhg(ih))
          endif
        endif
#endif
      enddo

#ifdef MPI
      if (Node .eq. 0) then
        call memory('D','D',size(buffer),'iohs')
        deallocate(buffer)
        deallocate(numhg)   
      endif
#endif

      if (node == 0) then
        call io_close( ius )
      endif

      call timer("write_HS_fmt",2)

    end subroutine write_hs_formatted

end module io_hsx_m
