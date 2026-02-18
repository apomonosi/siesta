! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
#include "mpi_macros.f"

module io_hsx_m

  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData1D
  use class_dSpData2D
  use kpoint_t_m

  implicit none

  public :: write_hsx
  public :: write_hs_formatted
  public :: HSX_version

  interface write_hsx
     module procedure write_hsx_1
     module procedure write_hsx_2
  end interface

  private

contains

  subroutine write_hsx_header_internal(iu, dit, sp, nspin, Ef, qtot, temp, kpoint, is_dp, gncol)

    !! Writes a common header for the HSX file
    !!
    !! This is shared usage of the write_hsx_1|2 methods for reduced code duplication.
    !! The header is generic for either version == 1 or 2.
    !! The header also contain basis set information which internally is a bit of a mess,
    !! because it infers a lot of module dependencies.
    !! Ideally some of these details should be stored in a type, so it can be passed
    !! succinctly.
    !! Note that not all details in the header is *that* important.
    !!
    !! This subroutine only writes binary code.

    use precision, only: dp
    use io_sparse_m, only: io_write
    use parallel, only : Node
    use atm_types, only : nspecies
    use atomlist, only : iphorb, iaorb, lasto
    use siesta_geom, only: na_u, ucell, xa, isa, nsc, isc_off
    use atmfuncs, only : nofis, labelfis, zvalfis
    use atmfuncs, only : cnfigfio, lofio, zetafio

    integer, intent(in) :: iu
    !! The file-unit to write to.
    type(OrbitalDistribution), intent(inout) :: dit
    !! The distribution of the sparsity pattern.
    !! Required for the internal writing of the sparsity pattern.
    !! I.e. for a distributed matrix, it requires communication to share
    !! the column information.
    type(Sparsity), intent(inout) :: sp
    !! The (possibly distributed) sparsity pattern.
    !! This content is also part of the header of the file.
    integer, intent(in) :: nspin
    !! Number of spin-components of the Hamiltonian.
    real(dp), intent(in) :: Ef, qtot, temp
    !! Quantities describing the Hamiltonian, fermi-level, the number of electrons
    !! used to determine the fermi-level, and the electronic smearing temperature used.
    type(kpoint_t), intent(in) :: kpoint
    !! K-point sampling of the Hamiltonian
    logical, intent(in) :: is_dp
    !! Whether to store the matrix data is double- or single-precision.
    integer, intent(inout) :: gncol(:)
    !! Global variable used to speed up the subsequent matrix writing.
    !! This prohibits the communication of the `ncol` variable in the sparsity pattern,
    !! to speed up the process of writing.

    integer :: is, io
    !! Internal variables.

    if ( Node == 0 ) then

      ! Write version specification (to easily distinguish between different versions)
      write(iu) 2
      ! And what precision
      write(iu) is_dp

      ! Write overall data
      write(iu) na_u, nrows_g(sp), nspin, nspecies, nsc
      write(iu) ucell, Ef, qtot, temp

      write(iu) isc_off, xa(:,1:na_u), isa(1:na_u), lasto(1:na_u)

      ! Write other useful info
      write(iu) (labelfis(is),zvalfis(is),nofis(is),is=1,nspecies)
      do is = 1, nspecies
        write(iu) (cnfigfio(is,io), lofio(is,io), zetafio(is,io),io=1,nofis(is))
      end do

      ! Write k-point information
      write(iu) kpoint%k_cell, kpoint%k_displ

    end if

    ! Here we signal that the gncol should be globalized
    ! to the 1st node
    gncol(1) = -1

    ! Write sparsity pattern...
    call io_write(iu, sp, dit=dit, gncol=gncol)

  end subroutine write_hsx_header_internal


  subroutine write_hsx_1(H, S, Ef, qtot, temp, kpoint, prec, extension)

    !! Writing a HSX file with a single Hamiltonian.
    !!
    !! Saves the Hamiltonian and overlap matrices together with
    !! additional data describing the electronic structure + geometry for
    !! the simulation.
    !! The file contains:
    !!
    !!  - Fermi-level, to be able to correct determine DoS etc.
    !!  - SCF k-point sampling, to know how fine the details where for the
    !!    SCF cycles.
    !!  - SCF temperature to know about the smearing of the states, although
    !!    not complete, it at least gives some details to the stored data.
    !!
    !! The data can be stored in either double- or single-precision.
    !! This is determined by the flag `prec` which should be the `kind`
    !! of the stored data-type.
    !! Optionally the user may opt for a specific extension of the file.

    use precision, only: dp
    use io_sparse_m, only: io_write, io_write_r
    use parallel, only : Node, Nodes
    use files, only : slabel

    type(dSpData2D), intent(inout) :: H
    !! The Hamiltonian, which also contains the sparsity pattern.
    type(dSpData1D), intent(inout) :: S
    !! The overlap matrix, which also contains the sparsity pattern.
    real(dp), intent(in) :: Ef, qtot, temp
    !! Quantities describing the Hamiltonian, fermi-level, the number of electrons
    !! used to determine the fermi-level, and the electronic smearing temperature used.
    type(kpoint_t), intent(in) :: kpoint
    !! K-point sampling of the Hamiltonian
    integer, intent(in), optional :: prec
    !! The `kind` value of the precision to store the data in.
    character(len=*), intent(in), optional :: extension
    !! The extension for the file.

    external :: io_assign, io_close

    ! Internal variables and arrays
    integer :: no_u
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    integer :: iu, nspin
    logical :: is_dp
    integer, allocatable, target :: gncol(:)
    character(len=32) :: lext

    call timer("writeHSX",1)

    dit => dist(H)
    sp => spar(H)
    ! total number of rows
    no_u = nrows_g(sp)
    ! total number of spin-components
    nspin = size(H, 2)

    ! Check which precision
    is_dp = .true.
    if ( present(prec) ) is_dp = prec == dp

    lext = "HSX"
    if ( present(extension) ) then
       lext = trim(extension)
    end if

    if ( Node == 0 ) then

      call io_assign( iu )
      open( iu, file=trim(slabel)//'.'//trim(lext), form='unformatted', status='unknown' )

    end if

    ! Allocate array to reduce communication
    allocate(gncol(no_u))
    call write_hsx_header_internal(iu, dit, sp, nspin, Ef, qtot, temp, kpoint, is_dp, gncol)

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

  end subroutine write_hsx_1

  subroutine write_hsx_2(H1, H2, S, Ef, qtot, temp, kpoint, prec, extension)

    !! See `write_hsx_1` for details.
    !!
    !! The only difference for this subroutine is that it combines two Hamiltonian
    !! objects into a single one. Hence, the reading utility does not know whether
    !! the data originates from 1 or 2 Hamiltonians.

    use precision, only: dp
    use io_sparse_m, only: io_write, io_write_r
    use parallel, only : Node, Nodes
    use files, only : slabel

    type(dSpData2D), intent(inout) :: H1, H2
    type(dSpData1D), intent(inout) :: S
    real(dp) :: Ef, qtot, temp
    type(kpoint_t), intent(in) :: kpoint
    integer, intent(in), optional :: prec
    character(len=*), intent(in), optional :: extension

    external :: io_assign, io_close

    ! Internal variables and arrays
    integer :: no_u
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    integer :: iu, nspin1, nspin2
    logical :: is_dp
    integer, allocatable, target :: gncol(:)
    character(len=32) :: lext

    call timer("writeHSX",1)

    dit => dist(H1)
    sp => spar(H1)
    ! total number of rows
    no_u = nrows_g(sp)
    ! total number of spin-components for each Hamiltonian
    nspin1 = size(H1, 2)
    nspin2 = size(H2, 2)

    ! Check which precision
    is_dp = .true.
    if ( present(prec) ) is_dp = prec == dp

    lext = "HSX"
    if ( present(extension) ) then
       lext = trim(extension)
    end if

    if ( Node == 0 ) then

      call io_assign( iu )
      open( iu, file=trim(slabel)//'.'//trim(lext), form='unformatted', status='unknown' )

    end if

    ! Allocate array to reduce communication
    allocate(gncol(no_u))
    call write_hsx_header_internal(iu, dit, sp, nspin1+nspin2, Ef, qtot, temp, kpoint, is_dp, gncol)

    ! Write H and overlap
    if ( is_dp ) then
      call io_write(iu, H1, gncol=gncol)
      call io_write(iu, H2, gncol=gncol)
      call io_write(iu, S, gncol=gncol)
    else
      call io_write_r(iu, H1, gncol=gncol)
      call io_write_r(iu, H2, gncol=gncol)
      call io_write_r(iu, S, gncol=gncol)
    end if

    deallocate(gncol)

    if ( Node == 0 ) then
      call io_close( iu )
    end if

    call timer("writeHSX",2)

  end subroutine write_hsx_2

  function HSX_version(fname) result(version)
    !! Determine the version number of a `HSX` file.
    !!
    !! The return value is an integer that specifies the version.
    !! It can be 0, 1, ... so long as the specification for the file.
    !! is implemented in Siesta.

    use parallel,     only : Node

    character(len=*), intent(in) :: fname
    !! The filename to check the version number of.

    integer :: version
    integer :: iu
    integer :: no_u, no_s, nspin, n_nzs, err

    external :: io_assign, io_close

    ! Initialize
    version = -1
    if ( Node /= 0 ) return

    ! Open file
    call io_assign( iu )
    open( iu, file=fname, form='unformatted', status='unknown' )

    ! old version = 0 files only had these 4 numbers in the first line
    read(iu, iostat=err) no_u, no_s, nspin, n_nzs
    if ( err == 0 ) then
       ! we can successfully read 4 integers
       version = 0
    else
       rewind(iu)
       read(iu,iostat=err) version
    end if

    call io_close(iu)

  end function HSX_version

!-----------------------------------------------------------------
  subroutine write_hs_formatted(no_u, nspin, &
        maxnh, numh, listhptr, listh, H, S)
! *********************************************************************
! Saves the Hamiltonian and overlap matrices in formatted form
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
      integer    MPIerror, BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
      MPI_REQUEST_TYPE :: Request
      MPI_STATUS_TYPE :: Status
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
