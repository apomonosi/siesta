!
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
#include "mpi_macros.f"

      module meshdscf
        !! Stores quantities that are connected with Dscf in mesh local
        !! form when data is distributed for parallel execution.
        !!
        !! This module first contains a routine to set up the arrays and
        !! communications needed for changing a matrix between an orbital
        !! (block-cyclic) distribution, and a local, mesh-related distributions.
        !!
        !! Then it contains other routines to handle the transformation itself
        !! to/from a mesh (M) distributed matrix from/to an orbital (O)
        !! distributed matrix.

      use precision, only: dp
      implicit none

      ! Variables here largely imitate the sparsity structure outside,
      ! but for the local copy of the density matrix, DscfL.

      integer, public            :: nrowsDscfL
        !! Number of rows of DscfL.
      integer,  pointer, public  :: listdl(:)
        !! List of non-zero elements in a row of DscfL
      integer,  pointer, public  :: listdlptr(:)
        !! Pointer to row in listdl.
      integer,  pointer, public  :: NeedDscfL(:)
        !! Pointer as to whether a row of Dscf is needed in DscfL.
      integer,  pointer, public  :: numdl(:)
        !! Number of non-zero elements in a row of DscfL.
      real(dp), pointer, public  :: DscfL(:,:)
        !! Local copy of Dscf elements needed for the local mesh.
      logical           :: first_time = .true.

      public :: matrixOtoM
      public :: matrixMtoO
      public :: matrixMtoOC
      public :: resetDscfPointers
      public :: CreateLocalDscfPointers

      private

      CONTAINS

      subroutine resetDscfPointers( )
        !! Resets all of this module's pointers.
      use alloc     , only : de_alloc
      use m_dscfComm, only : resetdscfComm
      implicit none
      call resetdscfComm( )
      call de_alloc( listdl,    'listdl',    'meshdscf' )
      call de_alloc( listdlptr, 'listdlptr', 'meshdscf' )
      call de_alloc( NeedDscfL, 'NeedDscfL', 'meshdscf' )
      call de_alloc( numdl,     'numdl',     'meshdscf' )
      call de_alloc( DscfL,     'DscfL',     'meshdscf' )
      end subroutine resetDscfPointers

      subroutine CreateLocalDscfPointers( nmpl, nuotot, numd,
     &                                    listdptr, listd )
        !! Computes the communications needed to move data ordered by orbitals
        !! to data ordered by mesh (function dscfComm), and sets up the
        !! sparsity arrays for the local copy of Dscf: numdl, listdl and
        !! listdlptr, which are part of this module.
        !!
        !! All-to-all communications has been substituted by point-to-point
        !! communications.
        !! Written by Rogeli Grima (BSC) Dec.2007
      use atomlist,     only : indxuo
      use meshphi,      only : endpht, lstpht
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use precision   , only : dp
      use alloc       , only : re_alloc, de_alloc
      use m_dscfComm,   only : dscfcomm, CommWaitRcv, CommWaitSnd,
     &                         DCtotal, DCmaxnd, DCncom, DCpid,
     &                         DCsnd, DCrcv, DCreq,
     &                         DCinvp, DCself, DCrce, DCsne

#ifdef MPI
      use mpi_siesta
#endif
      implicit none

      ! Input variables
      integer, intent(in) :: nmpl
        !! Number of mesh points in unit cell for this node (local).
      integer, intent(in) :: nuotot
        !! Total number of basis orbitals in unit cell.
      integer, intent(in) :: numd(*)
        !!  Number of nonzero density-matrix elements for each matrix row.
      integer, intent(in) :: listdptr(*)
        !! Pointer to start of rows of density-matrix.
      integer, intent(in) :: listd(*)
        !! Nonzero-density-matrix-element column indexes for each matrix row.

      ! Local variables
      integer          :: i, ii, j, io, iio, ip, imp, iu, numdele,
     &                    nsize, mm, oNode
      integer, pointer :: ibuffer(:)
#ifdef MPI
      integer          :: MPIerror
      MPI_STATUS_TYPE  :: Status
#endif

#ifdef DEBUG
      call write_debug( '      PRE CreateLocalDscfPointers' )
#endif

      if (first_time) then
        first_time = .false.
      else
        call resetDscfPointers( )
      endif

      nullify( NeedDscfL, listdl, listdlptr, numdl )

      ! Create pointer as to whether a given row of DscfL is needed in NeedDscfL
      ! NOTE: This pointer is never deallocated...
      call re_alloc( NeedDscfL, 1, nuotot, 'NeedDscfL', 'meshdscf' )

      !-------------------------------------------------------------------------
      ! STEP 1
      !   Calculate which orbitals are needed by the current node, according to
      !   its local grid points. Then, setup a communication schedule.

      ! This essentially checks if any of the orbitals is needed by any of
      ! the grid points assigned to this node.
      do io= 1, nuotot
        NeedDscfL(io) = 0
      enddo

      do ip = 1,nmpl
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          iu = indxuo(i)
          NeedDscfL(iu) = 1
        enddo
      enddo

      ! Then we calculate the amount of rows needed by the "local" copy of Dscf.
      nrowsDscfL = 0
      do i = 1,nuotot
        if (NeedDscfL(i).eq.1) then
          nrowsDscfL = nrowsDscfL + 1
          NeedDscfL(i) = nrowsDscfL
        endif
      enddo

      ! Computes the communications needed to move data ordered
      ! by orbitals to data ordered by mesh.
      call dscfComm( nuotot, nrowsDscfL, NeedDscfL )

      ! Allocate/reallocate memory for numdl and listdlptr
      call re_alloc( numdl, 1, max(1,nrowsDscfL), 'numdl', 'meshdscf' )

      !-------------------------------------------------------------------------
      ! STEP 2
      !   Build numdl, communicating the numh values from other nodes as much
      !   as required.

      ! This data is in the current process. No communications needed
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb( io, Node, Nodes, iio )
        numdele              = numd(iio)
        numdl(NeedDscfL(io)) = numdele
      enddo

      ! Distribute information about numd globally
      ! Use communications scheduling precomputed in dscfComm
      nullify( ibuffer )
      call re_alloc( ibuffer, 1, DCtotal, 'ibuffer', 'meshdscf' )

      i = DCself     ! index of ibuffer and DCinvp
      numdele = 0
#ifdef MPI
      ! First we do all asynchronous communications, then we process all
      ! of the data received.

      do ii= 1, DCncom ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Asynchronous receive: Process data in the next loop
          call MPI_Irecv( ibuffer(i+1), DCrcv(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(1,ii), MPIerror )
          i = i + DCrcv(ii)             ! Update ibuffer/DCinvp index
        endif
        if (DCsnd(ii)/=0) then
          ! Process data into a buffer and do an asynchronous send
          do j= 1, DCsnd(ii)            ! Loop orbitals needed by oNode
            io = DCinvp(i+j)
            call GlobalToLocalOrb( io, Node, Nodes, iio )
            ibuffer(i+j) = numd(iio)
            numdele = numdele + ibuffer(i+j)
          enddo
          call MPI_ISEND( ibuffer(i+1), DCsnd(ii), MPI_INTEGER,
     &      oNode, 0, MPI_Comm_World, DCreq(2,ii), MPIerror )
          i = i + DCsnd(ii)             ! Update ibuffer/DCinvp index
        endif
      enddo

      ! Process asynchronous received data
      i = DCself       ! Reset buffer index
      do ii= 1, DCncom ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Wait until Receive communication is complete
          call MPI_Wait( DCreq(1,ii), status, MPIerror )
          do j= 1, DCrcv(ii)  ! Loop orbitals received from oNode
            i = i + 1
            io = DCinvp(i)
            numdl(NeedDscfL(io)) = ibuffer(i)
            numdele = numdele + ibuffer(i)
          enddo
        endif
        if (DCsnd(ii)/=0) i = i + DCsnd(ii)   ! Update buffer index
      enddo
      ! Check that all asynchronous Sent communication has finished
      call CommWaitSnd( )
#endif
      call de_alloc( ibuffer, 'ibuffer', 'meshdscf' )
      DCmaxnd = numdele

      !-------------------------------------------------------------------------
      ! STEP 3
      !   Build listdlptr, the array containing pointers to the corresponding
      !   elements in listh, using numdl.
      call re_alloc( listdlptr, 1, max(1,nrowsDscfL),
     &               'listdlptr', 'meshdscf' )
      listdlptr(1) = 0
      do io = 2,nrowsDscfL
        listdlptr(io) = listdlptr(io-1) + numdl(io-1)
      enddo


      !-------------------------------------------------------------------------
      ! STEP 4
      !   Build listdl, the local copy of listd as required by the mesh points.
      !   Similarly to numdl, we first set up the elements present in the local
      !   node, and then we communicate the components of listd that we need
      !   from other nodes.

      ! Allocate/reallocate listdl
      if (nrowsDscfL.gt.0) then
        nsize = listdlptr(nrowsDscfL)+numdl(nrowsDscfL)
      else
        nsize = 1
      endif
      call re_alloc( listdl, 1, nsize, 'listdl', 'meshdscf' )

      ! This data is in the current process. No communications needed
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        do i = 1,numd(iio)
          listdl(listdlptr(iu)+i) = listd(listdptr(iio)+i)
        enddo
      enddo

      ! Distribute information about listd globally, using communications
      ! scheduling precomputed in dscfComm
      nullify( ibuffer )
      call re_alloc( ibuffer, 1, DCmaxnd, 'ibuffer', 'meshdscf' )

      i       = DCself                  ! Reset index of DCinvp
      numdele = 0                       ! Reset index of ibuffer
#ifdef MPI
      ! First we do all asynchronous communications, then we process all
      ! of the data received.

      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Count the size of the matrix that we will receive
          DCrce(ii) = 0
          do j= 1, DCrcv(ii)
            i = i + 1
            DCrce(ii) = DCrce(ii) + numdl(NeedDscfL(DCinvp(i)))
          enddo
          ! Asynchronous receive: Process data in the next loop
          call MPI_Irecv( ibuffer(numdele+1), DCrce(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(1,ii), MPIerror )
          numdele = numdele + DCrce(ii) ! Update ibuffer index
        endif
        if (DCsnd(ii)/=0) then
          ! Process data into a buffer and do an asynchronous send
          mm = numdele
          do j= 1, DCsnd(ii)            ! Loop orbitals needed by oNode
            ! Save matrix row in ibuffer
            i = i + 1
            io = DCinvp(i)
            call GlobalToLocalOrb( io, Node, Nodes, iio )
            ibuffer(mm+1:mm+numd(iio)) = listd(listdptr(iio)+1:
     &                                         listdptr(iio)+numd(iio))
            mm = mm + numd(iio)
          enddo
          DCsne(ii) = mm - numdele      ! Size of the communication
          call MPI_ISend( ibuffer(numdele+1), DCsne(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(2,ii), MPIerror )
          numdele = mm                  ! Update ibuffer index
        endif
      enddo

      ! Process asynchronous received data
      i       = DCself                  ! Reset index of DCinvp
      numdele = 0                       ! Reset index of ibuffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Wait until Receive communication is complete
          call MPI_Wait( DCreq(1,ii), status, MPIerror )
          do j= 1, DCrcv(ii)  ! Loop orbitals received from oNode
            i = i + 1
            io = DCinvp(i)
            iu = NeedDscfL(io)
            ! Update local matrix structure
            listdl(listdlptr(iu)+1:listdlptr(iu)+numdl(iu)) =
     &              ibuffer(numdele+1:numdele+numdl(iu))
            numdele = numdele + numdl(iu) ! Update buffer index
          enddo
        endif
        if (DCsnd(ii)/=0) then
          numdele = numdele + DCsne(ii) ! Update buffer index
          i = i + DCsnd(ii)             ! Update DCinvp index
        endif
      enddo
      ! Check that all asynchronous Sent comm has finished
      call CommWaitSnd( )
#endif
      call de_alloc( ibuffer, 'ibuffer', 'meshdscf' )

#ifdef DEBUG
      call write_debug( '      POS CreateLocalDscfPointers' )
#endif
      end subroutine CreateLocalDscfPointers

      subroutine matrixOtoM( maxnd, numd, listdptr, maxndl, nuo,
     &                       nspin, Dscf, Dscf_L )
        !! Transforms a matrix which is distributed by block cyclic distribution
        !! of orbitals to a matrix that contains all the orbital rows needed for
        !! a mesh point distribution over the nodes.
        !! Created by J.D.Gale, February 2000
        !!
        !! Update: All-to-all communications has been substituted by
        !! point-to-point communications. These have been precomputed in
        !! dscfComm.
        !! Written by Rogeli Grima (BSC) Dec.2007
      use precision   , only : dp
      use alloc       , only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel    , only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
      use m_dscfComm  , only : DCself, DCPid, DCinvp, DCncom, DCmaxnd,
     &                         DCrcv, DCsnd, DCrce, DCsne, DCreq,
     &                         CommWaitSnd
#endif
      implicit none

      ! I/O variables.
      integer, intent(in)   :: maxnd
        !! First dimension of Dscf, i.e. the total number of non-zero elements.
      integer, intent(in)   :: maxndl
        !! Same as maxnd but for the local copy of Dscf.
      integer, intent(in)   :: nspin
        !! Number of spin components.
      integer, intent(in)   :: nuo
        !! Local number of orbitals in unit cell.
      integer, intent(in)   :: numd(nuo)
        !! Number of non-zero elements in row of Dscf.
      integer, intent(in)   :: listdptr(nuo)
        !! Pointer to start of rows in Dscf.
      real(dp), intent(in)  :: Dscf(maxnd,nspin)
        !! Matrix in orbital distributed form.
      real(dp), intent(out) :: Dscf_L(maxndl,nspin)
        !! Matrix in mesh distributed form.

      ! Internal variables and arrays
      integer           :: io, ispin
#ifdef MPI
      integer           :: i, j, ii, iu, nn, iio, ini, oNode, MPIerror
      MPI_STATUS_TYPE ::  Status
      real(dp), pointer :: buffer(:)
#else
      integer           :: il
#endif

      call timer( "OtoM", 1 )

#ifdef MPI
      !---------------------
      ! PARALLEL EXCECUTION.
      !---------------------
      ! Here we first copy the data from the local node, and then we move
      ! onto the communications needed between nodes.

      ! Copy the data that is in the current process
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        do ispin = 1,nspin
          Dscf_L(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),ispin) =
     &    Dscf(listdptr(iio)+1:listdptr(iio)+numd(iio),ispin)
        enddo
      enddo
      nn = DCself

      ! Allocate local Dscf storage array
      nullify(buffer)
      call re_alloc( buffer, 1, DCmaxnd*nspin, 'buffer', 'meshdscf' )

      i  = DCself                       ! Reset index of DCinvp
      nn = 0                            ! Reset index of buffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Asynchronous receive: Process data in the next loop
          call MPI_Irecv( buffer(nn+1), DCrce(ii)*nspin,
     &      MPI_double_precision, oNode, 0, MPI_COMM_WORLD,
     &      DCreq(1,ii), MPIerror )
          nn = nn + DCrce(ii)*nspin ! Update buffer index
          i  = i + DCrcv(ii)        ! Update DCinvp index
        endif
        if (DCsnd(ii)/=0) then
          ! Process data into a buffer and do an asynchronous send
          ini = nn + 1
          do ispin = 1,nspin
            do j= 1, DCsnd(ii)
              io = DCinvp(i+j)
              call GlobalToLocalOrb( io, Node, Nodes, iio )
              buffer(nn+1:nn+numd(iio)) = Dscf(listdptr(iio)+1:
     &                           listdptr(iio)+numd(iio),ispin)
              nn = nn + numd(iio)       ! Update buffer index
            enddo
          enddo
          call MPI_ISend( buffer(ini), DCsne(ii)*nspin,
     &      MPI_double_precision, oNode, 0, MPI_COMM_WORLD,
     &      DCreq(2,ii), MPIerror )
          i  = i + DCsnd(ii)            ! Update DCinvp index
        endif
      enddo

      ! Process asynchronous received data
      i  = DCself                       ! Reset index of DCinvp
      nn = 0                            ! Reset index of ibuffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Wait until Receive communication is complete
          call MPI_Wait( DCreq(1,ii), status, MPIerror )

          do ispin = 1,nspin
            do j= 1, DCrcv(ii)          ! Loop orbitals received from oNode
              io = DCinvp(i+j)
              iu = NeedDscfL(io)
              ! Save received values in the local matrix
              Dscf_L(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),ispin) =
     &              buffer(nn+1:nn+numdl(iu))
              nn = nn + numdl(iu)
            enddo
          enddo
          i  = i + DCrcv(ii)        ! Update DCinvp index
        endif
        if (DCsnd(ii)/=0) then
          nn = nn + DCsne(ii)*nspin     ! Update buffer index
          i = i + DCsnd(ii)             ! Update DCinvp index
        endif
      enddo
      ! Check that all asynchronous Sent comm has finished
      call CommWaitSnd( )

      call de_alloc( buffer, 'buffer', 'meshdscf' )

#else
      !-------------------
      ! SERIAL EXCECUTION.
      !-------------------

      ! Loop over rows of Dscf checking to see if they are in Dscf_L
      do ispin = 1,nspin
        do io = 1,nuo

          ! Get pointer for this row of Dscf and see if it is needed for Dscf_L
          il = NeedDscfL(io)
          if (il.gt.0) then
            Dscf_L(listdlptr(il)+1:listdlptr(il)+numdl(il),ispin) =
     &        Dscf(listdptr(io)+1:listdptr(io)+numdl(il),ispin)
          endif

        enddo
      enddo
#endif
      call timer( "OtoM", 2 )
      end subroutine matrixOtoM


      subroutine matrixMtoO( maxnvl, maxnv, numVs, listVsptr, nuo,
     &                       nspin, VsL, Vs )
        !! Transforms a matrix which is distributed by mesh points to a matrix
        !! that is distributed by a block cyclic distribution over the orbitals.
        !! It is important to note that the values will be ADDED to the input
        !! matrix, this subroutine does not reset the orbital-distributed matrix.
        !! Created by J.D.Gale, February 2000
        !!
        !! Update: All-to-all communications has been substituted by
        !! point-to-point communications. These have been precomputed in
        !! dscfComm.
        !! Written by Rogeli Grima (BSC) Dec.2007
      use precision   , only : dp
      use alloc       , only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel    , only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb,
     &                         LocalToGlobalOrb
      use m_dscfComm  , only : DCself, DCPid, DCinvp, DCncom, DCmaxnd,
     &                         DCrcv, DCsnd, DCrce, DCsne, DCreq,
     &                         CommWaitRcv
#endif

      implicit none

      ! I/O variables.
      integer , intent(in)    :: maxnv
        !! First dimension of Vs, i.e. the maximum number of non-zero elements.
      integer , intent(in)    :: maxnvl
        !! First dimension of VsL, i.e. the maximum number of non-zero elements.
      integer , intent(in)    :: nspin
        !! Number of spin components.
      integer , intent(in)    :: nuo
        !! Local number of orbitals in unit cell
      integer , intent(in)    :: numVs(nuo)
        !! Number of non-zero elements in each row of Vs.
      integer , intent(in)    :: listVsptr(nuo)
        !! Pointer to the start of each row in Vs.
      real(dp), intent(in)    :: VsL(maxnvl,nspin)
        !! Mesh contributions to be added to Vs.
      real(dp), intent(inout) :: Vs(maxnv,nspin)
        !! Orbital-distributed matrix, to which the mesh contributions
        !! are summed up.

      ! Internal variables and arrays
      integer :: i, iu, ispin

#ifdef MPI
      integer           :: j, ii, ini, MPIerror, nn, io, iio, oNode
      MPI_STATUS_TYPE   :: Status
      real(dp), pointer :: buffer(:)
#endif
      call timer( "MtoO", 1 )

#ifdef MPI
      !---------------------
      ! PARALLEL EXCECUTION.
      !---------------------
      ! Here we first copy the data from the local node, and then we move
      ! onto the communications needed between nodes.

      ! Copy the data that is in the current process from VsL to Vs
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        do ispin = 1,nspin
          Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) =
     &    Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) +
     &    VsL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),ispin)
        enddo
      enddo
      nn = DCself

      ! Allocate a buffer to Send/Receive data
      nullify(buffer)
      call re_alloc( buffer, 1, DCmaxnd*nspin, 'buffer', 'meshdscf' )

      i  = DCself                       ! Reset index of DCinvp
      nn = 0                            ! Reset index of buffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          ! Process data into a buffer and do an asynchronous send
          ini = nn + 1
          do ispin = 1,nspin
            do j= 1, DCrcv(ii)
              io = DCinvp(i+j)
              iio = NeedDscfL(io)
              buffer(nn+1:nn+numdl(iio)) =
     &        VsL(listdlptr(iio)+1:listdlptr(iio)+numdl(iio),ispin)
              nn = nn + numdl(iio)       ! Update buffer index
            enddo
          enddo
          call MPI_ISend( buffer(ini), DCrce(ii)*nspin,
     &      MPI_double_precision, oNode, 0, MPI_COMM_WORLD,
     &      DCreq(1,ii), MPIerror )
          i  = i + DCrcv(ii)            ! Update DCinvp index
        endif
        if (DCsnd(ii)/=0) then
C         Asynchronous receive: Process data in the next loop
          call MPI_Irecv( buffer(nn+1), DCsne(ii)*nspin,
     &      MPI_double_precision, oNode, 0, MPI_COMM_WORLD,
     &      DCreq(2,ii), MPIerror )
          nn = nn + DCsne(ii)*nspin ! Update buffer index
          i  = i + DCsnd(ii)        ! Update DCinvp index
        endif
      enddo

      ! Process asynchronous received data
      i  = DCself                       ! Reset index of DCinvp
      nn = 0                            ! Reset index of ibuffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if (DCrcv(ii)/=0) then
          nn = nn + DCrce(ii)*nspin     ! Update buffer index
          i = i + DCrcv(ii)             ! Update DCinvp index
        endif

        if (DCsnd(ii)/=0) then
          ! Wait until Receive communication is complete.
          call MPI_Wait( DCreq(2,ii), status, MPIerror )
          do ispin = 1,nspin
            do j= 1, DCsnd(ii)          ! Loop orbitals received from oNode
              io = DCinvp(i+j)
              call GlobalToLocalOrb( io, Node, Nodes, iio )
              ! Accumulate received data into local VS
              Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) =
     &        Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) +
     &        buffer(nn+1:nn+numVs(iio))
              nn = nn + numVs(iio)      ! Update buffer index
            enddo
          enddo
          i  = i + DCsnd(ii)            ! Update DCinvp index
        endif
      enddo

      ! Check that all asynchronous Sent comm has finished
      call CommWaitRcv( )

      call de_alloc( buffer, 'buffer', 'meshdscf' )
#else
      !-------------------
      ! SERIAL EXCECUTION.
      !-------------------

      ! Add those elements that are needed locally to the values already
      ! stored in the orbital oriented array.
      do ispin = 1,nspin
        do i = 1,nuo
          iu = NeedDscfL(i)
          if (iu.gt.0) then
            Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) =
     &        Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) +
     &        VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),ispin)
          endif
        enddo
      enddo
#endif
      call timer( "MtoO", 2 )

      end subroutine matrixMtoO

      subroutine matrixMtoOC( maxnvl, maxnv, numVs, listVsptr, nuo,
     &                        VsL, Vs )
        !! Transforms a matrix which is distributed by mesh points to a matrix
        !! that is distributed by a block cyclic distribution over the orbitals.
        !! It is important to note that the values will be ADDED to the input
        !! matrix, this subroutine does not reset the orbital-distributed matrix.
        !!
        !! The main differnce with the real case is that this routine assumes
        !! that the first two spin components of the mesh-distributed matrix VsL
        !! are then collapsed in the complex matrix Vs; i.e., VsL(:,1) is the
        !! real part and VsL(:,2) is the imaginary part.
        !!
        !! Created by J.D.Gale, February 2000
        !! Modified for the complex case by J. Junquera, May 2012.
        !!
        !! Update: All-to-all communications has been substituted by
        !! point-to-point communications. These have been precomputed in
        !! dscfComm.
        !! Written by Rogeli Grima (BSC) Dec.2007

      use precision   , only : dp
      use alloc       , only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel    , only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb,
     &                         LocalToGlobalOrb
      use m_dscfComm  , only : DCself, DCPid, DCinvp, DCncom, DCmaxnd,
     &                         DCrcv, DCsnd, DCrce, DCsne, DCreq,
     &                         CommWaitRcv
#endif

      implicit none

      ! I/O variables.
      integer    , intent(in)    :: maxnv
        !! First dimension of Vs, i.e. the maximum number of non-zero elements.
      integer    , intent(in)    :: maxnvl
        !! First dimension of VsL, i.e. the maximum number of non-zero elements.
      integer    , intent(in)    :: nuo
        !! Local number of orbitals in unit cell
      integer    , intent(in)    :: numVs(nuo)
        !! Number of non-zero elements in each row of Vs.
      integer    , intent(in)    :: listVsptr(nuo)
        !! Pointer to the start of each row in Vs.
      real(dp)   , intent(in)    :: VsL(maxnvl,2)
        !! Mesh contributions to be added to Vs. The number of spin components
        !! is assumed to be 2.
      complex(dp), intent(inout) :: Vs(maxnv)
        !! Orbital-distributed matrix, to which the mesh contributions
        !! are summed up.

      integer :: i, iu, irelim
#ifdef MPI
      integer           :: jj, ii, ini, MPIerror, nn, io, iio, oNode
      MPI_STATUS_TYPE   :: Status
      real(dp), pointer :: buffer(:)
#endif
      external :: timer

      call timer( 'MtoOC', 1 )

#ifdef MPI
      !---------------------
      ! PARALLEL EXCECUTION.
      !---------------------
      ! Here we first copy the data from the local node, and then we move
      ! onto the communications needed between nodes.

      !  Copy the data that is in the current process from VsL to Vs
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) =
     &    Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) +
     &    cmplx(VsL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),1),
     &          VsL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),2), kind=dp)
      enddo

      ! Allocate a buffer to Send/Receive data
      nullify(buffer)
      call re_alloc( buffer, 1, DCmaxnd*2, 'buffer', 'meshdscf' )

      i  = DCself                       ! Reset index of DCinvp
      nn = 0                            ! Reset index of buffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if ( DCrcv(ii) /= 0 ) then
          ! Process data into a buffer and do an asynchronous send
          ini = nn + 1
          do irelim = 1, 2
            do jj = 1, DCrcv(ii)
              io  = DCinvp(i+jj)
              iio = NeedDscfL(io)
              buffer(nn+1:nn+numdl(iio)) =
     &        VsL(listdlptr(iio)+1:listdlptr(iio)+numdl(iio),irelim)
              nn = nn + numdl(iio)       ! Update buffer index
            enddo
          enddo

          call MPI_ISend( buffer(ini), DCrce(ii)*2,
     &      MPI_double_precision, oNode, 0, MPI_COMM_WORLD,
     &      DCreq(1,ii), MPIerror )
          i  = i + DCrcv(ii)            ! Update DCinvp index
        endif
        if ( DCsnd(ii) /= 0 ) then
          ! Asynchronous receive: Process data in the next loop
          call MPI_Irecv( buffer(nn+1), DCsne(ii)*2,
     &      MPI_double_precision, oNode, 0, MPI_COMM_WORLD,
     &      DCreq(2,ii), MPIerror )
          nn = nn + DCsne(ii)*2 ! Update buffer index
          i  = i + DCsnd(ii)    ! Update DCinvp index
        endif
      enddo

      ! Process asynchronous received data
      i  = DCself                       ! Reset index of DCinvp
      nn = 0                            ! Reset index of ibuffer
      do ii= 1, DCncom                  ! For all the needed communications
        oNode = DCpid(ii)               ! Process Id of the other NODE
        if ( DCrcv(ii) /= 0 ) then
          nn = nn + DCrce(ii)*2     ! Update buffer index
          i = i + DCrcv(ii)         ! Update DCinvp index
        endif

        if ( DCsnd(ii) /= 0 ) then
          ! Wait until Receive communication is complete
          call MPI_Wait( DCreq(2,ii), status, MPIerror )

          do irelim = 1, 2
            do jj = 1, DCsnd(ii)     ! Loop orbitals received from oNode
              io = DCinvp(i+jj)
              call GlobalToLocalOrb( io, Node, Nodes, iio )

              if( irelim == 1 ) then
                ! Accumulate received data into local VS
                Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) =
     &            Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) +
     &            buffer(nn+1:nn+numVs(iio))

              elseif ( irelim == 2 ) then
                Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) =
     &            Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) +
     &            buffer(nn+1:nn+numVs(iio)) *
     &            cmplx(0.0_dp,1.0_dp,kind=dp)
              endif
              nn = nn + numVs(iio)      ! Update buffer index
            enddo
          enddo
          i  = i + DCsnd(ii)
        endif
      enddo

      ! Check that all asynchronous Sent comm has finished
      call CommWaitRcv( )

      call de_alloc( buffer, 'buffer', 'meshdscf' )
#else
      !-------------------
      ! SERIAL EXCECUTION.
      !-------------------

      ! Add those elements that are needed locally to the values already
      ! stored in the orbital oriented array.
      do i = 1,nuo
        iu = NeedDscfL(i)
        if (iu.gt.0) then
          Vs(listVsptr(i)+1:listVsptr(i)+numVs(i)) =
     &      Vs(listVsptr(i)+1:listVsptr(i)+numVs(i)) +
     &      cmplx(VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),1),
     &            VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),2),kind=dp)
        endif
      enddo
#endif
      call timer( 'MtoOC', 2 )

      end subroutine matrixMtoOC

      end module meshdscf

