
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
#include "mpi_macros.f"

      module m_dscfComm
#ifdef MPI
        use MPI_SIESTA
#endif
        !! This module stores arrays and sizes used when communicating matrices
        !! between an orbital-based distribution and a meshpoint-based
        !! distribution.
        !!
        !! Written by R. Grima (BSC) Dec.2007
      implicit none
      private

      integer, public :: DCsize
        !! Total number of communications.
      integer, public :: DCncom
        !! Total number of communications for this node.
      integer, public :: DCself
        !! Number of orbitals that the current process should keep, i.e.
        !! orbitals which are in the same node for both distributions.
      integer, public :: DCmax
        !! Maximum number of orbitals that we should send/receive
        !! from/to another process.
      integer, public :: DCmaxnd
        !! Maximum size of matrix that we should send/receive from/to
        !! another process.
      integer, public :: DCtotal
        !! Total number of orbitals of the local matrix.
      integer, public, pointer :: DCsrc(:)
        !! For a given communication, the source node (size: DCsize).
      integer, public, pointer :: DCdst(:)
        !! For a given communication, the destination node (size: DCsize).
      integer, public, pointer :: DCsic(:)
        !! For a given communication, the number of orbitals exchanged
        !! (size: DCsize).
      integer, public, pointer :: DCBuffer(:)
        !! Buffer that contains all the data transferred.


      integer, public, pointer :: DCsnd(:)
        !! The sender for each communication involving this node (size: DCncom).
      integer, public, pointer :: DCrcv(:)
        !! The receiver for each communication involving this node
        !! (size: DCncom).
      integer, public, pointer :: DCinvp(:)
        !! Permutation vector used to transform the received data to the local
        !! ordering (size: DCtotal).
      integer, public, pointer :: DCpid(:)
        !! Unique ID for each communication involving this node. This is
        !! required for the non-blocking send/recv operations (size: DCncom).
      MPI_REQUEST_TYPE, public, pointer :: DCreq(:,:)
        !! MPI handles (requests) for non-blocking communications. This ID is
        !! stored to check whether the communications are done or not. The first
        !! index is 1 if the node is the sender, and 2 if it is the receiver
        !! (size: 2xDCncom).
      integer, public, pointer :: DCsne(:)
        !! Size of buffer for a given *send* communication.
      integer, public, pointer :: DCrce(:)
        !! Size of buffer for a given *recv* communication.

      public :: dscfComm
      public :: resetdscfComm
      public :: CommWaitSnd
      public :: CommWaitRcv

      CONTAINS

      ! Main subroutine.
      subroutine dscfComm( nuotot, nrowsDscfL, NeedDscfL )
        !! Precomputes the needed communications to transform a matrix that is
        !! distributed by orbitals to one distributed by mesh. All outputs are
        !! stored within this same module.
        !!
        !! We count the number of orbitals that we should receive from every
        !! process. Then we create a list of communications that every process
        !! should do.
      use parallel    , only : NODE, NODES
      use parallelsubs, only : WhichNodeOrb
      use alloc       , only : de_alloc, re_alloc
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
      integer, intent(in) :: nuotot
        !! Total number of basis orbitals in unit cell.
      integer            :: nrowsDscfL
        !! Local number of rows of matrix DSCF.
      integer            :: NeedDscfL(nuotot)
        !!  Permutation vector of orbitals from global to local.

      integer            :: ii, ind, io, oNode, ncom
#ifdef MPI
      integer              :: MPIerr
#endif
      integer, pointer   :: send(:), recv(:), tinv(:), tsrc(:), xnei(:),
     &                      invp(:)

      call timer( "dscfComm", 1 )

      ! Allocate memory for temporary buffers.
      nullify(send,recv,tinv,tsrc)
      call re_alloc( send, 0, NODES-1, 'send', 'dscfComm' )
      call re_alloc( recv, 0, NODES-1, 'recv', 'dscfComm' )
      call re_alloc( tinv, 1, nrowsDscfL, 'tinv', 'dscfComm' )
      call re_alloc( tsrc, 1, nrowsDscfL, 'tsrc', 'dscfComm' )
      do ii = 0, NODES-1
        send(ii) = 0
        recv(ii) = 0
      enddo

      ! Check which node has the orbitals we need. For each node, we store
      ! the amount of orbitals that we need from it.
      ii = 0
      do io = 1, nuotot
        if ( NeedDscfL( io ).ne.0 ) then
          call WhichNodeOrb( io, Nodes, oNode )
          recv(oNode) = recv(oNode) + 1
          ii = ii + 1
          tinv(ii) = io
          tsrc(ii) = oNode
        endif
      enddo

      ! We get "starting index" for orbitals that we need from other processes.
      ! The first node starts from 0, and then the other indexes correspond to
      ! the total mount of orbitals from the previous node.
      nullify(xnei)
      call re_alloc( xnei, 0, NODES, 'xnei', 'dscfComm' )
      xnei(0) = 0
      do ii= 1, NODES
        xnei(ii) = xnei(ii-1) + recv(ii-1)
      enddo

      ! We fill invp, which is the "reverse search" array, i.e. for the orbital
      ! ii that we need, that is contained in node tsrc(ii), we store its
      ! global index io. Note that in invp the orbitals that we need end up
      ! sorted by node, that's why we rely on xnei.
      nullify(invp)
      call re_alloc( invp, 1, nrowsDscfL, 'invp', 'dscfComm' )
      do ii= 1, nrowsDscfL
        io    = tinv(ii)
        oNode = tsrc(ii)
        xnei(oNode) = xnei(oNode) + 1
        invp(xnei(oNode)) = io
      enddo
      call de_alloc( tinv, 'tinv', 'dscfComm' )
      call de_alloc( tsrc, 'tsrc', 'dscfComm' )

      ! We rebuild xnei since we overwrote it in the previous step.
      xnei(0) = 0
      do ii= 1, NODES
        xnei(ii) = xnei(ii-1) + recv(ii-1)
      enddo


      ! We communicate all recvs and store them in "send". Then we calculate
      ! the amount of communications needed *by this node*.
#ifdef MPI
      call MPI_Alltoall( recv, 1, MPI_Integer, send, 1, MPI_Integer,
     &                   MPI_COMM_WORLD, MPIerr )
#endif

      ncom = 0
      do ii = 0, NODES-2
        oNode = MOD((ii/2+1)*
     &             MERGE(1,-1,IAND(NODE,ishft(1,ii/2))==0)*
     &             MERGE(1,-1,IAND(ii,1)==0)+NODES+NODE,NODES)
        if (send(oNode)/=0 .or. recv(oNode)/=0) ncom = ncom + 1
      enddo
      DCncom = ncom

      ! We create the communication-related arrays.
      call re_alloc( DCpid, 1, DCncom, 'DCpid', 'dscfComm' )
      call re_alloc( DCsnd, 1, DCncom, 'DCsnd', 'dscfComm' )
      call re_alloc( DCrcv, 1, DCncom, 'DCrcv', 'dscfComm' )
      call re_alloc( DCsne, 1, DCncom, 'DCsne', 'dscfComm' )
      call re_alloc( DCrce, 1, DCncom, 'DCrce', 'dscfComm' )
      allocate( DCreq(1:2,1:DCncom) )

      ! Here onode becomes a sort of "communication ID" for the current node.
      ! We store then for each communication: its ID, receiving node, and
      ! sending node.
      ncom = 0
      do ii = 0, NODES-2
        oNode = MOD((ii/2+1)*
     &             MERGE(1,-1,IAND(NODE,ishft(1,ii/2))==0)*
     &             MERGE(1,-1,IAND(ii,1)==0)+NODES+NODE,NODES)
        if (send(oNode)/=0 .or. recv(oNode)/=0) then
          ncom = ncom + 1
          DCpid(ncom) = oNode
          DCsnd(ncom) = send(oNode)
          DCrcv(ncom) = recv(oNode)
        endif
      enddo

      ! We store the total amount of data transferred in communications and the
      ! amount data that are internal for this node.
      DCself = send(Node)
      DCtotal = DCself + SUM(DCsnd) + SUM(DCrcv)

      call de_alloc( send, 'send', 'dscfComm' )
      call de_alloc( recv, 'recv', 'dscfComm' )


      ! Finally, we store DCinvp. This array ends up being sorted by
      ! communication now. For each communication, it contains the global io
      ! index for each orbital involved in said communication.
      nullify(DCinvp)
      call re_alloc( DCinvp, 1, DCtotal, 'DCinvp', 'dscfComm' )
      DCinvp(:DCself) = invp(xnei(node)+1:xnei(node+1))
      ind = DCself
#ifdef MPI
      do ii= 1, DCncom
        oNode = DCpid(ii)
        if (DCrcv(ii)/=0) then
          DCinvp(ind+1:ind+DCrcv(ii)) =
     &      invp(xnei(oNode)+1:xnei(oNode+1))
          call MPI_ISEND( DCinvp(ind+1:), DCrcv(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(1,ii), MPIerr )
          ind = ind + DCrcv(ii)
        endif
        if (DCsnd(ii)/=0) then
          call MPI_Irecv( DCinvp(ind+1:),  DCsnd(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(2,ii), MPIerr )
          ind = ind + DCsnd(ii)
        endif
      enddo

      call CommWaitAll( )
#endif
      call de_alloc( invp, 'invp', 'dscfComm' )
      call de_alloc( xnei, 'xnei', 'dscfComm' )

      call timer( "dscfComm", 2 )
      end subroutine dscfComm

      ! Helpers
      subroutine CommWaitAll( )
        !! Waits for all communications done by this node to finish.
        !! Uses the handles stored in DCreq.
#ifdef MPI
      use mpi_siesta

      implicit none
      integer :: ii, MPIerr
      MPI_STATUS_TYPE :: status

      do ii= 1, DCncom
        if (DCrcv(ii)/=0) call MPI_Wait( DCreq(1,ii), status, MPIerr )
        if (DCsnd(ii)/=0) call MPI_Wait( DCreq(2,ii), status, MPIerr )
      enddo
#endif
      end subroutine CommWaitAll

      subroutine CommWaitSnd( )
        !! Waits for all SEND communications done by this node to finish.
        !! Uses the handles stored in DCreq.
#ifdef MPI
      use mpi_siesta

      implicit none
      integer :: ii, MPIerr
      MPI_STATUS_TYPE :: status

      do ii= 1, DCncom
        if (DCsnd(ii)/=0) call MPI_Wait( DCreq(2,ii), status, MPIerr )
      enddo
#endif
      end subroutine CommWaitSnd

      subroutine CommWaitRcv( )
        !! Waits for all RECV communications done by this node to finish.
        !! Uses the handles stored in DCreq.
#ifdef MPI
      use mpi_siesta

      implicit none
      integer :: ii, MPIerr
      MPI_STATUS_TYPE :: status

      do ii= 1, DCncom
        if (DCrcv(ii)/=0) call MPI_Wait( DCreq(1,ii), status, MPIerr )
      enddo
#endif
      end subroutine CommWaitRcv

      subroutine resetdscfComm( )
        !! Resets the buffer and inverse lookup arrays (misleading name).
      use alloc, only : de_alloc
      implicit none
      call de_alloc( DCBuffer, 'DCBuffer', 'dscfComm' )
      call de_alloc( DCinvp, 'DCinvp', 'dscfComm' )
      call de_alloc( DCpid, 'DCpid', 'dscfComm' )
      call de_alloc( DCsnd, 'DCsnd', 'dscfComm' )
      call de_alloc( DCrcv, 'DCrcv', 'dscfComm' )
      call de_alloc( DCsne, 'DCsne', 'dscfComm' )
      call de_alloc( DCrce, 'DCrce', 'dscfComm' )

      if ( associated(DCreq) ) then
        deallocate(DCreq)
        nullify(DCreq)
      endif
      end subroutine resetdscfComm

      end module m_dscfComm
