!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
#include "mpi_macros.f"

      MODULE moreMeshSubs
        !! This module defines and handles different parallel mesh
        !! distributions. It contains routines to pre-compute, select and move
        !! data between distributions: uniform, quadratic and linear.
        !!
        !! Written by Rogeli Grima (BSC) Dec.2007.
        !!
        !! The current implementation allows different distributions to use
        !! different mesh-sub-divisions. This may be handy for external
        !! libraries. However, it is not fully implemented. Currently one cannot
        !! change between a SEQUENTIAL and CLUSTER distribution if they have
        !! different nsm values. It also won't work for KEEP with different nsm
        !! in clustered mode since there may be overlaps where the clustering of
        !! values is not complete.
        !!
        !! Nick Papior, June 2020.
        !!
        !! The  mesh-subdivisions work by relying on the nsm property of a
        !! distribution. This allows to distribute an (nsm>1, %box, %nMesh)
        !! distribution to another (nsm == 1, %box, %nMesh) distribution.
        !!
        !! Effectively this simply propagates to the boxintersection  and
        !! other related routines which uses the full grid. There are some
        !! restrictions on the kind of redistributions this code can handle.
        !!
        !! if in_distr%nsm == out_distr%nsm, *any* redistribution may be used
        !! if in_distr%nsm /= out_distr%nsm: *only* SEQUENTIAL data, using KEEP
        !!
        !! If data in clustered form needs to change nsm it needs to do 2
        !! redistributions:
        !!    1. TO_SEQUENTIAL
        !!    2. KEEP

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! IMPORTANT NOTICE                                            [[MBGN]] !!
      !!  To ease navigation, some tags have been implemented as comments:    !!
      !!   [[MINI]] - Distribution initializations.                           !!
      !!   [[MDST]] - Data exchange routines between distributions.           !!
      !!   [[MEXT]] - Mesh extencil routines.                                 !!
      !!   [[WKLD]] - Workload partition related routines.                    !!
      !!   [[MCMM]] - Mesh communication setup.                               !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use precision, only : grid_p, dp, i8b
      use parallel,  only : node, Nodes, ProcessorY
      use sys,       only : die
      use alloc,     only : re_alloc, de_alloc
      use fdf
#ifdef MPI
      use mpi_siesta
#endif
      implicit none

      ! initMeshDistr: Precompute a new data distribution for the grid mesh.
      !                to be used in Hamiltonian construction.
      ! setMeshDistr: Select a mesh distribution and set the grid sizes.
      ! distMeshData: Move data from one data distribution to another.
      ! allocASynBuffer: Allocate buffer for asynchronous communications if
      !                  necessary
      public :: initMeshDistr
      public :: setMeshDistr
      public :: allocExtMeshDistr
      public :: allocIpaDistr
      public :: distMeshData
      public :: resetMeshDistr
#ifdef MPI
      public :: initMeshExtencil
      public :: distExtMeshData
      public :: gathExtMeshData
#endif
      PUBLIC :: allocASynBuffer
      PUBLIC :: getMeshBox
!     Symbolic names for parallel mesh distributions

      ! Symbolic names for parallel mesh distributions
      integer, parameter, public :: UNIFORM   = 1
        !! Identifier for uniform distribution.
      integer, parameter, public :: QUADRATIC = 2
        !! Identifier for quadratic distribution.
      integer, parameter, public :: LINEAR    = 3
        !! Identifier for linear distribution.


      ! Symbolic names for "reord"-type operations.
      integer, parameter, public :: TO_SEQUENTIAL = +1
        !! Convert current distribution to sequential.
      integer, parameter, public :: TO_CLUSTER = -1
        !! Convert current distribution to clustered.
      integer, parameter, public :: KEEP       =  0
        !! Keep current distribution as it is.

      private
      interface initMeshDistr
        !! Interface for mesh distribution initialization, using either a
        !! workload array or by explicitely setting the limits for the
        !! distribution.
        module procedure initMeshDistr_wload
        module procedure initMeshDistr_explicit
      end interface initMeshDistr

      interface distMeshData
        !! Interface for mesh data distribution considering double precision
        !! reals and integers.
        module procedure distMeshData_rea
        module procedure distMeshData_int
      end interface distMeshData

      type meshDisType
        !! Private type to hold mesh distribution data
        integer :: nMesh(3)
          !! Number of mesh divisions across each axis.
        integer,  pointer :: box(:,:,:)
          !! Mesh box bounds of each node: box(1,iAxis,iNode) are the lower
          !! bounds and box(2,iAxis,iNode) are the upper bounds.
        integer :: nsm
          !! Number of mesh subdivisions for this particular distribution.
        integer,  pointer :: indexp(:)
          !! Contains the indices for the translation from the extended to the
          !! normal mesh indices (size: number of extended box points).
        integer,  pointer :: idop(:)
          !! Extended-mesh-index displacement of points of points within a
          !! a sphere of radius rmax (size: number of points where orbitals
          !! are non-zero).
        real(dp), pointer :: xdop(:,:)
          !! Relative distance of the point within the mesh cell.
        integer,  pointer :: ipa(:)
          !! Mesh cell where a given atom is.
      end type meshDisType

      type meshCommType
        !! Private type to hold communications to move data from one
        !! distribution to another.
        integer          :: ncom
          !! Number of needed communications.
        integer, pointer :: src(:) => null()
          !! Sources of communications
        integer, pointer :: dst(:) => null()
          !! Destination of communications
      end type meshCommType

      character(len=*),      parameter :: moduName = 'moreMeshSubs'
        !! Name of the module, used in message/printing routines.
      integer,               parameter :: maxDistr = 5
        !! Maximum number of data distributions that can be handled.
      integer,               parameter :: gp = grid_p
        !! Alias of the grid precision.
      type(meshDisType),  target, save :: meshDistr(maxDistr)
        !! Contains information of the several data distributions.
      type(meshCommType), target, save ::
     &                   meshCommu((maxDistr*(maxDistr-1))/2)
        !! Contains all the communications to move among  the several data
        !! distributions.
      type(meshCommType), target, save :: exteCommu(maxDistr,3)
        !! Contains all the needed communications to compute the extencil.

      real(grid_p),            pointer :: tSBuff(:)
        !! Buffer to send information.
      real(grid_p),            pointer :: tRBuff(:)
        !! Buffer to receive information.
      integer                          :: tBuffSize = 0
        !! Size of the buffer used in snd/recv.
      MPI_REQUEST_TYPE,        pointer :: send_rq(:)
        !! Send requests, a handle generated by MPI to identify a communication.
      MPI_REQUEST_TYPE,        pointer :: recv_rq(:)
        !! Recv requests, a handle generated by MPI to identify a communication.

#ifdef MPI
      type wloadBOX
        !! Internal type for workload-weighted recursive boxes exclusive to
        !! the current process.
        integer                 :: box(2,3)
          !! Box dimensions: box(1,:) has the starting coordinates and
          !! box(2,:) has the final coordinates for current process.
        integer                 :: rbox(2,3)
          !! Box dimensions for the "parent" box.
        integer(i8b)            :: total
          !! Total number of mesh points in the current box.
        integer(i8b),   pointer :: x(:)
          !! Total workload of mesh points in the X direction. (size: nmx)
        integer(i8b),   pointer :: y(:)
          !! Total workload of mesh points in the Y direction. (size: nmy)
        integer(i8b),   pointer :: z(:)
          !! Total workload of mesh points in the Z direction. (size: nmz)
        type(wloadBOX), pointer :: sb1 => null()
          !! Sub-box n°1
        type(wloadBOX), pointer :: sb2 => null()
          !! Sub-box n°2
        integer,        pointer :: wload(:)
          !! Workload for each point on the mesh (size nmx*nmy*nmz)
      contains
        procedure :: INIT  => InitWload
        procedure :: FREE  => FreeWload
        procedure :: CUT   => CutWload
        procedure :: SPLIT => SplitWload
        procedure :: getTotal
        procedure :: getWload
      end type wloadBOX

      type recBox
        !! Internal type for recursive boxes that comprise the whole domain.
        integer               :: box(2,3)
          !! Box dimensions: box(1,:) has the starting coordinates and
          !! box(2,:) has the final coordinates for current process.
        integer               :: np
          !! Total number of processes/nodes.
        integer(i8b)          :: gwload
          !! Total workload for this box.
        type(recBox), pointer :: sb1 => null()
          !! Sub-box n°1
        type(recBox), pointer :: sb2 => null()
          !! Sub-box n°2
      contains
        procedure :: INIT => InitRecBox
        procedure :: FREE => FreeRecBox
        procedure :: recSplit
        procedure :: ExtractBoxes
        procedure :: PrintBoxes
      end type recBox
#endif

      CONTAINS


      !< Initialization of the meshDistr (global)
      !<
      !< Ensures that we are correctly nullifying contained data
      subroutine init_globalmeshDistr()
         implicit none

         integer :: i, j
         logical, save :: skip = .false.

         if ( skip ) return
         skip = .true.
         do i= 1, maxDistr
           nullify(meshDistr(i)%box)
           nullify(meshDistr(i)%indexp)
           nullify(meshDistr(i)%idop)
           nullify(meshDistr(i)%xdop)
           nullify(meshDistr(i)%ipa)
         enddo
         do i= 1, (maxDistr*(maxDistr-1))/2
           nullify(meshCommu(i)%src)
           nullify(meshCommu(i)%dst)
         enddo
         do i= 1, maxDistr
           do j= 1, 3
             nullify(exteCommu(i,j)%src)
             nullify(exteCommu(i,j)%dst)
           enddo
         enddo

         nullify( tSBuff, tRBuff, send_rq, recv_rq )
         tBuffSize = 0
      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Mesh distribution creations and initializations.            [[MINI]] !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initMeshDistr_wload( iDistr, oDistr, nm, nsm, wload )
        !! Computes a new data distribution and the communications needed to
        !! move data from/to the current distribution to the existing ones.
        !!
        !! The limits of the new distributions are stored in this module:
        !!   * meshDistr(oDistr)
        !!   * meshCommu(((oDistr-2)*(oDistr-1))/2+1:(oDistr-1)*oDistr/2)
        !!
        !! If this is the first distribution, we split the mesh uniformly among,
        !! the several processes (we only split it in dimensions Y and Z). For
        !! the other data distributions we should split the vector wload.
        !!
        !! * The  subroutine splitwload will return the limits of the new data
        !!   distribution.
        !! * The subroutine compMeshComm will return the communications needed
        !!   to move data from/to the current distribution to/from the previous
        !!   ones.

      implicit none
      integer, optional,  intent(in) :: iDistr
        !! Distribution index of the input vector.
      integer,            intent(in) :: oDistr
        !! Index for the new data distribution.
      integer,            intent(in) :: nm(3)
        !! Number of Mesh divisions of each cell vector.
      integer,            intent(in) :: nsm
        !! Number of mesh-sub-divisions for this particular distribution.
      integer, optional,  intent(in) :: wload(:)
        !! Weights of every point of the mesh using the input distribution.

      character(len=*), parameter :: myName = moduName//'initMeshDistr '
      character(len=*), parameter :: errMsg = myName//'ERROR: '
      integer                     :: ii, jj, PY, PZ, PP, ProcessorZ,
     &                               blocY, blocZ, nremY, nremZ,
     &                               iniY, iniZ, dimY, dimZ
      type(meshDisType),  pointer :: distr
!------------------------------------------------------------------------- BEGIN
      call timer( 'INITMESH', 1 )

      ! Reset data if necessary
      call init_globalmeshDistr()

      ! Check the number of mesh distribution
      if (oDistr.gt.maxDistr)
     &  call die( errMsg // 'oDistr.gt.maxDistr' )

      distr => meshDistr(oDistr)
      ! Allocate memory for the current distribution
      call re_alloc( distr%box, 1, 2, 1, 3, 1, Nodes,
     &               'distr%box', moduName )

      ! Store distribution mesh-sub-divisions.
      distr%nsm = nsm

      ! The first distribution should always be the uniform distribution.
      if ( oDistr == UNIFORM ) then
        ProcessorZ = Nodes/ProcessorY
        blocY = (nm(2)/ProcessorY)
        blocZ = (nm(3)/ProcessorZ)
        nremY = nm(2) - blocY*ProcessorY
        nremZ = nm(3) - blocZ*ProcessorZ

        PP   = 1
        iniZ = 1
        do PZ = 1, ProcessorZ
          dimZ = blocZ
          if (PZ.LE.nremZ) dimZ = dimZ + 1
          iniY = 1
          do PY = 1, ProcessorY
            dimY = blocY
            if (PY.LE.nremY) dimY = dimY + 1

            distr%box(1,1,PP) = 1
            distr%box(2,1,PP) = nm(1)
            distr%box(1,2,PP) = iniY
            distr%box(2,2,PP) = iniY + dimY - 1
            distr%box(1,3,PP) = iniZ
            distr%box(2,3,PP) = iniZ + dimZ - 1

            iniY = iniY + dimY
            PP   = PP + 1
          enddo
          iniZ = iniZ + dimZ
        enddo
      else
        ! In order to compute the other data distributions, we should split
        ! the vector "wload" among the several processes
#ifdef MPI
        if (.NOT. present(iDistr) .OR.
     &      .NOT. present(wload) ) then
          call die( errMsg // 'Wrong parameters' )
        endif

        call buildDataDistr( oDistr, nm, wload )
#endif
      endif
      if (Node == 0) then
        write(6,"(a,'[',i0,']: sub = ',i0)")
     &      "New grid distribution ", oDistr, distr%nsm
#ifdef SIESTA__DEV_INFO
         do PP= 1, Nodes
            write(6,"(i12,3x,3(i5,a1,i5))") PP,
     $       (distr%box(1,jj,PP), ":", distr%box(2,jj,PP), jj=1,3)
         enddo
#endif
      endif

      call timer( 'INITMESH', 2 )
      end subroutine initMeshDistr_wload

      subroutine initMeshDistr_explicit(oDistr, start_nml, end_nml, nsm)
        !! Creates a new data distribution and the communications needed to
        !! move data from/to the current distribution to the existing ones.
        !! Instead of using the workload, the distribution is explicitly set
        !! by start_nml and end_nml, the starting and ending local mesh
        !! subdivisions for this node.
        !!
        !! initMeshDistr_wload MUST be called prior to this routine, since we
        !! rely on existing distributions to do the new one.
        !!
        !! The limits of the new distributions are stored in this module:
        !!   * meshDistr(oDistr)
        !!   * meshCommu(((oDistr-2)*(oDistr-1))/2+1:(oDistr-1)*oDistr/2)

#ifdef MPI
      use mpi_siesta
#endif
      implicit none
      integer, intent(in) :: oDistr
        !! Index for the new data distribution.
      integer, intent(in) :: start_nml(3)
        !! Start of local Mesh divisions for this node.
      integer, intent(in) :: end_nml(3)
        !! End of local Mesh divisions for this node.
      integer, intent(in) :: nsm
        !! Number of mesh-sub-divisions for this particular distribution.

      character(len=*), parameter :: myName =
     &    moduName//'initMeshDistr_explicit '
      character(len=*), parameter :: errMsg = myName//'ERROR: '
      integer                     :: ii, jj
#ifdef MPI
      integer :: ierr, se_nml(2,3)
#endif

      type(meshDisType),  pointer :: distr

      call timer( 'INITMESH_EX', 1 )

      ! Reset data if necessary
      call init_globalmeshDistr()

      ! Check the number of mesh distribution
      if (oDistr.gt.maxDistr)
     &  call die( errMsg // 'oDistr.gt.maxDistr' )
      distr => meshDistr(oDistr)

      ! Allocate memory for the current distribution
      call re_alloc( distr%box, 1, 2, 1, 3, 1, Nodes,
     &    'distr%box', moduName )

      ! Store distribution mesh-sub-divisions
      distr%nsm = nsm

#ifdef MPI
      ! Turn start and end nml into a proper box, and gather it for all nodes.
      se_nml(1,:) = start_nml(:)
      se_nml(2,:) = end_nml(:)

      call MPI_AllGather(se_nml(1,1), 6, MPI_Integer,
     &    distr%box(1,1,1), 6, MPI_Integer, MPI_Comm_World, ierr)

      ! Precompute the communications needed to move data between the new
      ! data distribution and the previous ones.
      jj = ((oDistr-2)*(oDistr-1))/2 + 1
      do ii=1, oDistr-1
        call compMeshComm( meshDistr(ii), distr, meshCommu(jj) )
        jj = jj + 1
      end do
#else
      ! Turn start and end nml into a proper box.
      distr%box(1,:,1) = start_nml(:)
      distr%box(2,:,1) = end_nml(:)
#endif

      if ( Node == 0 ) then
        write(6,"(a,'[',i0,']: sub = ',i0)")
     &      "New grid distribution ", oDistr, distr%nsm
#ifdef SIESTA__DEV_INFO
        do PP= 1, Nodes
          write(6,"(i12,3x,3(i5,a1,i5))") PP,
     $        (distr%box(1,jj,PP), ":", distr%box(2,jj,PP), jj=1,3)
        end do
#endif
      end if

      call timer( 'INITMESH_EX', 2 )
      end subroutine initMeshDistr_explicit

      subroutine allocASynBuffer( ndistr )
        !! Allocate memory buffers for asynchronous communications.
        !!
        !! The output values are stored in the current module:
        !!  * tBuff1(:):  Buffer for distribution 1
        !!  * tBuff2(:):  Buffer for other distributions

      implicit none
      integer, intent(in) :: ndistr
         !! Total number of distributions

      integer          :: bsize, maxncom
      integer, pointer :: box(:,:)
      integer          :: ii, nsp

      tBuffSize = 0
      do ii= 1, ndistr
        box => meshDistr(ii)%box(:,:,node+1)
        bsize = ( box(2,1) - box(1,1) + 1)*
     &          ( box(2,2) - box(1,2) + 1)*
     &          ( box(2,3) - box(1,3) + 1)
        nsp = meshDistr(ii)%nsm ** 3
        tBuffSize = max(tBuffSize,bsize*nsp)
      enddo

      maxncom   = 0
      do ii= 1, ndistr*(ndistr-1)
        maxncom   = max( maxncom, meshCommu(ii)%ncom )
      enddo

      ! Allocate memory for asynchronous communications
      call re_alloc( tSBuff, 1, tBuffSize, 'tSBuff', moduName )
      call re_alloc( tRBuff, 1, tBuffSize, 'tRBuff', moduName )
      allocate( send_rq(1:maxncom) )
      allocate( recv_rq(1:maxncom) )
      end subroutine allocASynBuffer


      subroutine allocExtMeshDistr( iDistr, nep, mop )
        !! Allocates the arrays for the current mesh distribution, and
        !! associates the pointers in the mesh module to this same distribution.
      use mesh, only: indexp, idop, xdop

      implicit none
      integer, intent(in) :: iDistr
        !! Index for the current mesh distribution.
      integer, intent(in) :: nep
        !! Number of extended box points.
      integer, intent(in) :: mop
        !! Number of points where orbitals are non-zero.

      type(meshDisType),  pointer :: distr

      distr => meshDistr(iDistr)
      call re_alloc( distr%indexp, 1, nep, 'distr%indexp', moduName )
      call re_alloc( distr%idop, 1, mop, 'distr%idop', moduName )
      call re_alloc( distr%xdop, 1, 3, 1, mop, 'distr%xdop', moduName )

      indexp => distr%indexp
      idop   => distr%idop
      xdop   => distr%xdop
      end subroutine allocExtMeshDistr

      subroutine allocIpaDistr( iDistr, na )
        !! Allocates the ipa array (indicating the mesh cell for a given atom)
        !! for current distribution. Associates the ipa pointer in mesh to the
        !! current distribution.

      use mesh, only: ipa
      implicit none
      integer,         intent(in) :: iDistr
        !! Index for current mesh distribution.
      integer,         intent(in) :: na
        !! Total number of atoms.

      type(meshDisType),  pointer :: distr

      distr => meshDistr(iDistr)
      call re_alloc( distr%ipa, 1, na, 'distr%ipa', moduName )
      ipa => meshDistr(iDistr)%ipa

      end subroutine allocIpaDistr

      subroutine setMeshDistr( iDistr, nml, nmpl, ntml, ntpl )
        !! Fixes the new data limits and dimensions of the mesh to those of
        !! the data distribution iDistr. Essentially this re-asocciates the
        !! pointers in the mesh module to those of the current data
        !! distribution, and resets the value for the number of points and
        !! divisions.
      use mesh, only: meshLim, indexp, ipa, idop, xdop

      implicit none
      integer,  intent(in) :: iDistr
        !! Index of the current distribution.
      integer, intent(out) :: nml(3)
        !! Local number of mesh divisions in each cell direction.
      integer, intent(out) :: nmpl
        !! Local number of mesh divisions (nmlx*nmly*nmlz).
      integer, intent(out) :: ntml(3)
        !! Local number of mesh points in each cell direction.
      integer, intent(out) :: ntpl
        !! Local number of mesh points (ntmlx*ntmly*ntmlz).

      type(meshDisType),  pointer :: distr

      distr => meshDistr(iDistr)

      meshLim = distr%box(1:2,1:3,node+1)
      nml(1) = (MeshLim(2,1)-MeshLim(1,1)) + 1
      nml(2) = (MeshLim(2,2)-MeshLim(1,2)) + 1
      nml(3) = (MeshLim(2,3)-MeshLim(1,3)) + 1
      nmpl   = nml(1)*nml(2)*nml(3)
      ntml   = nml*distr%nsm
      ntpl   = nmpl*distr%nsm**3

      indexp => distr%indexp
      idop   => distr%idop
      xdop   => distr%xdop
      ipa    => distr%ipa
      end subroutine setMeshDistr

      subroutine resetMeshDistr( iDistr )
        !! Reset the data of the distribution iDistr, deallocating its arrays.
        !! If no index is passed, it will reset ALL distributions.

      implicit none
      integer, optional, intent(in) :: iDistr
        !! Index of the distribution to be reset.

      integer                        :: idis, ini, fin, icom
      type(meshDisType),  pointer    :: distr
      type(meshCommType), pointer    :: mcomm

      if (present(iDistr)) then
        ini = iDistr
        fin = iDistr
      else
        ini = 1
        fin = maxDistr
      endif

      do idis= ini, fin
        distr => meshDistr(idis)

        distr%nMesh = 0
        distr%nsm = 0

        if (associated(distr%box)) then
          call de_alloc( distr%box, 'distr%box', moduName )
        endif

        if (associated(distr%indexp)) then
          call de_alloc( distr%indexp, 'distr%indexp', moduName )
        endif

        if (associated(distr%idop)) then
          call de_alloc( distr%idop, 'distr%idop', moduName )
        endif

        if (associated(distr%xdop)) then
          call de_alloc( distr%xdop, 'distr%xdop', moduName )
        endif

        if (associated(distr%ipa)) then
          call de_alloc( distr%ipa, 'distr%ipa', moduName )
        endif

        do icom=1, 3
          mcomm => exteCommu(idis,icom)
          if (associated(mcomm%src)) then
            call de_alloc( mcomm%src, 'mcomm%src', moduName )
          endif
          if (associated(mcomm%dst)) then
            call de_alloc( mcomm%dst, 'mcomm%dst', moduName )
          endif
          mcomm%ncom = 0
        enddo

        do icom= ((idis-2)*(idis-1))/2 + 1, ((idis-1)*idis)/2
          mcomm => meshCommu(icom)
          if (associated(mcomm%src)) then
            call de_alloc( mcomm%src, 'mcomm%src', moduName )
          endif
          if (associated(mcomm%dst)) then
            call de_alloc( mcomm%dst, 'mcomm%dst', moduName )
          endif
          mcomm%ncom = 0
        enddo
      enddo

      if (associated(tSBuff)) then
        call de_alloc( tSBuff, 'tSBuff', moduName )
      endif
      if (associated(tRBuff)) then
        call de_alloc( tRBuff, 'tRBuff', moduName )
      endif
      if (associated(send_rq)) then
        deallocate( send_rq )
      endif
      if (associated(recv_rq)) then
        deallocate( recv_rq )
      endif
      end subroutine resetMeshDistr


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Routines that handle the re-distribution of integer and     [[MDST]] !!
      !! real-type arrays between different parallel distributions.           !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine init_js( nsm, nsize, js )
        !! Initializes the JS array, which contains the distance between points
        !! of the same cluster.

      implicit none
      integer :: nsm
        !! Number of mesh subdivisions.
      integer :: nsize(3)
        !! Number of mesh divisions in each direction.
      integer :: js(*)
        !! Distance between points of the same cluster
      integer :: NN, I1, I2, I3, N1, N2, N3

      ! From clustered to sequential
      NN = 1
      I3 = 0
      DO N3=0, NSM-1
        I2 = 0
        DO N2= 0, NSM-1
          I1 = I2 + I3
          DO N1= 0, NSM-1
            JS(NN) = I1
            NN     = NN + 1
            I1     = I1 + 1
          ENDDO
          I2 = I2 + nsize(1)
        ENDDO
        I3 = I3 + nsize(1)*nsize(2)
      ENDDO
      end subroutine init_js


      subroutine localCopy( itr, insm, onsm, nsp, sbox, nsrc, dbox,
     &                      ndst, js, fsrc, fdst )
        !! Makes a copy of the data contained in an array with a given
        !! distribution into an array with another distribution, but with the
        !! elements that the node shares in both distributions (i.e., that
        !! the elements copied are strictly local for this node).
        !!
        !! Note: When changing between sequential and clustered, the number
        !! of mesh subdivisions (nsm) must be the same.
      use mesh, only : nsm

      implicit none
      integer     , intent(in)  :: itr
        !! Transform from cluster/sequential to sequential/cluster, or keep
        !! current distribution.
      integer     , intent(in)  :: insm
        !! Number of mesh subdivisions for source distribution.
      integer     , intent(in)  :: onsm
        !! Number of mesh subdivisions for output distribution.
      integer     , intent(in)  :: nsp
        !! Number of sub-points.
      integer     , intent(in)  :: sbox(2,3)
        !! Box dimensions for the source distribution.
      integer     , intent(in)  :: nsrc(3)
        !! Number of mesh divisions in each direction for source distribution.
      integer     , intent(in)  :: dbox(2,3)
        !! Box dimensions for the output distribution.
      integer     , intent(in)  :: ndst(3)
        !! Number of mesh divisions in each direction for output distribution.
      integer     , intent(in)  :: js(*)
        !! Array containing the distance between neighbouring points in the
        !! same cluster.
      real(grid_p), intent(in)  :: fsrc(*)
        !! Data array in the source distribution.
      real(grid_p), intent(out) :: fdst(*)
        !! Data array in the output distribution.

      integer :: I1, I2, I3, J1, J2, J3, K1, K2, K3
      integer :: NN
      integer :: lbox(2,3), lsize(3)
      logical :: inters

      if (itr==1) then
        ! From clustered to sequential
        call boxIntersection( Sbox, Dbox, Lbox, inters )

        J3 = (Lbox(1,3) - Sbox(1,3))*NSRC(1)*NSRC(2)*nsm
        K3 = (Lbox(1,3) - Dbox(1,3))*NDST(1)*NDST(2)*nsm
        do I3 = Lbox(1,3), Lbox(2,3)
          J2 = (Lbox(1,2) - Sbox(1,2))*NSRC(1)*nsm*nsm
          K2 = (Lbox(1,2) - Dbox(1,2))*NDST(1)*nsm
          do I2 = Lbox(1,2), Lbox(2,2)
            J1 = (Lbox(1,1) - Sbox(1,1))*NSP + 1 + J2 + J3
            K1 = (Lbox(1,1) - Dbox(1,1))*nsm + 1 + K2 + K3
            do I1 = Lbox(1,1), Lbox(2,1)
              DO NN= 1, NSP
                fdst(K1+JS(NN)) = fsrc(J1)
                J1 = J1 + 1
              ENDDO
              K1 = K1 + nsm
            enddo
            J2 = J2 + NSRC(1)*nsm*nsm
            K2 = K2 + NDST(1)*nsm
          enddo
          J3 = J3 + NSRC(1)*NSRC(2)*nsm
          K3 = K3 + NDST(1)*NDST(2)*nsm
        enddo
      else if (itr==-1) then
        ! From sequential to clustered
        call boxIntersection( Sbox, Dbox, Lbox, inters )

        J3 = (Lbox(1,3) - Sbox(1,3))*NSRC(1)*NSRC(2)*nsm
        K3 = (Lbox(1,3) - Dbox(1,3))*NDST(1)*NDST(2)*nsm
        do I3 = Lbox(1,3), Lbox(2,3)
          J2 = (Lbox(1,2) - Sbox(1,2))*NSRC(1)*nsm
          K2 = (Lbox(1,2) - Dbox(1,2))*NDST(1)*nsm*nsm
          do I2 = Lbox(1,2), Lbox(2,2)
            J1 = (Lbox(1,1) - Sbox(1,1))*nsm + 1 + J2 + J3
            K1 = (Lbox(1,1) - Dbox(1,1))*NSP + 1 + K2 + K3
            do I1 = Lbox(1,1), Lbox(2,1)
              DO NN= 1, NSP
                fdst(K1) = fsrc(J1+JS(NN))
                K1 = K1 + 1
              ENDDO
              J1 = J1 + nsm
            enddo
            J2 = J2 + NSRC(1)*nsm
            K2 = K2 + NDST(1)*nsm*nsm
          enddo
          J3 = J3 + NSRC(1)*NSRC(2)*nsm
          K3 = K3 + NDST(1)*NDST(2)*nsm
        enddo
      else
        ! From sequential to sequential or from clustered to clustered
        call boxIntersection( Sbox, Dbox, Lbox, inters, insm, onsm )

        Lsize = Lbox(2,:) - Lbox(1,:) + 1
        J3 = (Lbox(1,3) - (Sbox(1,3)-1)*insm-1)*NSRC(1)*NSRC(2)
        K3 = (Lbox(1,3) - (Dbox(1,3)-1)*onsm-1)*NDST(1)*NDST(2)
        do I3 = 1, Lsize(3)
          J2 = (Lbox(1,2) - (Sbox(1,2)-1)*insm-1)*NSRC(1)
          K2 = (Lbox(1,2) - (Dbox(1,2)-1)*onsm-1)*NDST(1)
          do I2 = 1, Lsize(2)
            J1 = Lbox(1,1) - (Sbox(1,1)-1)*insm-1 + 1 + J2 + J3
            K1 = Lbox(1,1) - (Dbox(1,1)-1)*onsm-1 + 1 + K2 + K3
            do I1 = 1, Lsize(1)
              fdst(K1) = fsrc(J1)
              K1 = K1 + 1
              J1 = J1 + 1
            enddo
            J2 = J2 + NSRC(1)
            K2 = K2 + NDST(1)
          enddo
          J3 = J3 + NSRC(1)*NSRC(2)
          K3 = K3 + NDST(1)*NDST(2)
        enddo
      endif
      end subroutine localCopy


      function packBuffer( itr, insm, onsm, nsp, sbox, nsrc, dbox,
     &                     js, fsrc, tbuf ) result(n)
        !! Generates the transmision buffer of the data array fsrc, in order
        !! to send the data to another node that needs it in a new distribution.
        !! As a function, this also returns the size of the buffer.
        !!
        !! Note: When changing between sequential and clustered, the number
        !! of mesh subdivisions (nsm) must be the same.
      use mesh, only : nsm

      implicit none
      integer     , intent(in)  :: itr
        !! Transform from cluster/sequential to sequential/cluster, or keep
        !! current distribution.
      integer     , intent(in)  :: insm
        !! Number of mesh subdivisions for source distribution.
      integer     , intent(in)  :: onsm
        !! Number of mesh subdivisions for output distribution.
      integer     , intent(in)  :: nsp
        !! Number of sub-points.
      integer     , intent(in)  :: sbox(2,3)
        !! Box dimensions for the source distribution.
      integer     , intent(in)  :: nsrc(3)
        !! Number of mesh divisions in each direction for source distribution.
      integer     , intent(in)  :: dbox(2,3)
        !! Box dimensions for the output distribution.
      integer     , intent(in)  :: js(*)
        !! Array containing the distance between neighbouring points in the
        !! same cluster.
      real(grid_p), intent(in)  :: fsrc(*)
        !! Data array in the source distribution.
      real(grid_p), intent(out) :: tbuf(*)
        !! Data buffer to be transferred.
      integer                   :: n
        !! Total size of the buffer.

      ! Local variables
      integer                   :: I1, I2, I3, J1, J2, J3, K1
      integer                   :: NN
      integer                   :: lbox(2,3), lsize(3)
      logical                   :: inters

      if (itr==1) then
        call boxIntersection( Sbox, Dbox, Lbox, inters )
        ! From clustered to sequential
        J3 = (Lbox(1,3) - Sbox(1,3))*NSRC(1)*NSRC(2)*nsm
        K1 = 1
        do I3 = Lbox(1,3), Lbox(2,3)
          J2 = (Lbox(1,2) - Sbox(1,2))*NSRC(1)*nsm*nsm
          do I2 = Lbox(1,2), Lbox(2,2)
            J1 = (Lbox(1,1) - Sbox(1,1))*NSP + 1 + J2 + J3
            do I1 = Lbox(1,1), Lbox(2,1)
              DO NN= 1, NSP
                TBUF(K1) = fsrc(J1)
                J1 = J1 + 1
                K1 = K1 + 1
              ENDDO
            enddo
            J2 = J2 + NSRC(1)*nsm*nsm
          enddo
          J3 = J3 + NSRC(1)*NSRC(2)*nsm
        enddo
      else if (itr==-1) then
        call boxIntersection( Sbox, Dbox, Lbox, inters )
        ! From sequential to clustered
        J3 = (Lbox(1,3) - Sbox(1,3))*NSRC(1)*NSRC(2)*nsm
        K1 = 1
        do I3 = Lbox(1,3), Lbox(2,3)
          J2 = (Lbox(1,2) - Sbox(1,2))*NSRC(1)*nsm
          do I2 = Lbox(1,2), Lbox(2,2)
            J1 = (Lbox(1,1) - Sbox(1,1))*nsm + 1 + J2 + J3
            do I1 = Lbox(1,1), Lbox(2,1)
              DO NN= 1, NSP
                TBUF(K1) = fsrc(J1+JS(NN))
                K1 = K1 + 1
              ENDDO
              J1 = J1 + nsm
            enddo
            J2 = J2 + NSRC(1)*nsm
          enddo
          J3 = J3 + NSRC(1)*NSRC(2)*nsm
        enddo
      else
        ! From sequential to sequential or from clustered to clustered
        call boxIntersection( Sbox, Dbox, Lbox, inters, insm, onsm )
        Lsize = Lbox(2,:) - Lbox(1,:) + 1
        J3 = (Lbox(1,3) - (Sbox(1,3)-1)*insm-1)*NSRC(1)*NSRC(2)
        K1 = 1
        do I3 = 1, Lsize(3)
          J2 = (Lbox(1,2) - (Sbox(1,2)-1)*insm-1)*NSRC(1)
          do I2 = 1, Lsize(2)
            J1 = Lbox(1,1) - (Sbox(1,1)-1)*insm-1 + 1 + J2 + J3
            do I1 = 1, Lsize(1)
              TBUF(K1) = fsrc(J1)
              K1 = K1 + 1
              J1 = J1 + 1
            enddo
            J2 = J2 + NSRC(1)
          enddo
          J3 = J3 + NSRC(1)*NSRC(2)
        enddo
      endif
      n = K1 - 1
      end function packBuffer


      subroutine unpackBuffer( itr, insm, onsm, nsp, sbox, dbox, ndst,
     &                         js, tbuf, fdst )
        !! Unpacks the transmision buffer sent by another node into the new
        !! data distribution fdst.
        !!
        !! Note: When changing between sequential and clustered, the number
        !! of mesh subdivisions (nsm) must be the same.
      use mesh, only : nsm

      implicit none
      ! Input/output variables

      integer     , intent(in)  :: itr
        !! Transform from cluster/sequential to sequential/cluster, or keep
        !! current distribution.
      integer     , intent(in)  :: insm
        !! Number of mesh subdivisions for source distribution.
      integer     , intent(in)  :: onsm
        !! Number of mesh subdivisions for output distribution.
      integer     , intent(in)  :: nsp
        !! Number of sub-points.
      integer     , intent(in)  :: sbox(2,3)
        !! Box dimensions for the source distribution.
      integer     , intent(in)  :: dbox(2,3)
        !! Box dimensions for the output distribution.
      integer     , intent(in)  :: ndst(3)
        !! Number of mesh divisions in each direction for output distribution.
      integer     , intent(in)  :: js(*)
        !! Array containing the distance between neighbouring points in the
        !! same cluster.
      real(grid_p), intent(in)  :: tbuf(*)
        !! Data buffer received.
      real(grid_p), intent(out) :: fdst(*)
        !! Data array in the output distribution.

      ! Local variables
      integer                   :: I1, I2, I3, J1, K1, K2, K3
      integer                   :: NN
      integer                   :: lbox(2,3), lsize(3)
      logical                   :: inters

      if (itr==1) then
        ! From clustered to sequential
        call boxIntersection( Sbox, Dbox, Lbox, inters )
        J1 = 1
        K3 = (Lbox(1,3) - Dbox(1,3))*NDST(1)*NDST(2)*nsm
        do I3 = Lbox(1,3), Lbox(2,3)
          K2 = (Lbox(1,2) - Dbox(1,2))*NDST(1)*nsm
          do I2 = Lbox(1,2), Lbox(2,2)
            K1 = (Lbox(1,1) - Dbox(1,1))*nsm + 1 + K2 + K3
            do I1 = Lbox(1,1), Lbox(2,1)
              DO NN= 1, NSP
                fdst(K1+JS(NN)) = TBUF(J1)
                J1 = J1 + 1
              ENDDO
              K1 = K1 + nsm
            enddo
            K2 = K2 + NDST(1)*nsm
          enddo
          K3 = K3 + NDST(1)*NDST(2)*nsm
        enddo
      else if (itr==-1) then
        ! From sequential to clustered
        call boxIntersection( Sbox, Dbox, Lbox, inters )
        J1 = 1
        K3 = (Lbox(1,3) - Dbox(1,3))*NDST(1)*NDST(2)*nsm
        do I3 = Lbox(1,3), Lbox(2,3)
          K2 = (Lbox(1,2) - Dbox(1,2))*NDST(1)*nsm*nsm
          do I2 = Lbox(1,2), Lbox(2,2)
            K1 = (Lbox(1,1) - Dbox(1,1))*NSP + 1 + K2 + K3
            do I1 = Lbox(1,1), Lbox(2,1)
              DO NN= 1, NSP
                fdst(K1) = TBUF(J1)
                K1 = K1 + 1
                J1 = J1 + 1
              ENDDO
            enddo
            K2 = K2 + NDST(1)*nsm*nsm
          enddo
          K3 = K3 + NDST(1)*NDST(2)*nsm
        enddo
      else
        ! From sequential to sequential or from clustered to clustered
        call boxIntersection( Sbox, Dbox, Lbox, inters, insm, onsm )
        Lsize = Lbox(2,:) - Lbox(1,:) + 1
        J1 = 1
        K3 = (Lbox(1,3) - (Dbox(1,3)-1)*onsm-1)*NDST(1)*NDST(2)
        do I3 = 1, Lsize(3)
          K2 = (Lbox(1,2) - (Dbox(1,2)-1)*onsm-1)*NDST(1)
          do I2 = 1, Lsize(2)
            K1 = Lbox(1,1) - (Dbox(1,1)-1)*onsm-1 + 1 + K2 + K3
            do I1 = 1, Lsize(1)
              fdst(K1) = TBUF(J1)
              K1 = K1 + 1
              J1 = J1 + 1
            enddo
            K2 = K2 + NDST(1)
          enddo
          K3 = K3 + NDST(1)*NDST(2)
        enddo
      endif
      end subroutine unpackBuffer

      subroutine distMeshData_rea( iDistr, fsrc, oDistr, fdst, itr )
        !! Move data from vector fsrc, that uses distribution iDistr, to vector
        !! fdst, that uses distribution oDistr. It also, re-orders a clustered
        !! data array into a sequential one and viceversa.
        !! If this is a sequential execution, it only reorders the data.
        !!
        !! We first check the communications that this process should do to move
        !! data from iDistr to odistr. We have 3 kind of communications (send,
        !! receive and keep on the same node). We also have 3 kind of
        !! reorderings (clustered to sequential, sequential to clustered and
        !! keep the same ordering). For the sequential code we call "reord".
        !!
        !! AG: Note that the integer version does NOT have the exact
        !! functionality of the real version. In particular, the integer version
        !! has no provision for a "serial fallback", and so this case has been
        !! trapped.

      use mesh,    only : nsm, nmeshg
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
      integer     , intent(in)  :: iDistr
        !! Distribution index of the input vector.
      integer     , intent(in)  :: oDistr
        !! Distribution index of the output vector.
      integer     , intent(in)  :: itr
        !! Exchange between sequential/clustered, or keep the current
        !! distribution.
      real(grid_p), intent(in)  :: fsrc(*)
        !! Input data vector.
      real(grid_p), intent(out) :: fdst(*)
        !! Output data vector.

      integer                     :: nm(3), nsize, ind, ncom,
     &                               nsp, me, nsrc(3), ndst(3), icom,
     &                               sb_index, csize, ierr,
     &                               nr, ns, index, insm, onsm
#ifdef MPI
      MPI_STATUS_TYPE             :: status
#endif
      integer,            pointer :: src(:), dst(:), sbox(:,:), js(:),
     &                               dbox(:,:), addr(:)
      type(meshDisType),  pointer :: idis, odis

      call timer( 'COMM_BSC', 1 )

      if (nodes == 1) then
        nm(1:3) = nmeshg(1:3) / nsm
        if (itr.gt.0) then
          ! Note that in reord the first argument is always clustered.
          call reord( fsrc, fdst, nm, nsm, TO_SEQUENTIAL )
        else if (itr .lt. 0) then
          call reord( fdst, fsrc, nm, nsm, TO_CLUSTER )
        else
          ! Copy source to destination. This will be executed only in serial
          ! mode,so we know that the size is the total number of (small) points,
          ! but maybe this information should be more explicit.
          nsize = product(nmeshg(1:3))
          fdst(1:nsize) = fsrc(1:nsize)
        endif
      else  ! nodes > 1
        ! The communications are stored in a triangular structure.
        if (iDistr.gt.oDistr) then
          ind  = ((iDistr-1)*(iDistr-2))/2 + oDistr
          ncom = meshCommu(ind)%ncom
          src => meshCommu(ind)%dst
          dst => meshCommu(ind)%src
        else
          ind = ((oDistr-1)*(oDistr-2))/2 + iDistr
          ncom = meshCommu(ind)%ncom
          src => meshCommu(ind)%src
          dst => meshCommu(ind)%dst
        endif

        idis => meshDistr(iDistr)
        insm = idis%nsm
        odis => meshDistr(oDistr)
        onsm = odis%nsm

        nsp = nsm*nsm*nsm
        me  = node + 1

        ! Get the boxes from send and destitnation distribution
        sbox => idis%box(:,:,me)
        nsrc = ( sbox(2,:) - sbox(1,:) + 1)*insm
        dbox => odis%box(:,:,me)
        ndst = ( dbox(2,:) - dbox(1,:) + 1)*onsm

        ! Initiallize JS array. JS contains the distance between
        ! points of the same cluster
        nullify( js )
        call re_alloc( js, 1, nsp, 'JS', moduName )
        if (itr==1) then
          if ( insm /= onsm )
     &        call die("insm /= onsm: should not be called: 1")
          call init_js( nsm, ndst, js )

        else if (itr==-1) then
          if ( insm /= onsm )
     &        call die("insm /= onsm: should not be called: -1")
          call init_js( nsm, nsrc, js )

        endif

        ! Allocate addr: It will store the position in the receive buffer
        nullify(addr)
        call re_alloc( addr, 1, ncom, 'addr', moduName )
        addr(1) = 1

        sb_index = 1
        nr       = 0
        ns       = 0
        do icom= 1, ncom
          ! For all the communications
          sbox => idis%box(:,:,src(icom))
          dbox => odis%box(:,:,dst(icom))

          if (src(icom).eq.ME) then
            if (dst(icom).eq.ME) then
              call localCopy( itr, insm, onsm, nsp, sbox, nsrc, dbox,
     &                        ndst, js, fsrc, fdst )
            else
              ns = ns + 1
              ! Copy data to local buffer before non blocking send
#ifdef MPI
              csize = packBuffer( itr, insm, onsm, nsp, sbox, nsrc,
     &                            dbox, js, fsrc, tSBuff(sb_index:) )
              call mpi_isend( tSBuff(sb_index), csize, MPI_grid_real,
     &                        dst(icom)-1, 0, MPI_COMM_WORLD,
     &                        send_rq(ns), ierr )
#endif
              sb_index = sb_index + csize
            endif
          else
            nr = nr + 1
            csize = nsp*
     &        ( min(sbox(2,1),dbox(2,1))- max(sbox(1,1),dbox(1,1)) + 1)*
     &        ( min(sbox(2,2),dbox(2,2))- max(sbox(1,2),dbox(1,2)) + 1)*
     &        ( min(sbox(2,3),dbox(2,3))- max(sbox(1,3),dbox(1,3)) + 1)
            ! Receive data into temporary buffer. We will process data later
#ifdef MPI
            call mpi_irecv( tRBuff(addr(nr)), csize, MPI_grid_real,
     &                        src(icom)-1, 0, MPI_COMM_WORLD,
     &                        recv_rq(nr), ierr )
#endif
            if (nr < ncom) addr(nr+1) = addr(nr) + csize
          endif
        enddo

#ifdef MPI
        dbox => odis%box(:,:,ME)
        do icom= 1, nr
          ! Check for any finished RECV
          call MPI_Waitany( nr, recv_rq, index, status, ierr )
          sbox => idis%box(:,:,MPI_STATUS_SOURCE(status)+1)
          ! Move data from temporary buffer to destination
          call unpackBuffer( itr, insm, onsm, nsp, sbox, dbox, ndst,
     &                       js, tRBuff(addr(index):), fdst )
        enddo
        ! Check that all non blocking send have finished
        call MPI_Waitall( ns, send_rq, MPI_STATUSES_IGNORE, ierr )
#endif

        call de_alloc( addr, 'addr', moduName )
        call de_alloc( js, 'js', moduName )
      endif
      call timer( 'COMM_BSC', 2 )
      end subroutine distMeshData_rea

      subroutine distMeshData_int( iDistr, fsrc, oDistr, fdst, itr )
        !! Move data from vector fsrc, that uses distribution iDistr, to vector
        !! fdst, that uses distribution oDistr. It also, re-orders a clustered
        !! data array into a sequential one and viceversa.
        !! If this is a sequential execution, it only reorders the data.
        !!
        !! We first check the communications that this process should do to move
        !! data from iDistr to odistr. We have 3 kind of communications (send,
        !! receive and keep on the same node). We also have 3 kind of
        !! reorderings (clustered to sequential, sequential to clustered and
        !! keep the same ordering). For the sequential code we call "reord".
        !!
        !! AG: Note that this integer version does NOT have the exact
        !! functionality of the real version. In particular, the integer version
        !! has no provision for a "serial fallback", and so this case has been
        !! trapped.
        !!
        !! This method is special given that it will never use mesh-subdivisions
        !! i.e. this is meant for distribution of numphi indices that correspond
        !! to orbitals located on each big mesh point.
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
      integer, intent(in)  :: iDistr
        !! Distribution index of the input vector.
      integer, intent(in)  :: oDistr
        !! Distribution index of the output vector.
      integer, intent(in)  :: itr
        !! Exchange between sequential/clustered, or keep the current
        !! distribution.
      integer, intent(in)  :: fsrc(*)
        !! Input data vector.
      integer, intent(out) :: fdst(*)
        !! Output data vector.

      character(len=*), parameter :: myName = moduName//'distMeshData '
      character(len=*), parameter :: errMsg = myName//'ERROR: '
      integer                     :: I1, I2, I3, J1, J2, J3, K1, K2, K3,
     &                               ind, ncom, icom, NSRC(3), NDST(3),
     &                               ME, MaxSize, Xsize, Ysize, Zsize,
     &                               Lbox(2,3), insm, onsm
      integer, pointer            :: src(:), dst(:), Sbox(:,:),
     &                               Dbox(:,:)
      type(meshDisType),  pointer :: idis, odis
      logical                     :: inters
      integer          ,  pointer :: TBUF(:)
#ifdef MPI
      integer                     :: MPIerror
      MPI_STATUS_TYPE             :: Status
#endif
!---------------------------------------------------------------------- BEGIN
      if (nodes == 1) then
         call die("Called _int version of distMeshData for n=1")
      else
        ! The communications are stored in a triangular structure.

        if (iDistr.gt.oDistr) then
          ind  = ((iDistr-1)*(iDistr-2))/2 + oDistr
          ncom = meshCommu(ind)%ncom
          src => meshCommu(ind)%dst
          dst => meshCommu(ind)%src
        else
          ind = ((oDistr-1)*(oDistr-2))/2 + iDistr
          ncom = meshCommu(ind)%ncom
          src => meshCommu(ind)%src
          dst => meshCommu(ind)%dst
        endif

        idis => meshDistr(iDistr)
        odis => meshDistr(oDistr)
        insm = idis%nsm
        onsm = idis%nsm

        ME  = Node + 1

        ! Compute the maximum size of the buffer needed to transfer data
        ! among the several processes
        maxSize = 0
        do icom= 1, ncom
          if (src(icom).ne.dst(icom)) then
            Sbox => idis%box(:,:,src(icom))
            Dbox => odis%box(:,:,dst(icom))
            call boxIntersection( Sbox, Dbox, Lbox, inters )
            Xsize = Lbox(2,1) - Lbox(1,1) + 1
            Ysize = Lbox(2,2) - Lbox(1,2) + 1
            Zsize = Lbox(2,3) - Lbox(1,3) + 1
            MaxSize = max(MaxSize,Xsize*Ysize*Zsize)
          endif
        enddo

        if ( MaxSize > 0 ) then
          nullify( TBUF )
          call re_alloc( TBUF, 1, MaxSize, 'TBUF', moduName )
        endif

        Sbox => idis%box(:,:,ME)
        NSRC(:) = Sbox(2,:) - Sbox(1,:) + 1
        Dbox => odis%box(:,:,ME)
        NDST(:) = Dbox(2,:) - Dbox(1,:) + 1

        if (itr == KEEP) then
          ! From sequential to sequential
          do icom= 1, ncom
            Sbox => idis%box(:,:,src(icom))
            Dbox => odis%box(:,:,dst(icom))
            call boxIntersection( Sbox, Dbox, Lbox, inters )
            Xsize = Lbox(2,1) - Lbox(1,1) + 1
            Ysize = Lbox(2,2) - Lbox(1,2) + 1
            Zsize = Lbox(2,3) - Lbox(1,3) + 1

            if (src(icom).eq.ME) then
              if (dst(icom).eq.ME) then
                ! SRC and DST are the current process
                J3 = (Lbox(1,3) - Sbox(1,3))*NSRC(1)*NSRC(2)
                K3 = (Lbox(1,3) - Dbox(1,3))*NDST(1)*NDST(2)
                do I3 = 1, Zsize
                  J2 = (Lbox(1,2) - Sbox(1,2))*NSRC(1)
                  K2 = (Lbox(1,2) - Dbox(1,2))*NDST(1)
                  do I2 = 1, Ysize
                    J1 = Lbox(1,1) - Sbox(1,1) + 1 + J2 + J3
                    K1 = Lbox(1,1) - Dbox(1,1) + 1 + K2 + K3
                    do I1 = 1, Xsize
                      fdst(K1) = fsrc(J1)
                      K1 = K1 + 1
                      J1 = J1 + 1
                    enddo
                    J2 = J2 + NSRC(1)
                    K2 = K2 + NDST(1)
                  enddo
                  J3 = J3 + NSRC(1)*NSRC(2)
                  K3 = K3 + NDST(1)*NDST(2)
                enddo
              else
                ! We should send data to process dst(icom)-1
                J3 = (Lbox(1,3) - Sbox(1,3))*NSRC(1)*NSRC(2)
                K1 = 1
                do I3 = 1, Zsize
                  J2 = (Lbox(1,2) - Sbox(1,2))*NSRC(1)
                  do I2 = 1, Ysize
                    J1 = Lbox(1,1) - Sbox(1,1) + 1 + J2 + J3
                    do I1 = 1, Xsize
                      TBUF(K1) = fsrc(J1)
                      K1 = K1 + 1
                      J1 = J1 + 1
                    enddo
                    J2 = J2 + NSRC(1)
                  enddo
                  J3 = J3 + NSRC(1)*NSRC(2)
                enddo
#ifdef MPI
                call MPI_Send( TBUF, Xsize*Ysize*Zsize,
     &                         MPI_Integer, dst(icom)-1, 1,
     &                         MPI_Comm_world, MPIerror )
#endif
              endif
            else
              ! We should receive data from process src(icom)-1
#ifdef MPI
              call mpi_recv( TBUF, Xsize*Ysize*Zsize,
     &                       MPI_Integer, src(icom)-1, 1,
     &                       MPI_Comm_world, Status, MPIerror )
#endif
              J1 = 1
              K3 = (Lbox(1,3) - Dbox(1,3))*NDST(1)*NDST(2)
              do I3 = 1, Zsize
                K2 = (Lbox(1,2) - Dbox(1,2))*NDST(1)
                do I2 = 1, Ysize
                  K1 = Lbox(1,1) - Dbox(1,1) + 1 + K2 + K3
                  do I1 = 1, Xsize
                    fdst(K1) = TBUF(J1)
                    K1 = K1 + 1
                    J1 = J1 + 1
                  enddo
                  K2 = K2 + NDST(1)
                enddo
                K3 = K3 + NDST(1)*NDST(2)
              enddo
            endif

          enddo
        else
          if (Node.eq.0) then
            write(*,*)'ERROR: Wrong parameter for function distMeshData'
          endif
          call die("stopping program")
        endif

        if (MaxSize.gt.0) then
          call de_alloc( TBUF, 'TBUF', moduName )
        endif
      endif
      end subroutine distMeshData_int

      subroutine boxIntersection(ibox1, ibox2, obox, inters, nsm1, nsm2)
        !! Checks if there is an intersection between 2 boxes and, if it exists,
        !! it returns the resulting intersected box.

      implicit none
      integer, intent(in)           :: ibox1(2,3)
        !! Input box 1.
      integer, intent(in)           :: ibox2(2,3)
        !! Input box 2.
      integer, intent(in), optional :: nsm1
        !! Number of mesh subdivisions for box1.
      integer, intent(in), optional :: nsm2
        !! Number of mesh subdivisions for box2.
      integer, intent(out)          :: obox(2,3)
        !! Intersection between ibox1 and ibox2.
      logical, intent(out)          :: inters
        !! Returns true if there was an intersection, and false if not.

      integer :: iaxis

      inters = .true.
      if ( present(nsm1) .and. present(nsm2) ) then
        do iaxis= 1, 3
          obox(1,iaxis) =
     &        max((ibox1(1,iaxis)-1)*nsm1+1,(ibox2(1,iaxis)-1)*nsm2+1)
          obox(2,iaxis) = min(ibox1(2,iaxis)*nsm1,ibox2(2,iaxis)*nsm2)
          if (obox(2,iaxis).lt.obox(1,iaxis)) inters = .false.
        enddo
      else
        do iaxis= 1, 3
          obox(1,iaxis) = max(ibox1(1,iaxis),ibox2(1,iaxis))
          obox(2,iaxis) = min(ibox1(2,iaxis),ibox2(2,iaxis))
          if (obox(2,iaxis).lt.obox(1,iaxis)) inters = .false.
        enddo
      endif

      end subroutine boxIntersection


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Mesh extencil related subroutines.                          [[MEXT]] !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! TO-DO: Clarify what "extencil" means.

#ifdef MPI
      subroutine initMeshExtencil( iDistr, nm )
        !! Compute the needed communications in order to send/receive the
        !! extencil (when the data is ordered in the distribution iDistr).
        !! The outputs are stored in this module under exteCommu(iDistr,1:3).
        !!
        !! For every dimension of the problem, search all the neightbours that
        !! we have. Given the current data distribution we compute the limits
        !! of our extencil and we check its intersection with all the other
        !! processes. Once we know all our neightbours we call subroutine
        !! scheduleComm in order to minimize the number of communications steps.
      use scheComm, only: scheduleComm, COMM_T

      implicit none
      integer,  intent(in) :: iDistr
        !! Distribution index to be used.
      integer,  intent(in) :: nm(3)
        !! Number of mesh divisions in each cell directions.

      integer                     :: Ubox(2,3), Lbox(2,3), Ibox(2,3),
     &                               ii, iaxis, ncom, Gcom, Lcom, P1, P2
      integer,            pointer :: src(:), dst(:), Dbox(:,:)
      type(meshDisType),  pointer :: idis
      type(meshCommType), pointer :: mcomm
      type(COMM_T)                :: comm
      logical                     :: inters

      idis => meshDistr(iDistr)

      do iaxis=1, 3
        ! One communication structure for every dimension
        mcomm => exteCommu(iDistr,iaxis)

        ! Count the number of communications needed to send/receive the extencil
        ncom = 0
        do P1= 1, Nodes
          ! Create the extencil boxes for both sides of the current partition.
          Ubox(1:2,1:3) = idis%box(1:2,1:3,P1)
          Ubox(1,iaxis) = Ubox(1,iaxis) - 1
          if (Ubox(1,iaxis).lt.1) Ubox(1,iaxis) = nm(iaxis)
          Ubox(2,iaxis)   = Ubox(1,iaxis)

          Lbox(1:2,1:3) = idis%box(1:2,1:3,P1)
          Lbox(2,iaxis) = Lbox(2,iaxis) + 1
          if (Lbox(2,iaxis).gt.nm(iaxis)) Lbox(2,iaxis) = 1
          Lbox(1,iaxis) = Lbox(2,iaxis)

          do P2= P1+1, Nodes
            Dbox => idis%box(:,:,P2)
            call boxIntersection( Dbox, Ubox, Ibox, inters )
            if (inters) then
              ncom = ncom + 1
            else
              call boxIntersection( Dbox, Lbox, Ibox, inters )
              if (inters) ncom = ncom + 1
            endif
          enddo
        enddo

        Gcom = ncom
        ! Create a list of communications needed to send/receive the extencil.
        if (Gcom.gt.0) then
          nullify( src, dst )
          call re_alloc( src, 1, Gcom, 'src', moduName )
          call re_alloc( dst, 1, Gcom, 'dst', moduName )

          ncom = 0
          do P1= 1, Nodes
            Ubox(1:2,1:3) = idis%box(1:2,1:3,P1)
            Ubox(1,iaxis) = Ubox(1,iaxis) - 1
            if (Ubox(1,iaxis).lt.1) Ubox(1,iaxis) = nm(iaxis)
            Ubox(2,iaxis)   = Ubox(1,iaxis)

            Lbox(1:2,1:3) = idis%box(1:2,1:3,P1)
            Lbox(2,iaxis) = Lbox(2,iaxis) + 1
            if (Lbox(2,iaxis).gt.nm(iaxis)) Lbox(2,iaxis) = 1
            Lbox(1,iaxis) = Lbox(2,iaxis)

            do P2= P1+1, Nodes
              Dbox => idis%box(:,:,P2)
              call boxIntersection( Dbox, Ubox, Ibox, inters )
              if (inters) then
                ncom = ncom + 1
                src(ncom) = P1
                dst(ncom) = P2
              else
                call boxIntersection( Dbox, Lbox, Ibox, inters )
                if (inters) then
                  ncom = ncom + 1
                  src(ncom) = P1
                  dst(ncom) = P2
                endif
              endif
            enddo
          enddo

          comm%np = Nodes
          ! Reschedule the communications in order to minimize the time.
          call scheduleComm( Gcom, src, dst, comm )

          ! Count the number of communications needed by the current process.
          ncom = 0
          do P1= 1, comm%ncol
            if (comm%ind(P1,Node+1).ne.0) ncom = ncom + 1
          enddo
          Lcom = ncom

          ! Store the ordered list of communications needed by the current
          ! process to send/receive the extencil.
          if (Lcom.gt.0) then
            nullify( mcomm%src, mcomm%dst )
            call re_alloc( mcomm%src, 1, Lcom, 'mcomm%src', moduName )
            call re_alloc( mcomm%dst, 1, Lcom, 'mcomm%dst', moduName )

            ncom = 0
            do P1= 1, comm%ncol
              ii = comm%ind(P1,Node+1)
              if (ii.ne.0) then
                ncom            = ncom + 1
                mcomm%src(ncom) = src(ii)
                mcomm%dst(ncom) = dst(ii)
              endif
            enddo
            mcomm%ncom = Lcom

            call de_alloc( comm%ind, 'comm%ind', 'scheComm' )
          endif

          call de_alloc( dst, 'dst', moduName )
          call de_alloc( src, 'src', moduName )
        endif
      enddo
      end subroutine initMeshExtencil

      subroutine distExtMeshData( iDistr, iaxis, BS, NSM, NN, NSPIN,
     &                            maxp, NMeshG, dens, BDENS )
        !! Send/receive the extencil information from the "dens" matrix to the
        !! temporary array "BDENS".
        !!
        !! We have a different code for every axis. We should find if we
        !! intersect with a neightbour node throught the upper, the lower,
        !! or both sides.

      use mpi_siesta

      implicit none
      integer , intent(in)  :: iDistr
        !! Distribution index to be used.
      integer , intent(in)  :: iaxis
        !! Axis to be partitioned.
      integer , intent(in)  :: BS
        !! Dimension of a plane in the current axis.
      integer , intent(in)  :: NSM
        !! Number of mesh sub-divisions in each direction.
      integer , intent(in)  :: NN
        !! Size of the extencil.
      integer , intent(in)  :: NSPIN
        !! Total number of spin components.
      integer , intent(in)  :: maxp
        !! Total number of points.
      integer , intent(in)  :: NMeshG(3)
        !! Number of mesh points in each cell direction.
      real(gp), intent(in)  :: DENS(maxp,NSPIN)
        !! Electron density matrix (or similar array).
      real(gp), intent(out) :: BDENS(BS,2*NN,NSPIN)
        !! Auxiliary array to store the extencil from other partitions.

      integer                     :: Ubox(2,3), Lbox(2,3), IUbox(2,3),
     &                               ILbox(2,3), nm(3), ispin, Cnode,
     &                               iniX, endX, iniY, endY, iniZ, endZ,
     &                               ix, iy, iz, tt, uu, dimB(3), ii, PP
      logical                     :: inter1, inter2
      integer,            pointer :: Dbox(:,:)
      real(gp),           pointer :: SBUF(:), RBUF(:)
      type(meshDisType),  pointer :: idis
      type(meshCommType), pointer :: mcomm
      integer                     :: MPIerror
      MPI_STATUS_TYPE             :: Status

      idis    => meshDistr(iDistr)
      mcomm   => exteCommu(iDistr,iaxis)
      nm      = NMeshG/NSM
      Cnode   = Node + 1
      dimB(1) = (idis%box(2,1,Cnode)-idis%box(1,1,Cnode)+1)*NSM
      dimB(2) = (idis%box(2,2,Cnode)-idis%box(1,2,Cnode)+1)*NSM
      dimB(3) = (idis%box(2,3,Cnode)-idis%box(1,3,Cnode)+1)*NSM

      if (.not.associated(mcomm%dst)) then
        write(6,*) 'ERROR: Trying to communicate extencil ',
     &             'with an uninitialized mesh distribution'
        call die("stopping program")
      endif
      if (.not.associated(mcomm%src)) then
        write(6,*) 'ERROR: Trying to communicate extencil ',
     &             'with an uninitialized mesh distribution'
        call die("stopping program")
      endif

      nullify(SBUF,RBUF)
      call re_alloc( SBUF, 1, BS*NN*nspin, 'SBUF', moduName )
      call re_alloc( RBUF, 1, BS*NN*nspin, 'RBUF', moduName )

      Ubox(1:2,1:3) = idis%box(1:2,1:3,Cnode)
      Ubox(1,iaxis) = Ubox(1,iaxis) - 1
      if (Ubox(1,iaxis).lt.1) Ubox(1,iaxis) = nm(iaxis)
      Ubox(2,iaxis)   = Ubox(1,iaxis)

      Lbox(1:2,1:3) = idis%box(1:2,1:3,Cnode)
      Lbox(2,iaxis) = Lbox(2,iaxis) + 1
      if (Lbox(2,iaxis).gt.nm(iaxis)) Lbox(2,iaxis) = 1
      Lbox(1,iaxis) = Lbox(2,iaxis)

      do ii= 1, mcomm%ncom
        if (Cnode.eq.mcomm%src(ii)) then
          PP = mcomm%dst(ii)
        else
          PP = mcomm%src(ii)
        endif
        Dbox => idis%box(:,:,PP)
        call boxIntersection( Dbox, Ubox, IUbox, inter1 )
        call boxIntersection( Dbox, Lbox, ILbox, inter2 )
        if (inter1) then
          if (iaxis.eq.1) then
            iniX = 1
            endX = NN
            iniY = (IUbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (IUbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            iniZ = (IUbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (IUbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= iniZ, endZ
                do iy= iniY, endY
                  uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = dens(uu,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter2) then
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(2)+iniY
                  do iy= iniY, endY
                    do ix= 1, NN
                      tt = tt + 1
                      BDENS(uu,NN+ix,ispin) = RBUF(tt)
                    enddo
                    uu = uu + 1
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(2)+iniY
                  do iy= iniY, endY
                    do ix= NN, 1, -1
                      tt = tt + 1
                      BDENS(uu,ix,ispin) = RBUF(tt)
                    enddo
                    uu = uu + 1
                  enddo
                enddo
              enddo
            endif
          else if (iaxis.eq.2) then
            iniX = (IUbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (IUbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniZ = (IUbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (IUbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= iniZ, endZ
                do iy= 1, NN
                  uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = dens(uu,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )

            if (inter2) then
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  do iy= 1, NN
                    uu = (iz-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,NN+iy,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  do iy= NN, 1, -1
                    uu = (iz-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,iy,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif
          else
            iniX = (IUbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (IUbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniY = (IUbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (IUbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= 1, NN
                do iy= iniY, endY
                  uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = dens(uu,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )

            if (inter2) then
              tt = 0
              do ispin= 1, nspin
                do iz= 1, NN
                  do iy= iniY, endY
                    uu = (iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,NN+iz,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= NN, 1, -1
                  do iy= iniY, endY
                    uu = (iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,iz,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif
          endif
        endif

        if (inter2) then
          if (iaxis.eq.1) then
            iniX = dimB(1)-NN+1
            endX = dimB(1)
            iniY = (ILbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (ILbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            iniZ = (ILbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (ILbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= iniZ, endZ
                do iy= iniY, endY
                  uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = dens(uu,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter1) then
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(2)+iniY
                  do iy= iniY, endY
                    do ix= NN, 1, -1
                      tt = tt + 1
                      BDENS(uu,ix,ispin) = RBUF(tt)
                    enddo
                    uu = uu + 1
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(2)+iniY
                  do iy= iniY, endY
                    do ix= 1, NN
                      tt = tt + 1
                      BDENS(uu,NN+ix,ispin) = RBUF(tt)
                    enddo
                    uu = uu + 1
                  enddo
                enddo
              enddo
            endif
          else if (iaxis.eq.2) then
            iniX = (ILbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (ILbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniZ = (ILbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (ILbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= iniZ, endZ
                do iy= dimB(2)-NN+1, dimB(2)
                  uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = dens(uu,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter1) then
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  do iy= NN, 1, -1
                    uu = (iz-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,iy,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= iniZ, endZ
                  do iy= 1, NN
                    uu = (iz-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,NN+iy,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif

          else
            iniX = (ILbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (ILbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniY = (ILbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (ILbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= dimB(3)-NN+1, dimB(3)
                do iy= iniY, endY
                  uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = dens(uu,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter1) then
              tt = 0
              do ispin= 1, nspin
                do iz= NN, 1, -1
                  do iy= iniY, endY
                    uu = (iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,iz,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= 1, NN
                  do iy= iniY, endY
                    uu = (iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      BDENS(uu,NN+iz,ispin) = RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif

          endif
        endif
      enddo

      call de_alloc( RBUF, 'RBUF', moduName )
      call de_alloc( SBUF, 'SBUF', moduName )
      end subroutine distExtMeshData

      subroutine gathExtMeshData( iDistr, iaxis, BS, NSM, NN, NSPIN,
     &                            maxp, NMeshG, BVXC, VXC )
        !! Send/receive the extencil information from the "BVXC" temporary array
        !! to the array "VXC". This is the inverse of the routine above.
        !!
        !! We have a different code for every axis. We should find if we
        !! intersect with a neightbour node throught the upper, the lower,
        !! or both sides.

      use mpi_siesta
      implicit none
      integer , intent(in)  :: iDistr
        !! Distribution index to be used.
      integer , intent(in)  :: iaxis
        !! Axis to be partitioned.
      integer , intent(in)  :: BS
        !! Dimension of a plane in the current axis.
      integer , intent(in)  :: NSM
        !! Number of mesh sub-divisions in each direction.
      integer , intent(in)  :: NN
        !! Size of the extencil.
      integer , intent(in)  :: NSPIN
        !! Total number of spin components.
      integer , intent(in)  :: maxp
        !! Total number of points.
      integer , intent(in)  :: NMeshG(3)
        !! Number of mesh points in each cell direction.
      real(gp), intent(in)  :: BVXC(BS,2*NN,NSPIN)
        !! Auxiliar array that contains the extencil of the
        !! exchange-correlation potential (or similar array).
      real(gp), intent(out) :: VXC(maxp,NSPIN)
        !! Exchange-correlation potential.

      integer                     :: Ubox(2,3), Lbox(2,3), IUbox(2,3),
     &                               ILbox(2,3), nm(3), ispin, Cnode,
     &                               iniX, endX, iniY, endY, iniZ, endZ,
     &                               ix, iy, iz, tt, uu, dimB(3), ii, PP
      logical                     :: inter1, inter2
      integer,            pointer :: Dbox(:,:)
      real(gp),           pointer :: SBUF(:), RBUF(:)
      type(meshDisType),  pointer :: idis
      type(meshCommType), pointer :: mcomm
      integer                     :: MPIerror
      MPI_STATUS_TYPE             :: Status

      idis    => meshDistr(iDistr)
      mcomm   => exteCommu(iDistr,iaxis)
      nm      = NMeshG/NSM
      Cnode   = Node + 1
      dimB(1) = (idis%box(2,1,Cnode)-idis%box(1,1,Cnode)+1)*NSM
      dimB(2) = (idis%box(2,2,Cnode)-idis%box(1,2,Cnode)+1)*NSM
      dimB(3) = (idis%box(2,3,Cnode)-idis%box(1,3,Cnode)+1)*NSM

      nullify(SBUF,RBUF)
      call re_alloc( SBUF, 1, BS*NN*nspin, 'SBUF', moduName )
      call re_alloc( RBUF, 1, BS*NN*nspin, 'RBUF', moduName )

      Ubox(1:2,1:3) = idis%box(1:2,1:3,Cnode)
      Ubox(1,iaxis) = Ubox(1,iaxis) - 1
      if (Ubox(1,iaxis).lt.1) Ubox(1,iaxis) = nm(iaxis)
      Ubox(2,iaxis)   = Ubox(1,iaxis)

      Lbox(1:2,1:3) = idis%box(1:2,1:3,Cnode)
      Lbox(2,iaxis) = Lbox(2,iaxis) + 1
      if (Lbox(2,iaxis).gt.nm(iaxis)) Lbox(2,iaxis) = 1
      Lbox(1,iaxis) = Lbox(2,iaxis)

      do ii= 1, mcomm%ncom
        if (Cnode.eq.mcomm%src(ii)) then
          PP = mcomm%dst(ii)
        else
          PP = mcomm%src(ii)
        endif
        Dbox => idis%box(:,:,PP)
        call boxIntersection( Dbox, Ubox, IUbox, inter1 )
        call boxIntersection( Dbox, Lbox, ILbox, inter2 )
        if (inter1) then
          if (iaxis.eq.1) then
            iniX = 1
            endX = NN
            iniY = (IUbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (IUbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            iniZ = (IUbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (IUbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do ix= 1, NN
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(2)+iniY
                  do iy= iniY, endY
                    tt = tt + 1
                    SBUF(tt) = BVXC(uu,ix,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )

            if (inter2) then
              tt = 0
              do ispin= 1, nspin
                do ix= dimB(1), dimB(1)-NN+1, -1
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iniY-1)*dimB(1)+ix
                    do iy= iniY, endY
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + dimB(1)
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do ix=  1, NN
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iniY-1)*dimB(1)+ix
                    do iy= iniY, endY
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + dimB(1)
                    enddo
                  enddo
                enddo
              enddo
            endif
          else if (iaxis.eq.2) then
            iniX = (IUbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (IUbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniZ = (IUbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (IUbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iy= 1, NN
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = BVXC(uu,iy,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter2) then
              tt = 0
              do ispin= 1, nspin
                do iy= dimB(2), dimB(2)-NN+1, -1
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iy=  1, NN
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif
          else
            iniX = (IUbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (IUbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniY = (IUbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (IUbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= 1, NN
                do iy= iniY, endY
                  uu = (iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = BVXC(uu,iz,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter2) then
              tt = 0
              do ispin= 1, nspin
                do iz=  dimB(3), dimB(3)-NN+1, -1
                  do iy= iniY, endY
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz=  1, NN
                  do iy= iniY, endY
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif
          endif
        endif

        if (inter2) then
          if (iaxis.eq.1) then
            iniX = dimB(1)-NN+1
            endX = dimB(1)
            iniY = (ILbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (ILbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            iniZ = (ILbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (ILbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do ix= NN+1, 2*NN
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(2)+iniY
                  do iy= iniY, endY
                    tt = tt + 1
                    SBUF(tt) = BVXC(uu,ix,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter1) then
              tt = 0
              do ispin= 1, nspin
                do ix=  1, NN
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iniY-1)*dimB(1)+ix
                    do iy= iniY, endY
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + dimB(1)
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do ix= dimB(1), dimB(1)-NN+1, -1
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iniY-1)*dimB(1)+ix
                    do iy= iniY, endY
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + dimB(1)
                    enddo
                  enddo
                enddo
              enddo
            endif
          else if (iaxis.eq.2) then
            iniX = (ILbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (ILbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniZ = (ILbox(1,3)-idis%box(1,3,Cnode))*NSM + 1
            endZ = (ILbox(2,3)-idis%box(1,3,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iy= NN+1, 2*NN
                do iz= iniZ, endZ
                  uu = (iz-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = BVXC(uu,iy,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter1) then
              tt = 0
              do ispin= 1, nspin
                do iy= 1, NN
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iy= dimB(2), dimB(2)-NN+1, -1
                  do iz= iniZ, endZ
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif

          else
            iniX = (ILbox(1,1)-idis%box(1,1,Cnode))*NSM + 1
            endX = (ILbox(2,1)-idis%box(1,1,Cnode)+1)*NSM
            iniY = (ILbox(1,2)-idis%box(1,2,Cnode))*NSM + 1
            endY = (ILbox(2,2)-idis%box(1,2,Cnode)+1)*NSM
            tt = 0
            do ispin= 1, nspin
              do iz= NN+1, 2*NN
                do iy= iniY, endY
                  uu = (iy-1)*dimB(1)+iniX
                  do ix= iniX, endX
                    tt = tt + 1
                    SBUF(tt) = BVXC(uu,iz,ispin)
                    uu = uu + 1
                  enddo
                enddo
              enddo
            enddo
            call MPI_SendRecv( SBUF, tt, MPI_grid_real, PP-1, 0,
     &                         RBUF, tt, MPI_grid_real, PP-1, 0,
     &                         MPI_Comm_world, Status, MPIerror )
            if (inter1) then
              tt = 0
              do ispin= 1, nspin
                do iz=  1, NN
                  do iy= iniY, endY
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            else
              tt = 0
              do ispin= 1, nspin
                do iz= dimB(3), dimB(3)-NN+1, -1
                  do iy= iniY, endY
                    uu = (iz-1)*dimB(1)*dimB(2)+(iy-1)*dimB(1)+iniX
                    do ix= iniX, endX
                      tt = tt + 1
                      VXC(uu,ispin) = VXC(uu,ispin) + RBUF(tt)
                      uu = uu + 1
                    enddo
                  enddo
                enddo
              enddo
            endif
          endif
        endif
      enddo

      call de_alloc( RBUF, 'RBUF', moduName )
      call de_alloc( SBUF, 'SBUF', moduName )
      end subroutine gathExtMeshData


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Subroutines and helpers for workload-partitions.            [[WKLD]] !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function inters( box, cBox, nsm1, nsm2 ) result(v)
        !! Returns true if there is an intersection between boxes box and cBox.
      implicit none
      integer          , intent(in) :: box(2,3)
        !! Box for a given distribution.
      integer          , intent(in) :: cBox(2,3)
        !! Box in another distribution.
      integer, optional, intent(in) :: nsm1
        !! Number of mesh subdivisions in distribution 1.
      integer, optional, intent(in) :: nsm2
        !! Number of mesh subdivisions in distribution 2.

      logical :: v
      integer :: boxSM(2,3), cBoxSM(2,3)

      if ( present(nsm1) .and. present(nsm2) ) then
        boxSM(:,:)  = (box(:,:)  - 1) * nsm1 + 1
        cBoxSM(:,:) = (cBox(:,:) - 1) * nsm2 + 1
        v = .not. ( ANY(boxSM(1,:) > cBoxSM(2,:)) .or.
     &              ANY(boxSM(2,:) < cBoxSM(1,:)) )
      else
        v = .not. ( ANY(box(1,:) > cBox(2,:)) .or.
     &              ANY(box(2,:) < cBox(1,:)) )
      endif
      end function inters

      function contained( box, cBox, nsm1, nsm2 ) result(v)
      implicit none
      integer          , intent(in) :: box(2,3), cBox(2,3)
      integer, optional, intent(in) :: nsm1, nsm2

      logical :: v
      integer :: boxSM(2,3), cBoxSM(2,3)

      if ( present(nsm1) .and. present(nsm2) ) then
        boxSM(:,:)  = (box(:,:)  - 1) * nsm1 + 1
        cBoxSM(:,:) = (cBox(:,:) - 1) * nsm2 + 1
        v = ( ALL(boxSM(1,:) >= cBoxSM(1,:)) .and.
     &        ALL(boxSM(2,:) <= cBoxSM(2,:)) )
      else
        v = ( ALL(box(1,:) >= cBox(1,:)) .and.
     &        ALL(box(2,:) <= cBox(2,:)) )
      endif
      end function contained

      subroutine InitWload( this, box, wload )
        !! Initializes a workload-partitioned box.

      implicit none
      class(wloadBOX), intent(inout) :: this
        !! Current workload-partitioned box.
      integer        , intent(in)  :: box(2,3)
        !! Global box to be partitioned.
      integer        , target      :: wload(:)
        !! Array containing the weights for each of the points in the box.
      integer         :: nm(3), i, j, k, ii

      nullify(this%x,this%y,this%z,this%sb1,this%sb2)
      this%rbox = box
      this%box  = box
      this%wload => wload
      nm = box(2,:) - box(1,:) + 1
      call re_alloc( this%x, 1, nm(1), "x", moduName )
      call re_alloc( this%y, 1, nm(2), "y", moduName )
      call re_alloc( this%z, 1, nm(3), "z", moduName )

      this%total = 0
      ii = 1
      do k=1, nm(3)
        do j=1, nm(2)
          do i=1, nm(1)
            this%x(i)  = this%x(i)  + this%wload(ii)
            this%y(j)  = this%y(j)  + this%wload(ii)
            this%z(k)  = this%z(k)  + this%wload(ii)
            ii = ii + 1
          enddo
        enddo
      enddo
      do i=2, nm(1)
        this%x(i)  = this%x(i) + this%x(i-1)
      enddo
      do i=2, nm(2)
        this%y(i)  = this%y(i) + this%y(i-1)
      enddo
      do i=2, nm(3)
        this%z(i)  = this%z(i) + this%z(i-1)
      enddo
      this%total = this%z(nm(3))

      end subroutine InitWload

      recursive subroutine FreeWload( this )
        !! Deallocates and resets a workload-distributed box.

      implicit none
      class(wloadBOX), intent(inout) :: this
        !! Current workload-partitioned box.

      call de_alloc( this%x, "x", moduName )
      call de_alloc( this%y, "y", moduName )
      call de_alloc( this%z, "z", moduName )
      if (associated(this%sb1)) then
         call this%sb1%FREE( )
         deallocate(this%sb1)
      end if
      if (associated(this%sb2)) then
         call this%sb2%FREE( )
         deallocate(this%sb2)
      end if
      end subroutine FreeWload

      subroutine SplitWload( this, rwBox, box )
        !! Generates a split from (parent) rwBox into the current wlbox, using
        !! another input box as the delimiter for the partition.

      implicit none
      class(wloadBOX), intent(inout) :: this
        !! Current workload-partitioned box.
      type(wloadBOX) , intent(in)    :: rwBox
        !! Parent workload-partitioned box.
      integer        , intent(in)    :: box(2,3)
        !! Box dimensions to slice the parent workload box.
      integer :: nm(3), nmr(3), i, j, k, ii, jj, kk

      nullify(this%x,this%y,this%z,this%sb1,this%sb2)
      this%box  = box
      this%rbox = rwBox%rbox
      this%wload => rwBox%wload
      nm = box(2,:) - box(1,:) + 1
      nmr = this%rbox(2,:) - this%rbox(1,:) + 1
      call re_alloc( this%x, 1, nm(1), "x", moduName )
      call re_alloc( this%y, 1, nm(2), "y", moduName )
      call re_alloc( this%z, 1, nm(3), "z", moduName )

      this%total = 0
      kk = (this%box(1,3)-this%rbox(1,3))*nmr(2)*nmr(1)
      do k=1, nm(3)
        jj = (this%box(1,2)-this%rbox(1,2))*nmr(1)
        do j=1, nm(2)
          ii = this%box(1,1)-this%rbox(1,1) + jj + kk
          do i=1, nm(1)
            ii = ii + 1
            this%x(i)  = this%x(i)  + this%wload(ii)
            this%y(j)  = this%y(j)  + this%wload(ii)
            this%z(k)  = this%z(k)  + this%wload(ii)
          enddo
          jj = jj + nmr(1)
        enddo
        kk = kk + nmr(2)*nmr(1)
      enddo
      do i=2, nm(1)
        this%x(i)  = this%x(i) + this%x(i-1)
      enddo
      do i=2, nm(2)
        this%y(i)  = this%y(i) + this%y(i-1)
      enddo
      do i=2, nm(3)
        this%z(i)  = this%z(i) + this%z(i-1)
      enddo
      this%total = this%z(nm(3))
      end subroutine SplitWload

      recursive subroutine CutWload( this, box, cdir )
        !! Cuts a workload-partitioned box into a given direction using another
        !! box as reference.
      implicit none
      class(wloadBOX), intent(inout) :: this
        !! Current workload-partitioned box.
      integer        , intent(in)    :: box(2,3)
        !! Box used to segment the workload-partitioned box.
      integer        , intent(in)    :: cdir
        !! Direction in which to split the workload-partitioned box.

      integer :: ibox(2,3)

      if (inters( this%box, box ).and.
     &    (.not. contained( this%box, box ))) then

        if ( associated(this%sb1) ) then
          call this%sb1%CUT( box, cdir )
          call this%sb2%CUT( box, cdir )
        else
          ibox = this%box
          ibox(2,cdir) = box(2,cdir)
          allocate(this%sb1)
          call this%sb1%SPLIT( this, ibox )
          ibox(1,cdir) = box(2,cdir)+1
          ibox(2,cdir) = this%box(2,cdir)
          allocate(this%sb2)
          call this%sb2%SPLIT( this, ibox )
        endif
      endif
      end subroutine CutWload

      recursive function getTotal( this ) result(v)
        !! Returns the total workload for a workload-partitioned box.
      use m_mpi_utils, only : globalize_sum

      implicit none
      class(wloadBOX), intent(in) :: this
        !! Current workload-partitioned box.
      integer(i8b) :: v, v_in

      v_in = this%total
      call globalize_sum( v_in, v )
      end function getTotal

      recursive function getWload( this, cBox, cdir ) result(v)
        !! Returns the total workload in a given direction for a
        !! workload-partitioned box.
      use m_mpi_utils, only : globalize_sum

      implicit none
      class(wloadBOX), intent(in) :: this
        !! Current workload-partitioned box.
      integer        , intent(in) :: cBox(2,3)
        !! Box used to slice the workload-partitioned box, "masking" it.
      integer        , intent(in) :: cdir
        !! Direction in which to calculate the workload.

      integer(i8b)    :: v, v_in
      integer         :: ind

      if (inters( this%box, cBox )) then

        if (.not. associated(this%sb1)) then

          ind = MIN(cBox(2,cdir),this%box(2,cdir))-this%box(1,cdir)+1

          if (cdir==1) then
            v_in = this%x(ind)
          else if (cdir==2) then
            v_in = this%y(ind)
          else
            v_in = this%z(ind)
          endif
          call globalize_sum( v_in, v )
        else

          if (inters( this%sb1%box, cBox )) then
            v = this%sb1%getWload( cBox, cdir )
          else
            v = this%sb2%getWload( cBox, cdir )
          endif
        endif
      else

        v_in = 0
        call globalize_sum( v_in, v )
      endif
      end function getWload

      subroutine InitRecBox( this, box, np )
        !! Initializes the full-domain-based recursive box.

      implicit none
      class(recBox), intent(inout) :: this
        !! Current recursive box.
      integer      , intent(in)  :: box(2,3)
        !! Full box dimensions.
      integer      , intent(in)  :: np
        !! Total number of points.

      nullify(this%sb1,this%sb2)
      this%box    = box
      this%np     = np
      this%gwload = 0
      end subroutine InitRecBox

      recursive subroutine FreeRecBox( this )
        !! Resets the full-domain-based recursive box.

      implicit none
      class(recBox), intent(inout) :: this
        !! Current recursive box.

      if (associated(this%sb1)) then
         call this%sb1%FREE( )
         deallocate(this%sb1)
      end if
      if (associated(this%sb2)) then
         call this%sb2%FREE( )
         deallocate(this%sb2)
      end if
      end subroutine FreeRecBox

      recursive subroutine recSplit( this, wBox, gwload )
        !! Recursively split a workload array between processors, using a
        !! workload array to calculate the partitions.

      implicit none
      class(recBox) , intent(inout) :: this
        !! Current global-domain-based recursive box.
      type(wloadBOX), intent(inout) :: wBox
        !! A parallel work load box structure.
      integer(i8b)  , intent(in)    :: gwload
        !! Size of the total workload.

      integer         :: nm(3), cdir, p1, p2, obj, box(2,3), frst,
     &                   last, i
      integer(i8b)    :: owload, fwload, lwload, cwload, eex1, eex2

      this%gwload = gwload
      if (this%np < 2) return
      ! Choose the cutting direction
      nm  = this%box(2,:) - this%box(1,:) + 1
      cdir = 3
      if (nm(2).gt.nm(cdir)) cdir = 2
      if (nm(1).gt.nm(cdir)) cdir = 1

      ! Split the number of parts and Objective wload.
      P1 = this%np/2
      P2 = this%np-P1

      owload = gwload*P1/this%np
      box = this%box

      ! Use a linear interpolation to find the cutting point.
      frst   = box(1,cdir)-1
      fwload = 0
      last   = box(2,cdir)
      lwload = gwload
      i = 0
      do while((last-frst).gt.1)
        i = i + 1
        eex1 = (last - frst) * (owload - fwload)
        eex2 = lwload - fwload

        ! Here, using NINT is not possible due to both integers being long,
        ! since we would need quadruple precision reals.
        ! We first perform an integer division
        obj = eex1 / eex2
        ! Then, if the remainder is closer to eex2 than to 0, we round up.
        if ( 2*( eex1 - obj * eex2 ) >= eex2 ) obj = obj+1
        obj = obj + frst

        if (obj==frst) then
          obj = frst + 1
        else if (obj==last) then
          obj = last -1
        endif
        box(2,cdir) = obj
        cwload = wBox%getWload( box, cdir )
        if (cwload.eq.owload) then
          frst   = obj
          fwload = cwload
          exit
        endif
        if (cwload.lt.owload) then
          frst   = obj
          fwload = cwload
        else
          last   = obj
          lwload = cwload
        endif
      enddo

      ! Calculate excess error
      eex1 = (owload-fwload)/P2
      eex2 = (lwload-owload)/P1
      if (eex1<=eex2) then
        obj    = frst
        owload = fwload
      else
        obj    = last
        owload = lwload
      endif
      box(2,cdir) = obj

      ! Cut the work load boxes.
      call wBox%CUT( box, cdir )

      allocate(this%sb1)
      call this%sb1%INIT( box, P1 )
      call this%sb1%recSplit( wBox, owload )

      box(1,cdir) = obj+1
      box(2,cdir) = this%box(2,cdir)
      allocate(this%sb2)
      call this%sb2%INIT( box, P2 )
      call this%sb2%recSplit( wBox, gwload-owload )
      end subroutine recSplit

      ! TO-DO: Clarify what on earth is a "leaf node".
      recursive function ExtractBoxes( this, Distr, ind ) result(v)
        !! Iterative method that traverse the 'recBox'  tree looking for leaf
        !! nodes. Once we find one, we save the box inside the mesh distribution.
        !!
        !! If this is indeed a leaf node, we save the limits within the mesh
        !! distribution at the position "ind".
      implicit none
      class(recBox)    , intent(in)           :: this
        !! Current global-domain-based recursive box.
      integer          , intent(in), optional :: ind
        !! Index of the next box. Only leaf boxes are counted.
      type(meshDisType), intent(inout)        :: Distr
        !! A mesh distribution.
      integer           :: v
        !!

      if ( present(ind) ) then
         v = ind
      else
         v = 0
      end if
      if (this%np==1) then
        v = v + 1
        Distr%box(:,:,v) = this%box
      else
        v = this%sb1%ExtractBoxes( Distr, v )
        v = this%sb2%ExtractBoxes( Distr, v )
      endif
      end function ExtractBoxes

      recursive function PrintBoxes( this, wBox, ind ) result(v)
        !! Print box and associated workload of everyprocess.
        !!
        !! Iterative method that traverse the 'recBox'  tree looking for leaf
        !! nodes. Once we find a leaf, we print its limits and its work load.

      implicit none
      class(recBox) , intent(in)           :: this
              !! Current global-domain-based recursive box.
      type(wloadBOX), intent(in)           :: wBox
        !! Work load box to calculate total workloads.
      integer       , intent(in), optional :: ind
        !! Index of the next box. Only leaf boxes are counted
      integer           :: v
        !!

      v = MERGE( ind, 0, present(ind) )
      if (this%np==1) then
        v = v + 1
        if (node==0) then
          write(*,"(A,I0.3,A,6I4,A,I12)") "(", v, ") -- ",
     &      this%box, ' - wload = ', this%gwload
        endif
      else
        v = this%sb1%PrintBoxes( wBox, v )
        v = this%sb2%PrintBoxes( wBox, v )
      endif
      end function PrintBoxes

      subroutine buildDataDistr( index, nm, wload )
        !! Build a data distribution using a work load array. We use a recursive
        !! parallel method to distribute the workload. At every step we try to
        !! cut a workload box equally between the involved processors.
        !!
        !! As outputs, this routine initializes meshDistr(index) and the
        !! communications between this and the previous distributions in
        !! meshCommu(*).

      implicit none
      integer, intent(in) :: index
        !! Index of the output distribution.
      integer, intent(in) :: nm(3)
        !! Number of Mesh divisions of each cell direction.
      integer, intent(in) :: wload(:)
        !! Weights of every point of the mesh using the input distribution.

      type(meshDisType), pointer :: iDistr, oDistr
      type(wloadBOX)             :: wlBox
      type(recBox)               :: reBox
      integer                    :: dbox(2,3), ii, jj
      integer,           pointer :: box(:,:)

      call timer( 'BuildDISTR', 1 )
      iDistr => meshDistr(1)
      oDistr => meshDistr(index)

      ! Initiallize a recursive box containing all the domain
      dbox(1,:) = 1
      dbox(2,:) = nm
      call reBox%INIT( dBox, NODES )

      ! Initiallize a recursive box containing local domain workload
      box => iDistr%box(:,:,NODE+1)
      call wlBox%INIT( box, wload )

      ! Split the global domain, using the local work load array
      call reBox%recSplit( wlBox, wlBox%getTotal( ) )

      ! Extract the new data distribution
      ii = reBox%ExtractBoxes( oDistr )

      ! Precompute the communications needed to move data between the new data
      ! distribution and the previous ones.
      jj = ((index-2)*(index-1))/2 + 1
      do ii=1, index-1
        iDistr => meshDistr(ii)
        call compMeshComm( iDistr, oDistr, meshCommu(jj) )
        jj = jj + 1
      enddo

      call reBox%FREE( )
      call wlBox%FREE( )
      call timer( 'BuildDISTR', 2 )
      if (fdf_get("mesh-comm-timings",.false.)) then
         call timer( 'BuildDISTR', 3 )
      endif
      end subroutine  buildDataDistr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute mesh distribution communications.                   [[MCMM]] !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine compMeshComm( distr1, distr2, mcomm )
        !! Find the communications needed to transform one array that uses
        !! distribution "distr1" to distribution "distr2". In order to do this,
        !! we count the number of intersections between the source distribution
        !! and the target distribution. Every intersection represents a
        !! communication. We then create a pattern of communications in order to
        !! create an scheduling:
        !!
        !!   * 1 to the left/right
        !!   * 1 to the right/left
        !!   * 2 to the left/right
        !!   * 2 to the right/left
        !!   * 3 ... and so on.

      implicit none
      type(meshDisType),   intent(in) :: distr1
        !! Source distribution
      type(meshDisType),   intent(in) :: distr2
        !! Targete distribution.
      type(meshCommType), intent(inout) :: mcomm
        !! Stores the communications needed.

      integer           :: ii, ib, i2, onode, fac1, fac2, ncom,
     &                     nsm1, nsm2
      integer,  pointer :: src(:), dst(:), sbox(:,:), rbox(:,:)

      call timer( "meshCom-Ca", 1 )
      nullify(src,dst)
      call re_alloc( src, 1, 2*NODES-1, 'src', moduName )
      call re_alloc( dst, 1, 2*NODES-1, 'dst', moduName )

      ncom = 0
      sbox => distr1%box(:,:,NODE+1)
      rbox => distr2%box(:,:,NODE+1)
      nsm1 = distr1%nsm
      nsm2 = distr2%nsm

      if ( inters( sbox, rbox, nsm1, nsm2 ) ) then
        ncom = ncom + 1
        src(ncom) = NODE+1
        dst(ncom) = NODE+1
      endif
      fac2 = 1

      do ii = 0, NODES-2
        i2 = ISHFT(ii,-1)+1
        IB = ISHFT(i2,1)
        fac1 = MERGE(1,-1,MOD(NODE,IB)/I2==0)

        oNode = NODE + (ii/2+1)*fac1*fac2
        if (oNode < 0)      oNode = oNode + NODES
        if (oNode >= NODES) oNode = oNode - NODES

        if ( NODE < oNODE ) then
          sbox => distr1%box(:,:,NODE+1)
          rbox => distr2%box(:,:,oNODE+1)
          if ( inters( sbox, rbox, nsm1, nsm2 ) ) then
            ncom = ncom + 1
            src(ncom) = NODE+1
            dst(ncom) = oNODE+1
          endif
          sbox => distr1%box(:,:,oNODE+1)
          rbox => distr2%box(:,:,NODE+1)
          if ( inters( sbox, rbox, nsm1, nsm2 ) ) then
            ncom = ncom + 1
            src(ncom) = oNODE+1
            dst(ncom) = NODE+1
          endif

        else
          sbox => distr1%box(:,:,oNODE+1)
          rbox => distr2%box(:,:,NODE+1)
          if ( inters( sbox, rbox, nsm1, nsm2 ) ) then
            ncom = ncom + 1
            src(ncom) = oNODE+1
            dst(ncom) = NODE+1
          endif
          sbox => distr1%box(:,:,NODE+1)
          rbox => distr2%box(:,:,oNODE+1)
          if ( inters( sbox, rbox, nsm1, nsm2 ) ) then
            ncom = ncom + 1
            src(ncom) = NODE+1
            dst(ncom) = oNODE+1
          endif

        endif
        fac2 = -fac2
      enddo


      ! Save data in the comm structure
      call re_alloc( mcomm%src, 1, ncom, 'mcomm%src', moduName )
      call re_alloc( mcomm%dst, 1, ncom, 'mcomm%dst', moduName )
      mcomm%src = src(:ncom)
      mcomm%dst = dst(:ncom)
      mcomm%ncom = ncom

      call de_alloc( src, 'src', moduName )
      call de_alloc( dst, 'dst', moduName )
      call timer( "meshCom-Ca", 2 )
      if (fdf_get("mesh-comm-timings",.false.)) then
         call timer( "meshCom-Ca", 3 )
      endif
      end subroutine compMeshComm

#endif

C ==================================================================
C Returns the box of selected distribution
C ==================================================================
C SUBROUTINE compMeshComm( distr, box )
C
C INPUT:
C meshDisType  distr     : Selected distribution
C
C OUTPUT:
C integer pointer box    : pointer to a box
C
C BEHAVIOR:
C Returns the box of selected distribution
C
C ==================================================================
      subroutine getMeshBox( distr, box )
      implicit none
      integer :: distr
      integer, pointer :: box(:,:,:)
      box => meshDistr(distr)%box
      end subroutine getMeshBox

      END MODULE moreMeshSubs
