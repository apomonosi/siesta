!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
#include "mpi_macros.f"
module m_mesh_node

! Module for retaining information about the mesh.
! It is used by the Hartree module and the bias module.

  use precision, only : dp

  implicit none

  public
  save

  ! The offset for the current node
  real(dp) :: offset_r(3) = 0._dp
  ! the voxel-vectors for each sub-mesh element
  real(dp) :: dL(3,3) = 0._dp
  ! the voxel length along each direction
  real(dp) :: dMesh(3) = 0._dp

  ! offsets for the local node
  integer :: meshl(3) = 0
  integer :: offset_i(3) = 0

contains

  subroutine init_mesh_node( ucell , meshG , meshLim , nsm)

    use intrinsic_missing, only : VNORM
    use parallel, only : Node, Nodes, IONode
    use units, only: Ang
#ifdef MPI
    use mpi_siesta
#endif

    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of mesh divisions of each lattice vector
    integer, intent(in) :: meshG(3), meshLim(2,3)
    ! Number of fine points per big point (see iogrid_netcdf)
    integer, intent(in) :: nsm 

    ! Number of big division points
    integer :: nm(3)
    ! Processor specifics
    integer :: ProcessorZ, blocY, blocZ, nremY, nremZ, idx(3)
    ! dimension tracking of the divisions
    integer :: iniX, iniY, iniZ, dimX, dimY, dimZ
    ! Local node dimensionality of the grid
    integer :: ldimX, ldimY, ldimZ
    ! Loop stuff
    integer :: node_Y, node_Z, cur_Node
    integer :: i
#ifdef MPI
    real(dp) :: ll(3)
    integer :: ix, iy, iz
    integer :: MPIerror
    MPI_STATUS_TYPE :: status
#endif
    
    ! We calculate the spacing in each direction
    ! Notice that we now have:
    ! dL(1,1) = dX for stepping in the x-direction
    ! dL(1,2) = dX for stepping in the y-direction
    ! dL(1,3) = dX for stepping in the z-direction
    do i = 1 , 3 
       ldimX = max(meshG(i),1)
       ! The dimension stepping in each direction.
       dL(:,i) = ucell(:,i) / ldimX
    end do
    ! The voxel box-size in Cartesian coordinates 
    ! is calculated by adding all three vectors
    dMesh = matmul(dL,(/1._dp,1._dp,1._dp/))

    ! For nodes == 1 we have no offset
    ! (also some of the arrays are not initialized, which
    !  could lead to errors)
    if ( Nodes == 1 ) then
       meshl = meshG ! same as meshLim(2,:) * nsm
       return
    end if

    ! Find quantities in mesh coordinates
    meshl(:) = (meshLim(2,:) - meshLim(1,:)+1)*nsm

    ! Calculate starting point for grid
    offset_i(:) = (meshLim(1,:) - 1)*nsm
    offset_r(:) = offset_i(1)*dL(:,1) + offset_i(2)*dL(:,2) + offset_i(3)*dL(:,3)

  end subroutine init_mesh_node

  elemental subroutine mesh_correct_idx(mesh,idx)
    integer, intent(in) :: mesh
    integer, intent(inout) :: idx
    ! negative "supercell"
    do while ( idx <= 0 )
       idx = idx + mesh
    end do
    ! positive "supercell"
    do while ( mesh < idx )
       idx = idx - mesh
    end do
    
  end subroutine mesh_correct_idx
  
end module m_mesh_node

