! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_tbt_tri_scat_spinor

  use precision, only : dp
  use units, only : Pi
  use m_region
  
  use class_zTriMat

  use m_ts_tri_scat_spinor, only : GF_Gamma_GF, dir_GF_Gamma_GF
  use m_ts_tri_common, only : GFGGF_needed_worksize

  use ts_electrode_m

  implicit none

  private

  public :: A_DOS   ! Spectral function density of states
  public :: GF_DOS  ! Green's function density of states
  public :: A_Gamma ! Calculate the transmission from spectral function . Gamma
  public :: A_Gamma_Block ! Calculate the transmission from spectral function . Gamma (in block form)
  public :: TT_eigen ! Eigenvalue calculation of the transmission eigenvalues
  public :: GF_Gamma ! Calculate the transmission from Green function . Gamma (same-lead contribution)
  public :: GF_T, GF_T_solve
  
  public :: insert_Self_Energy
  public :: insert_Self_Energy_Dev

  ! From ts_tri_scat
  public :: GF_Gamma_GF
  public :: dir_GF_Gamma_GF
  public :: GFGGF_needed_worksize
#ifdef NCDF_4
  public :: GF_COOP, A_COOP
  public :: GF_COHP, A_COHP
  public :: GF_COHP_add_dH, A_COHP_add_dH
  public :: orb_current
  public :: orb_current_add_dH
  public :: GF_DM, A_DM
#endif

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = cmplx( 0._dp, 0._dp, dp)
  complex(dp), parameter :: z1  = cmplx( 1._dp, 0._dp, dp)
  complex(dp), parameter :: zm1 = cmplx(-1._dp, 0._dp, dp)
  complex(dp), parameter :: zi  = cmplx( 0._dp, 1._dp, dp)
#ifdef USE_GEMM3M
# define GEMM zgemm3m
#else
# define GEMM zgemm
#endif

contains

  ! Calculate the DOS from a non-fully calculated Green function.
  ! We assume that the diagonal Green function matrices are already calculated
  ! and the remaining \tilde X and \tilde Y matrices are present.
  !    all GF_nn are in Gf_tri
  !    all \tilde Yn and \tilde Xn are in Gf_tri
  ! After this routine, all off-diagonal Gf blocks are in work_tri (correctly
  ! positioned).
  ! I.e. the full Gf in the blocks can be extracted from Gf_tri and work_tri.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine GF_DOS(r,Gfd_tri,Gfo_tri,S_1D,pvt,DOS)
    use class_Sparsity
    use class_zSpData1D
    
    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(zSpData1D), intent(inout) :: S_1D ! (transposed S(k))
    type(tRgn), intent(in) :: pvt
    real(dp), intent(out) :: DOS(:)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd(2,2)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: np, n, no_o, no_i
    integer :: br, io, ind, bc
    real(dp) :: lDOS(4)

#ifdef TBTRANS_TIMING
    call timer('Gf-DOS',1)
#endif

    np = parts(Gfd_tri)
    
    ! First calculate all off-diagonal green-function elements
    no_o = nrows_g(Gfd_tri,1)
    no_i = nrows_g(Gfd_tri,2)
    call calc(2,1)
    do n = 2, np - 1
      no_o = nrows_g(Gfd_tri,n)
      no_i = nrows_g(Gfd_tri,n + 1)
      call calc(n+1,n)
      no_i = nrows_g(Gfd_tri,n - 1)
      call calc(n-1,n)
    end do
    no_o = nrows_g(Gfd_tri,np)
    no_i = nrows_g(Gfd_tri,np-1)
    call calc(np-1,np)

    ! At this point we have calculated all Green function matrices
    ! All diagonal elements are in Gfd_tri,
    ! all off-diagonal elements are in Gfo_tri

    ! The DOS per orbital is calculated like this (.=matrix multiplication);
    !  (x)=tensor product):
    !
    !   DOS(io,l) = - \sum_{is} Im[ (Gf-Gf^\dagger) . \sigma_l(x)S ](io,is,io,is) / Pi
    !             = - \sum_{jo,is,js} Im[ {(Gf-Gf^\dagger)(io,is, jo,js)}
    !                                     * \sigma_l(js, is) * S(jo, io)] / Pi
    ! With the spinor matrices
    !  \sigma_1 = (1  0)  \sigma_2 = (0  0)
    !             (0  0)             (0  1)
    !  \sigma_3 = (0  1)  \sigma_4 = (0 -i)
    !             (1  0)             (i  0)
    !
    ! DOS(io,1) = DOS
    ! DOS(io,2) = Magnetization of orbital along x
    ! DOS(io,3) = Magnetization of orbital along y
    ! DOS(io,4) = Magnetization of orbital along z
    !
    ! The fact that we need Gf - Gf^\dagger can be checked
    ! by a simple tight-binding calculation with large overlap matrices
    ! and *very* small eta values (and using sum_elec ADOS == DOS).
    ! In this case the k-resolved DOS is only correct if one uses
    ! the above equation. It should however be noted that the full
    ! DOS is independent on Gf or Gf - Gf^\dagger choice!

    sp => spar(S_1D)
    S => val(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

!$OMP parallel do default(shared), private(br,io,lDOS,ind,bc,GfGfd)
    do br = 1, r%n
      io = r%r(br)

      ! Loop columns in S(k)^T (actually the rows)
      lDOS = 0._dp
      do ind = l_ptr(io) + 1, l_ptr(io) + ncol(io)
        bc = pvt%r(l_col(ind))
        if ( bc > 0 ) then
          call calc_GfGfd(br, bc, GfGfd)
          lDOS(1) = lDOS(1) + aimag( GfGfd(1,1) * S(ind) )
          lDOS(2) = lDOS(2) + aimag( GfGfd(2,2) * S(ind) )
          lDOS(3) = lDOS(3) + aimag( (GfGfd(1,2) + GfGfd(2,1)) * S(ind) )
          lDOS(4) = lDOS(4) +  real( (GfGfd(1,2) - GfGfd(2,1)) * S(ind),dp )
        end if
      end do
       
      ! DOS(br,:) = - lDOS(:) / (2._dp * Pi)
      DOS(br) = - lDOS(1) / (2._dp * Pi)
       
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('Gf-DOS',2)
#endif

  contains

    subroutine calc(m,n)
      integer, intent(in) :: m,n
      complex(dp), pointer :: Gf(:), Mnn(:), XY(:)

      XY => val(Gfd_tri,m,n)
      Mnn => val(Gfd_tri,n,n)
      Gf => val(Gfo_tri,m,n)
      
      ! We need to calculate the 
      ! Mnm1n/Mnp1n Green's function
      call GEMM ('N','N',no_i,no_o,no_o, &
          zm1,XY(1),no_i,Mnn(1),no_o,z0,Gf(1),no_i)
      
    end subroutine calc

    ! Instead of evaluating a single element, we opearate on 2x2 blocks
    ! GfGfd(br,bc) --> ((br,UP  , bc,UP  ) (br,UP  , bc,DOWN))
    !                  ((br,DOWN, bc,UP  ) (br,DOWN, bc,DOWN))
    pure subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G(2,2)
      integer :: p_r, i_r, j_r, N_r, p_c, i_c, j_c, N_c, i

      call part_index(Gfo_tri, 2*br-1, p_r, i_r)
      call part_index(Gfo_tri, 2*bc-1, p_c, i_c)
      N_r = Gfo_tri%data%tri_nrows(p_r)
      N_c = Gfo_tri%data%tri_nrows(p_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            G(j_r,j_c) = Gfd(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfd(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      else
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            i = index_block(Gfo_tri, p_r, p_c)
            G(j_r,j_c) = Gfo(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            i = index_block(Gfo_tri, p_c, p_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfo(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      end if

    end subroutine calc_GfGfd

  end subroutine GF_DOS


  ! Calculate the DOS from a fully calculated spectral function.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine A_DOS(r,A_tri,S_1D,pvt,DOS)
    use class_Sparsity
    use class_zSpData1D

    type(tRgn), intent(in) :: r ! BTD matrix elements
    type(zTriMat), intent(inout) :: A_tri
    type(zSpData1D), intent(inout) :: S_1D ! (transposed S(k))
    type(tRgn), intent(in) :: pvt ! from sparse matrix to BTD
    real(dp), intent(out) :: DOS(:)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:), A(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: io, ind, idx(2,2), br, bc, ii, jj
    real(dp) :: lDOS(4)

#ifdef TBTRANS_TIMING
    call timer('A-DOS',1)
#endif

    ! Get data arrays
    A => val(A_tri)

    sp => spar(S_1D)
    S => val(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    ! The DOS per orbital is calculated like this:
    !   ADOS(io,l) = - \sum_{is} Re[ A . \sigma_l(x)S ](io,is,io,is) / 2Pi
    !              = - \sum_{jo,is,js} Re[A(io,is, jo,js) * \sigma_l(js, is) * S(jo, io)] / 2Pi
    !
    ! With the spinor matrices
    !  \sigma_1 = (1  0)  \sigma_2 = (0  0)
    !             (0  0)             (0  1)
    !  \sigma_3 = (0  1)  \sigma_4 = (0 -i)
    !             (1  0)             (i  0)
    !
    ! DOS(io,1) : DOS spin up
    ! DOS(io,2) = DOS spin down
    ! DOS(io,3) = Magnetization of orbital along x
    ! DOS(io,4) = Magnetization of orbital along y

!$OMP parallel do default(shared), private(br,io,lDOS,ind,bc,idx,ii,jj)
    do br = 1, r%n
      io = r%r(br)

      ! Loop columns in S(k)^T (actually the rows)
      lDOS = 0._dp
      do ind = l_ptr(io) + 1, l_ptr(io) + ncol(io)
        bc = pvt%r(l_col(ind))
        if ( bc > 0 ) then
          do ii = 1, 2
            do jj = 1, 2
              idx(ii,jj) = index(A_tri, 2*(br-1)+ii, 2*(bc-1)+jj)
            end do
          end do
          lDOS(1) = lDOS(1) +  real( A(idx(1,1)) * S(ind),dp )
          lDOS(2) = lDOS(2) +  real( A(idx(2,2)) * S(ind),dp )
          lDOS(3) = lDOS(3) +  real( (A(idx(1,2)) + A(idx(2,1))) * S(ind),dp )
          lDOS(4) = lDOS(4) - aimag( (A(idx(1,2)) - A(idx(2,1))) * S(ind) )
        end if
      end do
      
      ! DOS(br,:) = lDOS(:) / (2._dp * Pi)
      DOS(br) = lDOS(1) / (2._dp * Pi)
      
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('A-DOS',2)
#endif

  end subroutine A_DOS
  
#ifdef NCDF_4

  ! Calculate the COOP contribution from a fully calculated Green function.
  ! We assume that the Green function distribution like this:
  !    all GF_nn are in Gfd_tri (diagonal)
  !    all GF_mn (m/=n) are in Gfo_tri (off-diagonal)
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine GF_COOP(r,Gfd_tri,Gfo_tri,pvt,sp,S,sc_off,k,ph,COOP)
    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(tRgn), intent(in) :: pvt
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: S(:)
    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(dSpData1D), intent(inout) :: COOP

    type(Sparsity), pointer :: c_sp
    real(dp), pointer :: C(:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd(2,2)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer :: no_u, br, io, ind, iind, bc
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('Gf-COOP',1)
#endif

#ifdef TBT_PHONON
    call die('Currently not implemented for PHtrans')
#endif

    ! Extract COOP/COHP by looping the sparse matrix
    ! The following discussion in concerning COOP, but
    ! there is no ambiguity in the two methods.
    
    ! The COOP calculation can be written as
    !
    !   COOP(io,jo) = - Im{ [Gf - Gf^\dagger](io,jo) * S(jo,io) * e^(-ik.R) } / 2Pi
    ! Here we want:
    !   DOS(io) = \sum_jo COOP(io,jo)
    ! since we know that COOP(io,jo) is the io -> jo DOS.
    ! As COOP is interesting in the supercell picture we have
    ! to calculate it with the daggered component (Gf - Gf^\dagger) (also why we need /2).
    ! Note that this is not necessary if S is S(k). I.e. it is because
    ! we want the cross-cell COOP curves as well.

    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    call attach(sp,nrows_g=no_u, n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    c_sp => spar(COOP)
    C => val(COOP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

    C(:) = 0._dp

!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss,GfGfd,is,js)
    do br = 1, r%n
      io = r%r(br)
      
      ! Get lookup columns for the COOP
      call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))
      
      ! Loop on overlap entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)
        
        ! Check if the orbital exists in the region
        iind = cptr(io) + ssearch_find(ss, l_col(ind))
        
        ! if zero the element does not exist
        ! This is the case on the elements connecting out
        ! of the device region
        if ( iind <= cptr(io) ) cycle
        
        ! In non-collinear case, we need to take the trace over of spin indices.
        ! Using : S(jo,js,io,is) = S(jo,io) delta_{js,is}
        !   COOP(iind) =
        !      - Im[ sum_{is} (G(io,is,jo,is) - G^\dagger(io,is,jo,is)) * S(jo,io)} ] / 2Pi
        bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
        call calc_GfGfd(br, bc, GfGfd)

        C(iind) = -aimag( (GfGfd(1,1)+GfGfd(2,2)) * S(ind) * ph( (l_col(ind)-1)/no_u ))

      end do
          
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('Gf-COOP',2)
#endif

  contains

  ! Instead of evaluating a single element, we opearate on 2x2 blocks
  ! GfGfd(br,bc) --> ((br,UP  , bc,UP  ) (br,UP  , bc,DOWN))
  !                  ((br,DOWN, bc,UP  ) (br,DOWN, bc,DOWN))
  pure subroutine calc_GfGfd(br, bc, G)
    integer, intent(in) :: br, bc
    complex(dp), intent(inout) :: G(2,2)
    integer :: p_r, i_r, j_r, N_r, p_c, i_c, j_c, N_c, i

    call part_index(Gfo_tri, 2*br-1, p_r, i_r)
    call part_index(Gfo_tri, 2*bc-1, p_c, i_c)
    N_r = Gfo_tri%data%tri_nrows(p_r)
    N_c = Gfo_tri%data%tri_nrows(p_c)
    
    if ( p_r == p_c ) then
      i = index_block(Gfo_tri, p_r, p_c)
      ! Loop over spin components
      do j_r = 1, 2
        do j_c = 1, 2
          G(j_r,j_c) = Gfd(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
          G(j_r,j_c) = G(j_r,j_c) - conjg(Gfd(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
        end do
      end do
    else
      ! Loop over spin components
      do j_r = 1, 2
        do j_c = 1, 2
          i = index_block(Gfo_tri, p_r, p_c)
          G(j_r,j_c) = Gfo(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
          i = index_block(Gfo_tri, p_c, p_r)
          G(j_r,j_c) = G(j_r,j_c) - conjg(Gfo(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
        end do
      end do
    end if

  end subroutine calc_GfGfd

  end subroutine GF_COOP

  ! Calculate the COOP contribution from a fully calculated Green function.
  ! We assume that the Green function distribution like this:
  !    all GF_nn are in Gfd_tri (diagonal)
  !    all GF_mn (m/=n) are in Gfo_tri (off-diagonal)
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine GF_COHP(r,Gfd_tri,Gfo_tri,pvt,sp,H,sc_off,k,ph,COHP)
    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    use m_spin, only: spin

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(tRgn), intent(in) :: pvt
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: H(:,:)
    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(dSpData1D), intent(inout) :: COHP

    type(Sparsity), pointer :: c_sp
    real(dp), pointer :: C(:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd(2,2), H2D(2,2)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer :: no_u, br, io, ind, iind, bc, is, js
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('Gf-COHP',1)
#endif

#ifdef TBT_PHONON
    call die('Currently not implemented for PHtrans')
#endif

    ! Extract COHP by looping the sparse matrix
    ! The following discussion in concerning COOP, but
    ! there is no ambiguity in the two methods.
    
    ! The COOP calculation can be written as
    !
    !   COOP(io,jo) = - Im{ [Gf - Gf^\dagger](io,jo) * H(jo,io) * e^(-ik.R) } / 2Pi

    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    call attach(sp,nrows_g=no_u, n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    c_sp => spar(COHP)
    C => val(COHP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

    C(:) = 0._dp

    if ( spin%Ncol) then 
!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss,GfGfd,H2D,is,js)
      do br = 1, r%n
        io = r%r(br)
        
        ! Get lookup columns for the COOP
        call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))
        
        ! Loop on overlap entries here...
        do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)
          
          ! Check if the orbital exists in the region
          iind = cptr(io) + ssearch_find(ss, l_col(ind))
          
          ! if zero the element does not exist
          ! This is the case on the elements connecting out
          ! of the device region
          if ( iind <= cptr(io) ) cycle
          
          ! In non-collinear case, we need to take the trace over of spin indices
          !   COHP(iind) =
          !      - Im[ sum_{is,js} (G(io,is,jo,js) - G^\dagger(io,is,jo,js)) * H(jo,js,io,is)} ] / 2Pi
          bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
          call calc_GfGfd(br, bc, GfGfd)
          H2D(1,1) = cmplx(H(ind,1), H(ind,5), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(1,2) = cmplx(H(ind,3),-H(ind,4), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(2,1) = cmplx(H(ind,3), H(ind,4), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(2,2) = cmplx(H(ind,2), H(ind,6), dp) * ph( (l_col(ind)-1)/no_u )

          do is = 1,2
            do js = 1,2
              C(iind) = C(iind) - aimag( GfGfd(is,js) * H2D(js,is) )
            end do
          end do

        end do
            
      end do
!$OMP end parallel do

    else ! spin%SO
!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss,GfGfd,H2D,is,js)
      do br = 1, r%n
        io = r%r(br)
        
        ! Get lookup columns for the COOP
        call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))
        
        ! Loop on overlap entries here...
        do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)
          
          ! Check if the orbital exists in the region
          iind = cptr(io) + ssearch_find(ss, l_col(ind))
          
          ! if zero the element does not exist
          ! This is the case on the elements connecting out
          ! of the device region
          if ( iind <= cptr(io) ) cycle
          
          ! In non-collinear case, we need to take the trace over of spin indices
          !   COHP(iind) =
          !      - Im[ sum_{is,js} (G(io,is,jo,js) - G^\dagger(io,is,jo,js)) * H(jo,js,io,is)} ] / 2Pi
          bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
          call calc_GfGfd(br, bc, GfGfd)
          H2D(1,1) = cmplx(H(ind,1), H(ind,5), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(1,2) = cmplx(H(ind,3),-H(ind,4), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(2,1) = cmplx(H(ind,7), H(ind,8), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(2,2) = cmplx(H(ind,2), H(ind,6), dp) * ph( (l_col(ind)-1)/no_u )

          do is = 1,2
            do js = 1,2
              C(iind) = C(iind) - aimag( GfGfd(is,js) * H2D(js,is) )
            end do
          end do

        end do
            
      end do
!$OMP end parallel do
    end if

#ifdef TBTRANS_TIMING
    call timer('Gf-COHP',2)
#endif

  contains

    ! Instead of evaluating a single element, we opearate on 2x2 blocks
    ! GfGfd(br,bc) --> ((br,UP  , bc,UP  ) (br,UP  , bc,DOWN))
    !                  ((br,DOWN, bc,UP  ) (br,DOWN, bc,DOWN))
    pure subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G(2,2)
      integer :: p_r, i_r, j_r, N_r, p_c, i_c, j_c, N_c, i

      call part_index(Gfo_tri, 2*br-1, p_r, i_r)
      call part_index(Gfo_tri, 2*bc-1, p_c, i_c)
      N_r = Gfo_tri%data%tri_nrows(p_r)
      N_c = Gfo_tri%data%tri_nrows(p_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            G(j_r,j_c) = Gfd(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfd(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      else
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            i = index_block(Gfo_tri, p_r, p_c)
            G(j_r,j_c) = Gfo(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            i = index_block(Gfo_tri, p_c, p_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfo(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      end if

    end subroutine calc_GfGfd

  end subroutine GF_COHP

  subroutine Gf_COHP_add_dH(dH_3D,sc_off,k,ph,Gfd_tri,Gfo_tri,r,COHP,pvt)

    use class_Sparsity
    use class_zSpData3D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    type(zSpData3D), intent(in) :: dH_3D
    real(dp), intent(in) :: sc_off(:,:), k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: COHP
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: dH(:,:,:)
    type(Sparsity), pointer :: c_sp
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd(2,2)
    real(dp), pointer :: C(:)
    integer :: no_u, br, io, jo, ind, iind, is, js

#ifdef TBTRANS_TIMING
    call timer('COHP-Gf-dH',1)
#endif

    ! Retrieve dH
    sp => spar(dH_3D)
    dH => val(dH_3D)
    call attach(sp, nrows_g=no_u, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    c_sp => spar(COHP)
    C => val(COHP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    ! Create the phases
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

!$OMP parallel do default(shared), private(br,io,iind,jo,ind,col,GfGfd,is,js)
    do br = 1, r%n
      io = r%r(br)
      
      ! Loop on the COHP indices
      do iind = cptr(io) + 1, cptr(io) + cncol(io)

        ! Here we will calculate the COHP contribution from dH
        !  COHP(iind) = -Im{ [Gf(io, jo) - Gf^\dagger(io,jo)] * dH(jo, io) } / 2pi

        ! Since we are looping the dH indices we have to 

        ! Get column Gf orbital
        jo = ucorb(ccol(iind), no_u)

        ! Check if the jo,io orbital exists in dH
        if ( l_ncol(jo) > 0 ) then
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))

          ! Note that we here find the dH(jo,io) value (in the supercell picture)
          ind = l_ptr(jo) + SFIND(col, TO(ccol(iind)) + io)

          if ( ind > l_ptr(jo) ) then

            call calc_GfGfd(br, pvt%r(jo), GfGfd)
            ! COHP(iind) += sum_{is,js} 
            !     - Im[ (G(io,is,jo,js) - G^\dagger(io,is,jo,js)) * dH(jo,js,io,is)] / 2Pi
            do is = 1,2
              do js = 1,2
                C(iind) = C(iind) - aimag( GfGfd(is,js) * dH(ind,js,is) * ph( (l_col(ind)-1)/no_u ) )
              end do
            end do

          end if

        end if
        
      end do
      
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('COHP-Gf-dH',2)
#endif

  contains

    ! Instead of evaluating a single element, we opearate on 2x2 blocks
    ! GfGfd(br,bc) --> ((br,UP  , bc,UP  ) (br,UP  , bc,DOWN))
    !                  ((br,DOWN, bc,UP  ) (br,DOWN, bc,DOWN))
    pure subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G(2,2)
      integer :: p_r, i_r, j_r, N_r, p_c, i_c, j_c, N_c, i

      call part_index(Gfo_tri, 2*br-1, p_r, i_r)
      call part_index(Gfo_tri, 2*bc-1, p_c, i_c)
      N_r = Gfo_tri%data%tri_nrows(p_r)
      N_c = Gfo_tri%data%tri_nrows(p_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            G(j_r,j_c) = Gfd(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfd(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      else
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            i = index_block(Gfo_tri, p_r, p_c)
            G(j_r,j_c) = Gfo(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            i = index_block(Gfo_tri, p_c, p_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfo(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      end if

    end subroutine calc_GfGfd

    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i

      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)
        
        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('Gf_COHP_add_dH: could not find transpose supercell index')

    end function TO
    
  end subroutine Gf_COHP_add_dH

  ! Calculate the COOP contribution from a fully calculated spectral function.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine A_COOP(r,A_tri,pvt,sp,S,sc_off,k,ph,COOP)
    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: A_tri
    type(tRgn), intent(in) :: pvt
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: S(:)
    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(dSpData1D), intent(inout) :: COOP

    type(Sparsity), pointer :: c_sp
    real(dp), pointer :: C(:)
    complex(dp), pointer :: A(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer :: no_u, br, io, ind, iind, bc
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('A-COOP',1)
#endif

#ifdef TBT_PHONON
    call die('Currently not implemented for PHtrans')
#endif

    ! Extract COOP/COHP by looping the sparse matrix
    ! The following disôcussion in concerning COOP, but
    ! there is no ambiguity in the two methods.

    ! The COOP calculation can be written as
    !
    !   COOP(io,jo) = Re{ A(io,jo) * S(jo,io) * e^(ik.R) } / 2Pi
    ! Here we want:
    !   ADOS(io) = \sum_jo COOP(io,jo)
    ! since we know that COOP(io,jo) is the io -> jo ADOS.
    ! Note that this is not necessary if S is S(k). I.e. it is because
    ! we want the cross-cell COOP curves as well.

    ! Create the phases
    ! Since we have to do A.S we simply
    ! create the S(-k) (which is S^T)
    ! and thus get the correct values.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    call attach(sp,nrows_g=no_u, n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    c_sp => spar(COOP)
    C => val(COOP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    A => val(A_tri)

    C(:) = 0._dp

!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss)
    do br = 1, r%n
      io = r%r(br)

      ! Get lookup columns for the COOP
      call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))

      ! Loop on overlap entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

        ! Check if the orbital exists in the region
        iind = cptr(io) + ssearch_find(ss, l_col(ind))

        ! if zero the element does not exist
        ! This is the case on the elements connecting out
        ! of the device region
        if ( iind <= cptr(io) ) cycle

        ! COOP(iind) = Re[ Tr{ A(io,jo) * S(jo,io) } ] / (2 pi)
        !            = Re [ sum_{is} A(io,is,jo,is) * S(jo,io) ] / (2 pi)
        bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix

        ! UP-UP
        bc = index(A_tri, 2*br-1, 2*bc-1)
        C(iind) = real(A(bc) * S(ind) * ph( (l_col(ind)-1)/no_u ), dp)

        ! DN-DN
        bc = index(A_tri, 2*br, 2*bc)
        C(iind) = C(iind) + real(A(bc) * S(ind) * ph( (l_col(ind)-1)/no_u ), dp)

      end do
          
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('A-COOP',2)
#endif

  end subroutine A_COOP

  

  ! Calculate the COOP contribution from a fully calculated spectral function.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine A_COHP(r,A_tri,pvt,sp,H,sc_off,k,ph,COHP)
    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    use m_spin, only: spin

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: A_tri
    type(tRgn), intent(in) :: pvt
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: H(:,:)
    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(dSpData1D), intent(inout) :: COHP

    type(Sparsity), pointer :: c_sp
    real(dp), pointer :: C(:)
    complex(dp), pointer :: A(:)
    complex(dp) :: H2D(2,2)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer :: no_u, br, io, ind, iind, bc, ii, jj
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('A-COHP',1)
#endif

#ifdef TBT_PHONON
    call die('Currently not implemented for PHtrans')
#endif

    ! Extract COOP/COHP by looping the sparse matrix
    ! The following disôcussion in concerning COOP, but
    ! there is no ambiguity in the two methods.

    ! The COOP calculation can be written as
    !
    !   COOP(io,jo) = Re{ A(io,jo) * S(jo,io) * e^(ik.R) } / 2Pi
    ! Here we want:
    !   ADOS(io) = \sum_jo COOP(io,jo)
    ! since we know that COOP(io,jo) is the io -> jo ADOS.
    ! Note that this is not necessary if S is S(k). I.e. it is because
    ! we want the cross-cell COOP curves as well.

    ! Create the phases
    ! Since we have to do A.S we simply
    ! create the S(-k) (which is S^T)
    ! and thus get the correct values.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    call attach(sp,nrows_g=no_u, n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    c_sp => spar(COHP)
    C => val(COHP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    A => val(A_tri)

    C(:) = 0._dp

    if ( spin%Ncol ) then 

!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss,ii,jj,H2D)
      do br = 1, r%n
        io = r%r(br)

        ! Get lookup columns for the COOP
        call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))

        ! Loop on overlap entries here...
        do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

          ! Check if the orbital exists in the region
          iind = cptr(io) + ssearch_find(ss, l_col(ind))

          ! if zero the element does not exist
          ! This is the case on the elements connecting out
          ! of the device region
          if ( iind <= cptr(io) ) cycle

          ! COOP(iind) = Re[ Tr{ A(io,jo) * H(jo,io) } ] / (2 pi)
          !            = Re[ sum_{is,js} A(io,is,jo,is) * S(jo,js,io,is) ] / (2 pi)
          bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
          
          H2D(1,1) = cmplx(H(ind,1), H(ind,5), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(1,2) = cmplx(H(ind,3),-H(ind,4), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(2,1) = cmplx(H(ind,3), H(ind,4), dp) * ph( (l_col(ind)-1)/no_u )
          H2D(2,2) = cmplx(H(ind,2), H(ind,6), dp) * ph( (l_col(ind)-1)/no_u )

          do ii = 1, 2
            do jj = 1, 2
              bc = index(A_tri, 2*(br-1)+ii, 2*(bc-1)+jj)
              C(iind) = C(iind) + real(A(bc) * H2D(jj,ii))
            end do 
          end do

        end do
            
      end do
!$OMP end parallel do

    else ! spin%SO

      !$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss,ii,jj,H2D)
            do br = 1, r%n
              io = r%r(br)
      
              ! Get lookup columns for the COOP
              call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))
      
              ! Loop on overlap entries here...
              do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)
      
                ! Check if the orbital exists in the region
                iind = cptr(io) + ssearch_find(ss, l_col(ind))
      
                ! if zero the element does not exist
                ! This is the case on the elements connecting out
                ! of the device region
                if ( iind <= cptr(io) ) cycle
      
                ! COOP(iind) = Re[ Tr{ A(io,jo) * H(jo,io) } ] / (2 pi)
                !            = Re[ sum_{is,js} A(io,is,jo,is) * S(jo,js,io,is) ] / (2 pi)
                bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
                
                H2D(1,1) = cmplx(H(ind,1), H(ind,5), dp) * ph( (l_col(ind)-1)/no_u )
                H2D(1,2) = cmplx(H(ind,3),-H(ind,4), dp) * ph( (l_col(ind)-1)/no_u )
                H2D(2,1) = cmplx(H(ind,7), H(ind,8), dp) * ph( (l_col(ind)-1)/no_u )
                H2D(2,2) = cmplx(H(ind,2), H(ind,6), dp) * ph( (l_col(ind)-1)/no_u )
      
                do ii = 1, 2
                  do jj = 1, 2
                    bc = index(A_tri, 2*(br-1)+ii, 2*(bc-1)+jj)
                    C(iind) = C(iind) + real(A(bc) * H2D(jj,ii))
                  end do 
                end do
      
              end do
                  
            end do
      !$OMP end parallel do
    end if

#ifdef TBTRANS_TIMING
    call timer('A-COHP',2)
#endif

  end subroutine A_COHP

  subroutine A_COHP_add_dH(dH_3D,sc_off,k,ph,A_tri,r,COHP,pvt)

    use class_Sparsity
    use class_zSpData3D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    type(zSpData3D), intent(in) :: dH_3D
    real(dp), intent(in) :: sc_off(:,:), k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: A_tri
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: COHP
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: dH(:,:,:)
    type(Sparsity), pointer :: c_sp
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp), pointer :: A(:)
    real(dp), pointer :: C(:)
    integer :: no_u, iu, io, i, ind, iind, jo, iA, ii, jj

#ifdef TBTRANS_TIMING
    call timer('COHP-A-dH',1)
#endif

    ! Retrieve dH
    sp => spar(dH_3D)
    dH => val(dH_3D)
    call attach(sp, nrows_g=no_u, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    c_sp => spar(COHP)
    C => val(COHP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)
    
    ! Create the phases
    do i = 1 , size(sc_off, dim=2)
      ph(i-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,i)), dp)) / (Pi * 2._dp)
    end do

    A => val(A_tri)
    
!$OMP parallel do default(shared), private(iu,io,iind,jo,ind,col,iA,ii,jj)
    do iu = 1, r%n
      io = r%r(iu)

      ! Loop on the COHP indices
      do iind = cptr(io) + 1, cptr(io) + cncol(io)

        ! Here we will calculate the COHP contribution from dH
        !  COHP(iind) == A(io, jo) * dH(jo, io) / 2pi

        ! Get column A orbital
        jo = ucorb(ccol(iind), no_u)

        ! Check if the jo,io orbital exists in dH
        if ( l_ncol(jo) > 0 ) then
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))

          ! Note that we here find the dH(jo,io) value (in the supercell picture)
          ind = l_ptr(jo) + SFIND(col, TO(ccol(iind)) + io)

          if ( ind > l_ptr(jo) ) then
            do ii = 1, 2
              do jj = 1, 2
                iA = index(A_tri, 2*(iu-1)+ii, 2*(pvt%r(jo)-1)+jj)
                C(iind) = C(iind) + real(A(iA) * dH(ind,jj,ii) * ph( (l_col(ind)-1)/no_u ), dp)
              end do 
            end do

          end if

        end if

      end do
      
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('COHP-A-dH',2)
#endif
    
  contains
    
    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i
      
      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)

        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('A_COHP_add_dH: could not find transpose supercell index')

    end function TO
    
  end subroutine A_COHP_add_dH



#endif

  ! The simplest routine to do the transport calculation
  ! It takes the spectral function and multiplies it with
  ! the scattering matrix of the down-projected self-energy
  ! and calculates the transmission.
  ! We do this by taking advantage of the transposed scattering
  ! matrix: \Gamma
  subroutine A_Gamma(A_tri,El,T)

    type(zTriMat), intent(inout) :: A_tri ! Spectral function
    type(electrode_t), intent(in) :: El ! contains Gamma == (Sigma - Sigma^\dagger)^T
    real(dp), intent(out) :: T

    ! Here we need a double loop
    integer :: no, noGF
    integer :: i_Elec, ii, isN, in, A_i
    integer :: j_Elec, jj, jsN, jn, A_j
    integer :: o
    integer, pointer :: crows(:)
    complex(dp), pointer :: A(:)

    ! External BLAS routine
    complex(dp), external :: zdotu
    
#ifdef TBTRANS_TIMING
    call timer('A-Gamma',1)
#endif

    ! Get data from tri-diagonal matrix
    crows => cum_rows(A_tri)

    no = El%inDpvt%n
    noGF = 2 * no

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region
    T = 0._dp

    ! Loop columns
    i_Elec = 1
    do while ( i_Elec <= no ) 

      ! We start by creating a region of consecutive memory.
      call consecutive_index(A_tri,El,i_Elec,in,ii)
      isN = nrows_g(A_tri,in)
      
      ! Get starting placement of column in the current block
      ! of the spectral function (zero based)
      if ( in == 1 ) then
        A_i = 2 * (El%inDpvt%r(i_Elec) - 1)
      else
        A_i = 2 * (El%inDpvt%r(i_Elec) - 1) - crows(in-1)
      end if

      if ( ii == noGF ) then

        ! The easy calculation, note that ii == no, only
        ! if the entire electrode sits in one block
        A => val(A_tri,in,in)
        do o = 0 , noGF - 1
          T = T + zdotu(noGF,A((A_i+o)*isN+A_i+1),1,El%Gamma(o*noGF+1),1)
        end do

        ! Quick break of loop
        exit

      end if
      
      ! Loop rows
      j_Elec = 1
      do while ( j_Elec <= no ) 

        ! We start by creating a region of consecutive memory.
        call consecutive_index(A_tri,El,j_Elec,jn,jj)
        jsN = nrows_g(A_tri,jn)

        ! Get the block with the spectral function
        A => val(A_tri,jn,in)

        if ( jn == 1 ) then
          A_j = 2 * El%inDpvt%r(j_Elec) - 1
        else
          A_j = 2 * El%inDpvt%r(j_Elec) - 1 - crows(jn-1)
        end if

        do o = 0 , ii - 1
          T = T + zdotu(jj,A((A_i+o)*jsN+A_j),1, &
          El%Gamma((2*(i_Elec-1)+o)*noGF + 2*j_Elec-1), 1)
        end do

        j_Elec = j_Elec + jj/2

      end do

      i_Elec = i_Elec + ii/2

    end do

#ifdef TBTRANS_TIMING
    call timer('A-Gamma',2)
#endif
    
  end subroutine A_Gamma

  
  ! On entry A_tri is the spectral function
  ! on return the first El%o_inD%n x El%o_inD%n will be the
  ! G.Gamma.Gf.El%Gamma matrix
  ! This will enable eigenvalue calculators and possibly
  ! speed up the calculation of the transmission.
  subroutine A_Gamma_Block(A_tri,El,T,nwork,work)

    use intrinsic_missing, only : transpose, trace
    
    type(zTriMat), intent(inout) :: A_tri ! Spectral function
    type(electrode_t), intent(inout) :: El
    real(dp), intent(out) :: T
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(:)

    ! Here we need a double loop
    integer :: no, no2
    integer :: i_Elec, ii, isN, in, A_i
    integer :: j_Elec, jj, jsN, jn, A_j
    integer, pointer :: crows(:)
    complex(dp), pointer :: A(:)
    complex(dp) :: z
    
#ifdef TBTRANS_TIMING
    call timer('A-Block-Gamma',1)
#endif

    ! Get data from tri-diagonal matrix
    crows => cum_rows(A_tri)

    no = El%inDpvt%n
    no2 = 2 * no
    if ( no ** 2 > nwork ) then
       call die('A_Gamma_Block: Insufficient work-size')
    end if

    ! "sadly" Gamma is saved in transposed form, hence
    ! we transpose, and return it to original form, when returning
    call transpose(no,El%Gamma)

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Loop columns
    i_Elec = 1
    ! The first column calculation initializes the result
    z = z0
    do while ( i_Elec <= no ) 
      
      ! We start by creating a region of consecutive memory.
      call consecutive_index(A_tri,El,i_Elec,in,ii)
      isN = nrows_g(A_tri,in)

      ! Get starting placement of column in the current block
      ! of the spectral function (zero based)
      if ( in == 1 ) then
        A_i = 2 * (El%inDpvt%r(i_Elec) - 1)
      else
        A_i = 2 * (El%inDpvt%r(i_Elec) - 1) - crows(in-1) 
      end if

      if ( ii == no2 ) then
        ! The easy calculation, note that ii == no, only
        ! if the entire electrode sits in one block
        A => val(A_tri,in,in)

        call GEMM ('N','N',no2,no2,no2, z1, A(A_i*(isN+1)+1), isN, &
            El%Gamma(1), no2, z0, work(1), no2)

        ! Quick break of loop
        exit

      end if

      ! Loop rows
      j_Elec = 1
      do while ( j_Elec <= no ) 

        ! We start by creating a region of consecutive memory.
        call consecutive_index(A_tri,El,j_Elec,jn,jj)
        jsN = nrows_g(A_tri,jn)

        ! Get the block with the spectral function
        A => val(A_tri,jn,in)

        if ( jn == 1 ) then
          A_j = 2 * El%inDpvt%r(j_Elec) - 1
        else
          A_j = 2 * El%inDpvt%r(j_Elec) - 1 - crows(jn-1)
        end if

        call GEMM ('N','N',2*jj,no2,2*ii, z1, A(A_i*jsN + A_j), jsN, &
            El%Gamma(2*i_Elec-1), no2, z, work(2*j_Elec-1), no2)

        j_Elec = j_Elec + jj

      end do
       
      i_Elec = i_Elec + ii
      ! Now we have already filled the first entries, sum...
      z = z1
      
    end do
    
    ! Calculate transmission
    T = real( trace(no, work(:)), dp)
    
    ! Now we have the square matrix product
    !   tt = G \Gamma_1 G^\dagger \Gamma_El
    
    call transpose(no,El%Gamma)
    
#ifdef TBTRANS_TIMING
    call timer('A-Block-Gamma',2)
#endif
    
  end subroutine A_Gamma_Block


  subroutine TT_eigen(n,tt,nwork,work,eig)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: tt(:)
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(:)
    complex(dp), intent(inout) :: eig(:)

    real(dp) :: rwork(n*2)
    complex(dp) :: z
    integer :: i, j

#ifdef TBTRANS_TIMING
    call timer('TT-eig',1)
#endif

    ! To remove any singular values we add a 1e-3 to the diagonal
    do i = 1 , n
      tt((i-1)*n+i) = tt((i-1)*n+i) + 1.e-3_dp
    end do

    call zgeev('N','N',n,tt(1),n,eig(1),work(1),1,work(1),1, &
        work(1),nwork,rwork(1),i)
    if ( i /= 0 ) then
      print *,i
      call die('TT_eigen: Could not calculate eigenvalues.')
    end if

    ! Sort the eigenvalues, and simultaneously shift them back
    eig(1) = eig(1) - 1.e-3_dp
    do i = 2 , n
      eig(i) = eig(i) - 1.e-3_dp
      do j = 1 , i - 1
        if ( real(eig(j),dp) < real(eig(i),dp) ) then
          z = eig(j)
          eig(j) = eig(i)
          eig(i) = z
        end if
      end do
    end do

#ifdef TBTRANS_TIMING
    call timer('TT-eig',2)
#endif
    
  end subroutine TT_eigen
  
  subroutine GF_Gamma(Gfcol,El,T)

    use m_ts_trimat_invert, only : TriMat_Bias_idxs

    type(zTriMat), intent(inout) :: Gfcol
    type(electrode_t), intent(inout) :: El
    real(dp), intent(out) :: T

    complex(dp), pointer :: Gf(:)
    complex(dp), pointer :: z(:)

    integer :: no, np
    integer :: i, ii, i_Elec
    integer, pointer :: crows(:)

    integer :: sN, n, nb
    ! BLAS routines
    complex(dp), external :: zdotu, zdotc

#ifdef TBTRANS_TIMING
    call timer('Gf-Gamma',1)
#endif

    no = El%inDpvt%n
    np = parts(Gfcol)
    crows => cum_rows(Gfcol)

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Point to the matrices
    z => val(Gfcol,all=.true.)

    T = 0._dp

    i_Elec = 1
    do while ( i_Elec <= no ) 

      ! We start by creating a region of consecutive memory.
      call consecutive_index(Gfcol,El,i_Elec,n,nb)
      sN = nrows_g(Gfcol,n)

      ! get placement of the diagonal block in the column
      call TriMat_Bias_idxs(Gfcol,no,n,i,ii)

      i = i + El%inDpvt%r(i_Elec) - (crows(n)-sN) - 1
      Gf => z(i:ii)

#ifdef TBT_T_G_GAMMA_OLD
       
      ! Number of columns that we want to do product of
      ii = 1
      do i = 1 , no
        T = T - aimag( zdotu(nb,Gf(ii),1,El%Gamma(i_Elec+(i-1)*no),1) ) ! G \Gamma
        ii = ii + sN
      end do
      ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
      !    Tr[(G \Gamma)^\dagger]
      ! Hence the below calculation shouldn't be necessary
      ii = (i_Elec - 1) * no + 1
      do i = 1 , nb
        T = T + aimag( zdotc(no,Gf(i),sN,El%Gamma(ii),1) )! G^\dagger \Gamma
        ii = ii + no
      end do

#else
      
      ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
      !    Tr[(G \Gamma)^\dagger]
      ! Hence we only calculate one of the contributions and double
      ! it after
      ! Indeed we actually need to calculate:
      !    i Tr[G \Gamma] and since Gamma is not having the i factor
      ! we may take the negative real part.
      ii = 1
      do i = 1 , no
        T = T - aimag( zdotu(nb,Gf(ii),1,El%Gamma(i_Elec+(i-1)*no),1) )! G \Gamma
        ii = ii + sN
      end do
       
#endif

      i_Elec = i_Elec + nb

    end do

    ! Now we have:
    !   T = Tr[G \Gamma - G^\dagger \Gamma]
#ifndef TBT_T_G_GAMMA_OLD
    T = T * 2._dp
#endif

#ifdef TBTRANS_TIMING
    call timer('Gf-Gamma',2)
#endif
    
  end subroutine GF_Gamma

  subroutine GF_T(Gfcol,El,T_Gf,T_self,nzwork,zwork)

    use intrinsic_missing, only: TRACE
    use m_ts_trimat_invert, only : TriMat_Bias_idxs

    type(zTriMat), intent(inout) :: Gfcol
    type(electrode_t), intent(inout) :: El
    real(dp), intent(out) :: T_Gf, T_self
    integer, intent(in) :: nzwork
    complex(dp), intent(inout) :: zwork(:)

    complex(dp), pointer :: z(:)

    integer :: no, np, no2
    integer :: i, ii, i_Elec, itest
    integer, pointer :: crows(:)
    integer :: sN, n
    integer :: nb
    ! BLAS routines
    complex(dp), external :: zdotu

#ifdef TBTRANS_TIMING
    call timer('Gf-T',1)
#endif

    no = El%inDpvt%n
    no2 = 2*El%inDpvt%n
    np = parts(Gfcol)
    crows => cum_rows(Gfcol)

#ifndef TS_NOCHECKS
    if ( no2**2 > nzwork ) call die('GF_T: no2**2 < nzwork')

    ! First we check that we can use the first elements
    ! of Gfcol as temporary storage
    call TriMat_Bias_idxs(Gfcol,no2,1,i,ii)
    if ( i < no2 ** 2 ) then
      write(*,'(a)') 'Remove TBT.T.Gf from your fdf file. &
          &It is not possible in your current setup.'
      call die('GF_T: Size of temporary array not possible.')
    end if

#endif

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Point to the matrices
    z => val(Gfcol,all=.true.)

    i_Elec = 1
    do while ( i_Elec <= no ) 

      ! We start by creating a region of consecutive memory.
      call consecutive_index(Gfcol,El,i_Elec,n,nb)
      sN = nrows_g(Gfcol,n)

      ! get placement of the diagonal block in the column
      call TriMat_Bias_idxs(Gfcol,no2,n,i,ii)

      i = i + 2 * (El%inDpvt%r(i_Elec) - 1) - (crows(n)-sN)

      ! Calculate the G \Gamma+
      call GEMM ('N','T',nb,no2,no2, z1, z(i),sN, &
          El%Gamma(1), no2, z0, zwork(2*(i_Elec-1)+1), no2)
      i_Elec = i_Elec + nb

    end do

    ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
    !    Tr[(G \Gamma)^\dagger]
    ! Now we have:
    !    = G \Gamma
    T_Gf = - aimag( TRACE(no2, zwork) ) * 2._dp

    ! Now we need to correct for the current electrode
    
    ! Now we can calculate the spectral function for this
    ! electrode where it lives
    i_Elec = 1
    do while ( i_Elec <= no ) 
       
      ! We start by creating a region of consecutive memory.
      call consecutive_index(Gfcol,El,i_Elec,n,nb)
      sN = nrows_g(Gfcol,n)

      ! get placement of the diagonal block in the column
      call TriMat_Bias_idxs(Gfcol,no2,n,i,ii)

      i = i + 2 * (El%inDpvt%r(i_Elec) - 1) - (crows(n)-sN)

      ii = 2 * ( i_Elec-1 ) * no2 + 1
      ! Calculate the G \Gamma G^\dagger
      call GEMM ('N','C',no2,nb,no2, z1, zwork(1),no2, &
          z(i), sN, z0, z(ii), no2)
       
      i_Elec = i_Elec + nb

    end do

    ! The remaining calculation is very easy as we simply
    ! need to do the sum of the trace
    T_self = - aimag( zdotu(no2*no2,z(1),1,El%Gamma(1),1) )

#ifdef TBTRANS_TIMING
    call timer('Gf-T',2)
#endif
    
  end subroutine GF_T

  subroutine GF_T_solve(N_Elec,T,has_all)
    integer, intent(in) :: N_Elec
    real(dp), intent(inout) :: T(N_Elec+1,N_Elec)
    ! Whether or not we have calculated all transmissions
    logical, intent(in) :: has_all

    real(dp) :: TT(3)

    select case ( N_Elec ) 
    case ( 1 )

      ! For one electrodes, we simply return immediately
      return

    case ( 2 )

      ! The simple case is when we have 2 electrodes
      
      T(2,1) = T(N_Elec+1,1) - T(1,1)
      if ( has_all ) then
        T(1,2) = T(N_Elec+1,2) - T(2,2)
      else
        T(1,2) = T(2,1)
      end if

      return

    case ( 3 )
      
      if ( .not. has_all ) then
        call die('GF_T_solve: Can not separate transmissions (need all bulk).')
      end if

    case default

      call die('Calculating transmission from underdetermined &
          &system is not allowed. Remove TBT.T.Gf.')
      
    end select

    ! RHS
    TT(1) = T(N_Elec+1,1) - T(1,1)
    TT(2) = T(N_Elec+1,2) - T(2,2)
    TT(3) = T(N_Elec+1,3) - T(3,3)

    ! Calculate them exactly
    ! We do not need LAPACK here as this can
    ! only be definitely solved for N_Elec == 3
    T(2,1) = (TT(1) - TT(3) + TT(2)) * 0.5_dp
    T(1,2) = T(2,1)
    T(3,1) = TT(1) - T(2,1)
    T(1,3) = T(3,1)
    T(3,2) = TT(3) - T(3,1)
    T(2,3) = T(3,2)

  end subroutine GF_T_solve

  subroutine consecutive_index(Tri,El,current,p,n)
    type(zTriMat), intent(inout) :: Tri
    type(electrode_t), intent(in) :: El
    integer, intent(in) :: current
    integer, intent(out) :: p, n

    ! Local variables
    integer :: idx_Elec, i, sIdx, eIdx
    integer, pointer :: crows(:)
    
    idx_Elec = El%inDpvt%r(current)
    p = which_part(Tri,2*idx_Elec)
    crows => cum_rows(Tri)
    eIdx = crows(p)
    sIdx = eIdx - nrows_g(Tri,p) + 1

    n = 1
    do while ( current + n <= El%inDpvt%n )
      i = El%inDpvt%r(current+n)
      ! In case it is not consecutive
      if ( i - idx_Elec /= n ) exit
      ! In case the block changes, then
      ! we cut the block size here.
      if ( 2 * i - 1 < sIdx .or. eIdx < 2 * i ) exit
      n = n + 1
    end do
    n = 2 * n

  end subroutine consecutive_index


#ifdef NCDF_4
  subroutine orb_current(sp,H,S,sc_off,k,ph,cE,A_tri,r,orb_J,pvt)

    use class_Sparsity
    use class_dSpData2D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    use m_ts_cctype, only: ts_c_idx

    type(Sparsity), intent(inout) :: sp
    ! We require that the input Hamiltonian is Hermitian
    real(dp), intent(in) :: H(:,:), S(:), sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(ts_c_idx) :: cE
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of orb_J
    type(tRgn), intent(in) :: r
    type(dSpData2D), intent(inout) :: orb_J
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: i_sp
    integer, pointer :: i_ncol(:), i_ptr(:), i_col(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp), pointer :: A(:)
    complex(dp) :: Hi(2,2)
    real(dp), pointer :: J(:,:)
    real(dp) :: E
    integer :: no_u, iu, io, ind, iind, ju, jo, iA(2,2), ii, jj, iph
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('orb-current',1)
#endif

    ! Retrieve energy
    E = real(cE%e,dp)

    call attach(sp, nrows_g=no_u, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    i_sp => spar(orb_J)
    J => val (orb_J)
    call attach(i_sp, n_col=i_ncol, list_ptr=i_ptr, list_col=i_col)

    ! Create the phases
    ! We are using the symmetric H(j, i) = H(i, j) relation.
    ! So since we are taking the complex part on the first entry we retrieve the H(j,i) (in k-space)
    ! component.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp))
    end do

    A => val(A_tri)

    ! we need this in case the device region gets enlarged due to dH
    J(:,:) = 0._dp

!$OMP parallel do default(shared), private(iu,io,ju,jo,iind,ind,Hi,iA,ii,jj,iph)
    do iu = 1, r%n
      io = r%r(iu)

#ifndef TS_NOCHECKS
      if ( i_ncol(io) == 0 ) call die('orb_current: J has zero columns &
          &for at least one row')
#endif

      ! Loop on the orbital current indices
      do iind = i_ptr(io) + 1, i_ptr(io) + i_ncol(io)

        ! Here we will calculate the orbital current from dH
        ! onto orbital:
        !  J(iind) == J(io, jo)

        ! Get jo orbital
        jo = ucorb(i_col(iind), no_u)
        ju = pvt%r(jo) ! pivoted orbital index in tri-diagonal matrix
        
        ! Check if the jo, io orbital exists in H
        if ( l_ncol(jo) > 0 ) then
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))
          ! Get transpose element
          jj = TO(i_col(iind)) + io
          ind = l_ptr(jo) + SFIND(col, jj)

          if ( ind > l_ptr(jo) ) then

            ! Check for the Hamiltonian element H_ji
            iph = (l_col(ind)-1)/no_u
            Hi(1,1) = (cmplx(H(ind,1), H(ind,5), dp) - E * S(ind)) * ph(iph) ! Hji^uu
            Hi(1,2) = (cmplx(H(ind,3),-H(ind,4), dp)             ) * ph(iph) ! Hji^ud
            Hi(2,1) = (cmplx(H(ind,7), H(ind,8), dp)             ) * ph(iph) ! Hji^du
            Hi(2,2) = (cmplx(H(ind,2), H(ind,6), dp) - E * S(ind)) * ph(iph) ! Hji^dd
            
            do ii = 1, 2
              do jj = 1, 2
                iA(ii,jj) = index(A_tri, 2*(iu-1)+ii, 2*(ju-1)+jj)  ! Aij^{ii,jj}
              end do 
            end do

            ! Jij                      Aij   * Hji
            J(iind,1) = J(iind,1) + aimag( A(iA(1,1)) * Hi(1,1) )
            J(iind,2) = J(iind,2) + aimag( A(iA(2,2)) * Hi(2,2) )
            J(iind,3) = J(iind,3) + aimag( A(iA(1,2)) * Hi(2,1) )
            J(iind,4) = J(iind,4) +  real( A(iA(1,2)) * Hi(2,1) )
            J(iind,5) = J(iind,5) -  real( A(iA(1,1)) * Hi(1,1) ,dp)
            J(iind,6) = J(iind,6) -  real( A(iA(2,2)) * Hi(2,2) , dp)
            J(iind,7) = J(iind,7) + aimag( A(iA(1,2)) * Hi(2,1) )
            J(iind,8) = J(iind,8) -  real( A(iA(1,2)) * Hi(2,1) )

          end if
        end if

        ! Check if the io, jo orbital exists in H
        if ( l_ncol(io) > 0 ) then
          col => l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io))
          ind = l_ptr(io) + SFIND(col, i_col(iind))
          if ( ind > l_ptr(io) ) then

            ! Check for the Hamiltonian element H_ji
            iph = (l_col(ind)-1)/no_u
            Hi(1,1) = (cmplx(H(ind,1), H(ind,5), dp) - E * S(ind)) * ph(iph) ! Hij^uu
            Hi(1,2) = (cmplx(H(ind,3),-H(ind,4), dp)             ) * ph(iph) ! Hij^ud
            Hi(2,1) = (cmplx(H(ind,7), H(ind,8), dp)             ) * ph(iph) ! Hij^du
            Hi(2,2) = (cmplx(H(ind,2), H(ind,6), dp) - E * S(ind)) * ph(iph) ! Hij^dd
            
            do ii = 1, 2
              do jj = 1, 2
                iA(ii,jj) = index(A_tri, 2*(ju-1)+jj, 2*(iu-1)+ii) ! Aji^{jj,ii}
              end do 
            end do

            ! Jij                      Aij   * Hji
            J(iind,1) = J(iind,1) - aimag( A(iA(1,1)) * Hi(1,1) )
            J(iind,2) = J(iind,2) - aimag( A(iA(2,2)) * Hi(2,2) )
            J(iind,3) = J(iind,3) - aimag( A(iA(1,2)) * Hi(2,1) )
            J(iind,4) = J(iind,4) -  real( A(iA(1,2)) * Hi(2,1) )
            J(iind,5) = J(iind,5) +  real( A(iA(1,1)) * Hi(1,1) ,dp)
            J(iind,6) = J(iind,6) +  real( A(iA(2,2)) * Hi(2,2) , dp)
            J(iind,7) = J(iind,7) - aimag( A(iA(1,2)) * Hi(2,1) )
            J(iind,8) = J(iind,8) +  real( A(iA(1,2)) * Hi(2,1) )


          end if
        end if

      end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('orb-current',2)
#endif

  contains

    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i

      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)

        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('orb_current: could not find transpose supercell index')
      
    end function TO

  end subroutine orb_current
  
  subroutine orb_current_add_dH(dH_1D,sc_off,k,ph,A_tri,r,orb_J,pvt)

    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    type(zSpData1D), intent(in) :: dH_1D
    real(dp), intent(in) :: sc_off(:,:), k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of orb_J
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: orb_J
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: dH(:)
    type(Sparsity), pointer :: i_sp
    integer, pointer :: i_ncol(:), i_ptr(:), i_col(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp), pointer :: A(:)
    real(dp), pointer :: J(:)
    integer :: no_u, iu, io, ind, iind, ju, jo, jj

#ifdef TBTRANS_TIMING
    call timer('orb-current-dH',1)
#endif

    ! Retrieve dH
    sp => spar(dH_1D)
    dH => val (dH_1D)
    call attach(sp, nrows_g=no_u, &
        n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    i_sp => spar(orb_J)
    J => val (orb_J)
    call attach(i_sp, n_col=i_ncol, list_ptr=i_ptr, list_col=i_col)

    ! Create the phases
    ! We are using the explicit H(j, i) and thus the phases are consistent with +
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp))
    end do

    A => val(A_tri)

!$OMP parallel do default(shared), &
!$OMP&private(iu,io,iind,jo,ju,ind,col,jj)
    do iu = 1, r%n
      io = r%r(iu)

      ! Loop on the orbital current indices
      do iind = i_ptr(io) + 1, i_ptr(io) + i_ncol(io)

        ! Here we will calculate the orbital current from dH
        ! onto orbital:
        !  J(iind) == J(io, jo)

        ! Get jo orbital
        jo = ucorb(i_col(iind), no_u)
        ju = pvt%r(jo) ! pivoted orbital index in tri-diagonal matrix
        
        ! Check if the jo, io orbital exists in dH
        if ( l_ncol(jo) > 0 ) then
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))
          ! Get transpose element
          jj = TO(i_col(iind)) + io
          ind = l_ptr(jo) + SFIND(col, jj)

          if ( ind > l_ptr(jo) ) then

            ! Check for the Hamiltonian element H_ji
            jj = index(A_tri,iu,ju) ! A_ij

            ! Jij                      Aij   * Hji
            J(iind) = J(iind) + aimag( A(jj) * dH(ind) * ph( (l_col(ind)-1)/no_u ) )

          end if
        end if

        ! Check if the io, jo orbital exists in dH
        if ( l_ncol(io) > 0 ) then
          col => l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io))
          ind = l_ptr(io) + SFIND(col, i_col(iind))
          if ( ind > l_ptr(io) ) then

            ! Check for the Hamiltonian element H_ij
            jj = index(A_tri,ju,iu) ! A_ji

            ! Jij -=                   Aji   * Hij
            J(iind) = J(iind) - aimag( A(jj) * dH(ind) * ph( (l_col(ind)-1)/no_u ) )

          end if
        end if

      end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('orb-current-dH',2)
#endif

  contains

    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i

      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)

        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('orb_current_add_dH: could not find transpose supercell index')
      
    end function TO
    
  end subroutine orb_current_add_dH


  subroutine GF_DM(sc_off,k,ph,Gfd_tri,Gfo_tri,r,pvt,spDM)

    use class_Sparsity
    use class_dSpData2D

    use m_spin, only: spin
    
    use geom_helper,       only : UCORB

    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    ! The region that specifies the size of spDM
    type(tRgn), intent(in) :: r
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt
    type(dSpData2D), intent(inout) :: spDM

    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    type(Sparsity), pointer :: sp
    real(dp), pointer :: DM(:,:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd(2,2), kph

    integer :: no_u, iu, io, ind, ju

#ifdef TBTRANS_TIMING
    call timer('Gf-DM',1)
#endif

    sp => spar(spDM)
    DM => val(spDM)
    call attach(sp, nrows_g=no_u, n_col=ncol, list_ptr=l_ptr, list_col=l_col)

    ! Create the phases
    ! Since we have to do Gf.exp(ikR) we simply
    ! create exp(-ikR) for the supercell connections.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

    ! we need this in case the device region gets enlarged due to dH
    DM(:,:) = 0._dp

    if ( spin%NCol ) then
    !$OMP parallel do default(shared), private(iu,io,ind,ju,GfGfd,kph)
    do iu = 1, r%n
      io = r%r(iu)

#ifndef TS_NOCHECKS
      if ( ncol(io) == 0 ) call die('Gf_DM: DM has zero columns &
          &for at least one row')
#endif

      ! Loop on DM entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

        kph = ph((l_col(ind) - 1) / no_u)
        ju = pvt%r(ucorb(l_col(ind), no_u))
        call calc_GfGfd(iu, ju, GfGfd)
        DM(ind, 1) = - aimag(           GfGfd(1,1)             * kph )
        DM(ind, 2) = - aimag(           GfGfd(2,2)             * kph )
        DM(ind, 3) = - aimag( 0.5_dp * (GfGfd(1,2)+GfGfd(2,1)) * kph )
        DM(ind, 4) = -  real( 0.5_dp * (GfGfd(1,2)+GfGfd(2,1)) * kph,dp )

      end do
    end do
!$OMP end parallel do
    else if ( spin%SO ) then

!$OMP parallel do default(shared), private(iu,io,ind,ju,GfGfd,kph)
      do iu = 1, r%n
        io = r%r(iu)

#ifndef TS_NOCHECKS
        if ( ncol(io) == 0 ) call die('Gf_DM: DM has zero columns &
            &for at least one row')
#endif

        ! Loop on DM entries here...
        do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)
          
          kph = ph((l_col(ind) - 1) / no_u)
          ju = pvt%r(ucorb(l_col(ind), no_u))
          call calc_GfGfd(iu, ju, GfGfd)
          DM(ind, 1) = - aimag( GfGfd(1,1) * kph )
          DM(ind, 2) = - aimag( GfGfd(2,2) * kph )
          DM(ind, 3) = - aimag( GfGfd(1,2) * kph )
          DM(ind, 4) = -  real( GfGfd(1,2) * kph,dp )
          DM(ind, 5) =    real( GfGfd(1,1) * kph,dp )
          DM(ind, 6) =    real( GfGfd(2,2) * kph,dp )
          DM(ind, 7) = - aimag( GfGfd(2,1) * kph )
          DM(ind, 8) =    real( GfGfd(2,1) * kph,dp )

        end do
      end do
!$OMP end parallel do
    else
      call die("tbt_tri_scat_spinor: Spin has to be non-collinear or spin-orbit")
    endif 

#ifdef TBTRANS_TIMING
    call timer('Gf-DM',2)
#endif
    
  contains
    
    pure subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G(2,2)
      integer :: p_r, i_r, j_r, N_r, p_c, i_c, j_c, N_c, i

      call part_index(Gfo_tri, 2*br-1, p_r, i_r)
      call part_index(Gfo_tri, 2*bc-1, p_c, i_c)
      N_r = Gfo_tri%data%tri_nrows(p_r)
      N_c = Gfo_tri%data%tri_nrows(p_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            G(j_r,j_c) = Gfd(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfd(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      else
        ! Loop over spin components
        do j_r = 1, 2
          do j_c = 1, 2
            i = index_block(Gfo_tri, p_r, p_c)
            G(j_r,j_c) = Gfo(i + (i_r+j_r-1) + (i_c+j_c-2)* N_r)
            i = index_block(Gfo_tri, p_c, p_r)
            G(j_r,j_c) = G(j_r,j_c) - conjg(Gfo(i + (i_c+j_c-1) + (i_r+j_r-2)* N_c))
          end do
        end do
      end if

    end subroutine calc_GfGfd
    
  end subroutine GF_DM

  subroutine A_DM(sc_off,k,ph,A_tri,r,pvt,spDM)

    use class_Sparsity
    use class_dSpData2D
    use geom_helper,       only : UCORB

    use m_spin, only: spin

    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of spDM
    type(tRgn), intent(in) :: r
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt
    type(dSpData2D), intent(inout) :: spDM

    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: A(:)
    real(dp), pointer :: DM(:,:)
    integer :: no_u, iu, io, ind, ju, ii, jj, idx(2,2)
    complex(dp) :: kph

#ifdef TBTRANS_TIMING
    call timer('A-DM',1)
#endif

    sp => spar(spDM)
    DM => val(spDM)
    call attach(sp, nrows_g=no_u, n_col=ncol, list_ptr=l_ptr, list_col=l_col)

    ! Create the phases
    ! Since we have to do Gf.exp(ikR) we simply
    ! create exp(-ikR) for the supercell connections.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,io)), dp)) / (Pi * 2._dp)
    end do

    A => val(A_tri)

    ! we need this in case the device region gets enlarged due to dH
    DM(:,:) = 0._dp

    if ( spin%NCol ) then

!$OMP parallel do default(shared), private(iu,io,ind,ju,ii,jj,idx,kph)
      do iu = 1, r%n
        io = r%r(iu)

#ifndef TS_NOCHECKS
        if ( ncol(io) == 0 ) call die('A_DM: DM has zero columns &
            &for at least one row')
#endif

        ! Loop on DM entries here...
        do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

          ju = pvt%r(ucorb(l_col(ind), no_u))
          do ii = 1, 2
            do jj = 1, 2
              idx(ii,jj) = index(A_tri, 2*(iu-1)+ii, 2*(ju-1)+jj)
            end do
          end do
          kph = ph((l_col(ind) - 1) / no_u)
          DM(ind,1) =   real(          A(idx(1,1))                * kph, dp)
          DM(ind,2) =   real(          A(idx(2,2))                * kph, dp)
          DM(ind,3) =   real(0.5_dp * (A(idx(1,2)) + A(idx(2,1))) * kph, dp)
          DM(ind,4) = -aimag(0.5_dp * (A(idx(1,2)) + A(idx(2,1))) * kph)

        end do
      end do
!$OMP end parallel do

    else ! spin%SO

!$OMP parallel do default(shared), private(iu,io,ind,ju,ii,jj,idx,kph)
      do iu = 1, r%n
        io = r%r(iu)

#ifndef TS_NOCHECKS
        if ( ncol(io) == 0 ) call die('A_DM: DM has zero columns &
            &for at least one row')
#endif

        ! Loop on DM entries here...
        do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

          ju = pvt%r(ucorb(l_col(ind), no_u))
          do ii = 1, 2
            do jj = 1, 2
              idx(ii,jj) = index(A_tri, 2*(iu-1)+ii, 2*(ju-1)+jj)
            end do
          end do
          kph = ph((l_col(ind) - 1) / no_u)
          DM(ind,1) =   real(A(idx(1,1)) * kph, dp)
          DM(ind,2) =   real(A(idx(2,2)) * kph, dp)
          DM(ind,3) =   real(A(idx(1,2)) * kph, dp)
          DM(ind,4) = -aimag(A(idx(1,2)) * kph)
          DM(ind,5) =  aimag(A(idx(1,1)) * kph)
          DM(ind,6) =  aimag(A(idx(2,2)) * kph)
          DM(ind,7) =   real(A(idx(2,1)) * kph, dp)
          DM(ind,8) =  aimag(A(idx(2,1)) * kph)

        end do
      end do
!$OMP end parallel do

    end if

#ifdef TBTRANS_TIMING
    call timer('A-DM',2)
#endif

  end subroutine A_DM

#endif

  subroutine insert_Self_energy(n1,n2,M,r,El,off1,off2)

    ! The sizes of the matrix
    integer, intent(in) :: n1, n2
    complex(dp), intent(inout) :: M(:,:,:,:)
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    ! Electrodes...
    type(electrode_t), intent(inout) :: El
    ! The offsets of the matrix
    integer, intent(in) :: off1, off2

    ! local variables
    integer :: j, je, i, ie, no, idx

    idx = El%idx_o - 1
    no = El%device_orbitals()

    ! We are dealing with the intrinsic electrode
    ! self energy
    ! Here we have two options,
    ! Bulk) We are dealing with a bulk electrode
    ! not bulk) A non-bulk electrode

    if ( El%Bulk ) then
!$OMP do private(j,je,i,ie)
      do j = 1 , n2
        je = r%r(off2+j) - idx
        if ( 1 <= je .and. je <= no ) then
          je = 2 * (je - 1)
          do i = 1 , n1
            ie = r%r(off1+i) - idx
            if ( 1 <= ie .and. ie <= no ) then
              ie = 2 * ie - 1
                
              M(1,i,1,j) = El%Sigma(2*no*(je  )+ ie  )
              M(2,i,1,j) = El%Sigma(2*no*(je  )+ ie+1)
              M(1,i,2,j) = El%Sigma(2*no*(je+1)+ ie  )
              M(2,i,2,j) = El%Sigma(2*no*(je+1)+ ie+1)
            end if
          end do
        end if
      end do
!$OMP end do
    else
!$OMP do private(j,je,i,ie)
      do j = 1 , n2
        je = r%r(off2+j) - idx
        if ( 1 <= je .and. je <= no ) then
          je = 2 * (je - 1)
          do i = 1 , n1
            ie = r%r(off1+i) - idx
            if ( 1 <= ie .and. ie <= no ) then
              ie = 2 * ie - 1
    
              M(1,i,1,j) = M(1,i,1,j) - El%Sigma(2*no*(je  )+ ie  )
              M(2,i,1,j) = M(2,i,1,j) - El%Sigma(2*no*(je  )+ ie+1)
              M(1,i,2,j) = M(1,i,2,j) - El%Sigma(2*no*(je+1)+ ie  )
              M(2,i,2,j) = M(2,i,2,j) - El%Sigma(2*no*(je+1)+ ie+1)
            end if
          end do
        end if
      end do
!$OMP end do
    end if

  end subroutine insert_Self_energy


  subroutine insert_Self_energy_Dev(Gfinv_tri,Gfinv,r,El)

    type(zTriMat), intent(inout) :: GFinv_tri
    complex(dp), intent(inout) :: Gfinv(:)
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    type(electrode_t), intent(in) :: El

    ! local variables
    integer :: j, jj, je, i, ie,ii, idx, idxe, no

    no = El%o_inD%n

    ! A down-folded self-energy, this
    ! is always considered to be "non-bulk" as
    ! we have it downfolded.

!$OMP do private(j,jj,je,i,ii,ie,idx,idxe)
    do j = 1 , no
      ! grab the index in the full tri-diagonal matrix
      je = El%inDpvt%r(j)
      do jj = 1, 2
        do i = 1 , no
          do ii = 1, 2
            ie = El%inDpvt%r(i)
            
            idx = index(GFinv_tri, 2*(ie-1) + ii, 2*(je-1) + jj )
            
            idxe = 2*no*( 2*(j-1) + jj-1 ) +  2*(i-1) + ii
            Gfinv(idx) = Gfinv(idx) - El%Sigma(idxe)
          end do
        end do
        
        
      end do
    end do
!$OMP end do

  end subroutine insert_Self_energy_Dev


#undef GEMM

end module m_tbt_tri_scat_spinor
