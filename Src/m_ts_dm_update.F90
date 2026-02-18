!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************


! In general many of the loops below can be followed by examining this loop:

!    ! This loop is across the local rows...
!    do lio = 1 , lnr
!
!       ! Quickly go past the empty regions... (we have nothing to update)
!       if ( l_ncol(lio) == 0 ) cycle
!
!       ! obtain the global index of the local orbital.
!       io = index_local_to_global(dit,lio,Node)
!
!       ! Quickly go past the empty regions... (we have nothing to update)
!       if ( up_ncol(io) == 0 ) cycle
!
!       ! Do a loop in the local sparsity pattern...
!       ! The local sparsity pattern is more "spread", hence
!       ! we do fewer operations by having this as an outer loop
!       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
!
!          jo = UCORB(l_col(lind),nr)
!
!          ! Now search the update region
!          ! This one must *per definition* have less elements.
!          ! Hence, we can exploit this, and find equivalent
!          ! super-cell orbitals.
!          rind = up_ptr(io)
!          ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
!          if ( ind <= rind ) cycle ! The element does not exist
!
!          ! Obtain the phase for the orbital ij
!          kx = k(1) * xij(1,pnt(lind)) + &
!               k(2) * xij(2,pnt(lind)) + &
!               k(3) * xij(3,pnt(lind))
!
!          ph = fact * exp(cmplx(0._dp,-kx,dp))
!
!          ! The integration is:
!          ! \rho = e^{-i.k.R} [ \int (Gf^R-Gf^A)/2 dE + \int Gf^R\Gamma Gf^A dE ]
!          ! the passed array for the equilibrium part is Gf^R - Gf^A
!
!          if ( non_Eq ) then
!
!             ! The integration is:
!             ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
!             kx = real( ph*zDu(ind) , dp)
!             dD(lind) = dD(lind) + kx
!             dE(lind) = dE(lind) + real( ph*zEu(ind) ,dp)
!
!          else
!
!             ! The integration is:
!             ! \rho = e^{-i.k.R} \int (Gf^R-Gf^A)/2 dE
!             ! Gf^R - Gf^A should be passed from outside
!             
!             dD(lind) = dD(lind) + aimag( ph*zDu(ind) )
!             dE(lind) = dE(lind) + aimag( ph*zEu(ind) )
!
!          end if
!             
!       end do
!    end do
    

module m_ts_dm_update

  use precision, only : dp
  use geom_helper, only : UCORB
  use intrinsic_missing, only : SFIND

  implicit none

  private

  public :: init_DM
  interface init_DM
    module procedure init_DM_1D
    module procedure init_DM_2D
  end interface init_DM

  public :: update_DM
  interface update_DM
    module procedure update_DM_sp2D
    module procedure update_DM_sp3D
  end interface update_DM

  public :: update_zDM
  interface update_zDM
    module procedure update_zDM_sp2D
    module procedure update_zDM_sp3D
  end interface update_zDM

  public :: add_Gamma_DM
  interface add_Gamma_DM
    module procedure add_Gamma_DM_2D
    module procedure add_Gamma_DM_3D
  end interface add_Gamma_DM

  public :: add_k_DM
  interface add_k_DM
    module procedure add_k_DM_2D
    module procedure add_k_DM_3D
  end interface add_k_DM

contains

  ! ***
  ! The following scheme should be followed:
  !   add_*_DM routines are constructed to be able to handle different schemes
  !   They are called after each k-point and thus the arguments are the following:
  !     1. the local update sparsity pattern
  !     2. the global update sparsity pattern (suffix 'u')
  ! The k-point routine is constructed to handle three different methods of doing
  ! the weighting which is performed here.
  ! The routines require the integral part of the input
  !  For eq. this is G^R-G^A
  !  For neq. this is G.Gamma.G^A

  subroutine add_k_DM_2D(spDM,spuDM,D_dim2, spEDM, spuEDM, E_dim2, &
       n_s,sc_off,k, non_Eq)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_zSpData2D

! *********************
! * INPUT variables   *
! *********************
    ! The local integrated sparsity arrays
    type(dSpData2D), intent(inout) :: spDM, spEDM
    ! The current k-point global sparsity arrays
    type(zSpData2D), intent(inout) :: spuDM, spuEDM
    ! current update region of last dimension
    integer, intent(in) :: D_dim2, E_dim2
    ! The k-point
    real(dp), intent(in) :: k(3)
    ! The supercell offsets
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    logical, intent(in) :: non_Eq

    ! Arrays needed for looping the sparsity
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:), upp_col(:)
    real(dp), pointer :: dD(:,:) , dE(:,:)
    complex(dp), pointer :: zDu(:,:), zEu(:,:)
    integer :: lnr, lio, lind, io, ind, nr, jo, rind
    logical :: hasEDM
    complex(dp) :: ph(0:n_s-1)

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spuDM)) ) return

    hasEDM = initialized(spEDM)

    ! get the distribution
    dit => dist(spDM)

    l_s => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD => val(spDM)
    if ( hasEDM ) dE => val(spEDM)

    up_s => spar(spuDM)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    zDu => val(spuDM)
    if ( hasEDM ) zEu => val(spuEDM)

    if ( size(zDu,2) < D_dim2 .or. size(dD,2) < D_dim2 ) then
       call die('add_k_DM: Error in code')
    end if

    if ( hasEDM ) then
       if ( size(zEu,2) < E_dim2 .or. size(dE,2) < E_dim2 ) then
          call die('add_k_DM: Error in code')
       end if
    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected k-DM.')

    do jo = 0 , n_s - 1
       ph(jo) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,jo)),dp))
    end do

    if ( non_Eq ) then

     if ( hasEDM ) then
          
! No data race will occur, sparsity pattern only tranversed once
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)
             
             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
             
               jo = (l_col(lind)-1) / nr
             
               ! The integration is:
               ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
               dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + real( ph(jo)*zDu(ind,1:D_dim2) ,dp)
               dE(lind,1:E_dim2) = dE(lind,1:E_dim2) + real( ph(jo)*zEu(ind,1:E_dim2) ,dp)
             end if
             
          end do

          end if
          end if

       end do
!$OMP end parallel do

     else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)
             
             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
             
               jo = (l_col(lind)-1) / nr
             
               dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + real( ph(jo)*zDu(ind,1:D_dim2) ,dp)
             end if

          end do

          end if
          end if

       end do
!$OMP end parallel do

      end if
    
    else ! non_eq == .false.

     if ( hasEDM ) then
          
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)

             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
             
               jo = (l_col(lind)-1) / nr

               ! The integration is
               ! \rho = -\Im e^{-i.k.R} \int (Gf^R-Gf^A) dE
               ! Gf^R-Gf^A from outside
               dD(lind,1:D_dim2) = dD(lind,1:D_dim2) - aimag( ph(jo)*zDu(ind,1:D_dim2) )
               dE(lind,1:E_dim2) = dE(lind,1:E_dim2) - aimag( ph(jo)*zEu(ind,1:E_dim2) )
             end if

          end do

          end if
          end if

       end do
!$OMP end parallel do

    else
          
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)

             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
             
               jo = (l_col(lind)-1) / nr

               dD(lind,1:D_dim2) = dD(lind,1:D_dim2) - aimag( ph(jo)*zDu(ind,1:D_dim2) )
             end if

          end do

          end if
          end if

       end do
!$OMP end parallel do

     end if
    
    end if

  end subroutine add_k_DM_2D

  subroutine add_k_DM_3D(spDM,spuDM,D_dim2, D_dim3, spEDM, spuEDM, E_dim2, E_dim3,&
       n_s,sc_off,k, non_Eq)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData3D
    use class_zSpData3D
    use m_spin, only: spin
#ifdef TS_SOC_DEBUG
    use parallel,only: Node
#endif

! *********************
! * INPUT variables   *
! *********************
    ! The local integrated sparsity arrays
    type(dSpData3D), intent(inout) :: spDM, spEDM
    ! The current k-point global sparsity arrays
    type(zSpData3D), intent(inout) :: spuDM, spuEDM
    ! current update region of last dimension
    integer, intent(in) :: D_dim2, D_dim3, E_dim2, E_dim3
    ! The k-point
    real(dp), intent(in) :: k(3)
    ! The supercell offsets
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    logical, intent(in) :: non_Eq

    ! Arrays needed for looping the sparsity
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:), upp_col(:)
    real(dp), pointer :: dD(:,:,:) , dE(:,:,:)
    complex(dp), pointer :: zDu(:,:,:), zEu(:,:,:)
    integer :: lnr, lio, lind, io, ind, nr, jo, rind
    logical :: hasEDM
    complex(dp) :: ph(0:n_s-1)
    complex(dp) :: zDtemp(1:D_dim3), zEtemp(1:E_dim3)

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spuDM)) ) return

    hasEDM = initialized(spEDM)

    ! get the distribution
    dit => dist(spDM)

    l_s => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD => val(spDM)
    if ( hasEDM ) dE => val(spEDM)

    up_s => spar(spuDM)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    zDu => val(spuDM)
    if ( hasEDM ) zEu => val(spuEDM)

    ! zDu is complex, threfore the first dimension should always be 4 (never 8)
    if ( size(zDu,1) < 4 .or. size(dD,1) < D_dim2 ) then
       call die('add_k_DM_3D: Error in code')
    end if
    if ( size(zDu,3) < D_dim3 .or. size(dD,3) < D_dim3 ) then
       call die('add_k_DM_3D: Error in code')
    end if

    if ( hasEDM ) then
       ! zEu is complex, threfore the first dimension should always be 4 (never 8)
       if ( size(zEu,1) < 4 .or. size(dE,1) < E_dim2 ) then
          call die('add_k_DM_3D: Error in code')
       end if
       if ( size(zEu,3) < E_dim3 .or. size(dE,3) < E_dim3 ) then
          call die('add_k_DM_3D: Error in code')
       end if
    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.

    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected k-DM.')

    do jo = 0 , n_s - 1
       ph(jo) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,jo)),dp))
    end do

    if (spin%NCol) then
      call add_k_DM_3D_nc
    else if (spin%SO) then
      call add_k_DM_3D_so
    else
      call die('add_k_DM_3D: Error in code. This function should only be &
                &called during non-collinear or spin-orbit calculations.')
    end if

    contains

    subroutine add_k_DM_3D_nc

      if ( non_Eq ) then

        if ( hasEDM ) then

! No data race will occur, sparsity pattern only tranversed once
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp,zEtemp)
          do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

                jo = UCORB(l_col(lind),nr)

                ind = rind + SFIND(upp_col,jo)
                if ( rind < ind ) then

                  jo = (l_col(lind)-1) / nr

                  ! The integration is:
                  ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
                  !
                  !
                  ! After adding the phase factor we have to separate the complex and
                  ! real parts of the DM elements and re-order them to the usual
                  ! convention. We discard the imaginary parts on the spin-box 
                  ! diagonal and one of the off-diagonal spin-box elements:
                  !
                  !        | dD(1,:)                dD(3,:) - i dD(4,:) |
                  !   DM = |                                            |
                  !        | dD(3,:) + i dD(4,:)    dD(2,:)             |
                  !
                  !        | Re(zDu(1,:)*ph)        zDu(2,:)*ph         |
                  !      = |                                            |
                  !        | zDu(3,:)*ph            Re(zDu(4,:)*ph)     |
                  !
                  ! Same for the EDM.

                  ! UP, UP
                  zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                  dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) + real( zDtemp ,dp)
                  ! DOWN, DOWN
                  zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                  dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) + real( zDtemp ,dp)
                  ! UP, DOWN
                  zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                  dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) + real( zDtemp ,dp)
                  dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - aimag( zDtemp)

                  ! UP, UP
                  zEtemp = ph(jo)*zEu(1,ind,1:E_dim3)
                  dE(1,lind,1:E_dim3) = dE(1,lind,1:E_dim3) + real( zEtemp ,dp)
                  ! DOWN, DOWN
                  zEtemp = ph(jo)*zEu(4,ind,1:E_dim3)
                  dE(2,lind,1:E_dim3) = dE(2,lind,1:E_dim3) + real( zEtemp ,dp)
                  ! UP, DOWN
                  zEtemp = ph(jo)*zEu(2,ind,1:E_dim3)
                  dE(3,lind,1:E_dim3) = dE(3,lind,1:E_dim3) + real( zEtemp ,dp)
                  dE(4,lind,1:E_dim3) = dE(4,lind,1:E_dim3) - aimag( zEtemp)
                end if
#ifdef TS_SOC_DEBUG
        if ((sc_off(1,jo) .eq. 0.) .and. (sc_off(2,jo) .eq. 0.) .and. &
            (sc_off(3,jo) .eq. 0.)) then
            write(1340+Node, "(4I10)") io, UCORB(lio, nr), UCORB(l_col(lind),nr)
            write(1440+Node, "(8(spE26.15))") dD(:,lind,1)
            write(1540+Node, "(8(spE26.15))") dD(:,lind,2)
        end if
#endif

            end do

            end if
            end if

          end do
!$OMP end parallel do

        else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp)
          do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

                jo = UCORB(l_col(lind),nr)

                ind = rind + SFIND(upp_col,jo)
                if ( rind < ind ) then

                  jo = (l_col(lind)-1) / nr
                  !
                  ! After adding the phase factor we have to separate the complex and
                  ! real parts of the DM elements and re-order them to the usual
                  ! convention. We discard the imaginary parts on the spin-box 
                  ! diagonal and one of the off-diagonal spin-box elements:
                  !
                  !        | dD(1,:)                dD(3,:) - i dD(4,:) |
                  !   DM = |                                            |
                  !        | dD(3,:) + i dD(4,:)    dD(2,:)             |
                  !
                  !        | Re(zDu(1,:)*ph)        zDu(2,:)*ph         |
                  !      = |                                            |
                  !        | zDu(3,:)*ph            Re(zDu(4,:)*ph)     |
                  !

                  ! UP, UP
                  zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                  dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) + real( zDtemp ,dp)
                  ! DOWN, DOWN
                  zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                  dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) + real( zDtemp ,dp)
                  ! UP, DOWN
                  zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                  dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) + real( zDtemp ,dp)
                  dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - aimag( zDtemp)
                end if

            end do

            end if
            end if

          end do
!$OMP end parallel do

        end if

      else ! non_eq == .false.

        if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp,zEtemp)
          do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

                jo = UCORB(l_col(lind),nr)

                ind = rind + SFIND(upp_col,jo)
                if ( rind < ind ) then

                  jo = (l_col(lind)-1) / nr

                  ! The integration is
                  ! \rho = i e^{-i.k.R} \int (Gf^R-Gf^A) dE
                  ! After adding the phase factor we have to separate the complex and
                  ! real parts of the DM elements and re-order them to the usual
                  ! convention. We discard the imaginary parts on the spin-box 
                  ! diagonal and one of the off-diagonal spin-box elements:
                  !
                  !        | dD(1,:)                dD(3,:) - i dD(4,:) |
                  !   DM = |                                            |
                  !        | dD(3,:) + i dD(4,:)    dD(2,:)             |
                  !
                  !         |-Im(zDu(1,:)*ph)        -Im(zDu(2,:)*ph) + i Re(zDu(2,:)*ph) |
                  !      =  |                                                             |
                  !         |-Im(zDu(2,:)*ph) - i Re(zDu(2,:)*ph)        -Im(zDu(4,:)*ph) |
                  !
                  ! Same for the EDM.

                  ! UP, UP
                  zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                  dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) - aimag( zDtemp)
                  ! DOWN, DOWN
                  zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                  dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) - aimag( zDtemp)
                  ! UP, DOWN
                  zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                  dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) - aimag( zDtemp)
                  dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - real( zDtemp, dp)

                  ! UP, UP
                  zEtemp = ph(jo)*zEu(1,ind,1:E_dim3)
                  dE(1,lind,1:E_dim3) = dE(1,lind,1:E_dim3) - aimag( zEtemp)
                  ! DOWN, DOWN
                  zEtemp = ph(jo)*zEu(4,ind,1:E_dim3)
                  dE(2,lind,1:E_dim3) = dE(2,lind,1:E_dim3) - aimag( zEtemp)
                  ! UP, DOWN
                  zEtemp = ph(jo)*zEu(2,ind,1:E_dim3)
                  dE(3,lind,1:E_dim3) = dE(3,lind,1:E_dim3) - aimag( zEtemp)
                  dE(4,lind,1:E_dim3) = dE(4,lind,1:E_dim3) - real( zEtemp ,dp)
                end if
#ifdef TS_SOC_DEBUG
        if ((sc_off(1,jo) .eq. 0.) .and. (sc_off(2,jo) .eq. 0.) .and. &
            (sc_off(3,jo) .eq. 0.)) then
            write(1340+Node, "(4I10)") io, UCORB(lio, nr), UCORB(l_col(lind),nr)
            write(1440+Node, "(8(spE26.15))") dD(:,lind,1)
            write(1540+Node, "(8(spE26.15))") dD(:,lind,2)
        end if
#endif

            end do

            end if
            end if

          end do
!$OMP end parallel do

      else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp)
          do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

                jo = UCORB(l_col(lind),nr)

                ind = rind + SFIND(upp_col,jo)
                if ( rind < ind ) then

                  jo = (l_col(lind)-1) / nr

                  ! The integration is
                  ! \rho = i e^{-i.k.R} \int (Gf^R-Gf^A) dE
                  !
                  !
                  ! After adding the phase factor we have to separate the complex and
                  ! real parts of the DM elements and re-order them to the usual
                  ! convention. We discard the imaginary parts on the spin-box 
                  ! diagonal and one of the off-diagonal spin-box elements:
                  !
                  !        | dD(1,:)                dD(3,:) - i dD(4,:) |
                  !   DM = |                                            |
                  !        | dD(3,:) + i dD(4,:)    dD(2,:)             |
                  !
                  !         |-Im(zDu(1,:)*ph)        -Im(zDu(2,:)*ph) + i Re(zDu(2,:)*ph) |
                  !      =  |                                                             |
                  !         |-Im(zDu(2,:)*ph) - i Re(zDu(2,:)*ph)        -Im(zDu(4,:)*ph) |
                  !

                  ! UP, UP
                  zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                  dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) - aimag( zDtemp)
                  ! DOWN, DOWN
                  zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                  dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) - aimag( zDtemp)
                  ! UP, DOWN
                  zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                  dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) - aimag( zDtemp)
                  dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - real( zDtemp, dp)
                end if

            end do

            end if
            end if

          end do
!$OMP end parallel do

        end if

      end if

    end subroutine add_k_DM_3D_nc

    subroutine add_k_DM_3D_so

      if ( non_Eq ) then

        if ( hasEDM ) then

! No data race will occur, sparsity pattern only tranversed once
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp,zEtemp)
        do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

              jo = UCORB(l_col(lind),nr)

              ind = rind + SFIND(upp_col,jo)
              if ( rind < ind ) then

                jo = (l_col(lind)-1) / nr

                ! The integration is:
                ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
                !
                !
                ! After adding the phase factor we have to separate the complex and
                ! real parts of the DM elements and re-order them to the usual
                ! convention:
                !
                !        | dD(1,:) + i dD(5,:)    dD(3,:) - i dD(4,:) |
                !   DM = |                                            |
                !        | dD(7,:) + i dD(8,:)    dD(2,:) + i dD(6,:) |
                !
                !        | zDu(1,:)*ph            zDu(2,:)*ph         |
                !      = |                                            |
                !        | zDu(3,:)*ph            zDu(4,:)*ph         |
                !
                !   Similar for the energy density matrix, but we discard the
                !   imaginary parts on the spin-box diagonal and one of 
                !   the off-diagonal spin-box elements:
                !
                !         | E(1,:)              E(3,:) - iE(4,:) | 
                !   EDM = |                                      |
                !         | E(3,:) + iE(4,:)    E(2,:)           |
                !
                !         | Re(zDu(1,:)*ph)      zDu(2,:)*ph         |
                !      =  |                                          |
                !         | (zDu(2,:)*ph)*       Re(zDu(4,:)*ph)     |

                ! UP, UP
                zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(5,lind,1:D_dim3) = dD(5,lind,1:D_dim3) + aimag( zDtemp)
                ! DOWN, DOWN
                zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(6,lind,1:D_dim3) = dD(6,lind,1:D_dim3) + aimag( zDtemp)
                ! UP, DOWN
                zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - aimag( zDtemp)
                ! DOWN, UP
                zDtemp =  ph(jo)*zDu(3,ind,1:D_dim3)
                dD(7,lind,1:D_dim3) = dD(7,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(8,lind,1:D_dim3) = dD(8,lind,1:D_dim3) + aimag( zDtemp)

                ! UP, UP
                zEtemp = ph(jo)*zEu(1,ind,1:E_dim3)
                dE(1,lind,1:E_dim3) = dE(1,lind,1:E_dim3) + real( zEtemp ,dp)
                ! DOWN, DOWN
                zEtemp = ph(jo)*zEu(4,ind,1:E_dim3)
                dE(2,lind,1:E_dim3) = dE(2,lind,1:E_dim3) + real( zEtemp ,dp)
                ! UP, DOWN
                zEtemp = ph(jo)*zEu(2,ind,1:E_dim3)
                dE(3,lind,1:E_dim3) = dE(3,lind,1:E_dim3) + real( zEtemp ,dp)
                dE(4,lind,1:E_dim3) = dE(4,lind,1:E_dim3) - aimag( zEtemp)
              end if
#ifdef TS_SOC_DEBUG
        if ((sc_off(1,jo) .eq. 0.) .and. (sc_off(2,jo) .eq. 0.) .and. &
            (sc_off(3,jo) .eq. 0.)) then
            write(1340+Node, "(4I10)") io, UCORB(lio, nr), UCORB(l_col(lind),nr)
            write(1440+Node, "(8(spE26.15))") dD(:,lind,1)
            write(1540+Node, "(8(spE26.15))") dD(:,lind,2)
        end if
#endif

            end do

            end if
            end if

        end do
!$OMP end parallel do

      else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp)
        do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

              jo = UCORB(l_col(lind),nr)

              ind = rind + SFIND(upp_col,jo)
              if ( rind < ind ) then

                jo = (l_col(lind)-1) / nr
                !
                !
                ! After adding the phase factor we have to separate the complex and
                ! real parts of the DM elements and re-order them to the usual
                ! convention:
                !
                !        | dD(1,:) + i dD(5,:)    dD(3,:) - i dD(4,:) |
                !   DM = |                                            |
                !        | dD(7,:) + i dD(8,:)    dD(2,:) + i dD(6,:) |
                !
                !        | zDu(1,:)*ph            zDu(2,:)*ph         |
                !      = |                                            |
                !        | zDu(3,:)*ph            zDu(4,:)*ph         |
                !

                ! UP, UP
                zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(5,lind,1:D_dim3) = dD(5,lind,1:D_dim3) + aimag( zDtemp)
                ! DOWN, DOWN
                zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(6,lind,1:D_dim3) = dD(6,lind,1:D_dim3) + aimag( zDtemp)
                ! UP, DOWN
                zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - aimag( zDtemp)
                ! DOWN, UP
                zDtemp =  ph(jo)*zDu(3,ind,1:D_dim3)
                dD(7,lind,1:D_dim3) = dD(7,lind,1:D_dim3) + real( zDtemp ,dp)
                dD(8,lind,1:D_dim3) = dD(8,lind,1:D_dim3) + aimag( zDtemp)
              end if

            end do

            end if
            end if

        end do
!$OMP end parallel do

        end if

      else ! non_eq == .false.

      if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp,zEtemp)
        do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

              jo = UCORB(l_col(lind),nr)

              ind = rind + SFIND(upp_col,jo)
              if ( rind < ind ) then

                jo = (l_col(lind)-1) / nr

                ! The integration is
                ! \rho = i e^{-i.k.R} \int (Gf^R-Gf^A) dE
                !
                !
                ! After adding the phase factor we have to separate the complex and
                ! real parts of the DM elements and re-order them to the usual
                ! convention:
                !
                !        | dD(1,:) + i dD(5,:)    dD(3,:) - i dD(4,:) |
                !   DM = |                                            |
                !        | dD(7,:) + i dD(8,:)    dD(2,:) + i dD(6,:) |
                !
                !        | i zDu(1,:)*ph          i zDu(2,:)*ph       |
                !      = |                                            |
                !        | i zDu(3,:)*ph          i zDu(4,:)*ph       |
                !
                !   Similar for the energy density matrix, but we discard the
                !   imaginary parts on the spin-box diagonal and one of 
                !   the off-diagonal spin-box elements:
                !
                !         | E(1,:)              E(3,:) - iE(4,:) | 
                !   EDM = |                                      |
                !         | E(3,:) + iE(4,:)    E(2,:)           |
                !
                !         | Re(zDu(1,:)*ph)      zDu(2,:)*ph         |
                !      =  |                                          |
                !         | (zDu(2,:)*ph)*       Re(zDu(4,:)*ph)     |
                !

                ! UP, UP
                zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) - aimag( zDtemp)
                dD(5,lind,1:D_dim3) = dD(5,lind,1:D_dim3) + real( zDtemp, dp)
                ! DOWN, DOWN
                zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) - aimag( zDtemp)
                dD(6,lind,1:D_dim3) = dD(6,lind,1:D_dim3) + real( zDtemp, dp)
                ! UP, DOWN
                zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) - aimag( zDtemp)
                dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - real( zDtemp, dp)
                ! DOWN, UP
                zDtemp =  ph(jo)*zDu(3,ind,1:D_dim3)
                dD(7,lind,1:D_dim3) = dD(7,lind,1:D_dim3) - aimag( zDtemp)
                dD(8,lind,1:D_dim3) = dD(8,lind,1:D_dim3) + real( zDtemp, dp)

                ! UP, UP
                zEtemp = ph(jo)*zEu(1,ind,1:E_dim3)
                dE(1,lind,1:E_dim3) = dE(1,lind,1:E_dim3) - aimag( zEtemp)
                ! DOWN, DOWN
                zEtemp = ph(jo)*zEu(4,ind,1:E_dim3)
                dE(2,lind,1:E_dim3) = dE(2,lind,1:E_dim3) - aimag( zEtemp)
                ! UP, DOWN
                zEtemp = ph(jo)*zEu(2,ind,1:E_dim3)
                dE(3,lind,1:E_dim3) = dE(3,lind,1:E_dim3) - aimag( zEtemp)
                dE(4,lind,1:E_dim3) = dE(4,lind,1:E_dim3) - real( zEtemp ,dp)
              end if
#ifdef TS_SOC_DEBUG
        if ((sc_off(1,jo) .eq. 0.) .and. (sc_off(2,jo) .eq. 0.) .and. &
            (sc_off(3,jo) .eq. 0.)) then
            write(1340+Node, "(4I10)") io, UCORB(lio, nr), UCORB(l_col(lind),nr)
            write(1440+Node, "(8(spE26.15))") dD(:,lind,1)
            write(1540+Node, "(8(spE26.15))") dD(:,lind,2)
        end if
#endif

            end do

            end if
            end if

        end do
!$OMP end parallel do

      else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind,zDtemp)
        do lio = 1 , lnr

            if ( l_ncol(lio) /= 0 ) then
            io = index_local_to_global(dit,lio)
            if ( up_ncol(io) /= 0 ) then
            rind = up_ptr(io)
            upp_col => up_col(rind+1:rind+up_ncol(io))

            do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

              jo = UCORB(l_col(lind),nr)

              ind = rind + SFIND(upp_col,jo)
              if ( rind < ind ) then

                jo = (l_col(lind)-1) / nr

                ! The integration is
                ! \rho = i e^{-i.k.R} \int (Gf^R-Gf^A) dE
                !
                !
                ! After adding the phase factor we have to separate the complex and
                ! real parts of the DM elements and re-order them to the usual
                ! convention:
                !
                !        | dD(1,:) + i dD(5,:)    dD(3,:) - i dD(4,:) |
                !   DM = |                                            |
                !        | dD(7,:) + i dD(8,:)    dD(2,:) + i dD(6,:) |
                !
                !        | i zDu(1,:)*ph          i zDu(2,:)*ph       |
                !      = |                                            |
                !        | i zDu(3,:)*ph          i zDu(4,:)*ph       |

                ! UP, UP
                zDtemp =  ph(jo)*zDu(1,ind,1:D_dim3)
                dD(1,lind,1:D_dim3) = dD(1,lind,1:D_dim3) - aimag( zDtemp)
                dD(5,lind,1:D_dim3) = dD(5,lind,1:D_dim3) + real( zDtemp, dp)
                ! DOWN, DOWN
                zDtemp =  ph(jo)*zDu(4,ind,1:D_dim3)
                dD(2,lind,1:D_dim3) = dD(2,lind,1:D_dim3) - aimag( zDtemp)
                dD(6,lind,1:D_dim3) = dD(6,lind,1:D_dim3) + real( zDtemp, dp)
                ! UP, DOWN
                zDtemp =  ph(jo)*zDu(2,ind,1:D_dim3)
                dD(3,lind,1:D_dim3) = dD(3,lind,1:D_dim3) - aimag( zDtemp)
                dD(4,lind,1:D_dim3) = dD(4,lind,1:D_dim3) - real( zDtemp, dp)
                ! DOWN, UP
                zDtemp =  ph(jo)*zDu(3,ind,1:D_dim3)
                dD(7,lind,1:D_dim3) = dD(7,lind,1:D_dim3) - aimag( zDtemp)
                dD(8,lind,1:D_dim3) = dD(8,lind,1:D_dim3) + real( zDtemp, dp)
              end if

            end do

            end if
            end if

        end do
!$OMP end parallel do

      end if

      end if
    end subroutine add_k_DM_3D_so

  end subroutine add_k_DM_3D

  subroutine add_Gamma_DM_2D(spDM,spuDM,D_dim2,spEDM,spuEDM,E_dim2)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D

! *********************
! * INPUT variables   *
! *********************
    ! The local integrated sparsity arrays
    type(dSpData2D), intent(inout) :: spDM
    ! The current Gamma-point global sparsity arrays
    type(dSpData2D), intent(inout) :: spuDM
    integer, intent(in) :: D_dim2
    ! The local integrated sparsity arrays
    type(dSpData2D), intent(inout) :: spEDM
    ! The current Gamma-point global sparsity arrays
    type(dSpData2D), intent(inout) :: spuEDM
    integer, intent(in) :: E_dim2

    ! Arrays needed for looping the sparsity
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:), upp_col(:)
    real(dp), pointer :: dD(:,:), dE(:,:)
    real(dp), pointer :: dDu(:,:), dEu(:,:)
    integer :: lnr, lio, lind, io, ind, nr, jo, rind
    logical :: hasEDM

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spuDM)) ) return

    hasEDM = initialized(spEDM)

    ! get distribution
    dit => dist(spDM)

    l_s => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD => val(spDM)
    if ( hasEDM ) dE  => val(spEDM)

    up_s => spar(spuDM)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    dDu => val(spuDM)
    if ( hasEDM ) dEu => val(spuEDM)

    if ( size(dDu,2) < D_dim2 .or. size(dD,2) < D_dim2 ) then
       call die('add_Gamma_DM: Error in code')
    end if

    if ( hasEDM ) then
       if ( size(dEu,2) < E_dim2 .or. size(dE,2) < E_dim2 ) then
          call die('add_Gamma_DM: Error in code')
       end if
    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected G-DM.')

    if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr
          
          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! we might still have SIESTA-non-Gamma
             jo = ucorb(l_col(lind),nr)

             ! This sparsity pattern is in UC
             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
               dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + dDu(ind,1:D_dim2)
               dE(lind,1:E_dim2) = dE(lind,1:E_dim2) + dEu(ind,1:E_dim2)
             end if
             
          end do

          end if
          end if

       end do
!$OMP end parallel do
       
    else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             jo = ucorb(l_col(lind),nr)

             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
               dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + dDu(ind,1:D_dim2)
             end if
             
          end do

          end if
          end if

       end do
!$OMP end parallel do

    end if

  end subroutine add_Gamma_DM_2D

  subroutine add_Gamma_DM_3D(spDM,spuDM,D_dim2,D_dim3,spEDM,spuEDM,E_dim2,E_dim3)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData3D

! *********************
! * INPUT variables   *
! *********************
    ! The local integrated sparsity arrays
    type(dSpData3D), intent(inout) :: spDM
    ! The current Gamma-point global sparsity arrays
    type(dSpData3D), intent(inout) :: spuDM
    integer, intent(in) :: D_dim2, D_dim3
    ! The local integrated sparsity arrays
    type(dSpData3D), intent(inout) :: spEDM
    ! The current Gamma-point global sparsity arrays
    type(dSpData3D), intent(inout) :: spuEDM
    integer, intent(in) :: E_dim2, E_dim3

    ! Arrays needed for looping the sparsity
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:), upp_col(:)
    real(dp), pointer :: dD(:,:,:), dE(:,:,:)
    real(dp), pointer :: dDu(:,:,:), dEu(:,:,:)
    integer :: lnr, lio, lind, io, ind, nr, jo, rind
    logical :: hasEDM

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spuDM)) ) return

    hasEDM = initialized(spEDM)

    ! get distribution
    dit  => dist(spDM)

    l_s  => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD   => val(spDM)
    if ( hasEDM ) dE  => val(spEDM)

    up_s => spar(spuDM)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    dDu  => val(spuDM)
    if ( hasEDM ) dEu => val(spuEDM)

    if ( size(dDu,1) < D_dim2 .or. size(dD,1) < D_dim2 ) then
       call die('add_Gamma_DM: Error in code')
    end if
    if ( size(dDu,3) < D_dim3 .or. size(dD,3) < D_dim3 ) then
       call die('add_Gamma_DM: Error in code')
    end if

    if ( hasEDM ) then
       if ( size(dEu,1) < E_dim2 .or. size(dE,1) < E_dim2 ) then
          call die('add_Gamma_DM: Error in code')
       end if
       if ( size(dEu,3) < E_dim3 .or. size(dE,3) < E_dim3 ) then
          call die('add_Gamma_DM: Error in code')
       end if
    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.

    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected G-DM.')

    if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! we might still have SIESTA-non-Gamma
             jo = ucorb(l_col(lind),nr)

             ! This sparsity pattern is in UC
             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
               dD(1:D_dim2,lind,1:D_dim3) = dD(1:D_dim2,lind,1:D_dim3) + dDu(1:D_dim2,ind,1:D_dim3)
               dE(1:E_dim2,lind,1:E_dim3) = dE(1:E_dim2,lind,1:E_dim3) + dEu(1:E_dim2,ind,1:E_dim3)
             end if

          end do

          end if
          end if

       end do
!$OMP end parallel do

    else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,upp_col,lind,jo,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          if ( up_ncol(io) /= 0 ) then
          rind = up_ptr(io)
          upp_col => up_col(rind+1:rind+up_ncol(io))

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             jo = ucorb(l_col(lind),nr)

             ind = rind + SFIND(upp_col,jo)
             if ( rind < ind ) then
               dD(1:D_dim2,lind,1:D_dim3) = dD(1:D_dim2,lind,1:D_dim3) + dDu(1:D_dim2,ind,1:D_dim3)
             end if

          end do

          end if
          end if

       end do
!$OMP end parallel do

    end if

  end subroutine add_Gamma_DM_3D

  subroutine update_DM_sp2D(dit,sp,n_nzs,DM, spDM, Ef, &
       EDM, spEDM, ipnt, UpSpGlobal)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_iSpData1D

! *********************
! * INPUT variables   *
! *********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData2D), intent(inout) :: spDM
    ! fermi-level, we shift the energy density matrix back
    real(dp), intent(in) :: Ef
    ! Sparse energy-DM-arrays (local)
    real(dp), intent(inout) :: EDM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData2D), intent(inout) :: spEDM
    ! I.e. a pointer from the local update sparsity to the local sparsity
    type(iSpData1D), intent(in), optional :: ipnt
    ! Whether the update sparsity pattern is a global update sparsity pattern
    logical, intent(in), optional :: UpSpGlobal

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:), lupp_col(:)
    integer, pointer :: pnt(:)
    real(dp), pointer :: dD(:,:), dE(:,:)
    integer :: lnr, nr, uind, lio, io, lind, ind, jo, rind
    logical :: hasipnt, hasEDM, lUpSpGlobal

#ifdef TRANSIESTA_TIMING
    call timer('TS_update_DM', 1)
#endif

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    dD => val(spDM)

    hasEDM = initialized(spEDM) 

    if ( hasEDM ) then

       dE => val(spEDM)

       ! the actual size of the shift
       jo = nnzs(spDM)

       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       call daxpy(jo,Ef,dD(1,1),1,dE(1,1),1)

    end if

    ! We have that the update sparsity pattern is in local
    ! form.
    ! this means that sp == s (besides the non-update objects)
    ! Hence we don't need to utilize index_local_to_global

    
    hasipnt = present(ipnt)
    if ( hasipnt ) hasipnt = initialized(ipnt)

    lUpSpGlobal = .false.
    if ( present(UpSpGlobal) ) lUpSpGlobal = UpSpGlobal
    
    if ( .not. lUpSpGlobal ) then

    if ( lnr /= nrows(s) ) &
         call die('The sparsity format is not as expected u-DM.')

    if ( hasipnt ) then

       ! The pointer
       pnt => val(ipnt)

       ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(io,uind,ind)
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
                
             ind = pnt(uind)
             
             DM(ind) = DM(ind) + dD(uind,1)
             if ( hasEDM ) EDM(ind) = EDM(ind) + dE(uind,1)
             
          end do

          end if
       end do
!$OMP end parallel do
       
    else 

       ! This loop is across the local rows...
! no data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(io,uind,jo,ind)
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
                
             ! We are dealing with a non-UC sparsity pattern
             jo = lup_col(uind)

             ! Now we loop across the local region
             ind = l_ptr(io) + &
                  minloc(abs(l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io))-jo),1)
             if ( l_col(ind) /= jo ) then
                do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                   if ( l_col(ind) == jo ) exit
                end do
                if ( l_col(ind) /= jo ) cycle
             end if
             
             ! We need to add in case of special weighting...
             DM(ind) = DM(ind) + dD(uind,1)
             if ( hasEDM ) EDM(ind) = EDM(ind) + dE(uind,1)
             
          end do
          end if
       end do
!$OMP end parallel do
    end if

    else
       ! We have a global update sparsity pattern

              ! This is the global sparsity pattern
       ! i.e. we require to call index_local_to_global
       ! The global sparsity pattern is not in supercell format

       ! This loop is across the local rows...
! We will never have a data race here (it is on local sparsity pattern)
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,lupp_col,lind,jo,ind)
       do lio = 1 , lnr

          ! obtain the global index of the local orbital.
          io = index_local_to_global(dit,lio)

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Retrieve pointer index
          rind = lup_ptr(io)
          lupp_col => lup_col(rind+1:rind+lup_ncol(io))

          ! Do a loop in the local sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! we need to compare with the global update sparsity
             jo = UCORB(l_col(lind),nr)

             ! Now we loop across the update region
             ! This one must *per definition* have less elements.
             ! Hence, we can exploit this, and find equivalent
             ! super-cell orbitals.
             ! Ok, this is Gamma (but to be consistent)
             ind = rind + SFIND(lupp_col,jo)
             if ( ind > rind ) then
               DM(lind) = DM(lind) + dD(ind,1)
               if ( hasEDM ) EDM(lind) = EDM(lind) + dE(ind,1)
             end if

          end do

          end if
       end do
!$OMP end parallel do

    end if

#ifdef TRANSIESTA_TIMING
    call timer('TS_update_DM', 2)
#endif

  end subroutine update_DM_sp2D

  subroutine update_DM_sp3D(dit,sp,n_nzs, DM, spDM, Ef, &
       EDM, spEDM, ipnt, UpSpGlobal)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData3D
    use class_iSpData1D
    use m_spin, only : spin
#ifdef TS_SOC_DEBUG
    use parallel,only: Node
#endif

! *********************
! * INPUT variables   *
! *********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs,spin%DM)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData3D), intent(inout) :: spDM
    ! fermi-level, we shift the energy density matrix back
    real(dp), intent(in) :: Ef
    ! Sparse energy-DM-arrays (local)
    real(dp), intent(inout) :: EDM(n_nzs,spin%EDM)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData3D), intent(inout) :: spEDM
    ! I.e. a pointer from the local update sparsity to the local sparsity
    type(iSpData1D), intent(in), optional :: ipnt
    ! Whether the update sparsity pattern is a global update sparsity pattern
    logical, intent(in), optional :: UpSpGlobal

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:), lupp_col(:)
    integer, pointer :: pnt(:)
    real(dp), pointer :: dD(:,:,:), dE(:,:,:)
    integer :: lnr, nr, uind, lio, io, lind, ind, jo, rind
    logical :: hasipnt, hasEDM, lUpSpGlobal

#ifdef TS_SOC_DEBUG
    character(len=1024) :: filename
#endif

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s  => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    dD => val(spDM)

    hasEDM = initialized(spEDM)

    if ( hasEDM ) then

       dE => val(spEDM)

       ! the actual size of the shift
       jo = nnzs(spDM)

       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       ! Because the sparsity dimension is in the 2nd place
       ! we need to only do it for the diagonal parts
       ! TODO NW This will be have to be changed for non-collinear case 8 -> spin%DM
       call daxpy(jo,Ef,dD(1,1,1),spin%DM,dE(1,1,1),spin%EDM)
       call daxpy(jo,Ef,dD(2,1,1),spin%DM,dE(2,1,1),spin%EDM)

    end if

#ifdef TS_SOC_DEBUG
      write (filename, "('DM-idx-',I0,'.dat')") Node
      open(unit=340+Node, file=filename)
      write (filename, "('DM-',I0,'.dat')") Node
      open(unit=440+Node, file=filename)
#endif
    ! We have that the update sparsity pattern is in local
    ! form.
    ! this means that sp == s (besides the non-update objects)
    ! Hence we don't need to utilize index_local_to_global


    hasipnt = present(ipnt)
    if ( hasipnt ) hasipnt = initialized(ipnt)

    lUpSpGlobal = .false.
    if ( present(UpSpGlobal) ) lUpSpGlobal = UpSpGlobal

    if ( .not. lUpSpGlobal ) then

    if ( lnr /= nrows(s) ) &
         call die('The sparsity format is not as expected u-DM.')

    if ( hasipnt ) then

       ! The pointer
       pnt => val(ipnt)

       ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(io,uind,ind)
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

             ind = pnt(uind)
             DM(ind,:) = DM(ind,:) + dD(:,uind,1)
             if ( hasEDM ) EDM(ind,:) = EDM(ind,:) + dE(:,uind,1)
#ifdef TS_SOC_DEBUG
             write(340+Node, "(4I10)") io, index_local_to_global(dit,io), l_col(ind), 1
             write(440+Node, "(8(spE26.15))") DM(ind,:)
#endif

          end do

          end if
       end do
!$OMP end parallel do

    else

       ! This loop is across the local rows...
! no data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(io,uind,jo,ind)
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

             ! We are dealing with a non-UC sparsity pattern
             jo = lup_col(uind)

             ! Now we loop across the local region
             ind = l_ptr(io) + &
                  minloc(abs(l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io))-jo),1)
             if ( l_col(ind) /= jo ) then
                do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                   if ( l_col(ind) == jo ) exit
                end do
                if ( l_col(ind) /= jo ) cycle
             end if

             ! We need to add in case of special weighting...
             DM(ind,:) = DM(ind,:) + dD(:,uind,1)
             if ( hasEDM ) EDM(ind,:) = EDM(ind,:) + dE(:,uind,1)
#ifdef TS_SOC_DEBUG
             write(340+Node, "(4I10)") ind, index_local_to_global(dit,io), UCORB(l_col(ind),nr), 2
             write(440+Node, "(8(spE26.15))") DM(ind,:)
#endif

          end do
          end if
       end do
!$OMP end parallel do
    end if

    else
       ! We have a global update sparsity pattern

              ! This is the global sparsity pattern
       ! i.e. we require to call index_local_to_global
       ! The global sparsity pattern is not in supercell format

       ! This loop is across the local rows...
! We will never have a data race here (it is on local sparsity pattern)
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,lupp_col,lind,jo,ind)
       do lio = 1 , lnr

          ! obtain the global index of the local orbital.
          io = index_local_to_global(dit,lio)

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Retrieve pointer index
          rind = lup_ptr(io)
          lupp_col => lup_col(rind+1:rind+lup_ncol(io))

          ! Do a loop in the local sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! we need to compare with the global update sparsity
             jo = UCORB(l_col(lind),nr)

             ! Now we loop across the update region
             ! This one must *per definition* have less elements.
             ! Hence, we can exploit this, and find equivalent
             ! super-cell orbitals.
             ! Ok, this is Gamma (but to be consistent)
             ind = rind + SFIND(lupp_col,jo)
             if ( ind > rind ) then
               DM(lind,:) = DM(lind,:) + dD(:,ind,1)
               if ( hasEDM ) EDM(lind,:) = EDM(lind,:) + dE(:,ind,1)
#ifdef TS_SOC_DEBUG
               write(340+Node, "(4I10)") ind, io, UCORB(l_col(lind),nr), 3
               write(440+Node, "(8(spE26.15))") DM(lind,:)
#endif

             end if
          end do

          end if
       end do
!$OMP end parallel do

    end if
#ifdef TS_SOC_DEBUG
    close (340+Node)
    close (440+Node)
#endif

  end subroutine update_DM_sp3D

  ! This routine will ONLY be called if .not. IsVolt,
  ! Hence we don't have any sparsity patterns with local sparsity patterns
  ! that is dealing with this routine (hence we do need the index_local_to_global)
  subroutine update_zDM_sp2D(dit,sp,n_nzs,DM,spDM, Ef, &
       EDM,spEDM, k, n_s, sc_off)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData2D

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(zSpData2D), intent(inout) :: spDM, spEDM
    ! The fermi level
    real(dp), intent(in) :: Ef
    ! The k-point...
    real(dp), intent(in) :: k(3)
    ! The supercell offset
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:), lupp_col(:)
    complex(dp), pointer :: zD(:,:), zE(:,:)
    complex(dp) :: ph(0:n_s-1)
    integer :: lio, io, jo, ind, nr
    integer :: lnr, lind, rind
    logical :: hasEDM
#ifdef TRANSIESTA_TIMING
    call timer('TS_update_zDM', 1)
#endif

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    zD => val(spDM)

    hasEDM = initialized(spEDM)

    if ( hasEDM ) then

       zE => val(spEDM)

       ! the actual size of the shift
       jo = nnzs(spDM)

       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       ph(0) = cmplx(Ef,0._dp,dp)
       call zaxpy(jo,ph(0),zD(1,1),1,zE(1,1),1)

    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern (but still the global sparsity pattern)

    if ( nr /= nrows(s) ) call die('The sparsity format is not as &
         &expected uz-DM.')

    do jo = 0 , n_s - 1
       ph(jo) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,jo)),dp))
    end do

    ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,lupp_col,lind,jo,ind)
    do lio = 1 , lnr

       ! obtain the global index of the local orbital.
       io = index_local_to_global(dit,lio)

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( lup_ncol(io) /= 0 ) then
       rind = lup_ptr(io)
       lupp_col => lup_col(rind+1:rind+lup_ncol(io))

       ! Do a loop in the local sparsity pattern...
       ! The local sparsity pattern is more "spread", hence
       ! we do fewer operations by having this as an outer loop
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          jo = UCORB(l_col(lind),nr)

          ! Now search the update region
          ! This one must *per definition* have less elements.
          ! Hence, we can exploit this, and find equivalent
          ! super-cell orbitals.
          ind = rind + SFIND(lupp_col,jo)

          if ( rind < ind ) then

            ! \rho = - \Im e^{-i.k.R} [ \int (Gf^R-Gf^A) dE + \int Gf^R\Gamma Gf^A dE ]
            ! Gf^R-Gf^A from outside
            jo = (l_col(lind)-1) / nr

            DM(lind) = DM(lind) - aimag( ph(jo)*zD(ind,1) )
            if ( hasEDM ) &
                EDM(lind) = EDM(lind) - aimag( ph(jo)*zE(ind,1) )
          end if

       end do

       end if
    end do
!$OMP end parallel do

#ifdef TRANSIESTA_TIMING
    call timer('TS_update_zDM', 2)
#endif
  end subroutine update_zDM_sp2D

  ! This routine will ONLY be called if .not. IsVolt,
  ! Hence we don't have any sparsity patterns with local sparsity patterns
  ! that is dealing with this routine (hence we do need the index_local_to_global)
  ! TODO NW remove ikpt argument
  subroutine update_zDM_sp3D(dit,sp,n_nzs,DM,spDM, Ef, &
       EDM,spEDM, k, n_s, sc_off, &
       ikpt)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData3D
    use m_spin, only : spin
#ifdef TS_SOC_DEBUG
    use parallel,only: Node
#endif

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs,spin%DM), EDM(n_nzs,spin%EDM)
    ! Updated sparsity arrays (they contain the current integration)
    type(zSpData3D), intent(inout) :: spDM, spEDM
    ! The fermi level
    real(dp), intent(in) :: Ef
    ! The k-point...
    real(dp), intent(in) :: k(3)
    ! The supercell offset
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)

    integer, intent(in), optional :: ikpt

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:), lupp_col(:)
    complex(dp), pointer :: zD(:,:,:), zE(:,:,:)
    complex(dp) :: ph(0:n_s-1), zDtemp, zEtemp
    integer :: lio, io, jo, ind, nr
    integer :: lnr, lind, rind, ispin
    logical :: hasEDM

#ifdef TS_SOC_DEBUG
    character(len=1024) :: filename
#endif

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    zD => val(spDM)

    hasEDM = initialized(spEDM)

    if ( hasEDM ) then

       zE => val(spEDM)

       ! the actual size of the shift
       jo = nnzs(spDM)

       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       ph(0) = cmplx(Ef,0._dp,dp)
       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       ! Because the sparsity dimension is in the 2nd place
       ! we need to only do it for the diagonal parts
       ! Note the fixed size of 4. In complex we always have a complex matrix
       ! with four entries per orbital no need for spin%DM/spin%EDM
       call zaxpy(jo,ph(0),zD(1,1,1),4,zE(1,1,1),4)
       call zaxpy(jo,ph(0),zD(4,1,1),4,zE(4,1,1),4)

    end if

#ifdef TS_SOC_DEBUG
    write (filename, "('DM-idx-',I0,'-',I0,'.dat')") Node, ikpt
    open(unit=1340+Node, file=filename)
    write (filename, "('DM-',I0,'-',I0,'.dat')") Node, ikpt
    open(unit=1440+Node, file=filename)
#endif
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern (but still the global sparsity pattern)

    if ( nr /= nrows(s) ) call die('The sparsity format is not as &
         &expected uz-DM.')

    do jo = 0 , n_s - 1
       ph(jo) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,jo)),dp))
    end do

    if (spin%NCol) then
      call update_zDM_sp3D_nc
    else if (spin%SO) then
      call update_zDM_sp3D_so
    else
      call die('update_zDM_sp3D: Error in code. This function should only be &
                &called during non-collinear or spin-orbit calculations.')
    end if

#ifdef TS_SOC_DEBUG
    close(1340+Node)
    close(1440+Node)
#endif

    contains

    subroutine update_zDM_sp3D_nc

    ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,lupp_col,lind,jo,ind,zDtemp,zEtemp)
      do lio = 1 , lnr

        ! obtain the global index of the local orbital.
        io = index_local_to_global(dit,lio)

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( lup_ncol(io) /= 0 ) then
        rind = lup_ptr(io)
        lupp_col => lup_col(rind+1:rind+lup_ncol(io))

        ! Do a loop in the local sparsity pattern...
        ! The local sparsity pattern is more "spread", hence
        ! we do fewer operations by having this as an outer loop
        do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

            jo = UCORB(l_col(lind),nr)

            ! Now search the update region
            ! This one must *per definition* have less elements.
            ! Hence, we can exploit this, and find equivalent
            ! super-cell orbitals.
            ind = rind + SFIND(lupp_col,jo)

            if ( rind < ind ) then

              ! \rho = - \Im e^{-i.k.R} [ \int (Gf^R-Gf^A) dE + \int Gf^R\Gamma Gf^A dE ]
              ! Gf^R-Gf^A from outside
              ! We have to separate the complex and
              ! real parts of the DM elements and re-order them to the usual
              ! convention:
              !
              !        | DM(1,:)                DM(3,:) - i DM(4,:) |
              !   DM = |                                            |
              !        | DM(3,:) + i DM(4,:)    DM(2,:)             |
              !
              !         |-Im(zD(1,:)*ph)        -Im(zD(2,:)*ph) + i Re(zD(2,:)*ph) |
              !      =  |                                                          |
              !         |-Im(zD(2,:)*ph) - i Re(zD(2,:)*ph)        -Im(zD(4,:)*ph) |
              !
              !   Similar for the energy density matrix, but we discard the
              !   imaginary parts on the spin-box diagonal and one of 
              !   the off-diagonal spin-box elements:
              !
              !         | EDM(1,:)              EDM(3,:) - iEDM(4,:) | 
              !   EDN = |                                            |
              !         | EDM(3,:) + iEDM(4,:)              EDM(2,:) |
              !
              !         |-Im(zE(1,:)*ph)        -Im(zE(2,:)*ph) + i Re(zE(2,:)*ph) |
              !      =  |                                                          |
              !         |-Im(zE(2,:)*ph) - i Re(zE(2,:)*ph)        -Im(zE(4,:)*ph) |

              jo = (l_col(lind)-1) / nr

              ! Up, Up
              zDtemp =  ph(jo)*zD(1,ind,1)
              DM(lind,1) = DM(lind,1) - aimag( zDtemp )
              ! Down, Down
              zDtemp =  ph(jo)*zD(4,ind,1)
              DM(lind,2) = DM(lind,2) - aimag( zDtemp )
              ! Up, Down
              zDtemp =  ph(jo)*zD(2,ind,1)
              DM(lind,3) = DM(lind,3) - aimag( zDtemp )
              DM(lind,4) = DM(lind,4) - real ( zDtemp, dp)

              if ( hasEDM ) then
                ! Up, Up
                zEtemp = ph(jo)*zE(1,ind,1)
                EDM(lind,1) = EDM(lind,1) - aimag( zEtemp)
                ! Down, Down
                zEtemp = ph(jo)*zE(4,ind,1)
                EDM(lind,2) = EDM(lind,2) - aimag( zEtemp)
                ! Up, Down
                zEtemp = ph(jo)*zE(2,ind,1)
                EDM(lind,3) = EDM(lind,3) - aimag( zEtemp)
                EDM(lind,4) = EDM(lind,4) - real ( zEtemp ,dp)
              end if
#ifdef TS_SOC_DEBUG
        if ((sc_off(1,jo) .eq. 0.) .and. (sc_off(2,jo) .eq. 0.) .and. &
            (sc_off(3,jo) .eq. 0.)) then
            write(1340+Node, "(4I10)") lind, ind, io, UCORB(l_col(lind),nr)
            write(1440+Node, "(8(spE26.15))") DM(lind,:)
        end if
#endif
            end if
        end do

        end if
      end do
!$OMP end parallel do

    end subroutine update_zDM_sp3D_nc

    subroutine update_zDM_sp3D_so

    ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,rind,lupp_col,lind,jo,ind,zDtemp,zEtemp)
      do lio = 1 , lnr

        ! obtain the global index of the local orbital.
        io = index_local_to_global(dit,lio)

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( lup_ncol(io) /= 0 ) then
        rind = lup_ptr(io)
        lupp_col => lup_col(rind+1:rind+lup_ncol(io))

        ! Do a loop in the local sparsity pattern...
        ! The local sparsity pattern is more "spread", hence
        ! we do fewer operations by having this as an outer loop
        do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

            jo = UCORB(l_col(lind),nr)

            ! Now search the update region
            ! This one must *per definition* have less elements.
            ! Hence, we can exploit this, and find equivalent
            ! super-cell orbitals.
            ind = rind + SFIND(lupp_col,jo)

            if ( rind < ind ) then

              ! \rho = - \Im e^{-i.k.R} [ \int (Gf^R-Gf^A) dE + \int Gf^R\Gamma Gf^A dE ]
              ! Gf^R-Gf^A from outside
              jo = (l_col(lind)-1) / nr

              ! We have to separate the complex and
              ! real parts of the DM elements and re-order them to the usual
              ! convention:
              !        | DM(1,:)                DM(3,:) - i DM(4,:) |
              !   DM = |                                            |
              !        | DM(3,:) + i DM(4,:)    DM(2,:)             |
              !
              !         |-Im(zD(1,:)*ph) + i Re(zD(5,:)*ph)  -Im(zD(2,:)*ph) + i Re(zD(2,:)*ph) |
              !      =  |                                                                       |
              !         |-Im(zD(2,:)*ph) - i Re(zD(2,:)*ph)  -Im(zD(4,:)*ph) + i Re(zD(6,:)*ph) |
              !
              !   Similar for the energy density matrix, but we discard the
              !   imaginary parts on the spin-box diagonal and one of 
              !   the off-diagonal spin-box elements:
              !
              !         | EDM(1,:)              EDM(3,:) - iEDM(4,:) | 
              !   EDN = |                                            |
              !         | EDM(3,:) + iEDM(4,:)              EDM(2,:) |
              !
              !         |-Im(zE(1,:)*ph)        -Im(zE(2,:)*ph) + i Re(zE(2,:)*ph) |
              !      =  |                                                          |
              !         |-Im(zE(2,:)*ph) - i Re(zE(2,:)*ph)        -Im(zE(4,:)*ph) |z

              ! Up, Up
              zDtemp =  ph(jo)*zD(1,ind,1)
              DM(lind,1) = DM(lind,1) - aimag( zDtemp )
              DM(lind,5) = DM(lind,5) + real ( zDtemp, dp)
              ! Down, Down
              zDtemp =  ph(jo)*zD(4,ind,1)
              DM(lind,2) = DM(lind,2) - aimag( zDtemp )
              DM(lind,6) = DM(lind,6) + real ( zDtemp, dp)
              ! Up, Down
              zDtemp =  ph(jo)*zD(2,ind,1)
              DM(lind,3) = DM(lind,3) - aimag( zDtemp )
              DM(lind,4) = DM(lind,4) - real ( zDtemp, dp)
              ! Down, Up
              zDtemp =  ph(jo)*zD(3,ind,1)
              DM(lind,7) = DM(lind,7) - aimag( zDtemp )
              DM(lind,8) = DM(lind,8) + real ( zDtemp, dp)

              if ( hasEDM ) then
                ! Up, Up
                zEtemp = ph(jo)*zE(1,ind,1)
                EDM(lind,1) = EDM(lind,1) - aimag( zEtemp)
                ! Down, Down
                zEtemp = ph(jo)*zE(4,ind,1)
                EDM(lind,2) = EDM(lind,2) - aimag( zEtemp)
                ! Up, Down
                zEtemp = ph(jo)*zE(2,ind,1)
                EDM(lind,3) = EDM(lind,3) - aimag( zEtemp)
                EDM(lind,4) = EDM(lind,4) - real ( zEtemp ,dp)
              end if
#ifdef TS_SOC_DEBUG
        if ((sc_off(1,jo) .eq. 0.) .and. (sc_off(2,jo) .eq. 0.) .and. &
            (sc_off(3,jo) .eq. 0.)) then
            write(1340+Node, "(4I10)") lind, ind, io, UCORB(l_col(lind),nr)
            write(1440+Node, "(8(spE26.15))") DM(lind,:)
        end if
#endif
            end if
        end do

        end if
      end do
!$OMP end parallel do

    end subroutine update_zDM_sp3D_so

  end subroutine update_zDM_sp3D

  subroutine init_DM_1D(dit,sp,n_nzs,DM,EDM, up_sp, Calc_Forces)
    ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! The updated sparsity arrays...
    type(Sparsity), intent(inout) :: up_sp
    logical, intent(in) :: Calc_Forces

    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:), lupp_col(:)
    integer :: lnr, lio, lind, io, jo, ind, nr

#ifdef TRANSIESTA_TIMING
    call timer('TS_init_DM', 1)
#endif

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(up_sp, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    if ( nr /= nrows(up_sp) ) call die('The sparsity format is not as &
         &expected i-DM.')

    if ( Calc_Forces ) then
     
       ! This loop is across the local rows...
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lupp_col,lind,jo,ind)
       do lio = 1 , lnr

          ! obtain the global index of the local orbital.
          io = index_local_to_global(dit,lio)

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ind = lup_ptr(io)
          lupp_col => lup_col(ind+1:ind+lup_ncol(io))

          ! Do a loop in the local sparsity pattern...
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             ! Search for unit-cell entry in update sparsity
             ! pattern (UC)
             jo = UCORB(l_col(lind),nr)
             ind = SFIND(lupp_col,jo)
             if ( ind > 0 ) then
                DM(lind) = 0._dp
                EDM(lind) = 0._dp
             end if

          end do
          end if
       end do
!$OMP end parallel do

    else
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lupp_col,lind,jo,ind)
       do lio = 1 , lnr
          io = index_local_to_global(dit,lio)
          if ( lup_ncol(io) /= 0 ) then
          ind = lup_ptr(io)
          lupp_col => lup_col(ind+1:ind+lup_ncol(io))
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             jo = UCORB(l_col(lind),nr)
             ind = SFIND(lupp_col,jo)
             if ( ind > 0 ) DM(lind) = 0._dp
          end do
          end if
       end do
!$OMP end parallel do
    end if

#ifdef TRANSIESTA_TIMING
    call timer('TS_init_DM', 2)
#endif

  end subroutine init_DM_1D

  subroutine init_DM_2D(dit,sp,n_nzs,spin,DM,EDM, up_sp, Calc_Forces)
    ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity
    use t_spin, only : tspin

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays and number of spin componenets
    integer, intent(in) :: n_nzs
    type(tspin), intent(in) :: spin
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs,spin%DM), EDM(n_nzs,spin%EDM)
    ! The updated sparsity arrays...
    type(Sparsity), intent(inout) :: up_sp
    logical, intent(in) :: Calc_Forces

    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:), lupp_col(:)
    integer :: lnr, lio, lind, io, jo, ind, nr

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(up_sp, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    if ( nr /= nrows(up_sp) ) call die('The sparsity format is not as &
         &expected i-DM.')

    if ( Calc_Forces ) then

       ! This loop is across the local rows...
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lupp_col,jo,ind,lind)
          do lio = 1 , lnr

             ! obtain the global index of the local orbital.
             io = index_local_to_global(dit,lio)

             ! Quickly go past the empty regions... (we have nothing to update)
             if ( lup_ncol(io) /= 0 ) then

             ind = lup_ptr(io)
             lupp_col => lup_col(ind+1:ind+lup_ncol(io))

             ! Do a loop in the local sparsity pattern...
             do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

                ! Search for unit-cell entry in update sparsity
                ! pattern (UC)
                jo = UCORB(l_col(lind),nr)
                ind = SFIND(lupp_col,jo)
                if ( ind > 0 ) then
                  DM(lind,:) = 0._dp
                  EDM(lind,:) = 0._dp
                end if

             end do
             end if
          end do
!$OMP end parallel do

    else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lupp_col,jo,ind,lind)
       do lio = 1 , lnr
          io = index_local_to_global(dit,lio)
          if ( lup_ncol(io) /= 0 ) then
          ind = lup_ptr(io)
          lupp_col => lup_col(ind+1:ind+lup_ncol(io))
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
            jo = UCORB(l_col(lind),nr)
            ind = SFIND(lupp_col,jo)
            if ( ind > 0 ) DM(lind,:) = 0._dp
          end do
          end if
       end do
!$OMP end parallel do
    end if

  end subroutine init_DM_2D

end module m_ts_dm_update
