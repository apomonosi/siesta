! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! Module to calculate the total charge.

module dm_charge_m

  implicit none

  private
  public :: dm_charge

  interface dm_charge
     module procedure dm_charge_class
     module procedure dm_charge_direct
     module procedure dm_charge_direct_nspin
  end interface

contains

  subroutine dm_charge_class(spin, DM_2D, S_1D, Q)
    use precision, only: dp

    use t_spin, only: tSpin
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D
    use class_OrbitalDistribution

    ! Containers for data
    type(tSpin), intent(in) :: spin
    type(dSpData2D), intent(inout) :: DM_2D
    type(dSpData1D), intent(inout) :: S_1D
    real(dp), intent(out) :: Q(spin%Grid)

    !  New interface data
    type(Sparsity), pointer :: sp

    integer :: lnr
    integer, pointer :: lptr(:), ncol(:)
    real(dp), pointer :: DM(:,:), S(:)

    ! get the distribution
    sp => spar(DM_2D)
    call attach(sp, n_col=ncol, list_ptr=lptr, nrows=lnr)

    S => val(S_1D)
    DM => val(DM_2D)

    call dm_charge(spin%Grid, lnr, ncol, lptr, S, DM, Q)

  end subroutine dm_charge_class


  subroutine dm_charge_direct(spin, no_l, ncol, list_ptr, S, DM, Q)

    use precision, only: dp
    use t_spin, only: tSpin

    ! Containers for data
    type(tSpin), intent(in) :: spin
    integer, intent(in) :: no_l, ncol(no_l), list_ptr(no_l)
    real(dp), intent(in) :: S(:), DM(:,:)
    real(dp), intent(out) :: Q(spin%Grid)

    call dm_charge(spin%Grid, no_l, ncol, list_ptr, S, DM, Q)

  end subroutine

  subroutine dm_charge_direct_nspin(nspin, no_l, ncol, list_ptr, S, DM, Q)

    use precision, only: dp

#ifdef MPI
    use m_mpi_utils,     only: globalize_sum
#endif

    ! Containers for data
    integer, intent(in) :: nspin
    integer, intent(in) :: no_l, ncol(no_l), list_ptr(no_l)
    real(dp), intent(in) :: S(:), DM(:,:)
    real(dp), intent(out) :: Q(nspin)

    integer :: ir, ind
#ifdef MPI
    real(dp) :: Qr(nspin)
#endif

    ! Initialize
    Q(:) = 0._dp

    if ( size(DM, 2) == 8 ) then

       do ir = 1, no_l
         do ind = list_ptr(ir) + 1, list_ptr(ir) + ncol(ir)
           ! In the SOC case, hermitify Dscf
           Q(1:2) = Q(1:2) + DM(ind,1:2) * S(ind)
           Q(3) = Q(3) + 0.5_dp*(DM(ind,3)+DM(ind,7)) * S(ind)
           Q(4) = Q(4) + 0.5_dp*(DM(ind,4)+DM(ind,8)) * S(ind)
         end do
       end do

    else

       ! NOTE: this also works for non-colinear
       do ir = 1, no_l
         do ind = list_ptr(ir) + 1, list_ptr(ir) + ncol(ir)
           Q(:) = Q(:) + DM(ind,:) * S(ind)
         end do
       end do

    endif

#ifdef MPI
    call globalize_sum(Q, Qr)
    Q = Qr
#endif

  end subroutine

end module dm_charge_m
