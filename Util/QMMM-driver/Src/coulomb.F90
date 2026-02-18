

module coulomb_m
  !! This module contains data related to Ewald summations.
  use precision, only : dp
  implicit none

  type, public :: ewald_data_t
    real(dp) :: alpha
    real(dp) :: sqr_alpha
    real(dp) :: sqr_alpha_pi
    real(dp) :: kcut
    real(dp) :: kcut_sq
  contains
    procedure :: set_coef => ewald_coeff
  end type ewald_data_t

  type(ewald_data_t), public :: ewald
    !! Ewald summation coefficients.
  character(len=80) , public :: coulombtype
    !! Ewald or cut-off.
contains

  subroutine ewald_coeff( self, rcut_ew )
    !! Computes Ewald's summation coefficients based on the cut-off radius
    !! passed in the fdf at the begining of the run.
    use precision, only : dp
    use units    , only : pi

    implicit none
    class(ewald_data_t), intent(inout) :: self
      !! Datatype containing information for Ewald summations.
    real(dp)           , intent(in)  :: rcut_ew
      !! Cut off radius for Ewald summations.

    real(dp) :: sfactor
    sfactor = 2.6_dp

    self%alpha = ( sfactor / rcut_ew ) **2
    self%kcut  = 2.0_dp * sfactor * sqrt( self%alpha )

    self%sqr_alpha    = sqrt( self%alpha )
    self%sqr_alpha_pi = sqrt( self%alpha / pi )
    self%kcut_sq      = self%kcut * self%kcut
  end subroutine ewald_coeff
end module coulomb_m