
module m_ncps_froyen_ps_t

  implicit none

  private

  integer, parameter  :: dp = selected_real_kind(14)

  public :: froyen_ps_t
  public :: pseudo_init_constant
  public :: pseudo_reset

  type froyen_ps_t
        character(len=2)        :: name
        integer                 :: nr
        integer                 :: nrval
        real(dp)                :: zval
        real(dp)                :: gen_zval  ! Generation valence charge
        character(len=40)       :: gen_config_string ! Packed config
        logical                 :: relativistic
        character(len=10)       :: correlation
        character(len=2)        :: icorr
        integer                 :: libxc_packed_code
        character(len=3)        :: irel
        character(len=4)        :: nicore
        real(dp)                :: a
        real(dp)                :: b
        character(len=10)       :: method(6)
        character(len=70)       :: text
        integer                 :: npotu
        integer                 :: npotd
        real(dp), pointer       :: r(:)        => null()
        real(dp), pointer       :: chcore(:)   => null()
        real(dp), pointer       :: chval(:)    => null()
        real(dp), pointer       :: vdown(:,:)  => null()
        real(dp), pointer       :: vup(:,:)    => null()
        integer, pointer        :: ldown(:)    => null()
        integer, pointer        :: lup(:)      => null()
     end type froyen_ps_t

      CONTAINS

      subroutine pseudo_init_constant(p)
      type(froyen_ps_t), intent(inout) :: p

      p%nr = 0
      p%nrval = 0
      p%zval = 0._dp
      p%gen_zval = 0._dp
      p%gen_config_string = ' '
      p%relativistic = .false.
      p%correlation = ' '
      p%icorr = ' '
      p%libxc_packed_code = 0
      p%irel = ' '
      p%nicore = ' '
      p%a = 0._dp
      p%b = 0._dp
      p%method(:) = ' '
      p%text = ' '
      p%npotu = 0
      p%npotd = 0

      end subroutine pseudo_init_constant

      subroutine pseudo_reset(p)
         implicit none
         type(froyen_ps_t), intent(inout) :: p

         call pseudo_init_constant(p)
         if ( associated(p%r)      ) deallocate(p%r)
         if ( associated(p%chcore) ) deallocate(p%chcore)
         if ( associated(p%chval)  ) deallocate(p%chval)
         if ( associated(p%vdown)  ) deallocate(p%vdown)
         if ( associated(p%vup)    ) deallocate(p%vup)
         if ( associated(p%ldown)  ) deallocate(p%ldown)
         if ( associated(p%lup)    ) deallocate(p%lup)

         nullify( p%r, p%chcore, p%chval, p%vdown, p%vup, p%ldown, p%lup)
      end subroutine pseudo_reset

   end module m_ncps_froyen_ps_t



