module m_ncps_utils

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)
  
  public :: get_n_semicore_shells
  public :: ncps_has_spin_orbit_potentials

CONTAINS

subroutine get_n_semicore_shells(p,nval_gs,nsemic)
  use m_ncps_froyen_ps_t,  only: froyen_ps_t
  
    type(froyen_ps_t), intent(in) :: p
    !> valence configuration, ('n' for each l)
    integer, intent(in)                 :: nval_gs(0:3)
    integer, intent(out)                :: nsemic(0:3)

  ! Returns an array with the number of semicore shells
  ! on channel l (0..3)
  ! It uses the information in Froyen-style pseudopotential files

  ! If no pseudo is available for a given l, a 0 is returned.
  ! If for some reason the generation step did not pseudize
  ! a proper valence state (e.g. Cu 3d), the number of semicore
  ! shells is returned as 0, and the fact is recorded in the output.

  character(len=1)     :: sym(0:4) = (/ 's','p','d','f','g' /)

  integer   :: lmax, inp_lun, l, n, i

  character(len=2), allocatable :: orb_arr(:)
  real(dp), allocatable         :: zdown_arr(:)
  real(dp), allocatable         :: zup_arr(:)
  real(dp), allocatable         :: rc_arr(:)
  integer,  allocatable         :: gen_n(:)

  real(dp) :: chgvps

  lmax = p%npotd-1
  allocate(orb_arr(0:lmax))
  allocate(zdown_arr(0:lmax))
  allocate(zup_arr(0:lmax))
  allocate(rc_arr(0:lmax))
  allocate(gen_n(0:lmax))

  write(6,"(/,a)") "---- Pseudopotential check for " // p%name
  call get_ps_conf(p%irel,lmax,p%text,chgvps, &
                   orb_arr,zdown_arr,zup_arr,rc_arr)
  if (len_trim(p%gen_config_string) /= 0) then
     write(6,fmt="(a,a,a,f10.6)") "Valence configuration for ps generation: ", &
          trim(p%gen_config_string), " Total charge: ", p%gen_zval
  else
     write(6,fmt="(a,a)") "Valence configuration for ps generation: ", &
          "(assumed as above)"
  endif
  
  nsemic(:) = 0
  do l = 0, lmax
     read(orb_arr(l),"(i1)") gen_n(l)
     nsemic(l) =  nval_gs(l) - gen_n(l)
  enddo

  if (any(nsemic > 0))  then
     write(6,fmt="(a,i1,2x,i1,1x,a)",advance="no") "Semicore shell(s):"
     do l = 0, lmax
        if (gen_n(l) < nval_gs(l)) then
           nsemic(l) =  nval_gs(l) - gen_n(l)
           do i = gen_n(l), nval_gs(l)-1
              write(6,fmt="(1x,i1,a1)",advance="no") i, sym(l)
           enddo
        endif
     enddo
     write(6,fmt=*)
  endif

  ! It could be (for example Cu 3d) that a 'valence' state is not pseudized
  ! In this case we reset nsemic to 0.

  if (any(nsemic < 0)) then
     write(6,fmt="(a,i1,1x,a)",advance="no") "Non-pseudized 'valence' shell(s):"
     do l = 0, lmax
        if (gen_n(l) > nval_gs(l)) then
           nsemic(l) = 0
           do i = nval_gs(l), gen_n(l)-1
              write(6,fmt="(1x,i1,a1)",advance="no") i, sym(l)
           enddo
        endif
     enddo
     write(6,fmt=*)
  endif

end subroutine get_n_semicore_shells

subroutine get_ps_conf(irel,lmax,text,chgvps, &
                       orb_arr,zdown_arr,zup_arr,rc_arr)
!
!     Attempt to decode the valence configuration used for
!     the generation of the pseudopotential
!     (At least, the valence charge)

      character(len=3), intent(in)  :: irel
      integer, intent(in)           :: lmax
      character(len=70), intent(in) :: text
      real(dp), intent(out)         :: chgvps
      character(len=2), intent(out) :: orb_arr(0:)
      real(dp), intent(out)         :: zdown_arr(0:)
      real(dp), intent(out)         :: zup_arr(0:)
      real(dp), intent(out)         :: rc_arr(0:)

      integer  :: l, itext
      real(dp) :: ztot, zup, zdown, rc_read
      character(len=2) :: orb

      chgvps=0.0_dp

            if(irel.eq.'isp') then

               write(6,'(/,2a)')  &
               'Pseudopotential generated from a ', &
                      'spin-dft atomic calculation'
               write(6,'(/,a)') 'Pseudized shells:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8080) &
                      orb, zdown, zup, rc_read
 8080             format(a2,f4.2,1x,f4.2,1x,f4.2)

                  orb_arr(l) = orb
                  zdown_arr(l) = zdown
                  zup_arr(l) = zup
                  rc_arr(l) = rc_read

                  chgvps = chgvps + zdown + zup
                  write(6,8085) orb, zdown, zup, rc_read
 8085             format(a2,'(',f4.2,',',f4.2,') rc: ',f4.2)
               enddo

            else
               if(irel.eq.'rel') then
                  write(6,'(/,2a)')  &
               'Pseudopotential generated from a ', &
                      'fully relativistic atomic calculation'
                  write(6,'(2a)')   &
               'There are spin-orbit semi-local pseudopotentials available'
               endif
               write(6,'(/,a)') 'Pseudized shells:'


               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8090) &
                      orb, ztot, rc_read
 8090             format(a2,f5.2,4x,f5.2)

                  orb_arr(l) = orb
                  zdown_arr(l) = ztot
                  zup_arr(l) = 0.0_dp
                  rc_arr(l) = rc_read

                  chgvps = chgvps + ztot
                  write(6,8095) orb, ztot, rc_read
 8095             format(a2,'(',f5.2,') rc: ',f4.2)
               enddo

           endif
           return

 5000    continue       ! Error return: set chgvps to zero
         call die("Error in get_ps_conf")

         end subroutine get_ps_conf

   function ncps_has_spin_orbit_potentials(p) result(valid)
      use m_ncps_froyen_ps_t,  only: froyen_ps_t

      type(froyen_ps_t), intent(in) :: p
      logical                :: valid

      valid = ((p%irel == "rel") .and. (p%npotu > 0))

   end function ncps_has_spin_orbit_potentials

end module m_ncps_utils

