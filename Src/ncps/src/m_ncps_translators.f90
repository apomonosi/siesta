!> An example of the use of the PSML processing library
!!
!! Siesta uses internally the old "Froyen" ps structure. This module provides
!! the functionality to fill that structure using calls to the PSML 
!! processing library.
!!
!> @author Alberto Garcia
!
module m_ncps_translators

  public :: ncps_psml2froyen

CONTAINS

  subroutine ncps_psml2froyen( ps, p, new_grid, a, b, rmax )

! Translate the more complete xml-adapted data structure 
! into the 'Froyen' ps type used in Atom and Siesta.
! Use accessors, and assume nothing about the grid in 
! the XML file. For this,
!
! p%nr, p%a, and p%b should be set on entry if we want
! a new grid. By default, a logarithmic grid with the
! vanilla ATOM parameters is used.
!

        use m_psml
        use m_ncps_psml_plugins
        use m_ncps_froyen_ps_t,  only: froyen_ps_t
        use m_libxc_compat, only: xc_id_t, get_xc_id_from_libxc
        use m_libxc_compat, only: xc_id_to_string

        implicit none 

        integer, parameter  :: dp = selected_real_kind(14)


        type(ps_t), intent(in)               :: ps
        type(froyen_ps_t), intent(inout)     :: p
        logical, intent(in), optional        :: new_grid 
        real(dp), intent(in), optional       :: a
        real(dp), intent(in), optional       :: b
        real(dp), intent(in), optional       :: rmax


        logical          :: want_new_grid
        integer          :: position, i, il, ir, l, n, nval_shells, ii
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld, zelu, occupation, jval
        real(dp)         :: r2, rc, rmax_grid, znuc
        character(len=64):: xc_string
        character(len=40):: method_string
        character(len=1), dimension(0:4) :: &
                         sym = (/ "s", "p", "d", "f", "g" /)
        integer          :: libxc_ids(2), n_xcfuncs
        integer          :: status

        ! These are the "current" ATOM parameters
        ! There is another set turned on by the UCB_COMPAT flag in ATOM
        real(dp), parameter :: aa_def = 6.0_dp
        real(dp), parameter :: bb_def = 80.0_dp    ! UCB_COMPAT: 40.0
        real(dp), parameter :: rmax_def = 120.0_dp ! UCB_COMPAT: 80.0

        type(xc_id_t)                        :: xc_id
        type(ps_annotation_t)                :: grid_annotation, annot
        character(len=40)                    :: strvalue
        logical                              :: log_grid_in_file

        integer, allocatable, dimension(:)   :: idxd, idxu, idxlj
        integer, allocatable, dimension(:)   :: nn, ll
        real(dp), allocatable, dimension(:)  :: rrc

        character(len=10) :: psml_version
        character(len=10) :: relativity
        logical           :: spin_polarized, core_corrections
        integer           :: l_shell, n_shell
        logical           :: shell_found
        real(dp)          :: zup_shell, zdown_shell, occ_shell
        
        integer :: iu, id, li, npotd, npotu, nscalar, nlj, lmax
        logical :: has_lj, has_sr, has_sr_so, has_nonrel
        logical :: has_up_down, has_spin_ave
        real(dp) :: v

        call ps_RootAttributes_Get(ps,version=psml_version)
        write(6,"(a)") "PSML file version: " // &
                        trim(psml_version)
        call ps_PseudoAtomSpec_Get(ps,atomic_symbol=p%name,z_pseudo=p%zval,&
             atomic_number=znuc, pseudo_flavor=method_string,&
             relativity=relativity,spin_dft=spin_polarized,&
             core_corrections=core_corrections)

        call setup_gen_valence_data(ps, nval_shells, p%gen_zval, &
                                        p%gen_config_string )

        !
        call ps_ExchangeCorrelation_Get(ps,annotation=annot,&
                                           n_libxc_functionals=n_xcfuncs)
        
!       Partial support for libxc functionals
!       (no 'cocktails')
!
        do i = 1, n_xcfuncs
           call ps_LibxcFunctional_Get(ps,i,code=libxc_ids(i))
        enddo
        ! pack these unconditionally
        call xcid_pack(n_xcfuncs,libxc_ids,p%libxc_packed_code)

        ! Try to convert to the legacy two-char code used by ATOM
        if (n_xcfuncs == 2) then
           call get_xc_id_from_libxc(libxc_ids,xc_id,status)
        else
           ! treat this below with the 'xc' code
           status = -1
        endif

        if (status == 0) then
              write(6,"(a,2i4)") "Using libxc ids: ", libxc_ids(:)
              write(6,"(a)") trim(xc_id_to_string(xc_id))
              p%icorr = xc_id%atom_id
        else
           ! Fall back to querying a possible XC annotation
           call get_annotation_value(annot,  &
                                      "atom-xc-code",p%icorr,status)
           if (status == 0) then
              write(6,"(a)") "Atom-xc-code from annotation: " // p%icorr
           else
              ! Fall back to using the "xc" pseudo-code
              ! The packed libxc ids have been computed above
              p%icorr = "xc"
           endif
        endif
!
!       Note that most (all?) formats include only "scalar-relativistic" plus
!       maybe spin-orbit components, but never "polarized" pseudos.
!       This is a shortcoming of the "Froyen" format: there is no support
!       for "scalar_relativistic" calculations in the "irel" label.

        if (trim(relativity) == "dirac" ) then
           p%irel    = 'rel'
           ispp      = 'r'
           polarized = .false.
        else
           if (spin_polarized) then
              p%irel    = 'isp'
              ispp      = 's'
              polarized = .true.
           else
              p%irel    = 'nrl'
              ispp      = ' '
              polarized = .false.
           end if
        endif

        if (core_corrections) then
            p%nicore = 'pcec'
         else
            p%nicore = 'nc'
         endif
!
!       Grid handling. We do not assume that the file
!       uses a logarithmic grid.

        want_new_grid = .false.
        if (present(new_grid)) then
           want_new_grid = new_grid
        endif

        ! We want to check for grid annotations, in case
        ! the grid is already of the "atom" type
        ! This is an example of a "PSML plugin"
        call check_atom_grid(ps,log_grid_in_file,p%nrval,p%a,p%b)

        if (log_grid_in_file) then
           
           if (p%nrval > 0) then  ! Normal ATOM grid
              print *, "Using ATOM log grid in PSML file ..."
              p%nr = p%nrval - 1  ! For backwards compatibility
           else    ! Sampled ATOM grid; nrval not available               
              print *, "Using ATOM log grid {a,b} parameters in PSML file ..."
              rmax_grid = rmax_def
              p%nr = nint(log(rmax_grid/(p%b)+1.0d0)/(p%a))
              p%nrval = p%nr + 1  ! Count also r=0
           endif

        else if (want_new_grid) then
           if (.not. present(a)) call die("new grid: a not present")
           if (.not. present(b)) call die("new grid: b not present")
           if (.not. present(rmax)) call die("new grid: rmax not present")
           rmax_grid = rmax
           if (rmax_grid < 1.0_dp) then
              ! If not specified, set it to the default
              rmax_grid = rmax_def
           endif
           p%a = a
           p%b = b
           p%nr = nint(log(rmax_grid/b+1.0d0)/a)
           p%nrval = p%nr + 1  ! Count also r=0

        else 
              print *, "Using ATOM defaults for log grid..."
              ! use the ATOM defaults 
              ! Note that a and b are interchanged in Siesta!
              p%a = 1.0_dp / bb_def
              p%b = exp(-aa_def)/znuc
              rmax_grid = rmax_def
              p%nr = nint(log(rmax_grid/(p%b)+1.0d0)/(p%a))
              p%nrval = p%nr + 1  ! Count also r=0
        endif

        ! Assume first action is ps generation...
        call ps_Provenance_Get(ps,level=1,creator=p%method(1),&
                               date=p%method(2))
        read(method_string,'(4a10)') (p%method(i),i=3,6) 

        has_nonrel = .false.
        has_sr = .false.
        has_sr_so = .false.
        has_up_down = .false.
        has_spin_ave = .false.
        has_lj = .false.

        select case (trim(relativity))
        case ("dirac")

          call ps_Potential_Filter(ps,set=SET_SREL,number=nscalar,indexes=idxd)

           if (nscalar == 0) then

              ! Will get the scalar-relativistic SL potentials
              ! from the lj set

              call ps_Potential_Filter(ps,set=SET_LJ,number=nlj,indexes=idxlj)
              if (nlj == 0) call die("Cannot find srel SL potentials for dirac case")
              has_lj = .true.
              npotd = 0
              npotu = 0
              do i = 1, nlj
                 call ps_Potential_Get(ps,idxlj(i),l=l,j=jval)
                 if ( (l==0) .or. (jval>l)) then
                    npotd = npotd + 1
                 else
                    npotu = npotu + 1
                 endif
              enddo

           else  ! nscalar > 0

              ! We have a scalar-relativistic set
              npotd = nscalar
              ! Get the information on the companion spin-orbit set
              call ps_Potential_Filter(ps,set=SET_SO,number=npotu,indexes=idxu)
              has_sr_so = .true.

           endif

        case ("scalar")

           call ps_Potential_Filter(ps,set=SET_SREL,number=npotd,indexes=idxd)
           if (npotd == 0) call die("Cannot find srel SL potentials for srel case")
           npotu = 0
           has_sr = .true.

           ! We assume that srel calculations are not polarized...
        case ("no")

           if (spin_polarized) then
              call ps_Potential_Filter(ps,set=SET_UP,number=npotu,indexes=idxu)
              call ps_Potential_Filter(ps,set=SET_DOWN,number=npotd,indexes=idxd)
              if (  (npotu > 0) .and. (npotd > 0) ) then

                 ! We have spin_up and spin_down potentials
                 ! Will get the average later

                 has_up_down = .true.

              else
                 ! We must have (at least) the spin_average
                 call ps_Potential_Filter(ps,set=SET_SPINAVE,number=npotd,indexes=idxd)
                 call ps_Potential_Filter(ps,set=SET_SPINDIFF,number=npotu,indexes=idxu)
                 if (npotd == 0) call die("Cannot get spin-averaged SL potentials")
                 has_spin_ave = .true.
              endif

           else ! not polarized

              call ps_Potential_Filter(ps,set=SET_NONREL,number=npotd,indexes=idxd)
              if (npotd == 0) call die("Cannot get non-relativistic SL potentials")
              npotu = 0
              has_nonrel = .true.

           endif

        case default

           call die("Wrong relativity scheme in PSML file")

        end select

! Allocate the radial variables and semilocal potentials

        p%npotd = npotd
        p%npotu = npotu

        allocate(p%r(1:p%nrval))
        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%nrval,1:p%npotd))
           allocate(p%ldown(1:p%npotd))
           allocate(nn(1:npotd),ll(1:npotd),rrc(1:npotd))
        endif

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%nrval,1:p%npotu))
           allocate(p%lup(1:p%npotu))
        endif
! ---

! Calculate the points of the logarithmic radial grid 
        ! The first point here (NOT in the classic files
        ! produced by ATOM (vps,psf) is r=0
        do ir = 1, p%nrval
           p%r(ir) = p%b * (exp(p%a*(ir-1))-1)
        enddo
! ---

! Translate the valence charge density and the pseudo-core charge density,
! and define the value at the first point of the logarithmic grid
        if (core_corrections) then
           do ir = 2, p%nrval
              p%chcore(ir) = ps_CoreCharge_Value(ps,p%r(ir))
              p%chcore(ir) = p%chcore(ir) * (p%r(ir))**2
           enddo
        else
           p%chcore(2:p%nrval) = 0.0_dp
        endif

        do ir = 2, p%nrval
           p%chval(ir) = ps_ValenceCharge_Value(ps,p%r(ir))
           p%chval(ir) = p%chval(ir) * (p%r(ir))**2
        enddo
        !
        ! Note this re-scaling to comply with the psf (Froyen)
        ! convention of writing a neutral-atom total valence charge
        !
        if (abs(p%zval-p%gen_zval) > 1.0e-3_dp) then
           print "(a)", "Rescaling valence charge in psml file to neutral-atom"
           print "(a,2f10.4)", "Zval, GenerationZval:", p%zval, p%gen_zval
           do ir = 2, p%nrval
              p%chval(ir) = (p%zval/p%gen_zval) * p%chval(ir)
           enddo
        endif
!
        ! This is no longer necessary
        ! We can directly evaluate the PSML radfuncs at r=0
        ! without an explicit extrapolation

        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---

        if ( (has_sr_so) .or. (has_spin_ave) .or. (has_nonrel) .or. &
             (has_sr)   ) then
         ! No need for any extra computations
         do il = 1, p%npotd
           call ps_Potential_Get(ps,idxd(il),&
                     l=p%ldown(il),n=nn(il),rc=rrc(il))
           ll(il) = p%ldown(il)
                     
           do ir = 2, p%nrval
             p%vdown(ir,il) = p%r(ir) * &
                           ps_Potential_Value(ps,idxd(il),p%r(ir))
             p%vdown(ir,il) = p%vdown(ir,il) * 2.0_dp   ! rydberg

           enddo
           p%vdown(1,il) = p%vdown(2,il) - r2*(p%vdown(3,il)-p%vdown(2,il))
         enddo

         do il = 1, p%npotu
           call ps_Potential_Get(ps,idxu(il),l=p%lup(il))
           do ir = 2, p%nrval
              p%vup(ir,il) = p%r(ir) * &
                           ps_Potential_Value(ps,idxu(il),p%r(ir))
              p%vup(ir,il) = p%vup(ir,il) * 2.0_dp   ! rydberg
           enddo
           p%vup(1,il) = p%vup(2,il) - r2*(p%vup(3,il)-p%vup(2,il))
         enddo

        ! if not, we need to get the right averages
        else if (has_lj) then

           id = 0
           iu = 0
           do i = 1, nlj
              call ps_Potential_Get(ps,idxlj(i),l=l,j=jval)
              if ( (l==0) .or. (jval>l)) then
                 id = id + 1
                 ! If the lj slpots are not ordered by l in the psml
                 ! file this array will not be monotonic
                 p%ldown(id) = l
                 ! get some extra info needed later
                 ll(id) = l
                 call ps_Potential_Get(ps,idxlj(i),n=nn(id),rc=rrc(id))
              else
                 iu = iu + 1   
                 p%lup(iu) = l
              endif
           enddo
           call assert((id == npotd),"Wrong check on number of down pots")
           call assert((iu == npotu),"Wrong check on number of up pots")

           lmax = maxval(p%ldown(1:npotd))

           p%vdown(:,:) = 0.0_dp
           p%vup(:,:) = 0.0_dp

           do l = 0, lmax   ! Loop on l

              ! determine corresponding indexes in the arrays

              id = 0
              do ii = 1, npotd
                 if (p%ldown(ii) == l) then
                    id = ii
                    exit
                 endif
              enddo
              call assert((id>0),"l mismatch in down array")

              iu = 0   ! This will be zero for l=0, and skipped below
              do ii = 1, npotu
                 if (p%lup(ii) == l) then
                    iu = ii
                    exit
                 endif
              enddo

              do i = 1, nlj
                 call ps_Potential_Get(ps,idxlj(i),l=li,j=jval)
                 if (li /= l) cycle

                 ! Process the two (except for l=0) j channels for this l

                 if ((l==0) .or. jval > l) then  ! j=l+1/2 or l=0,j=0
                    !print *, "l,j+, i, id, iu ", l, jval, i, id, iu
                    do ir = 2, p%nrval
                       v = ps_Potential_Value(ps,idxlj(i),p%r(ir))
                       p%vdown(ir,id) =  p%vdown(ir,id) + (l+1)*v / dble(2*l+1)
                       if (iu>0) p%vup(ir,iu) = p%vup(ir,iu) + 2*v / dble(2*l+1)
                    enddo
                 else   ! j=l-1/2
                    !print *, "l,j=-, i, id, iu ", l, jval, i, id, iu
                    do ir = 2, p%nrval
                       v = ps_Potential_Value(ps,idxlj(i),p%r(ir))
                       p%vdown(ir,id) =  p%vdown(ir,id) + l*v/dble(2*l+1)
                       if (iu>0) p%vup(ir,iu) = p%vup(ir,iu) - 2*v/dble(2*l+1)
                    enddo
                 endif  ! j+ or j-
              enddo   ! over lj set

              ! rydberg and rV
              do ir = 2, p%nrval
                 p%vdown(ir,id) = 2.0_dp * p%r(ir)*p%vdown(ir,id) 
                 if (iu>0) p%vup(ir,iu) = 2.0_dp * p%r(ir)* p%vup(ir,iu)
              enddo
              ! extrapolate to r=0
              p%vdown(1,id) = p%vdown(2,id) - r2*(p%vdown(3,id)-p%vdown(2,id))
              if (iu>0) p%vup(1,iu) = p%vup(2,iu) - r2*(p%vup(3,iu)-p%vup(2,iu))

           enddo        ! over l values

           ! Do we need to re-sort the arrays if l is not monotonic?

        else if (has_up_down) then
           ! to be implemented
           call die("up/down case not implemented yet")
           !
        endif

        ! Encode generation configuration and cutoffs
        ! Since the Froyen record cannot hold multiple shells per l,
        ! the information here is limited:
        !
        ! n values are taken from the 'slps' elements in the PSML
        !     files, and the matching occupations are filled in.
        ! But note that the total valence charge for pseudopotential
        ! generation is not correctly represented here when semicore
        ! states are present. It is (correctly) stored in the
        ! p%gen_zval field, which should be taken into account (and
        ! output to .psf files).
        
        p%text = ' '
        position = 1
        ! We deal with "down" potentials only, as they are enough
        ! for this particular bookeeping
        do il = 1, p%npotd
           n = nn(il)  
           l = ll(il)  
           rc = rrc(il)

           if ( .not. polarized) then
              occupation = 0.0_dp
              shell_found = .false.
              do i = 1, nval_shells
                 call ps_ValenceShell_Get(ps,i,l=l_shell,n=n_shell, &
                                          occupation=occ_shell)
                 if ((l_shell == l) .and. (n_shell == n)) then
                    occupation = occ_shell
                    shell_found = .true.
                    exit
                 endif
              enddo
              if (.not. shell_found) then
                 ! Note that the psml file might not list empty
                 ! valence shells, so we do not raise an error
                 !! call die("Cannot find valence shell for legacy conf info")
              endif
              write(p%text(position:),9070) n, sym(l), occupation, ispp, rc
 9070         format(i1,a1,f5.2,a1,' r=',f5.2,'/')
              position = position + 17
           else
              zeld = 0.0_dp
              zelu = 0.0_dp
              shell_found = .false.
              do i = 1, nval_shells
                 call ps_ValenceShell_Get(ps,i,l=l_shell,n=n_shell, &
                               occ_up=zup_shell,occ_down=zdown_shell)
                 if ((l_shell == l) .and. (n_shell == n)) then
                    zeld = zdown_shell
                    zelu = zup_shell
                    shell_found = .true.
                    exit
                 endif
              enddo
              if (.not. shell_found) then
                 ! Note that the psml file might not list empty
                 ! valence shells, so we do not raise an error
                 !! call die("Cannot find valence shell for legacy conf info")
              endif
              write(p%text(position:),9090)    &
                   n, sym(l), zeld,zelu,ispp, rc
 9090         format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
              position = position + 17
           end if
        enddo

      end subroutine ncps_psml2froyen

      subroutine assert(cond,message)
        logical, intent(in) :: cond
        character(len=*) message

        if (.not. cond) call die(message)
      end subroutine assert

!> Packs libxc codes into a single integer
!> of the form XXXXCCCC      
subroutine xcid_pack(nfuncs,id,code)
  integer, intent(in)  :: nfuncs
  integer, intent(in)  :: id(nfuncs)
  integer, intent(out) :: code

  if (nfuncs == 2) then
     ! Separate X and C functionals
     code = 10000 * id(1) + id(2)
  else if (nfuncs == 1) then
     ! Just a single XC functional
     code = id(1)
  else
     code = 0
  endif

end subroutine xcid_pack

!> For PSML files, gets information about the valence charge configuration
!  used at the time of pseudopotential generation.

subroutine setup_gen_valence_data(ps,gen_nshells,gen_zval,gen_config_string)
  use m_psml, only: ps_t, ps_ValenceConfiguration_get
  use m_psml, only: ps_ValenceShell_get

  integer, parameter :: dp = selected_real_kind(10,100)
  character(len=1), dimension(0:4) :: &
                         sym = (/ "s", "p", "d", "f", "g" /)
  
  type(ps_t), intent(in) :: ps

  !> Total number of valence shells
  integer, intent(out)   :: gen_nshells

  !> The total valence charge density (this is important to avoid
  !  assuming that an ionic configuration was used when semicore
  !  states are present
  real(dp), intent(out)  :: gen_zval

  !> A string holding the configuration in compact form, for use
  !  only in reporting. No occupations are recorded.
  !  Note that complementary information is available about the
  !  ("minimum n" pseudized shells in the p%text record.
  character(len=*), intent(out) :: gen_config_string

  integer  :: i, l_shell, n_shell, str_pos
  real(dp) :: occupation  ! not used for now
    
  call ps_ValenceConfiguration_Get(ps,nshells=gen_nshells, &
                                   charge=gen_zval)
  gen_config_string = ""
  str_pos = 0
  do i = 1, gen_nshells
     call ps_ValenceShell_Get(ps,i,l=l_shell,n=n_shell, &
          occupation=occupation)
     write(gen_config_string(str_pos+1:str_pos+3),fmt="(i1,a1,a1)") &
          n_shell, sym(l_shell) , ":"
     str_pos = str_pos + 3
  enddo

end subroutine setup_gen_valence_data
  
  
end module m_ncps_translators
