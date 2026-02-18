! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module spin_subs_m
  
  use precision, only: dp
  use m_spin

  implicit none
  
  public :: init_spin

  public :: print_spin_options
  public :: init_spiral

  public :: fname_spin

contains
  
  subroutine init_spin(default_nspin)
    
    use sys, only: die
    use fdf, only : fdf_get, leqi, fdf_deprecated
    use alloc, only: re_alloc

    use files, only: slabel
    use parallel, only: IONode

    !< Externally option for a specific spin configuration for the Hamiltonian.
    !! Upon entry this value is used to setup the default settings, but
    !! fdf-flags still has precedence.
    !! Upon exit, it contains the requested spin-number for the Hamiltonian
    integer, intent(inout), optional :: default_nspin

    character(len=32) :: opt, opt_old

    ! Create pointer assignments...
    call int_pointer(spinor_dim, spin%spinor)
    call int_pointer(nspin     , spin%grid)
    call int_pointer(h_spin_dim, spin%H)
    call int_pointer(e_spin_dim, spin%EDM)

    ! Create pointer assignments...
    call log_pointer(NoMagn, spin%none)
    call log_pointer(SPpol , spin%Col)
    call log_pointer(NonCol, spin%NCol)
    call log_pointer(SpOrb , spin%SO)

    ! Time reversal symmetry
    TrSym  = .true.

    ! Initialize all components of the 'spin' variable
    spin%none = .false.
    spin%Col = .false.
    spin%NCol = .false.
    spin%SO = .false.
    spin%SO_offsite = .false.
    spin%SO_onsite  = .false.

    ! Read in old flags (discouraged)
    if ( present(default_nspin) ) then
      select case ( default_nspin )
      case ( 1 )
        ! default to non-polarized (the default)
      case ( 2 )
        spin%Col = .true.
      case ( 4 )
        spin%NCol = .true.
      case ( 8 )
        spin%SO = .true.
      end select
    end if
    spin%Col  = fdf_get('SpinPolarized', spin%Col)
    spin%NCol = fdf_get('NonCollinearSpin', spin%NCol)
    spin%SO   = fdf_get('SpinOrbit', spin%SO)

    ! Announce the deprecated flags (if used)...
    call fdf_deprecated('SpinPolarized','Spin')
    call fdf_deprecated('NonCollinearSpin','Spin')
    call fdf_deprecated('SpinOrbit','Spin')

    ! Set default option from "old" options
    if ( spin%SO ) then
       opt_old = 'spin-orbit'
    else if ( spin%NCol ) then
       opt_old = 'non-collinear'
    else if ( spin%Col ) then
       opt_old = 'polarized'
    else
       opt_old = 'none'
    end if

    ! In order to enable text input (and obsolete the
    ! 4 different options we use a single value now)
    opt = fdf_get('Spin', opt_old)
    
    if ( leqi(opt, 'none') .or. &
         leqi(opt, 'non-polarized') .or. &
         leqi(opt, 'non-polarised') .or. &
         leqi(opt, 'NP') .or. leqi(opt,'N-P') ) then

       spin%none = .true.
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%NCol)
       call warn_and_set_to_false(Spin%SO)
       
    else if ( leqi(opt, 'polarized') .or. &
         leqi(opt, 'collinear') .or. leqi(opt, 'colinear') .or. &
         leqi(opt, 'polarised') .or. leqi(opt, 'P') ) then
       
       spin%Col = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%NCol)
       call warn_and_set_to_false(Spin%SO)
       
    else if ( leqi(opt, 'non-collinear') .or.  leqi(opt, 'non-colinear') .or. &
         leqi(opt, 'NC') .or. leqi(opt, 'N-C') ) then
       
       spin%NCol = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%SO)
       
    else if ( leqi(opt, 'spin-orbit') .or. leqi(opt, 'S-O') .or. &
         leqi(opt, 'SOC') .or. leqi(opt, 'SO') .or. &
         leqi(opt, 'spin-orbit+offsite') .or. leqi(opt, 'S-O+offsite') .or. &
         leqi(opt, 'SOC+offsite') .or. leqi(opt, 'SO+offsite') ) then
       ! Spin-orbit defaults to the off-site (full) implementation
       
       spin%SO = .true.
       spin%SO_offsite = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%NCol)

    else if ( leqi(opt, 'spin-orbit+onsite') .or. leqi(opt, 'S-O+onsite') .or. &
         leqi(opt, 'SOC+onsite') .or. leqi(opt, 'SO+onsite') ) then
       
       spin%SO = .true.
       spin%SO_onsite = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%NCol)

    else
       write(*,*) 'Unknown spin flag: ', trim(opt)
       call die('Spin: unknown flag, please assert the correct input.')
    end if

    ! These are useful for fine control beyond the old "nspin". Some routines expect
    ! an argument 'nspin' which might really mean 'spin%spinor' (like diagon),
    ! 'spin%Grid' (like dhscf), etc.

    if ( spin%SO ) then
       ! Spin-orbit case

       ! Dimensions
       spin%H      = 8
       spin%DM     = 8
       spin%EDM    = 4
       spin%Grid   = 4
       spin%spinor = 2

       ! Flags
       spin%none = .false.
       spin%Col  = .false.
       spin%NCol = .false.
       spin%SO   = .true.
       ! off/on MUST already be set!

       ! should be moved...
       TRSym = .false.    ! Cannot be changed !

    else if ( spin%NCol ) then
       ! Non-collinear case

       ! Dimensions
       spin%H      = 4
       spin%DM     = 4
       spin%EDM    = 4
       spin%Grid   = 4
       spin%spinor = 2

       ! Flags
       spin%none       = .false.
       spin%Col        = .false.
       spin%NCol       = .true.
       spin%SO         = .false.
       spin%SO_onsite  = .false.
       spin%SO_offsite = .false.

       ! should be moved...
       TRSym = .false.    ! Cannot be changed !

    else if ( spin%Col ) then
       ! Collinear case

       ! Dimensions
       spin%H      = 2
       spin%DM     = 2
       spin%EDM    = 2
       spin%Grid   = 2
       spin%spinor = 2

       ! Flags
       spin%none       = .false.
       spin%Col        = .true.
       spin%NCol       = .false.
       spin%SO         = .false.
       spin%SO_onsite  = .false.
       spin%SO_offsite = .false.

       ! should be moved...
       TRSym = .true.
       ! Allow the user to choose not to have TRS
       TRSym  = fdf_get('TimeReversalSymmetry',TrSym)

    else if ( spin%none ) then
       ! No spin configuration...

       ! Dimensions
       spin%H      = 1
       spin%DM     = 1
       spin%EDM    = 1
       spin%Grid   = 1
       spin%spinor = 1

       ! Flags
       spin%none       = .true.
       spin%Col        = .false.
       spin%NCol       = .false.
       spin%SO         = .false.
       spin%SO_onsite  = .false.
       spin%SO_offsite = .false.

       ! should be moved...
       TRSym = .true.
       ! Allow the user to choose not to have TRS
       TRSym  = fdf_get('TimeReversalSymmetry',TrSym)

    end if

    if ( present(default_nspin) ) then
      default_nspin = spin%H
    end if

    if (.not. spin%none .and. IOnode ) then
       write(6,'(/a)') 'NOTE on units for spin output: '
       write(6,'(a)')  '{S} is the spin magnetic moment, in units of the Bohr magneton,'
       write(6,'(a)')  'assuming a g-factor for the electron spin of exactly 2.'
       write(6,'(a)')  'The spin magnetic moment is then numerically equal to the charge inbalance,'
       write(6,'(a/)') 'in units of electrons, between the spin-up and -down channels.'
    endif

  contains

    subroutine int_pointer(from, to)
      integer, pointer :: from
      integer, intent(in), target :: to

      from => to

    end subroutine int_pointer

    subroutine log_pointer(from, to)
      logical, pointer :: from
      logical, intent(in), target :: to

      from => to

    end subroutine log_pointer
    
    subroutine warn_and_set_to_false(var)
      logical, intent(inout) :: var
      if (var) then
         call message("WARNING","Deprecated spin keyword overridden by new-style 'Spin' input")
         call message("WARNING","Option from deprecated keyword: "//trim(opt_old))
         call message("WARNING","Option from 'Spin' keyword input: "//trim(opt))
         var = .false.
      endif
    end subroutine warn_and_set_to_false

end subroutine init_spin


  ! Print out spin-configuration options
  subroutine print_spin_options( )
    use parallel, only: IONode
    use m_cite, only : add_citation

    character(len=32) :: opt

    if ( .not. IONode ) return

    if ( spin%SO ) then
      if ( spin%SO_offsite ) then
        opt = 'spin-orbit+offsite'
        call add_citation("10.1088/0953-8984/24/8/086005")
      else if (spin%SO_onsite) then
        opt = 'spin-orbit+onsite'
        call add_citation("10.1088/0953-8984/19/19/489001")
      else
        call die("Inconsistent SO option")
      end if
    else if ( spin%NCol ) then
       opt = 'non-collinear'
    else if ( spin%Col ) then
       opt = 'collinear'
    else 
       opt = 'none'
    end if

    write(*,'(a,t53,''= '',a)') 'redata: Spin configuration',trim(opt)
    write(*,'(a,t53,''= '',i0)')'redata: Number of spin components',spin%H
    write(*,'(a,t53,''= '',l1)')'redata: Time-Reversal Symmetry',TRSym
    write(*,'(a,t53,''= '',l1)')'redata: Spin spiral', Spiral
    if ( Spiral ) then
      write(*,'(a,t53,''='',3(tr1,e12.6))') &
          'redata: Spin spiral pitch wave vector', qSpiral
      if ( .not. spin%NCol ) then
        write(*,'(a)') 'redata: WARNING: spin spiral requires non-collinear spin'
        call die("Spin spiral requires non-collinear spin")
      end if
    end if

    if ( spin%SO .and. spin%SO_onsite ) then
      write(*,'(a)') repeat('#',60)
      call message("WARNING","Using a local approximation for spin-orbit coupling.")
      write(*,'(a)') repeat('#',60)
    end if

  end subroutine print_spin_options
  

  subroutine init_spiral( ucell )
    use fdf, only: fdf_get
    use fdf, only: block_fdf, parsed_line
    use fdf, only: fdf_block, fdf_bline, fdf_bclose
    use fdf, only: fdf_bvalues
    use m_get_kpoints_scale, only: get_kpoints_scale

    ! Unit cell lattice vectors
    real(dp), intent(in) :: ucell(3,3)

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    ! Reciprocal cell vectors
    real(dp) :: rcell(3,3)
    integer :: ierr

    Spiral = fdf_block('Spin.Spiral', bfdf)

    if ( .not. Spiral ) return
    ! We cannot use TRS for spirals, regardless of user!
    TRSym = .false.

    ! Retrieve scale
    call get_kpoints_scale('Spin.Spiral.Scale', rcell, ierr)

    if ( ierr /= 0 ) then
      call die('init_spiral: ERROR in Spin.Spiral.Scale, could &
          &not find scale for spiral wave vector.')
    end if

    if (.not. fdf_bline(bfdf,pline)) &
        call die('init_spiral: ERROR in Spin.Spiral block, could &
        &not find spiral wave vector.')

    ! Read pitch wave-vector
    qSpiral(1) = fdf_bvalues(pline,1)
    qSpiral(2) = fdf_bvalues(pline,2)
    qSpiral(3) = fdf_bvalues(pline,3)
    qSpiral(:) = matmul(rcell, qSpiral)

    call fdf_bclose(bfdf)

  end subroutine init_spiral

  function fname_spin(nspin,ispin,slabel,suffix,basename) result(fname)
    integer, intent(in) :: nspin, ispin
    character(len=*), intent(in), optional :: slabel, suffix, basename
    character(len=200) :: fname
    
    if ( present(basename) ) then
       if ( nspin == 1 ) then
          fname = trim(basename)
       else
          if ( ispin == 1 ) fname = trim(basename)//'_UP'
          if ( ispin == 2 ) fname = trim(basename)//'_DN'
       end if
    else
       if ( .not. &
            ( present(slabel) .and. present(suffix) ) ) &
            call die('Error in filename input')
       if ( nspin == 1 ) then
          fname = trim(slabel)//'.'//trim(suffix)
       else
          if ( ispin == 1 ) fname = trim(slabel)//'.'//trim(suffix)//'_UP'
          if ( ispin == 2 ) fname = trim(slabel)//'.'//trim(suffix)//'_DN'
       end if
    end if
    
  end function fname_spin


end module spin_subs_m
