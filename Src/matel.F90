!
! Copyright (C) 1996-2025 The SIESTA group
! This file is distributed under the terms of the
! GNU General Public License: see COPYING in the top directory
! or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
!> Re-designed module to allow a pre-computation of the needed
!> matrix elements in parallel, followed by a globalization of the
!> data among all the MPI processes. Once the interpolation
!> tables are setup, further calls to the matrix-element evaluator
!> are cheap. This has a dramatic effect in some routines (such as
!> nlefsm) that had to perform the table-building operations under
!> conditions that did not scale in parallel.

!> Concept: Rogeli Grima (BSC) and Alberto Garcia (ICMAB)
!> Initial implementation: Rogeli Grima (BSC)
!> Further development: Alberto Garcia (ICMAB)

module matel_m

  use matel_ylm_m, only: spher_harm_t
  use matel_table_m, only: matel_t
  use matel_table_m, only: MODE_S, MODE_T, MODE_XYZ
  
  use precision, only : dp
  use alloc, only : alloc_default, allocDefaults

  use m_radfft, only : reset_radfft
  use m_matel_registry, only : EVALUATE, &
                              EVALUATE_X, EVALUATE_Y, EVALUATE_Z

  private

  !> Initialize main tables: S, T, TA
  public :: init_matel_main_tables

  !> Initialize tables X, Y, Z, for use in "molecule" optical calculations and
  !> Berry-phase-based polarization calculations
  public :: init_matel_orb_XYZ_orb   

  !> Initialize tables S_opt, X_opt, Y_opt, Z_opt for use in solid-state optical calculations
  public :: init_matel_optical_P

  !> Initialize table S_wann for overlaps of Wannier90 trial functions
  public :: init_matel_wannier       

  ! Main routine that evaluates matrix elements, keeping the traditional interface
  public :: matel

  ! Specific routine for wannier data, as the structures are more dynamic
  ! and might need to be reset at specific points of the program
  public :: reset_wannier_tables

  ! Global matel cleanup
  public :: cleanup_matel_tables
  
  ! These variables allow to call the initialization routines repeatedly 
  logical :: main_tables_initialized = .false.
  logical :: optical_P_initialized = .false.
  logical :: orb_XYZ_orb_initialized = .false.
  logical :: wannier_initialized = .false.
  
  !> These relate to the functions recorded in the 'matel registry';
  !> the number of each kind serve as markers for the different
  !> sections needed in the tables. 'init_matel_main_tables' fills
  !> the essential first four values, so that routine must be
  !> called first. This is enforced by 'check_main_tables'

  !> A prerequisite for the correct setup is that the
  !> atomic-information tables (in the 'species' data structure) has
  !> been correctly filled up. In future, this module might be more
  !> integrated with that data structure.
  
  integer :: num_orb  ! Number of different orbitals
  integer :: num_kb   ! Number of different KB projectors
  integer :: num_dftu ! Number of different DFT+U functions
  integer :: num_va   ! Number of different Vna functions
  integer :: total_num_wannier_projs = 0  ! Number of different Wannier trial orbitals

  !> k-space spherical-harmonic decomposition of registered
  !> functions and composite functions (such as x*Orb, y*KB, etc)
  
  type(spher_harm_t), target :: ylmk_val_all   ! orbs, kb, dftuprojs, vna
  type(spher_harm_t), target :: ylmk_x_orbs
  type(spher_harm_t), target :: ylmk_y_orbs
  type(spher_harm_t), target :: ylmk_z_orbs
  type(spher_harm_t), target :: ylmk_x_kbs
  type(spher_harm_t), target :: ylmk_y_kbs
  type(spher_harm_t), target :: ylmk_z_kbs
  type(spher_harm_t), target :: ylmk_val_wannier_projs

  !> There are different tables, each appropriate to a given operation and kind of function.
  !> The indexing is a bit cumbersome due to the one-dimensional nature of the matel_registry,
  !> in which all functions are stored in the same section: first orbitals, then KB projectors,
  !> then DFT+U projectors, and finally Vna.
  !> In addition, when using the Wannier interface, "trial orbitals" (called also "projectors"
  !> in the code (numproj of them) need to be dealt with.

  !> The tables are created once and for all for the duration of the program. There is currently
  !> no support for "coexisting qualities" or "on-the-fly" changes (such as could conceivably be
  !> done one day via the Lua interface)
  
  !> Overlaps of PAOs, KBs, or DFTUprojs with PAOs
  type(matel_t) :: tab_S           ! Unity (overlap). 

  !> Laplacian among PAOs
  type(matel_t) :: tab_T           ! -Laplacian.      

  !> Laplacian among Vna functions
  type(matel_t) :: tab_TA          ! -Laplacian.
  
  !> X, Y, Z among orbitals (for dielectric polarization calculations
  !> and "r" optical calculations
  type(matel_t) :: tab_X           ! X projection.   
  type(matel_t) :: tab_Y           ! Y projection.   
  type(matel_t) :: tab_Z           ! Z projection.   

  !> Overlaps of PAOs with KBs for "momentum" optical calculations (different from table "S" because
  !> the KBs are "function 2")
  type(matel_t) :: tab_S_opt       ! Unity (overlap)

  !> X,Y,Z PAOs with KBs for "momentum" optical calculations
  type(matel_t) :: tab_X_opt       ! X projection
  type(matel_t) :: tab_Y_opt       ! Y projection
  type(matel_t) :: tab_Z_opt       ! Z projection

  !> Overlaps of PAOs with wannier projectors
  type(matel_t) :: tab_S_wann      ! Unity (overlap).

  ! convenience variable
  type(allocDefaults) :: OLDEFS

contains

  !> Initialize the main MATEL tables S, T, TA
  subroutine init_matel_main_tables()
    use atm_types, only : nspecies, species

    implicit none
    ! Local Variables
    integer :: is, top, top_vna

    if (main_tables_initialized) return

    ! Set allocation defaults
    call alloc_default(old=oldefs, copy=.true., shrink=.false.)

    call reset_radfft()

    ! Get the number of orbitals, KB projs, DFTU projs and Neutral-atom potentials
    NUM_ORB = 0
    NUM_KB = 0
    NUM_DFTU = 0
    NUM_VA = 0
    do is = 1, nspecies
      NUM_ORB = NUM_ORB + species(is)%norbs
      NUM_KB = NUM_KB + species(is)%nprojs
      NUM_DFTU = NUM_DFTU + species(is)%nprojsdftu
      if (species(is)%z > 0) then ! Not a floating orbital
        NUM_VA = NUM_VA + 1
      endif
    enddo

    ! Find the spherical harmonic decomposition and
    ! Fourier transform of all the registry objects
    ! except the wannier trial functions
    top = num_orb + num_kb + num_dftu 
    top_vna = top + num_va
    call ylmk_val_all%compute_spha(1, top_vna, EVALUATE, 0)

    ! Overlap table from PAOs, KBs, and DFTU projs to PAOs.
    call tab_S%init(MODE_S, ylmk_val_all, 1, top, &
                   ylmk_val_all, 1, num_orb)

    ! Laplacian table from PAOs to PAOs.
    call tab_T%init(MODE_T, ylmk_val_all, 1, num_orb, &
                   ylmk_val_all, 1, num_orb)
    
    ! Laplacian table from Vnas to Vnas
    call tab_TA%init(MODE_T, ylmk_val_all, top+1, top_vna, &
                    ylmk_val_all, top+1, top_vna)

    main_tables_initialized = .true.
    
    ! Restore allocation defaults
    call alloc_default(restore=oldefs)
    
  end subroutine init_matel_main_tables

  subroutine init_matel_orb_XYZ_orb()
    if (orb_XYZ_orb_initialized) return

    call check_for_main_tables()
    
    ! Set allocation defaults
    call alloc_default(old=oldefs, copy=.true., shrink=.false.)

    call reset_radfft()

    ! Find the spherical harmonic decomposition and
    ! Fourier transform of {x,y,z}*\phi, where \phi is a PAO

    ! The X, Y, and Z tables are then overlaps of PAOs and xPAOs
    call ylmk_x_orbs%compute_spha(1, num_orb, EVALUATE_X, 1)
    call tab_X%init(MODE_XYZ, ylmk_val_all, 1, num_orb, &
                   ylmk_x_orbs, 1, num_orb)

    call ylmk_y_orbs%compute_spha(1, num_orb, EVALUATE_Y, 1)
    call tab_Y%init(MODE_XYZ, ylmk_val_all, 1, num_orb, &
                   ylmk_y_orbs, 1, num_orb)

    call ylmk_z_orbs%compute_spha(1, num_orb, EVALUATE_Z, 1)
    call tab_Z%init(MODE_XYZ, ylmk_val_all, 1, num_orb, &
                   ylmk_z_orbs, 1, num_orb)
    
    orb_XYZ_orb_initialized = .true.

    ! Restore allocation defaults
    call alloc_default(restore=oldefs)
  end subroutine init_matel_orb_XYZ_orb

  subroutine init_matel_optical_P()
    if (optical_P_initialized) return

    call check_for_main_tables()

    ! Set allocation defaults
    call alloc_default(old=oldefs, copy=.true., shrink=.false.)

    call reset_radfft()

    ! Overlap table from PAOs to KBs.
    call tab_S_opt%init(MODE_S, ylmk_val_all, 1, num_orb, &
                      ylmk_val_all, num_orb+1, num_orb + num_kb)

    ! Find the spherical harmonic decomposition and
    ! Fourier transform of {x,y,z}*\chi, where \chi is a KB proj.

    ! The X_opt, Y_opt, and Z_opt tables are then overlaps of PAOs and xKBs
    call ylmk_x_kbs%compute_spha(num_orb+1, num_orb+num_kb, &
                                EVALUATE_X, 1)

    call tab_X_opt%init(MODE_XYZ, ylmk_val_all, 1, num_orb, &
                       ylmk_x_kbs, 1, num_kb)

    call ylmk_y_kbs%compute_spha(num_orb+1, num_orb+num_kb, &
                                EVALUATE_Y, 1)

    call tab_Y_opt%init(MODE_XYZ, ylmk_val_all, 1, num_orb, &
                       ylmk_y_kbs, 1, num_kb)

    call ylmk_z_kbs%compute_spha(num_orb+1, num_orb+num_kb, &
                                EVALUATE_Z, 1)

    call tab_Z_opt%init(MODE_XYZ, ylmk_val_all, 1, num_orb, &
                       ylmk_z_kbs, 1, num_kb)
       
    optical_P_initialized = .true.

    ! Restore allocation defaults
    call alloc_default(restore=oldefs)
  end subroutine init_matel_optical_P
      
  subroutine init_matel_wannier(numproj)
    ! Initialize the table of overlaps between orbitals and
    ! Wannier trial functions, once the number of functions is known

    implicit none
    integer, intent(in) :: numproj

    ! Local Variables
    integer :: top_vna

    ! Since there can be several manifolds with their
    ! own wannier projectors, we re-create the table
    ! each time, enlarging it

    ! First clean up previous table if it exists
    if (wannier_initialized) then
      call tab_S_wann%delete()
 
      ! The ylmk_val_wannier_projs structure might be referenced by
      ! the table we just deleted, but we also need to clean it up
      ! directly since we're about to recompute it
      if (associated(ylmk_val_wannier_projs%F)) then
        call ylmk_val_wannier_projs%delete()
      endif
    endif

    ! This variable is an accumulator
    ! It is set to zero at the start of the program, at
    ! declaration time, and can be reset using the
    ! 'reset_wannier_tables' routine
    total_num_wannier_projs = total_num_wannier_projs + numproj
    
    call check_for_main_tables()

    ! Set allocation defaults
    call alloc_default(old=oldefs, copy=.true., shrink=.false.)
    call reset_radfft()

    ! This table refers to overlaps between orbitals and
    ! Wannier trial functions, which are registered in the
    ! pool after every other function
    top_vna = num_orb + num_kb + num_dftu + num_va

    call ylmk_val_wannier_projs%compute_spha(top_vna+1, &
                           top_vna+total_num_wannier_projs, &
                           EVALUATE, 0)

    call tab_S_wann%init(MODE_S, ylmk_val_all, 1, num_orb, &
                     ylmk_val_wannier_projs, &
                     1, total_num_wannier_projs)

    ! This is not really used to avoid re-computation, since
    ! we *do* recreate the table each time (see above)
    wannier_initialized = .true.

    ! Restore allocation defaults
    call alloc_default(restore=oldefs)
  end subroutine init_matel_wannier

  ! Convenience subroutine to reset wannier
  ! structures when required
  subroutine reset_wannier_tables()
    if (wannier_initialized) then
      call tab_S_wann%delete()
      if (associated(ylmk_val_wannier_projs%F)) then
        call ylmk_val_wannier_projs%delete()
      endif
      total_num_wannier_projs = 0
      wannier_initialized = .false.
    endif
  end subroutine reset_wannier_tables
      
  subroutine check_for_main_tables()
    if (.not. main_tables_initialized) then
      call init_matel_main_tables()
    endif
  end subroutine check_for_main_tables
      
  subroutine matel(OPERAT, IG1, IG2, R12, S12, DSDR)
    ! Finds two-center matrix elements between 'atomic orbitals' 
    ! with finite radial and angular momentum cutoffs.

    ! The routine dispatches execution
    ! to the new table objects (Alberto Garcia, July 2019, after Rogeli Grima)      

    ! The tables must be initialized previously by calls to the init_matel_XXXX
    ! routines. 
          
    ! INPUT 
    ! CHARACTER OPERAT : Operator to be used. The valid options are:
    !   'S' => Unity (overlap). Uppercase required for all values.
    !   'T' => -Laplacian
    !   'U' => 1/|r'-r| (with evaluate returning charge distributions) *NOT IMPLEMENTED
    !   'X' => x, returning <phi1(r-R12)|x|phi2(r)> (**origin on second atom**)
    !   'Y' => y, returning <phi1(r-R12)|y|phi2(r)>
    !   'Z' => z, returning <phi1(r-R12)|z|phi2(r)>
    ! INTEGER IG1    : Global index of 1st function (must be positive)
    ! INTEGER IG2    : Global index of 2nd function
    !                    Indexes IG1, IG2 are used only to call 
    !                    routines LCUT, RCUT and EVALUATE (in matel_registry), and 
    !                    may have other meanings within those routines
    ! REAL*8  R12(3) : Vector from first to second atom
    ! ************************* OUTPUT **********************************
    ! REAL*8 S12      : Matrix element between orbitals.
    ! REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.

    ! Length units are arbitrary, but must be consistent in MATEL, RCUT
    ! and EVALUATE. The laplacian unit is (length unit)**(-2).
    ! ************************* BEHAVIOUR *******************************
    ! If |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2), returns exactly zero.

    !> Operation code
    character(len=1), intent(in) :: OPERAT
    !> Global index of 1st function
    integer, intent(in)      :: IG1
    !> Global index of 2nd function
    integer, intent(in)      :: IG2

    !> Vector from first to second atom
    real(dp), intent(in)    :: R12(3)

    !> Matrix element
    real(dp), intent(out)    :: S12
    !> Derivative (gradient) of S12 with respect to R12.
    real(dp), intent(out)    :: DSDR(3)

    integer :: top
    
    top = num_orb + num_kb + num_dftu 

    select case (operat)
    case ('S')
      if ((ig1 <= top) .and. (is_orb(ig2))) then
        call tab_S%get_matel(IG1, IG2, R12, S12, DSDR)
      else if (is_orb(ig1) .and. (is_kb(ig2))) then
        call tab_S_opt%get_matel(IG1, IG2, R12, S12, DSDR)
      else if (is_orb(ig1) .and. (is_wannier_proj(ig2))) then
        call tab_S_wann%get_matel(IG1, IG2, R12, S12, DSDR)
      else
        call die("Cannot process 'S' in matel")
      endif

    case ('T')
      if (is_orb(ig1) .and. (is_orb(ig2))) then
        call tab_T%get_matel(IG1, IG2, R12, S12, DSDR)
      else if (is_vna(ig1) .and. (is_vna(ig2))) then
        call tab_TA%get_matel(IG1, IG2, R12, S12, DSDR)
      else
        call die("Cannot process 'T' in matel")
      endif
       
    case ('X')
      if (is_orb(ig1) .and. (is_orb(ig2))) then
        call tab_X%get_matel(IG1, IG2, R12, S12, DSDR)
      else if (is_orb(ig1) .and. (is_kb(ig2))) then
        call tab_X_opt%get_matel(IG1, IG2, R12, S12, DSDR)
      else
        call die("Cannot process 'X' in matel")
      endif
       
    case ('Y')
      if (is_orb(ig1) .and. (is_orb(ig2))) then
        call tab_Y%get_matel(IG1, IG2, R12, S12, DSDR)
      else if (is_orb(ig1) .and. (is_kb(ig2))) then
        call tab_Y_opt%get_matel(IG1, IG2, R12, S12, DSDR)
      else
        call die("Cannot process 'Y' in matel")
      endif
       
    case ('Z')
      if (is_orb(ig1) .and. (is_orb(ig2))) then
        call tab_Z%get_matel(IG1, IG2, R12, S12, DSDR)
      else if (is_orb(ig1) .and. (is_kb(ig2))) then
        call tab_Z_opt%get_matel(IG1, IG2, R12, S12, DSDR)
      else
        call die("Cannot process 'Z' in matel")
      endif
       
    case ('U')
      call die("Cannot process 'U' in matel")

    case default
      call die("Unrecognized 'operat' in matel")
    end select
  end subroutine matel

  ! Convenience functions
  function is_orb(ig) result(res)
    integer, intent(in) :: ig
    logical :: res

    res = (ig > 0 .and. ig <= num_orb)
  end function is_orb

  function is_kb(ig) result(res)
    integer, intent(in) :: ig
    logical :: res

    res = (ig > num_orb .and. ig <= (num_orb+num_kb))
  end function is_kb

  function is_vna(ig) result(res)
    integer, intent(in) :: ig
    logical :: res

    integer :: top

    top = num_orb + num_kb + num_dftu 
    res = (ig > top .and. ig <= (top + num_va))
  end function is_vna

  function is_wannier_proj(ig) result(res)
    integer, intent(in) :: ig
    logical :: res

    integer :: top_vna

    top_vna = num_orb + num_kb + num_dftu + num_va
    res = (ig > top_vna .and. &
           ig <= (top_vna + total_num_wannier_projs))
  end function is_wannier_proj

  subroutine cleanup_matel_tables()
    ! This routine cleans up all matel tables and associated data structures
    ! Should be called during program termination in siesta_end
    
    if (orb_XYZ_orb_initialized) then
      ! Clean up orbital position operator tables
      call tab_X%delete()
      call tab_Y%delete()
      call tab_Z%delete()
    
      ! Clean up associated spherical harmonic expansions
      if (associated(ylmk_x_orbs%F)) call ylmk_x_orbs%delete()
      if (associated(ylmk_y_orbs%F)) call ylmk_y_orbs%delete()
      if (associated(ylmk_z_orbs%F)) call ylmk_z_orbs%delete()
      orb_XYZ_orb_initialized = .false.
    endif

    if (optical_P_initialized) then
      ! Clean up optical tables
      call tab_S_opt%delete()
      call tab_X_opt%delete()
      call tab_Y_opt%delete()
      call tab_Z_opt%delete()
    
      ! Clean up associated spherical harmonic expansions
      if (associated(ylmk_x_kbs%F)) call ylmk_x_kbs%delete()
      if (associated(ylmk_y_kbs%F)) call ylmk_y_kbs%delete()
      if (associated(ylmk_z_kbs%F)) call ylmk_z_kbs%delete()
      optical_P_initialized = .false.
    endif

    call reset_wannier_tables()

    if (main_tables_initialized) then
      ! Clean up main overlap and kinetic energy tables
      call tab_S%delete()
      call tab_T%delete()
      call tab_TA%delete()
    
      ! Clean up the associated spherical harmonic expansion
      if (associated(ylmk_val_all%F)) call ylmk_val_all%delete()
      main_tables_initialized = .false.
    endif
  end subroutine cleanup_matel_tables
      
end module matel_m
