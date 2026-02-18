module restr_subs
  !! This module contains r0utines that deal with harmonic restraints.
  !! There are eight types of restraints available:
  !!   1) Distances
  !!   2) Angles
  !!   3) Dihedrals
  !!   4) Average distance from atom 1 to other 4 atoms (essentially a distance
  !!      to plane).
  !!   5) Weighted combination of distances: c1*r1 + c2*r2 + c3*r3 + ...
  !!
  !! The way this works, it involves a double loop in the main driver.
  !! We first select the amount of "steps" that we need for the restraints
  !! optimization, i.e. the different values of the restraints that we will be
  !! exploring. For example, take a distance restraint that goes from 2.0 Ang to
  !! 4.0 Ang. If the number of steps we choose is 10, then we will go from 2.0
  !! Ang to 4.0 Ang every 0.2 Ang.
  !!
  !! Then within each "restraint step", we run the regular MD. So if we set
  !! MD.FinalTimeStep 300, then, for the example above, we will be running a
  !! total of 300 * 10 = 3000 MD steps.
  !!
  !! Input in fdf file:
  !! The block Harmonic Restraints must first start by indicating the total
  !! amount of restraints and the number of restraint "steps". Then, for each
  !! restraint:
  !!  * First line indicates the type of restraint, the force constant (in
  !!    kcal/mol*Ang, the initial restraint equilibrium position and the final
  !!    restraint equilibrium position (in Ang for distance, degrees for angle
  !!    and dihedral).
  !!  * Second line indicates the GLOBAL atom indexes involved in each
  !!    restraint (two for distances, three for angles, etc.)
  !!
  !! An exception is made for the "combination of distances":
  !!  * First line: same as before.
  !!  * Second line: Number of distances involved in the combination.
  !!  * Third line: Weight for each of the distances (one per distance).
  !!  * Fourth line: Atom indexes for each of the distances (two per distance).
  !!
  !! See the example below:
  !!
  !! %block HarmonicRestraints
  !!   5 20
  !!
  !!   5 300.0 2.0 4.0
  !!   2
  !!   -2.0 3.0
  !!   474 475 543 545
  !!
  !!   1 400.0 1.5 1.5
  !!   232 233
  !!
  !!   2 1.0 105.0 105.0
  !!   239 238 240
  !!
  !!   3 0.1 90.0 90.0
  !!   405 403 448 450
  !!
  !!   4 250.0 3.0 6.0
  !!   350 332 333 334 335
  !! %endblock HarmonicRestraints
  !!
  !! Blank lines are optional, they are added to improve readability. You can
  !! have as many restraints as you want of any type.
  !!
  !! In this case, the first line indicates that we have 5 different restraints
  !! which will explore 20 different values ("steps"). Only the first and last
  !! restraints actually change value, the others are constant.
  !!
  !! The first restraint is a distance combination restraint, with a force
  !! constant of 300 kcal/mol*Ang, changing from 2.0 Ang to 4.0 Ang. It involves
  !! TWO distances, with weights -2.0 and 3.0, using the four atoms written
  !! just below. This essentially means:
  !!        r = -2.0 * (r474 - r475) + 3.0 (r543 - r545)
  !!
  !! The second is a simple distance restraint, with a force constant of 400
  !! kcal/mol*Ang, and an equilibrium distance of 1.5 Ang which is constant
  !! throughout the simulation.
  !!
  !! The third is an angle restraint, with a force constant of 1 kcal/mol*rad,
  !! with a constant equilibrium position of 105°. Note that this indicates the
  !! angle from atom 239 to atom 240 with 238 in the middle.
  !!
  !! The fourth is an angle restraint, with a force constant of 0.1kcal/mol*rad,
  !! with a constant equilibrium position of 90°. Note that the atoms at both
  !! ends of the dihedral are 405 and 450.
  !!
  !! The last one is the average distance from atom 350 to atoms 332, 333, 334
  !! and 335, wit a force constant of 250 kcal/mol*Ang. In this case, the
  !! restraint equilibrium position will change from 3.0 to 6.0 across the 20
  !! restraint steps.
  !!
  !!
  !! Data on restraints is printed to the systemlabel.ctr.ene file. A restart
  !! for restrained simulations is available via the systemlabel.wrt file.
  use precision, only : dp

  implicit none
  public :: restr_read
  public :: restr_calc
  public :: restr_update
  public :: restr_deinit

  private
  logical, public :: hasRestraints = .false.
    !! Whether we have harmonic restraints.
  integer, public :: restr_steps   = 0
    !! Amount of restraint steps.
  integer, public :: nrestr        = 0
    !! Number of restraints.

  type :: restraint_t
    integer  :: type = 1
      !! Type of restraint (distance, angle, etc).
    real(dp) :: r0   = 0.0_dp
      !! Initial position for the restraint. It later stores the current
      !! equilibrium distance for the restraint.
    real(dp) :: rfin = 0.0_dp
      !! Final value for the restraint.
    real(dp) :: dr   = 0.0_dp
      !! R increment for each restrained step.

    real(dp) :: kf   = 0.0_dp
      !! Force constant for the restraint in kcal/mol.
    integer  :: natoms
      !! Number of atoms for restraint.
    integer, allocatable :: atmidx(:)
      !! Indexes of the atoms involved in the restraint.

    integer  :: ndists
      !! Number of distances involved for restraint type 8.
    real(dp), allocatable :: weights(:)
      !! Weights for each of the components of the restraint. Used only in
      !! coupled distances.
  end type restraint_t
  type(restraint_t), allocatable :: harmrestr(:)
  real(dp)         , allocatable, public :: restr_r(:)

contains

  subroutine restr_read( )
    !! This subroutine reads the system harmonic restraints.
    use fdf         , only : block_fdf, parsed_line
    use fdf         , only : fdf_block, fdf_bline, fdf_bintegers,&
                             fdf_bvalues, fdf_bclose, fdf_get, fdf_defined
    use precision   , only : dp
    use qmmm_files_m, only : qmmm_files
    use sys         , only : die

    implicit none
    integer            :: irs, junit, irestr, ios, iline
    logical            :: found
    real(dp)           :: dummyr

    ! Externals
    external           :: io_assign, io_close

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    if ( fdf_block('ConstrainedOpt',bfdf) ) &
      call die('restr_read: ConstrainedOpt is no longer compatible.'//&
               'use the %block HarmonicRestraints instead.')

    hasRestraints = fdf_defined( 'HarmonicRestraints' )

    if ( .not. hasRestraints ) return

    ! Read input variables
    if ( .not. fdf_block('HarmonicRestraints',bfdf) ) &
      call die('restr_read: You must specify the HarmonicRestraints block')

    if ( .not. fdf_bline(bfdf,pline) ) then
      write(6, *) "WARNING: HarmonicRestraints block found, but no restraints."
      hasRestraints = .false.
      restr_steps   = 0
      call fdf_bclose( bfdf )
      return
    endif

    nrestr      = fdf_bintegers(pline, 1)
    restr_steps = fdf_bintegers(pline, 2)

    if ( restr_steps == 0 ) then
      call die('restr_read: restr_steps must be larger than 0')
    elseif( restr_steps > 100 ) then
      call die('restr_read: restr_steps must be lower than 100')
    endif

    allocate( harmrestr(nrestr), restr_r(nrestr) )

    do irestr = 1, nrestr
      if ( .not. fdf_bline(bfdf,pline) ) &
        write( 6, * ) "restr_read: Problem reading HarmonicRestraints block."
      harmrestr(irestr)%type = fdf_bintegers(pline, 1)
      harmrestr(irestr)%kf   = fdf_bvalues(pline, 2)
      harmrestr(irestr)%r0   = fdf_bvalues(pline, 3)
      harmrestr(irestr)%rfin = fdf_bvalues(pline, 4)

      irs = 0
      if ( .not. fdf_bline(bfdf,pline) ) &
        write( 6, * ) "restr_read: Problem reading HarmonicRestraints block."
      select case( harmrestr(irestr)%type )
      case (1)
        harmrestr(irestr)%natoms = 2
        allocate( harmrestr(irestr)%atmidx(2) )
        harmrestr(irestr)%atmidx(1) = fdf_bintegers(pline, 1)
        harmrestr(irestr)%atmidx(2) = fdf_bintegers(pline, 2)

      case (2)
        harmrestr(irestr)%natoms = 3
        allocate( harmrestr(irestr)%atmidx(3) )
        harmrestr(irestr)%atmidx(1) = fdf_bintegers(pline, 1)
        harmrestr(irestr)%atmidx(2) = fdf_bintegers(pline, 2)
        harmrestr(irestr)%atmidx(3) = fdf_bintegers(pline, 3)

      case (3)
        harmrestr(irestr)%natoms = 4
        allocate( harmrestr(irestr)%atmidx(4) )
        harmrestr(irestr)%atmidx(1) = fdf_bintegers(pline, 1)
        harmrestr(irestr)%atmidx(2) = fdf_bintegers(pline, 2)
        harmrestr(irestr)%atmidx(3) = fdf_bintegers(pline, 3)
        harmrestr(irestr)%atmidx(4) = fdf_bintegers(pline, 4)

      case (4)
        harmrestr(irestr)%natoms = 5
        allocate( harmrestr(irestr)%atmidx(5) )
        harmrestr(irestr)%atmidx(1) = fdf_bintegers(pline, 1)
        harmrestr(irestr)%atmidx(2) = fdf_bintegers(pline, 2)
        harmrestr(irestr)%atmidx(3) = fdf_bintegers(pline, 3)
        harmrestr(irestr)%atmidx(4) = fdf_bintegers(pline, 4)
        harmrestr(irestr)%atmidx(5) = fdf_bintegers(pline, 5)

      case (5)
        harmrestr(irestr)%ndists = fdf_bintegers(pline, 1)
        allocate( harmrestr(irestr)%weights( harmrestr(irestr)%ndists ) )

        if ( .not. fdf_bline(bfdf,pline) ) &
          write( 6, * ) "restr_read: Problem while reading Contrained weights."
        do irs = 1, harmrestr(irestr)%ndists
          harmrestr(irestr)%weights(irs) = fdf_bvalues(pline, irs)
        enddo

        if ( .not. fdf_bline(bfdf,pline) ) &
          write( 6, * ) "restr_read: Problem while reading Contrained atoms."

        harmrestr(irestr)%natoms = 2*harmrestr(irestr)%ndists
        allocate( harmrestr(irestr)%atmidx( harmrestr(irestr)%natoms ) )
        do irs = 1, harmrestr(irestr)%natoms
          harmrestr(irestr)%atmidx(irs) = fdf_bintegers(pline, irs)
        enddo

      case default
        call die( 'restr_read: Wrong restraint type.' )
      end select
    enddo
    call fdf_bclose( bfdf )

    write(6,'(/,a)') 'restr_read: Doing a restrained optimization run.'

    ! Calculates the initial equilibrium r0 for all types of restraints.
    harmrestr(:)%dr = ( harmrestr(:)%rfin - harmrestr(:)%r0 ) / restr_steps

    if ( restr_steps == 1 ) restr_steps = 0
    if ( all( abs(harmrestr(:)%dr) < 1.0e-14_dp ) ) restr_steps = 0


    ! Reads from .rce of a former run.
    found = .false.
    inquire( file = qmmm_files%wrt, exist = found )
    if ( found ) then
      call io_assign( junit )
      open( junit, file = qmmm_files%wrt )

      iline = 0
      ios   = 0
      do while ( ios == 0 )
        iline = iline +1
        read( junit, *, iostat = ios ) dummyr
        if ( ios < 0 ) exit
        if ( ios > 0 ) &
          call die('restr_read: Problem while reading form rce file.')
      enddo
      call io_close( junit )

      iline = iline -1

      if ( restr_steps > iline ) then
        restr_steps = restr_steps - iline
        harmrestr(:)%r0 = harmrestr(:)%r0 + harmrestr(:)%dr * iline

        write(6,'(/,a)') &
          'restr_read: Re-starting a restrained optimization run.'
      endif
    endif !found
  end subroutine restr_read

  subroutine restr_calc( natot, r_inp, fdummy, istep, &
                         istep_restr, ucell, lattice_type )
    !! Calculates restrained energies and forces contributions.
    use alloc    , only : re_alloc, de_alloc
    use functions, only : norm_v2, scalar_v2, angle_v2, dihedral_rescaled, &
                          dihedral_gradient
    use mm_units , only : rads
    use precision, only : dp
    use qmmm_pbc , only : pbc_displ_vector, reccel
    use sys      , only : die
    use units    , only : Ang

    implicit none
    integer  , intent(in)    :: natot
      !! Total number of atoms.
    integer  , intent(in)    :: istep
      !! Current MD step.
    integer  , intent(in)    :: istep_restr
      !! Current restrained optimization step.
    real(dp) , intent(in)    :: r_inp(3,natot)
      !! Atomic cartesian coordinates.
    real(dp) , intent(in)    :: ucell(3,3)
      !! Unit cell vectors (in a.u.)
    character, intent(in)    :: lattice_type
      !! Type of periodic lattice used.
    real(dp) , intent(inout) :: fdummy(3,natot)
      !! Dummy forces array over which we add restrained forces.


    integer  :: npi, irestr, iat, ratom(8), idst
    real(dp) :: fce, fnew(3,10), rp(3), r12, r32, rdst(3), rdst12(3), &
                rdst32(3), rdst34(3), fdihe(3), Fterm, rtot, req, kf, &
                scal, dscalar(3), dr12r32(3), amber_cell(3,3), &
                amber_kcell(3,3), rdst31(3), rdst41(3), rdst42(3)

    real(dp), pointer :: rclas(:,:)

    ! Convert units from atomic to angstrom.
    nullify( rclas )
    call re_alloc( rclas, 1, 3, 1, natot, 'rclas', 'restr_calc' )

    rclas(1:3,1:natot)  = r_inp(1:3,1:natot) / Ang
    amber_cell(1:3,1:3) = ucell(1:3,1:3)     / Ang
    call reccel( 3, amber_cell, amber_kcell, 0 )

    ! Loop over all restraints.
    do irestr = 1, nrestr
      fnew = 0.0_dp

      select case( harmrestr(irestr)%type )
      case (1) ! Simple distance restraint.
        ratom(1:2) = harmrestr(irestr)%atmidx(1:2)
        rdst(1:3)  = rclas(1:3,ratom(1)) - rclas(1:3,ratom(2))

        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst )
        rtot       = norm_v2( rdst )
        restr_r(irestr) = rtot
        req        = harmrestr(irestr)%r0
        kf         = harmrestr(irestr)%kf

        rdst(1:3) =  rdst(1:3) * 2.0_dp * kf * (rtot - req) / rtot
        fnew(:,1) = -rdst(1:3)
        fnew(:,2) =  rdst(1:3)

        ! Adding fnew to fdummy.
        fdummy(1:3,ratom(1)) = fdummy(1:3,ratom(1)) + fnew(1:3,1)
        fdummy(1:3,ratom(2)) = fdummy(1:3,ratom(2)) + fnew(1:3,2)

      case (2) ! Angle restraint.
        ratom(1:3)  = harmrestr(irestr)%atmidx(1:3)
        rdst12(1:3) = rclas(1:3,ratom(1)) - rclas(1:3,ratom(2))

        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst12 )

        rdst32(1:3) = rclas(1:3,ratom(3)) - rclas(1:3,ratom(2))
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst32 )

        rtot       = angle_v2( rdst12, rdst32 )
        restr_r(irestr) = rtot
        req        = harmrestr(irestr)%r0
        kf         = harmrestr(irestr)%kf

        scal = scalar_v2( rdst12, rdst32 )
        r12  = norm_v2( rdst12 )
        r32  = norm_v2( rdst32 )
        fce  = 2.0_dp * kf * (rtot - req) / rads

        scal = scal / ( r12 * r12 * r32 * r32 )

        ! Forces over atom 1:
        rdst(1:3) = scal * r32 * rdst12(:) / r12 - rdst32(:) / (r12 * r32)
        rdst(1:3) = rdst(1:3) * fce / sqrt( 1.0_dp - ( scal * r12 * r32 ) ** 2 )

        fnew(1:3,1) = rdst(1:3)

        ! Forces over atom 2:
        dscalar(1:3) = - ( rdst12(1:3) + rdst32(1:3) )
        dr12r32(1:3) = - ( r32 * rdst12(1:3) / r12 + r12 * rdst32(1:3) / r32 )

        rdst(1:3) = dscalar(1:3) / ( r12 * r32 ) - scal * dr12r32(1:3)
        rdst(1:3) = rdst(1:3)*fce / (sqrt( 1.0_dp - ( scal * r12 * r32 ) ** 2 ))

        fnew(1:3,2) = rdst(1:3)

        ! Forces over atom 3:
        fnew(1:3,3) = -fnew(1:3,1)

        ! Adding fnew to fdummy.
        do iat = 1, 3
          fdummy(1:3,ratom(iat))= fdummy(1:3,ratom(iat)) + fnew(1:3,iat)
        enddo

      case (3) ! Dihedral restraint.
        ratom(1:4) = harmrestr(irestr)%atmidx(1:4)
        rdst12(1:3) = rclas(1:3,ratom(1)) - rclas(1:3,ratom(2))
        rdst31(1:3) = rclas(1:3,ratom(3)) - rclas(1:3,ratom(1))
        rdst32(1:3) = rclas(1:3,ratom(3)) - rclas(1:3,ratom(2))
        rdst41(1:3) = rclas(1:3,ratom(4)) - rclas(1:3,ratom(1))
        rdst42(1:3) = rclas(1:3,ratom(4)) - rclas(1:3,ratom(2))
        rdst34(1:3) = rclas(1:3,ratom(4)) - rclas(1:3,ratom(3))

        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst12 )
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst32 )
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst34 )
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst31 )
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst41 )
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst42 )

        rtot = dihedral_rescaled( rdst12, rdst32, rdst34, rdst31, rdst41 )

        req = harmrestr(irestr)%r0
        kf  = harmrestr(irestr)%kf
        if ( (req < 90.0_dp ) .and. (rtot > 180.0_dp) ) rtot = rtot - 360.0_dp
        if ( (req > 270.0_dp) .and. (rtot < 180.0_dp) ) rtot = rtot + 360.0_dp

        Fterm = 2.0_dp * kf * (rtot - req) / rads

        do iat = 1, 4
          call dihedral_gradient( rdst12, rdst32, rdst34, rdst31, rdst41, &
                                  IAT, Fterm, fdihe )
          fnew(:,iat) = fdihe(:)
        enddo

        ! Adding fnew to fdummy.
        if ( (.not. ((rtot < 0.0_dp) .or. (rtot > 180.0_dp))) .or. &
             (rtot > 360_dp) ) then
          fnew(1:3,1:4) = (-1.0_dp) * fnew(1:3,1:4)
        elseif ( ((rtot > 180.0_dp) .and. (rtot < 360.0_dp)) .or. &
                 (rtot < 0.0_dp) ) then
          fnew(1:3,1:4) = fnew(1:3,1:4)
        else
          call die('restr_calc: Wrong dihedral angle value.')
        endif
        do iat = 1, 4
          fdummy(1:3,ratom(iat))= fdummy(1:3,ratom(iat)) + fnew(1:3,iat)
        enddo
        restr_r(irestr) = rtot

      case (4) ! Average distance to 2-4 atoms.
        npi = 0
        do iat = 2, 5
          if ( harmrestr(irestr)%atmidx(iat) /= 0 ) npi = npi +1
        enddo
        if ( npi == 0 ) &
          call die( 'restr_calc: atoms in average can not be zero.' )

        ratom(1:5)=harmrestr(irestr)%atmidx(1:5)
        rp(1:3) = ( rclas(1:3,ratom(2)) + rclas(1:3,ratom(3)) + &
                    rclas(1:3,ratom(4)) + rclas(1:3,ratom(5)) ) / npi

        rdst12(1:3) = rclas(1:3,ratom(1)) - rp(1:3)
        call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst12 )
        rtot = norm_v2( rdst12 )

        restr_r(irestr) = rtot
        req        = harmrestr(irestr)%r0
        kf         = harmrestr(irestr)%kf

        ! Force is calculated only over atom n°1.
        ! Is this physically correct though???
        fnew(1:3,1) = -rdst12(1:3) * 2.0_dp * kf * (rtot - req) / rtot

        fdummy(1:3,ratom(1)) = fdummy(1:3,ratom(1)) + fnew(1:3,1)

      case (5) ! Coefficients * positions.

        rtot = 0.0_dp
        do idst = 1, harmrestr(irestr)%ndists
          ratom(1) = harmrestr(irestr)%atmidx(2*(idst-1)+1)
          ratom(2) = harmrestr(irestr)%atmidx(2*(idst-1)+2)

          rdst12(1:3) = rclas(1:3,ratom(1)) - rclas(1:3,ratom(2))
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst12 )
          r12 = norm_v2( rdst12 )

          rtot = rtot + harmrestr(irestr)%weights(idst) * r12
        enddo

        restr_r(irestr) = rtot
        req        = harmrestr(irestr)%r0
        kf         = harmrestr(irestr)%kf

        ! We calculate forces by pairs.
        do idst = 1, harmrestr(irestr)%ndists
          ratom(1) = harmrestr(irestr)%atmidx(2*(idst-1)+1)
          ratom(2) = harmrestr(irestr)%atmidx(2*(idst-1)+2)

          rdst12(1:3) = rclas(1:3,ratom(1)) - rclas(1:3,ratom(2))
          call pbc_displ_vector( lattice_type, amber_cell, amber_kcell, rdst12 )
          r12 = norm_v2( rdst12 )

          rdst(1:3)   =  rdst12(1:3) * 2.0_dp * kf * (rtot - req) / r12
          fnew(1:3,1) = -rdst(1:3) * harmrestr(irestr)%weights(idst)
          fnew(1:3,2) =  rdst(1:3) * harmrestr(irestr)%weights(idst)

          ! Adding fnew to fdummy
          fdummy(1:3,ratom(1))= fdummy(1:3,ratom(1)) + fnew(1:3,1)
          fdummy(1:3,ratom(2))= fdummy(1:3,ratom(2)) + fnew(1:3,2)
        enddo
      case default

      end select ! Restraint type.
    enddo ! Loop ver restraints.
    call de_alloc( rclas, 'rclas', 'restr_calc' )

    ! Write outputs only in the first step of each restr_step.
    if ( istep > 1 ) return

    write(*,'(/,A)') 'restr_calc: harmonic restraints'
    write(*,'(/,A)') '-------------------------------'
    write(*,'(/,A)') ' '

    if ( istep_restr == 1 ) write(*,'(A,i4)')    'steps    :', restr_steps

    do irestr = 1, nrestr
      write(*,'(/,A)')     '--------------------'
      write(*,'(A,i4)')    'irestr   :', irestr
      write(*,'(A,F8.3)')  'targeted :', harmrestr(irestr)%r0
      write(*,'(A,F8.3)')  'current  :', restr_r(irestr)

      if ( istep_restr == 1 ) then
        write(*,'(A,F8.3)')  'final    :', harmrestr(irestr)%rfin
        write(*,'(A,F8.3)')  'constant :', harmrestr(irestr)%kf
        write(*,'(A,i4)')    'type     :', harmrestr(irestr)%type

        select case( harmrestr(irestr)%type )
        case (1)
          write(*,'(A,2i4)')   'atoms    :', &
            harmrestr(irestr)%atmidx(1), harmrestr(irestr)%atmidx(2)

        case (2)
          write(*,'(A,3i4)')   'atoms    :', &
            harmrestr(irestr)%atmidx(1), harmrestr(irestr)%atmidx(2), &
            harmrestr(irestr)%atmidx(3)

        case (3)
          write(*,'(A,4i4)')   'atoms    :', &
            harmrestr(irestr)%atmidx(1), harmrestr(irestr)%atmidx(2), &
            harmrestr(irestr)%atmidx(3), harmrestr(irestr)%atmidx(4)

        case (4)
          write(*,'(A,5i4)')   'atoms    :', &
            harmrestr(irestr)%atmidx(1), harmrestr(irestr)%atmidx(2), &
            harmrestr(irestr)%atmidx(3), harmrestr(irestr)%atmidx(4), &
            harmrestr(irestr)%atmidx(5)

        case (5)
          write(*,'(A,i4)')    'ndists   :', harmrestr(irestr)%ndists
          do idst = 1, harmrestr(irestr)%ndists
            write(*,'(A,F5.2)')  'weights  :', harmrestr(irestr)%weights(idst)
          enddo
          do idst = 1, 2*harmrestr(irestr)%ndists
            write(*,'(A,i4)')    'atoms    :', harmrestr(irestr)%atmidx(idst)
          enddo

        case default

        end select
      endif
    enddo ! Loop ver restraints.

    write(*,'(/,A)') ' '

  end subroutine restr_calc

  subroutine restr_update( )
    !! Updates r0 if necessary and prints energy.
    use precision, only: dp

    implicit none
    harmrestr(:)%r0 = harmrestr(:)%r0 + harmrestr(:)%dr

  end subroutine restr_update

  subroutine restr_deinit( )
    !! Destroys the restraint data structure and resets variables.
    implicit none
    integer :: irestr

    if ( allocated(harmrestr) ) then
      do irestr = 1, nrestr
        if ( allocated( harmrestr(irestr)%atmidx ) ) &
          deallocate( harmrestr(irestr)%atmidx )
        if ( allocated( harmrestr(irestr)%weights ) ) &
          deallocate( harmrestr(irestr)%weights )
      enddo

      deallocate( harmrestr )
    endif

    nrestr        = 0
    restr_steps   = 0
    hasRestraints = .false.
  end subroutine restr_deinit

end module restr_subs