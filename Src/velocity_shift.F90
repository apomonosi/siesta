!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!< Module implementing an eigenspectrum shift according to the group velocity of an
!< eigenstate.
!<
!< This is developed in two parts:
!< 1. The diagkp_velocity code implements the diagonalization routines + the
!<    velocity calculations.
!<    The velocities are calculated using:
!<      < psi_i | dH_kx/ky/kz - epsilon_i dS_kx/ky/kz | psi_i >
!< 2. A shift of the eigenvalue epsilon_i depending on the applied bias and direction.
!<    The group velocity is projected onto a "field" direction and if the scalar projection
!<    is positive (> 0) the eigenspectrum is shifted V / 2 down, to fill it, while negative
!<    projections are shifted V / 2 up, to empty it.
module velocity_shift_m

  use precision, only: dp

  implicit none

  save
  private

  !< Logical to control whether the band-shifts are being performed.
  logical :: use_velocity_shift = .false.

  !< Logical to control whether the current is calculated
  logical :: calc_velocity_current = .true.

  !< The halve bias applied to the bands (we shift V/2 up and V/2 down)
  real(dp) :: velocity_h_bias = 0._dp

  !< The direction along which the bands will be displaced.
  real(dp) :: velocity_dir(3) = 0._dp

  !< The projection factor required before the velocity-shift is used
  real(dp) :: velocity_tolerance = 0._dp

  public :: use_velocity_shift
  public :: calc_velocity_current
  public :: read_velocity_shift
  public :: velocity_dir
  public :: velocity_shift
  public :: velocity_results
  interface velocity_results
    module procedure velocity_results_pol
    module procedure velocity_results_nc
  end interface velocity_results
  public :: velocity_results_print
  interface velocity_results_print
    module procedure velocity_results_print_pol
    module procedure velocity_results_print_nc
  end interface velocity_results_print

  ! Parameters for conversion
  real(dp), parameter :: Coulomb = 1.6021766208e-19_dp
  real(dp), parameter :: hbar_Rys = 4.8377647940592375e-17_dp

contains

  subroutine read_velocity_shift()
    use fdf
    use parallel, only: IONode
    use units, only: eV, Ang
    use siesta_cml
    use intrinsic_missing, only: VNORM
    use m_spin, only: TrSym
    use m_cite, only: add_citation

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    real(dp) :: velocity_bias

    character(len=*), parameter :: form_afa = '("redata: ",a,t53,"= ",f10.4,a)'
    character(len=*), parameter :: form_ava = '("redata: ",a,t53,"= ",e10.4,a)'
    character(len=*), parameter :: form_avvv = '("redata: ",a,t53,"=",3(tr1,e10.4))'
    character(len=*), parameter :: form_al = '("redata: ",a,t53,"= ",l1)'

    ! Read in the velocity related options
    velocity_bias = fdf_get('BulkBias.V', 0._dp, 'Ry')
    velocity_bias = fdf_get('BulkBias.Voltage', velocity_bias, 'Ry')

    ! The velocity tolerance in atomic units
    velocity_tolerance = fdf_get('BulkBias.Tolerance', 1.e-15_dp)

    ! TODO libfdf
    ! Once libfdf 0.1 is merged with trunk we can use fdf_islist for reals as well
!    if ( fdf_islist('BulkBias.Direction') ) then
!      ! In this case the BulkBias.V is determining the size of the field
!      call fdf_list('BulkBias.Direction', 3, velocity_dir)
    if ( fdf_block('BulkBias.Direction', bfdf) ) then
      ! Read first line
      if ( .not. fdf_bline(bfdf, pline) ) then
        call die('velocity_shift: could not read BulkBias.Direction block')
      end if

      ! Read in the field
      velocity_dir(1) = fdf_bvalues(pline, 1)
      velocity_dir(2) = fdf_bvalues(pline, 2)
      velocity_dir(3) = fdf_bvalues(pline, 3)

      call fdf_bclose(bfdf)
    else
      velocity_dir = 0._dp
    end if

    ! Determine whether we should use the bulkbias
    ! Only for applied bias above 1e-6 eV
    use_velocity_shift = abs(velocity_bias) > 1.e-6_dp * eV

    calc_velocity_current = fdf_get('BulkBias.Current', .true.)

    ! Immediately return if not used
    if ( .not. use_velocity_shift ) return

    if (IOnode) call add_citation("10.1103/PhysRevB.106.155401")

    ! Ensure that all k-points are not using TRS
    call fdf_overwrite("TimeReversalSymmetryForKpoints false")
    TrSym = .false.

    if ( VNORM(velocity_dir) == 0._dp ) then
      call die('velocity_shift: BulkBias.Direction direction has not been specified!')
    end if

    ! Normalize the direction
    velocity_dir = velocity_dir / VNORM(velocity_dir)

    if ( IONode ) then
      write(*,form_afa) 'Bulk velocity bias ', velocity_bias / eV, ' Volts'
      write(*,form_avvv) 'Bulk velocity bias direction', velocity_dir
      write(*,form_ava) 'Bulk velocity bias projection tolerance', velocity_tolerance * 1.e-12_dp / Ang / hbar_Rys, ' Ang/ps'
      write(*,form_al) 'Bulk velocity bias current calculation', calc_velocity_current
    end if

    if ( cml_p ) then
      call cmlStartPropertyList(xf=mainXML, title='BulkBias', &
          dictRef='siesta:BulkBias')
      call cmlAddProperty( xf=mainXML, value=velocity_bias/eV, &
          dictRef='siesta:BulkBias.V', units='siestaUnits:eV')
      call cmlAddProperty( xf=mainXML, value=velocity_tolerance * 1.e-12_dp / Ang / hbar_Rys, &
          dictRef='siesta:BulkBias.Tolerance', units='siestaUnits:Ang/ps')
      call cmlAddProperty( xf=mainXML, value=velocity_dir, &
          dictRef='siesta:BulkBias.Direction')
      call cmlEndPropertyList(mainXML)
    end if

    ! Only take half the bias (so we dont have to divide by 2 all the time)
    velocity_h_bias = velocity_bias * 0.5_dp

  end subroutine read_velocity_shift

  !< Shift the eigenspectrum according to the velocities
  subroutine velocity_shift(sign, ne, e, v)
    integer, intent(in) :: sign
    integer, intent(in) :: ne
    real(dp), intent(inout) :: e(ne)
    real(dp), intent(in) :: v(3, ne)

    integer :: ie
    real(dp) :: p

    if ( sign > 0 ) then

      do ie = 1, ne

        p = dot_product(velocity_dir, v(:,ie))
        if ( p > velocity_tolerance ) then
          ! The velocity is along the direction of the potential drop
          ! This means that the electron will be pushed in the same direction and filled
          ! up faster.
          e(ie) = e(ie) - velocity_h_bias
        else if ( p < - velocity_tolerance ) then
          e(ie) = e(ie) + velocity_h_bias
        end if

      end do

    else

      ! Shift back the eigenvalues

      do ie = 1, ne
        p = dot_product(velocity_dir, v(:,ie))
        if ( p > velocity_tolerance ) then
          e(ie) = e(ie) + velocity_h_bias
        else if ( p < - velocity_tolerance ) then
          e(ie) = e(ie) - velocity_h_bias
        end if
      end do

    end if

  end subroutine velocity_shift

  !< Calculate current from an eigenspectrum according to the velocities and calculate current
  subroutine velocity_results_pol(ne, e, o, v, Ef, w, Temp, results)
    use m_fermid, only: stepf
    integer, intent(in) :: ne
    ! The eigenvalues `e` *must* be un-shifted
    real(dp), intent(in) :: e(ne), o(ne), v(3,ne), Ef, w, Temp
    real(dp), intent(inout) :: results(5)

    integer :: ie
    real(dp) :: p, oc, lresults(3), EEf

    lresults(:) = 0._dp
    do ie = 1, ne
      eEf = e(ie) - Ef

      ! velocity_dir is normalized, so this gives the correct
      ! scalar value.
      p = dot_product(velocity_dir, v(:,ie))
      if ( p > velocity_tolerance ) then
        oc = stepf((eEf - velocity_h_bias)/Temp) - stepf(eEf/Temp)

        ! Calculate number of electrons for positive/negative direction
        results(4) = results(4) + o(ie)

      else if ( p < - velocity_tolerance ) then
        oc = stepf((eEf + velocity_h_bias)/Temp) - stepf(eEf/Temp)

        ! Calculate number of electrons for positive/negative direction
        results(4) = results(4) - o(ie)

      else
        oc = 0._dp
        results(5) = results(5) + o(ie)
      end if

      ! this is velocity projected onto bias direction
      lresults(1) = lresults(1) + v(1, ie) * oc
      lresults(2) = lresults(2) + v(2, ie) * oc
      lresults(3) = lresults(3) + v(3, ie) * oc

    end do

    results(1) = results(1) + lresults(1) * w
    results(2) = results(2) + lresults(2) * w
    results(3) = results(3) + lresults(3) * w

  end subroutine velocity_results_pol

  !< Calculate current from an eigenspectrum according to the velocities and calculate current
  subroutine velocity_results_nc(ne, e, o, v, S, Ef, w, Temp, results)
    use m_fermid, only: stepf
    integer, intent(in) :: ne
    ! The eigenvalues `e` *must* be un-shifted
    real(dp), intent(in) :: e(ne), o(ne), v(3,ne), S(3,ne), Ef, w, Temp
    real(dp), intent(inout) :: results(3,0:4) ! never reference last entry

    integer :: ie
    real(dp) :: p, oc, lresults(3,0:3), lv(3), eEf

    lresults(:,:) = 0._dp
    do ie = 1, ne
      eEf = e(ie) - Ef

      ! velocity_dir is normalized, so this gives the correct
      ! scalar value.
      p = dot_product(velocity_dir, v(:,ie))
      if ( p > velocity_tolerance ) then
        oc = stepf((eEf - velocity_h_bias)/Temp) - stepf(eEf/Temp)

        ! Calculate number of electrons for positive/negative direction
        results(1,4) = results(1,4) + o(ie)

      else if ( p < - velocity_tolerance ) then
        oc = stepf((eEf + velocity_h_bias)/Temp) - stepf(eEf/Temp)

        ! Calculate number of electrons for positive/negative direction
        results(1,4) = results(1,4) - o(ie)

      else
        oc = 0._dp
        results(2,4) = results(2,4) + o(ie)
      end if

      ! this is velocity projected onto bias direction
      lv(:) = v(:, ie) * oc
      ! lresults =
      !  (direction, spin)
      lresults(:,0) = lresults(:,0) + lv(:)
      lresults(:,1) = lresults(:,1) + S(1,ie) * lv(:)
      lresults(:,2) = lresults(:,2) + S(2,ie) * lv(:)
      lresults(:,3) = lresults(:,3) + S(3,ie) * lv(:)

    end do

    do ie = 0, 3
      results(:,ie) = results(:,ie) + lresults(:,ie) * w
    end do

  end subroutine velocity_results_nc

  subroutine velocity_results_print_pol(spin, ucell, cell_periodic, results)
    use parallel, only: IONode
    use intrinsic_missing, only: VNORM
    use units, only: eV, Ang
    use t_spin, only: tSpin
    use m_energies, only: E_bulk_bias

    type(tSpin), intent(in) :: spin
    real(dp), intent(in) :: ucell(3,3)
    logical, intent(in) :: cell_periodic(3)
    real(dp), intent(in) :: results(5,spin%spinor)

    real(dp), external :: volcel
    real(dp) :: vcross(3), fac
    real(dp) :: I(4,spin%spinor)
    character(len=8) :: suffix

    ! Calculate the correction to the total energy due to the
    ! different populations.
    ! This will only work in case the bias is always +- V/2
    ! results(4) == dq
    E_bulk_bias = - velocity_h_bias * sum(results(4,:))

    ! Quick escape
    if ( .not. IONode ) return

    ! Current I is in [e Bohr Ry], then convert to [e Bohr/s]
    ! Note this equation also holds for NC/SOC since there spinor will
    ! be 2.
    I(2:4,:) = results(1:3,:) * 2._dp / spin%spinor
    I(1,:) = matmul(velocity_dir, I(2:4,:))

    select case ( count( cell_periodic(:) ) )
    case ( 3 )
      ! We are dealing with the volume
      fac = 1._dp / volcel(ucell) * Ang ** 2
      suffix = 'uA/Ang^2'
    case ( 2 )
      if ( .not. cell_periodic(1) ) then
        call cross(ucell(:, 2), ucell(:, 3), vcross)
      else if ( .not. cell_periodic(2) ) then
        call cross(ucell(:, 1), ucell(:, 3), vcross)
      else if ( .not. cell_periodic(3) ) then
        call cross(ucell(:, 1), ucell(:, 2), vcross)
      end if
      fac = 1._dp / VNORM(vcross) * Ang
      suffix = 'uA/Ang'
    case ( 1 )
      if ( cell_periodic(1) ) then
        fac = 1._dp / VNORM(ucell(:,1))
      else if ( cell_periodic(2) ) then
        fac = 1._dp / VNORM(ucell(:,2))
      else if ( cell_periodic(3) ) then
        fac = 1._dp / VNORM(ucell(:,3))
      end if
      suffix = 'uA'
    end select
    ! Current I is in [e Bohr Ry], then convert to [e Bohr/s]
    fac = fac * Coulomb / hbar_Rys * 1.e6_dp
    I(:,:) = I(:,:) * fac

    ! Write out the current along the bulk-bias direction and each of the other ones
    select case ( spin%spinor )
    case ( 1 )
      write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
          'bulk-bias: {v}    , |v|    [',trim(suffix),'] {', I(2:4,1), " }   ", I(1,1)

      write(*,'(tr5,a,2(tr1,e15.7))') 'bulk-bias: dq , q0 = ', results(4:5,1)
    case ( 2 )
      write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
          'bulk-bias: {v}    , |v|    [',trim(suffix),'] {', I(2:4,1)+I(2:4,2), " }   ", I(1,1)+I(1,2)
      write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
          'bulk-bias: {v}_Sz , |v|_Sz [',trim(suffix),'] {', I(2:4,1)-I(2:4,2), " }   ", I(1,1)-I(1,2)

      write(*,'(tr5,a,2(tr1,e15.7))') 'bulk-bias: dq , q0 = ', results(4:5,1)+results(4:5,2)
    end select

  end subroutine velocity_results_print_pol

  subroutine velocity_results_print_nc(ucell, cell_periodic, results)
    use parallel, only: IONode
    use intrinsic_missing, only: VNORM
    use units, only: eV, Ang
    use m_energies, only: E_bulk_bias

    real(dp), intent(in) :: ucell(3,3)
    logical, intent(in) :: cell_periodic(3)
    real(dp), intent(in) :: results(3,0:4)

    real(dp), external :: volcel
    real(dp) :: vcross(3), fac
    real(dp) :: I, IS(3,0:3)
    character(len=8) :: suffix

    ! Calculate the correction to the total energy due to the
    ! different populations.
    ! This will only work in case the bias is always +- V/2
    ! results(4) == dq
    E_bulk_bias = - velocity_h_bias * results(1,4)

    ! Quick escape
    if ( .not. IONode ) return

    select case ( count( cell_periodic(:) ) )
    case ( 3 )
      ! We are dealing with the volume
      fac = 1._dp / volcel(ucell) * Ang ** 2
      suffix = 'uA/Ang^2'
    case ( 2 )
      if ( .not. cell_periodic(1) ) then
        call cross(ucell(:, 2), ucell(:, 3), vcross)
      else if ( .not. cell_periodic(2) ) then
        call cross(ucell(:, 1), ucell(:, 3), vcross)
      else if ( .not. cell_periodic(3) ) then
        call cross(ucell(:, 1), ucell(:, 2), vcross)
      end if
      fac = 1._dp / VNORM(vcross) * Ang
      suffix = 'uA/Ang'
    case ( 1 )
      if ( cell_periodic(1) ) then
        fac = 1._dp / VNORM(ucell(:,1))
      else if ( cell_periodic(2) ) then
        fac = 1._dp / VNORM(ucell(:,2))
      else if ( cell_periodic(3) ) then
        fac = 1._dp / VNORM(ucell(:,3))
      end if
      suffix = 'uA'
    end select

    ! Current I is in [e Bohr Ry], then convert to [e Bohr/s]
    fac = fac * Coulomb / hbar_Rys * 1.e6_dp
    IS(:,:) = results(:,0:3) * fac

    ! Write out the current along the bulk-bias direction and each of the other ones
    I = dot_product(velocity_dir, IS(:,0))
    write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
        'bulk-bias: {v}    , |v|    [',trim(suffix),'] {', IS(:,0), " }   ", I
    I = dot_product(velocity_dir, IS(:,1))
    write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
        'bulk-bias: {v}_Sx , |v_Sx| [',trim(suffix),'] {', IS(:,1), " }   ", I
    I = dot_product(velocity_dir, IS(:,2))
    write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
        'bulk-bias: {v}_Sy , |v_Sy| [',trim(suffix),'] {', IS(:,2), " }   ", I
    I = dot_product(velocity_dir, IS(:,3))
    write(*,'(tr5,3a,3(tr1,e15.7),a,e15.7)') &
        'bulk-bias: {v}_Sz , |v_Sz| [',trim(suffix),'] {', IS(:,3), " }   ", I

    write(*,'(tr5,a,2(tr1,e15.7))') 'bulk-bias: dq , q0 = ', results(1:2,4)

  end subroutine velocity_results_print_nc

end module velocity_shift_m
