module m_get_kpoints_scale
  public :: get_kpoints_scale
CONTAINS
  !
  ! Determines the reciprocal-lattice vectors to use
  ! as basis for the k-point coordinates entered in
  ! the bands and/or wavefunction computation blocks
  !
  subroutine get_kpoints_scale(block_string,rcell,ierr)

    use fdf,            only: leqi, fdf_get
    use parallel,       only: IOnode
    use siesta_geom,    only: ucell
    use units,          only: Pi, Ang

    implicit none

    integer, parameter :: dp = selected_real_kind(10,100)

    character(len=*), intent(in) :: block_string  ! BandLinesScale, etc
    real(dp), intent(out)        :: rcell(3,3)    ! Reciprocal unit cell
    integer,  intent(out)        :: ierr          ! Error code: 0: success

    character(len=30)  :: scale
    real(dp)           :: alat, sum, norm
    integer            :: i

    ierr = 0

    scale = fdf_get( trim(block_string), 'pi/a' )

    if (leqi(scale,'pi/a')) then
       ! Check if the cell is orthogonal.
       do i = 1, 3
         sum  = abs(ucell(i,1) + ucell(i,2) + ucell(i,3))
         norm = sqrt( ucell(i,1)**2 + ucell(i,2)**2 + ucell(i,3)**2 )

         if ( abs( sum - norm ) > 1e-9_dp ) then
           ierr = -2
           if (IOnode) &
             write(6,'(a,/,a)') 'ERROR: pi/a option should only be used for '//&
                                'orthorhombic cells. No calculation performed.'
           return
         endif
       enddo

       alat = fdf_get( 'LatticeConstant', Ang, 'Bohr' )
       if ( abs(alat) < 1e-16_dp) then
          if (IOnode ) then
             write(6,'(a)') 'ERROR: Lattice constant required for ' // &
                  trim(block_string) // ". No calculation performed."
          endif
          ierr = -1
       else
          if (IOnode ) then
            write(6,'(a,f12.6,a)') 'Using LatticeConstant from fdf file for '//&
                  trim(block_string) // ":", alat, " Bohr"
            write(6,'(a,f12.6,a)') &
                  'Beware any cell changes by the end of the run.'
          endif
          rcell(:,:) = 0.0_dp
          do i = 1, 3
             rcell(i,i) = pi/alat
          enddo
       endif

    elseif (leqi(scale,'ReciprocalLatticeVectors')) then
       !
       ! Bands/wfns initialization takes place at the beginning
       ! of the run, in 'siesta_init'. The ucell stored in the
       ! siesta_geom module at this point is exactly what would
       ! be read using any of the input options in effect.
       !
       call reclat( ucell, rcell, 1 )
       if (IOnode ) then
         write(6,'(a,f12.6,a)') 'Using current reciprocal lattice vectors for '//&
            trim(block_string)
         write(6,'(a,f12.6,a)') 'Beware any cell changes by the end of the run.'
       endif

    else
       ierr = -10
       if (IOnode) then
          write(6,'(a,/,2a,/,a)')                                      &
             'WARNING: Invalid value for ' // trim(block_string),   &
             'Allowed values are pi/a and',                         &
             ' ReciprocalLatticeVectors',                     &
             'No calculation performed'
       endif
    endif

  end subroutine get_kpoints_scale
end module m_get_kpoints_scale
