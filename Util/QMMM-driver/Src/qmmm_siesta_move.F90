! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module qmmm_move_m

  private
  public :: qmmm_move
contains

  subroutine qmmm_move( istep, relaxd, na_qm, na_mm, qm_atoms, mm_atoms )
    use atomlist            , only : iza, amass
    use mm_topology         , only : mm_atom_t, qm_atom_t
    use m_energies          , only : Ekinion
    use m_mpi_utils         , only : barrier
    use m_check_walltime    , only : check_walltime
    use parallel            , only : IOnode
    use precision           , only : dp
    use siesta_geom         , only : na_u, na_s, isa, xa, xa_last, va, ucell, &
                                     ucell_last, scell, scell_last, vcell
    use siesta_options      , only : broyden_optim, fire_optim, idyn, nmove, &
                                     RelaxCellOnly, tp, strtol, varcel, ftol,&
                                     dxmax, usesavecg, iquench, ianneal,     &
                                     taurelax, bulkm, tt, dt, dx, ia1, mn, mpr
    use sys                 , only : die, message
    use units               , only : Kelvin
    use zmatrix             , only : lUseZmatrix

    ! MD-related modules
    use m_forces            , only : cfa, ntcon
    use m_kinetic           , only : vn, vpr, kn, kpr, tempion
    use m_steps             , only : inicoor, fincoor, istp
    use m_stress            , only : cstress
    use m_target_stress     , only : subtract_target_stress
    use m_dynamics          , only : nose, verlet2, npr, anneal, prdyn=>pr
    use m_broyden_optim     , only : broyden_optimizer
    use m_fire_optim        , only : fire_optimizer
    use m_zm_broyden_optim  , only : zm_broyden_optimizer
    use m_zm_fire_optim     , only : zm_fire_optimizer
    use m_cell_broyden_optim, only : cell_broyden_optimizer
    use m_cell_fire_optim   , only : cell_fire_optimizer

    ! Outputs
    use m_ioxv              , only : ioxv
    use siesta_cml          , only : cml_p, cmlEndStep, mainXML
    use qmmm_write          , only : qmmm_write_positions

    implicit none
    integer, intent(inout) :: istep
      !! Current MD step.
    logical, intent(out)   :: relaxd
      !! Whether the current structure is relaxed or not.
    integer, intent(in) :: na_qm
      !! Number of QM atoms.
    integer, intent(in) :: na_mm
      !! Number of MM atoms.
    type(qm_atom_t), intent(inout) :: qm_atoms(na_qm)
      !! Information for QM atoms.
    type(mm_atom_t), intent(inout) :: mm_atoms(na_mm)
      !! Information for MM atoms.

    integer            :: ix, iadispl, ixdispl
    real(dp)           :: P_int, tp_pr, eff_stress(3,3)
    logical            :: foundxv, foundzm, time_is_up
    character(len=40)  :: tmp_str
    integer, parameter :: iunit = 2
      !! Physical-units option for MD: 1=>(eV,Ang), 2=>(Ry,Bohr)

    real(dp), external :: volcel
    external           :: cgvc, cgvc_zmatrix, iozm, pxfflush, write_md_record

    call timer( 'qmmm_siesta_move', 1 )

    ! Save the last geometry.
    xa_last(1:3,1:na_s) = xa(1:3,1:na_s)
    ucell_last(1:3,1:3) = ucell(1:3,1:3)
    scell_last(1:3,1:3) = scell(1:3,1:3)

    ! Move atoms. This routine should only record checkpointing information.
    ! There is a pletora of output routines:
    !      -- ioxv
    !      -- iozm
    !      -- write_positions  (which calls write_struct and zm_canonical...)
    !      -- write_md_record

    ! In the following, it is better to put all the logic in the individual
    ! blocks, avoiding branching and re-checks of the 'idyn' variable (which,
    ! by the way, should have symbolic values instead of hardwired
    ! numerical ones)
    relaxd = .false.
    select case( idyn )
    case(0)
      ! The original "relaxation" case, but note that it also covers
      ! pure single-point calculations (with MD.NumCGsteps = 0)

      ! Pure single-point calculations will not enter this block.
      if ( nmove /= 0 ) then    ! That is, if requesting "CG" steps
        ! Here we want checkpointing, except if the structure is already relaxed
        ! Note that these routines do not update the coordinates if
        ! the force/stress criterion for relaxation is satisfied

        if ( RelaxCellOnly ) then
          if ( broyden_optim ) then
            call cell_broyden_optimizer( na_u, xa, ucell, cstress, tp, strtol,&
                                         varcel, relaxd )
          elseif ( fire_optim ) then
            call cell_fire_optimizer( na_u, xa, ucell, cstress, tp, strtol,&
                                      varcel, relaxd )
          else
            call die( "Cell-only optimization must be either Broyden or FIRE." )
          endif

        else ! Coordinate relaxation (and maybe cell)
          if ( lUseZmatrix ) then
            if ( broyden_optim ) then
              call zm_broyden_optimizer( na_u, xa, ucell, cstress, tp, &
                                         strtol, varcel, relaxd )
            elseif (fire_optim) then
              call zm_fire_optimizer( na_u, xa, cfa, ucell, cstress, dxmax, &
                                      tp, ftol, strtol, varcel, relaxd )
            else
              call cgvc_zmatrix( na_u, xa, cfa, ucell, cstress, dxmax, tp, &
                                 ftol, strtol, varcel, relaxd, usesavecg )
            endif

          else
            if ( broyden_optim ) then
              call broyden_optimizer( na_u, xa, cfa, ucell, cstress, tp, &
                                      ftol, strtol, varcel, relaxd )
            elseif ( fire_optim ) then
              call fire_optimizer( na_u, xa, cfa, ucell, cstress, dxmax, &
                                   tp, ftol, strtol, varcel, relaxd )
            else
              call cgvc( na_u, xa, cfa, ucell, cstress, dxmax, tp, ftol, &
                         strtol, varcel, relaxd, usesavecg )
            endif
          endif
        endif ! RelaxCellOnly

        if ( relaxd ) then
          ! Will not call ioxv et al in the block below, so that:
          !  - The XV file will contain the xa coords written in the previous step.
          !  - The .STRUCT_NEXT_ITER file would be that written in the previous step
          !    (if any), and it will contain the same coords as the .STRUCT_OUT.
        else
          call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )

          if ( lUseZmatrix ) call iozm( 'write', ucell, vcell, xa, foundzm )
            ! This writes the "next_iter" STRUCT and canonical zmatrix files
            call qmmm_write_positions( .true. )

            ! Save atomic positions and velocities accumulatively and
            ! accumulate coor in Xmol file for animation
            call write_md_record( istep )
        endif
      endif

    case(1)  ! Micro-canonical MD with fixed cell (Verlet)
      ! Check convergence for quenching runs (which are really relaxations)
      ! Avoid calling the routine if relaxed, so that the coordinates
      ! are not changed.
      if ( iquench /= 0 ) then
        relaxd = all( abs(cfa(1:3,1:na_u)) < ftol )
        ! We might want to fall back on a different relaxation scheme when we
        ! reach the slow-moving part (check temp_ion?;  check cfa behavior?)
      endif
      if ( .not. relaxd ) then
        call verlet2( istp, iunit, iquench, na_u, cfa, dt, amass, ntcon, &
                      va, xa, Ekinion, tempion )
        if ( IOnode ) &
          write(6,'(/,a,f12.3,a)') 'siestaqmmm: Temp_ion =', tempion, ' K'

        call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )

        call qmmm_write_positions( .true. )
        call write_md_record( istep )
      endif


    case (2) ! Canonical MD with fixed cell (Nose-Hoover)
      call nose( istp, iunit, na_u, cfa, tt, dt, amass, mn, ntcon, &
                 va, xa, Ekinion, kn, vn, tempion )
      if ( IOnode ) &
        write(6,'(/,a,f12.3,a)') 'siestaqmmm: Temp_ion =', tempion, ' K'

      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )
      call qmmm_write_positions( .true. )
      call write_md_record( istep )

    case (3) ! Micro-canonical variable-cell MD (Parrinello-Rahman)
      ! Check convergence for quenching runs (which are really relaxations)
      ! Avoid calling the routine if relaxed, so that the coordinates
      ! are not changed.
      if ( iquench /= 0 ) then
        ! Allow a general target stress in this relaxation mode
        call subtract_target_stress( cstress, eff_stress )

        relaxd  = all( abs(cfa(1:3,1:na_u)) < ftol   ) .and. &
                  all( abs(eff_stress)      < strtol )
        cstress = eff_stress
        tp_pr   = 0.0_dp  ! As we are already passing the modified stress
      else
        ! Standard MD mode
        relaxd = .false.
        tp_pr = tp      ! Keep the original target pressure
      endif

      if ( .not. relaxd ) then
        call prdyn( istp, iunit, iquench, na_u, cfa, cstress, tp_pr, dt,  &
                    amass, mpr, ntcon, va, xa, vcell, ucell, Ekinion, kpr,&
                    vpr, tempion, P_int )

        if ( IOnode ) then
          write(6,'(/,a,f12.3,a)') 'siestaqmmm: E_kin PR =', kpr/Kelvin, ' K'
          write(6,'(/,a,f12.3,a)') 'siestaqmmm: Temp_ion =', tempion   , ' K'
        endif

        call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )
        call qmmm_write_positions( .true. )
        call write_md_record( istep )
      endif

    case (4)  ! Canonical variable-cell MD (Nose-Parrinello-Rahman)
      call npr( istp, iunit, na_u, cfa, cstress, tp, tt, dt, amass, mn, mpr,&
                ntcon, va, xa, vcell, ucell, Ekinion, kn, kpr, vn, vpr,     &
                tempion, P_int )

      if (IOnode) &
        write(6,'(/,a,f12.3,a)') 'siestaqmmm: Temp_ion =', tempion, ' K'

      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )
      call qmmm_write_positions( .true. )
      call write_md_record( istep )

    case (5)  ! Annealings
      call anneal( istp, iunit, ianneal, taurelax, bulkm, na_u, cfa, &
                   cstress, tp, tt, dt, amass, ntcon, va, xa, ucell, &
                   Ekinion, tempion, P_int )
      if (IOnode) &
        write(6,'(/,a,f12.3,a)') 'siestaqmmm: Temp_ion =', tempion, ' K'

      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )
      call qmmm_write_positions( .true. )
      call write_md_record( istep )

    case (6) ! Force-constant-matrix calculation
      ! Output of coordinates is meaningless here in general, except maybe to
      ! checkpoint a calculation to be restarted later, by reading the XV file
      ! and choosing inicoor and fincoor appropriately.

      ! Save current (displaced) atomic positions
      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv )

      ! Undo the last atom displacement
      if ( istep > inicoor ) then ! inicoor is always 0 for FC.
        ix      = mod( istep-1, 6 ) +1
        iadispl = ia1 + ( istep - ix +1 ) / 6
        ixdispl = ( ix - mod(ix-1,2) +1 ) / 2
        xa(ixdispl,iadispl) = xa(ixdispl,iadispl) - dx
      endif

      ! Displace atom by dx. NOTE: MOVE is before istep increment, hence
      ! fincoor == istep will be the last step and the initial geometry
      ! should be retained.
      if ( istep < fincoor ) then ! fincoor = (ia2-ia1+1)*6
        ix      = mod( istep-1, 6 ) +1
        iadispl = ia1 + ( istep +2 - ix ) / 6
        ixdispl = ( ix - mod(ix-1,2) +1 ) / 2
        dx      = -dx
        xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
      endif

    case default
      call die( 'siestaqmmm: Wrong MD type of run.' )
    end select

    do ix = 1, na_qm
      qm_atoms(ix)%r(:) = xa(:,ix)
    enddo
    do ix = 1, na_mm
      mm_atoms(ix)%r(:) = xa(:,ix+na_qm)
    enddo

    ! Output memory use at the end of this geometry step
    if ( cml_p ) call cmlEndStep( mainXML )

    call timer( 'qmmm_siesta_move', 2 )

    ! End of one MD step - flush stdout
    if ( ionode ) call pxfflush( 6 )

    ! Check whether we are short of time to continue
    call check_walltime( time_is_up )
    if ( time_is_up ) then
      ! Do any other bookeeping not done by "die".
      call timer( 'all', 2 )
      call timer( 'all', 3 )

      if ( .not. relaxd ) then
        call message( 'WARNING', 'Coordinate relaxation not converged' )
        write( tmp_str, "(a,i5)" ) 'Geometry step: ', istep
        call message( ' (info)', trim(tmp_str) )
      endif

      call barrier( )
      call die( "OUT_OF_TIME: Time is up at end of geometry step." )
    endif

  end subroutine qmmm_move

end module qmmm_move_m
