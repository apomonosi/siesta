!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_transiesta

  use precision, only : dp

  use files, only : slabel

  use m_ts_sparse, only : ts_sparse_init
  use m_ts_method, only : ts_method
  use m_ts_method, only : ts_A_method, TS_BTD_A_COLUMN

  use m_ts_method, only : TS_FULL, TS_BTD
#ifdef SIESTA__MUMPS
  use m_ts_method, only : TS_MUMPS
#endif

  use m_ts_tri_init, only : ts_tri_init

  use m_ts_fullg
  use m_ts_fullg_spinor
  use m_ts_fullk
  use m_ts_fullk_spinor

  use m_ts_trig
  use m_ts_trig_spinor
  use m_ts_trik
  use m_ts_trik_spinor

#ifdef SIESTA__MUMPS
  use m_ts_mumpsg
  use m_ts_mumpsk
#endif

  implicit none

  public :: transiesta
  private

contains

  subroutine transiesta(TSiscf,nspin, &
       sp_dist, sparse_pattern, &
       no_aux_cell, ucell, nsc, isc_off, no_u, na_u, lasto, xa, n_nzs, &
       H, S, DM, EDM, Ef, &
       Qtot)

    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use units, only : eV
    use alloc, only : re_alloc, de_alloc

    use parallel, only : IONode
    use m_os, only: file_delete

    use class_OrbitalDistribution
    use class_Sparsity

    use ts_kpoint_scf_m, only : ts_kpoint_scf, ts_gamma_scf

    use ts_electrode_m

    use m_ts_options, only : N_Elec, Elecs
    use m_ts_options, only : IsVolt, Calc_Forces

    use m_ts_options, only : BTD_method

    use m_ts_contour_eq , only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use ts_charge_m
    use ts_dq_m
    
    use m_ts_gf, only : read_Green
    use m_interpolate

    use m_energies, only: NEGF_DE

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in)  :: TSiscf
    integer, intent(in)  :: nspin
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    logical, intent(in)  :: no_aux_cell
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: nsc(3), no_u, na_u
    integer, intent(in) :: isc_off(3,product(nsc))
    integer, intent(in)  :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: H(n_nzs,nspin), S(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,min(nspin,4))
    real(dp), intent(in) :: Qtot
    real(dp), intent(inout) :: Ef

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
! ************************************************************

    ! * local variables
    integer :: iEl, NEn, no_used, no_used2
    logical :: converged, dq_run, truncated
    ! In case of Fermi-correction, we save the previous steps
    ! and do a spline interpolation... :)
    integer :: N_F, i_F, ioerr
    real(dp) :: dEf
    real(dp), pointer :: Q_Ef(:,:) => null()

    ! Open GF files...
    ! Read-in header of Green functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    ! Initialize the DE_NEGF energy
    ! Note that this will *only* be updated in case V /= 0.
    ! This energy corresponds to the non-equilibrium energies:
    !   e \sum_i N_i * \mu_i
    NEGF_DE = 0._dp

    if ( TSiscf == 1 ) then
      ! We need to initialize TRANSIESTA
      
      call timer('TS_init',1)

      if ( IONode ) then
        if ( file_delete('TS_FERMI') ) then
          write(*,'(/a/)') 'ts: Deleted TS_FERMI'
        end if
      end if
      
      call ts_sparse_init(slabel, IsVolt, N_Elec, Elecs, &
          ucell, nsc, na_u, xa, lasto, sp_dist, sparse_pattern, no_aux_cell, &
          isc_off)
      
      if ( ts_method == TS_BTD ) then
        ! initialize the tri-diagonal partition
        call ts_tri_init( sp_dist, sparse_pattern , N_Elec, &
            Elecs, IsVolt, ucell, na_u, xa, lasto ,nsc, isc_off, &
            BTD_method )
      end if

      ! print out estimated memory usage...
      call ts_print_memory(ts_gamma_scf, nspin)

      call ts_charge_print(N_Elec,Elecs, Qtot, sp_dist, sparse_pattern, &
          nspin, n_nzs, DM, S, TS_Q_INFO_FULL)
      
      if ( .not. Calc_Forces ) then
        if ( IONode ) then
          write(*,'(a)') 'transiesta: *** The forces are NOT updated ***'
          write(*,'(a)') 'transiesta: ***  Will set the energy-density matrix to 0!  ***'
        end if
        EDM(:,:) = 0._dp
      end if

      call timer('TS_init',2)
      
    end if

    
    call timer('TS',1)

    
    ! Total number of energy-points...
    NEn = N_Eq_E() + N_nEq_E()

    ! in case the file-descriptor is negative it basically 
    ! means "out-of-core" calculation.
    allocate(uGF(N_Elec),nq(N_Elec))
    uGF(:) = -1
    do iEl = 1 , N_Elec

       ! Calculate number of Bloch expansion k-points
       nq(iEl) = Elecs(iEl)%Bloch%size()

       ! Allocate the electrode quantities
       nullify(Elecs(iEl)%HA,Elecs(iEl)%SA,Elecs(iEl)%Gamma)

       ! We allocate for once as much space as needed,

       ! Allocate the non-repeated hamiltonian and overlaps...
       no_used = Elecs(iEl)%no_used
       if ( Elecs(iEl)%pre_expand > 1 ) then ! > 1 also expand H, S before writing
          no_used = Elecs(iEl)%device_orbitals()
          nq(iEl) = 1
       end if
       ! Adjust number of orbitals for calculation with non-collinear spin
       if (nspin > 2) no_used = 2 * no_used

       if ( IsVolt .or. .not. Elecs(iEl)%Bulk ) then
          ! If we using bulk electrodes, we need not the Hamiltonian, 
          ! nor the overlap...
          call re_alloc(Elecs(iEl)%HA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
          call re_alloc(Elecs(iEl)%SA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
       end if

       no_used = Elecs(iEl)%device_orbitals()
       if ( IsVolt ) then
          ! We need Gamma's with voltages (now they are both GAA and GammaT)
          no_used2 = no_used
       else 
          ! This is only for having space for GA
          if ( Elecs(iEl)%pre_expand > 0 ) then
             no_used2 = no_used
          else
             no_used2 = Elecs(iEl)%no_used
          end if
       end if
       ! Adjust number of orbitals for calculation with non-collinear spin
       if (nspin > 2) then
           no_used  = 2 * no_used
           no_used2 = 2 * no_used2
       end if
       call re_alloc(Elecs(iEl)%Gamma,1,no_used*no_used2,routine='transiesta')

       ! This seems stupid, however, we never use the expansion array and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       if ( Elecs(iEl)%pre_expand == 0 ) then
           no_used2 = Elecs(iEl)%no_used
           if (nspin > 2) no_used2 = 2 * no_used2
       end if
       Elecs(iEl)%GA => Elecs(iEl)%Gamma(1:no_used*no_used2)

    end do

    ! start calculation
    converged = .false.
    i_F = -1
    dq_run = ts_dq%run
    if ( dq_run ) then

      ! Allocate for interpolation, this is typically
      ! converging in 3-5 steps
      N_F = 10
      i_F = 0
      call re_alloc(Q_Ef,1,N_F,1,2)
    end if

    do while ( .not. converged )
#ifdef TRANSIESTA_TIMING
      call timer('TS_IterdQ', 1)
#endif
      call open_GF(N_Elec,Elecs,uGF,NEn)

      select case ( ts_method )
      case ( TS_FULL )
        select case ( nspin )
        case ( 1, 2 )
          if ( ts_gamma_scf ) then
            call ts_fullg(N_Elec,Elecs, &
                nq, uGF, nspin, na_u, lasto, &
                sp_dist, sparse_pattern, &
                no_u, n_nzs, &
                H, S, DM, EDM, Ef, NEGF_DE)
          else
            call ts_fullk(N_Elec,Elecs, &
                nq, uGF, &
                ucell, nspin, na_u, lasto, &
                sp_dist, sparse_pattern, &
                no_u, n_nzs, &
                H, S, DM, EDM, Ef, NEGF_DE)
          end if
        case ( 4, 8 )
          if ( ts_gamma_scf ) then
            call ts_fullg_spinor(N_Elec,Elecs, &
                nq, uGF, nspin, na_u, lasto, &
                sp_dist, sparse_pattern, &
                no_u, n_nzs, &
                H, S, DM, EDM, Ef, NEGF_DE)
          else
            call ts_fullk_spinor(N_Elec,Elecs, &
                nq, uGF, &
                ucell, nspin, na_u, lasto, &
                sp_dist, sparse_pattern, &
                no_u, n_nzs, &
                H, S, DM, EDM, Ef, NEGF_DE)
          end if
        case default
          call die('Programming error, unrecognized number &
              &of spin components.')
        end select
      case ( TS_BTD )
        select case ( nspin )
        case ( 1, 2 )
          if ( ts_gamma_scf ) then
            call ts_trig(N_Elec,Elecs, &
                nq, uGF, nspin, na_u, lasto, &
                sp_dist, sparse_pattern, &
                no_u, n_nzs, &
                H, S, DM, EDM, Ef, NEGF_DE)
          else
            call ts_trik(N_Elec,Elecs, &
                nq, uGF, &
                ucell, nspin, na_u, lasto, &
                sp_dist, sparse_pattern, &
                no_u, n_nzs, &
                H, S, DM, EDM, Ef, NEGF_DE)
          end if
      case ( 4, 8 )
        if ( ts_gamma_scf ) then
          call ts_trig_spinor(N_Elec,Elecs, &
              nq, uGF, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        else
          call ts_trik_spinor(N_Elec,Elecs, &
              nq, uGF, &
              ucell, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        end if
      case default
        call die('Programming error, unrecognized number &
            &of spin components.')
      end select
#ifdef SIESTA__MUMPS
      case ( TS_MUMPS )
        if ( ts_gamma_scf ) then
          call ts_mumpsg(N_Elec,Elecs, &
              nq, uGF, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        else
          call ts_mumpsk(N_Elec,Elecs, &
              nq, uGF, &
              ucell, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        end if
#endif
      case default
        call die('Error in code')
      end select

      ! Close files
      do iEl = 1 , N_Elec
        if ( IONode .and. Elecs(iEl)%out_of_core ) then
          call io_close(uGF(iEl))
        end if
      end do

      if ( dq_run ) then

        ! This is the 1. part of the Fermi correction
        !
        ! In this section we estimate the new Fermi level by
        ! calculating the charge Q(E_F) and then correct
        ! E_F such that dQ = Q(TS) - Qtot -> 0.
        i_F = i_F + 1
        if ( N_F < i_F ) then
          N_F = N_F + 10
          call re_alloc(Q_Ef, 1, N_F, 1, 2, copy=.true.)
        end if

        ! Save current fermi level and charge
        call ts_charge_get(N_Elec, sp_dist, sparse_pattern, &
            nspin, n_nzs, DM, S, Qtot = Q_Ef(i_F,1) )
        Q_Ef(i_F,2) = Ef

        ! Determine whether the charge *has* converged
        converged = abs(Q_Ef(i_F,1) - Qtot) < TS_DQ_FERMI_TOLERANCE

        if ( i_F < 2 ) then

          call ts_dq%calculate_dEf(Qtot, sp_dist, nspin, n_nzs, DM, S, dEf)
          ! In case the charge *is* converged
          ! Then we further reduce dEf since
          ! it causes SCF fluctuations
          if ( converged ) then
            Ef = Ef + dEf * min(abs(Q_Ef(i_F,1) - Qtot), 0.05_dp)
          else
            Ef = Ef + dEf
          end if

        else

          ! Instead of doing Q(E_F) for every change
          ! we will perform spline interpolation ones we have 2 estimated
          ! Q(E_F).
          ! This tends to drastically speed up the convergence of the dQ -> 0.
          ! In case we have accumulated 2 or more points
          dEf = Ef
          call interp_spline(i_F,Q_Ef(1:i_F,1),Q_Ef(1:i_F,2),Qtot,Ef)
          if ( ieee_is_nan(Ef) ) then
            ! There was an error in the interpolation
            ! So we'll stick with the current one and continue
            ! It probably means we are close anyway.
            Ef = dEf
            converged = .true.

          else
            dEf = Ef - dEf

            ! Truncate to the maximum allowed change in Fermi-level
            call ts_dq_truncate(Q_Ef(i_F,2), &
                TS_DQ_FERMI_MAX, Ef, truncated=truncated)
            converged = converged .and. (.not. truncated)

            if ( IONode ) then
              write(*,'(a,es13.6)') 'ts-dq-iscf: cubic spline dq = ', &
                  Q_Ef(i_F,1) - qtot
            end if
          end if


        end if

        ! Signal only doing this for the first run
        ! for subsequent dq-scf we interpolate
        ts_dq%run = .false.

      else

        ! If no Fermi-correction, we are converged
        converged = .true.

      end if

#ifdef TRANSIESTA_TIMING
      call timer('TS_IterdQ', 2)
#endif

    end do

    ! We will only do major extrapolation when we are
    ! far from the charge neutrality point.
    ! Otherwise we will iterate forever...
    if ( dq_run ) then

      ! Store the actual number of charge iterations runned
      N_F = i_F

      if ( N_F > 1 ) then

        if ( IONode ) then
          ! After converge we write out the convergence
          call io_assign(iEl)
          inquire(file='TS_FERMI', exist=converged)
          if ( converged ) then
            open(unit=iEl,file='TS_FERMI',position='append',form='formatted', &
                status='old',iostat=ioerr)
            write(iEl,'(/,a,i0)') '# TSiscf = ', TSiscf
          else
            open(unit=iEl,file='TS_FERMI',form='formatted', &
                status='new')
            write(iEl,'(a,i0)') '# TSiscf = ', TSiscf
          end if
          write(iEl,'(a,i0)')'# ', N_F ! Number of iterations
          do i_F = 1 , N_F
            write(iEl,'(2(tr1,e24.17))') Q_Ef(i_F,2)/eV, Q_Ef(i_F,1) - Qtot
          end do

          call io_close(iEl)

        end if

        ! When N_F > 1 we *know* that the initial run did
        ! not converge the charges

        ! This is the 2nd step of dEF correction.
        ! At this point we have corrected E_f for the current
        ! iteration. But generally the Hartree potential will counter
        ! the change in E_F. So to speed up convergence
        ! we do a spline interpolation of the dE_F by doing a spline
        ! interpolation of the ISCF corrections.
        ! Say TS corrects EF at iterations 50 and 80
        ! which means TS_FERMI may look like this:
        !#####
        ! # TSiscf = 50
        ! # 3
        !   -0.204782E+01    0.172948E-01
        !   -0.205297E+01    0.834017E-03
        !   -0.205323E+01   -0.250016E-05
        !
        ! # TSiscf = 80
        !# 3
        !   -0.207423E+01    0.967930E-02
        !   -0.207710E+01    0.490200E-03
        !   -0.207726E+01   -0.817259E-06
        !#####
        !

        ! Guess-stimate the actual Fermi-shift
        ! the above will typically be "too" little.
        ! So we interpolate between all previous
        ! estimations for this geometry...
        call ts_dq_Fermi_file(Q_Ef(1,2), dEf, dq=Q_Ef(1,1) - Qtot)
        ! If the fermi-file extrapolation does not work
        ! we get back a dEf == 0.
        ! In this case we rely on the simpler dQ/dQ@Ef
        if ( abs(dEf) > 1.e-11_dp ) then
          Ef = Q_Ef(1,2) + dEf
        end if

      end if

      ! We have now calculated the new Ef
      if ( IONode ) then
        write(*,'(a,es13.6,a)') 'ts-dq: dEf = ', (Ef - Q_Ef(1,2))/eV, ' eV'
      end if

      call de_alloc(Q_Ef)

    end if

    !***********************
    !       Clean up
    !***********************
    do iEl = 1 , N_Elec
      if ( .not. Elecs(iEl)%out_of_core ) then
        call Elecs(iEl)%delete()
      end if
    end do

    !***********************
    !  Clean up electrodes
    !***********************
    do iEl = 1 , N_Elec
      if ( associated(Elecs(iEl)%HA) ) then
        call de_alloc(Elecs(iEl)%HA,routine='transiesta')
        call de_alloc(Elecs(iEl)%SA,routine='transiesta')
      end if
      call de_alloc(Elecs(iEl)%Gamma,routine='transiesta')
      nullify(Elecs(iEl)%GA)
    end do

    deallocate(uGF,nq)

    ! We do the charge correction of the transiesta
    ! computation here (notice that the routine will automatically
    ! return if no charge-correction is requested)
    if ( TS_DQ_METHOD == TS_DQ_METHOD_BUFFER ) then
      call ts_dq_buffer(N_Elec, sp_dist, &
          sparse_pattern, nspin, n_nzs, DM, EDM, S, Qtot)
    end if

    call ts_charge_print(N_Elec, Elecs, Qtot, sp_dist, sparse_pattern, &
        nspin, n_nzs, DM, S, &
        method = TS_Q_INFO_SCF)

    call timer('TS',2)

  contains

    subroutine open_GF(N_Elec,Elecs,uGF,NEn)
      integer, intent(in) :: N_Elec
      type(electrode_t), intent(inout) :: Elecs(N_Elec)
      integer, intent(out) :: uGF(N_Elec)
      integer, intent(in) :: NEn

      ! Local variables
      integer :: iEl
      
      do iEl = 1 , N_Elec

        ! Initialize k-points (never seen k-point)
        Elecs(iEl)%bkpt_cur(:) = 2352345._dp

        if ( Elecs(iEl)%out_of_core ) then

          if ( IONode ) then
            call io_assign(uGF(iEl))
            open(file=Elecs(iEl)%GFfile,unit=uGF(iEl),form='unformatted')
          end if

        else

          ! prepare the electrode to create the surface self-energy
          call Elecs(iEl)%prepare_SE(IO=.false.)

        end if

        if ( Elecs(iEl)%out_of_core ) then
          call read_Green(uGF(iEl),Elecs(iEl), ts_kpoint_scf%N, NEn )
        end if
         
      end do
      
    end subroutine open_GF

  end subroutine transiesta

  subroutine ts_print_memory(ts_gamma_scf, nspin)

    use parallel, only : IONode
    use precision, only : i8b

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Max
    use mpi_siesta, only : MPI_Integer
#endif 

    use class_Sparsity
    use m_ts_options, only : IsVolt, Calc_Forces
    use m_ts_options, only : N_mu, N_Elec, Elecs
    use m_ts_contour_neq, only : N_nEq_id
    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc, ltsup_sp_sc
    use ts_electrode_m

    use m_ts_tri_init, only : c_Tri
    use m_ts_tri_common, only : GFGGF_needed_worksize
    use m_ts_tri_common, only : nnzs_tri
    use m_ts_method, only : no_Buf
    use byte_count_m, only: byte_count_t

    logical, intent(in) :: ts_gamma_scf ! transiesta Gamma
    integer, intent(in) :: nspin
    integer :: i, no_E, no_used, spin_factor
    integer(i8b) :: nel
    integer :: padding, worksize
    character(len=32) :: c_tmp
    type(byte_count_t) :: m_elec, m_glob, m_local, m_inv
#ifdef MPI
    integer :: MPIerror
#endif
    spin_factor = 1
    ! For NC/SO calculations memory requirements are four times higher
    if ( nspin > 2 ) spin_factor = 2

    ! Estimate electrode sizes
    do i = 1 , N_Elec

      no_used = Elecs(i)%no_used * spin_factor
      no_E = Elecs(i)%device_orbitals() * spin_factor

      if ( IsVolt .or. .not. Elecs(i)%Bulk ) then
        ! Hamiltonian and overlap
        if ( Elecs(i)%pre_expand > 1 ) then
          call m_elec%add_type(cmplx(0, 0, dp), no_E, no_E, 2)
        else
          call m_elec%add_type(cmplx(0, 0, dp), no_E, no_used, 2)
        end if
      end if

      if ( IsVolt ) then
        call m_elec%add_type(cmplx(0, 0, dp), no_E, no_E) ! GS/Gamma
      else
        if ( Elecs(i)%pre_expand > 0 ) then
          call m_elec%add_type(cmplx(0, 0, dp), no_E, no_E)
        else
          call m_elec%add_type(cmplx(0, 0, dp), no_E, no_used)
        end if
      end if
       
    end do

    call m_elec%get_string(c_tmp)
    if ( IONode ) then
      write(*,'(/,a,t55,a)') 'transiesta: [memory] electrodes (static): ', trim(c_tmp)
    end if

    ! Global arrays
    nel = nnzs(ts_sp_uc) * spin_factor**2
    if ( ts_gamma_scf ) then
      call m_glob%add_type(0._dp, nel, 2_i8b)
    else
      call m_glob%add_type(cmplx(0, 0, dp), nel, 2_i8b)
    end if

    ! global sparsity update
    nel = nnzs(tsup_sp_uc) * spin_factor**2
    if ( Calc_Forces ) then
      i = max(N_mu,N_nEq_id) + N_mu
    else
      i = max(N_mu,N_nEq_id)
    end if
    if ( ts_gamma_scf ) then
      call m_glob%add_type(0._dp, nel, int(i, i8b))
    else
      call m_glob%add_type(cmplx(0, 0, dp), nel, int(i, i8b))
    end if

    call m_glob%get_string(c_tmp)
    if ( IONode ) then
      write(*,'(a,t55,a)') 'transiesta: [memory] global update arrays (static): ', trim(c_tmp)
    end if

    ! Local sparsity update
    if ( IsVolt ) then
      no_E = nnzs(ltsup_sp_sc) * spin_factor**2
#ifdef MPI
      call MPI_Reduce(no_E,i,1,MPI_INTEGER, &
          MPI_MAX, 0, MPI_Comm_World, MPIerror)
      no_E = i
#endif
      if ( Calc_Forces ) then
        call m_local%add_type(0._dp, no_E, 2 * N_mu + N_nEq_id)
      else
        call m_local%add_type(0._dp, no_E, N_mu + N_nEq_id)
      end if

      ! Bias local sparsity pattern is always in double precision
      call m_local%get_string(c_tmp)
      if ( IONode ) then
        write(*,'(a,t55,a)') 'transiesta: [memory] bias sparse arrays (max): ', trim(c_tmp)
      end if
    end if

    if ( ts_method == TS_BTD ) then

      ! initialize padding and work-size query
      padding = 0
      worksize = 0

      if ( ts_A_method == TS_BTD_A_COLUMN .and. IsVolt ) then
        ! Calculate size of the tri-diagonal matrix
        call GFGGF_needed_worksize(c_Tri%n,c_Tri%r, &
            N_Elec, Elecs, padding, worksize)
      end if

      nel = nnzs_tri(c_Tri%n,c_Tri%r) * spin_factor**2
      ! Check we do not cross allowed indexing routines
      if ( nel + max(padding, worksize) > huge(1) ) then
        call die('transiesta: Memory consumption is too large!')
      end if
      call m_inv%add_type(cmplx(0, 0, dp), nel, 2_i8b)
      call m_inv%add_type(cmplx(0, 0, dp), padding)
      call m_inv%add_type(cmplx(0, 0, dp), worksize)
      call m_inv%get_string(c_tmp)
      if ( IONode ) then
        write(*,'(a,t55,a)') 'transiesta: [memory] tri-diagonal matrices: ', trim(c_tmp)
      end if

    else if ( ts_method == TS_FULL ) then

      ! Calculate size of the full matrices
      ! Here we calculate number of electrodes not needed to update the cross-terms
      no_E = sum(Elecs(:)%device_orbitals(),Elecs(:)%DM_update==0) * spin_factor
      i = (nrows_g(ts_sp_uc) - no_Buf) * spin_factor
      ! LHS
      call m_inv%add_type(cmplx(0, 0, dp), i, i)
      ! RHS
      if ( IsVolt ) then
        call m_inv%add_type(cmplx(0, 0, dp), i, max(i-no_E,sum(Elecs%device_orbitals()) * spin_factor))
      else
        call m_inv%add_type(cmplx(0, 0, dp), i, i-no_E)
      end if
      call m_inv%get_string(c_tmp)
      if ( IONode ) then
        write(*,'(a,t55,a)') 'transiesta: [memory] full matrices: ', trim(c_tmp)
      end if
#ifdef SIESTA__MUMPS
    else if ( ts_method == TS_MUMPS ) then
      if ( IONode ) then
        write(*,'(a)')'transiesta: [memory] search TS_MUMPS_<Node>.dat for: ### Minimum memory.'
      end if
#endif
    end if

    ! Sum different parts
    m_glob = m_elec + m_glob + m_local + m_inv

    call m_glob%get_string(c_tmp)
    if ( IONode ) then
      write(*,'(a,t55,a)') 'transiesta: [memory] total: ', trim(c_tmp)
    end if

  end subroutine ts_print_memory

end module m_transiesta

