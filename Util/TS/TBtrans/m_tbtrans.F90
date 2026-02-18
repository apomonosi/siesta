! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_tbtrans

  use precision, only : dp

  implicit none

  public :: tbt
  private

contains

  subroutine tbt(TSHS, kT)

#ifdef MPI
    use mpi_siesta, only : MPI_Barrier, MPI_Comm_World
#endif
    use units, only : eV
    use alloc, only : re_alloc, de_alloc

    use fdf, only: fdf_get

    use parallel, only : IONode

    use class_OrbitalDistribution
    use class_Sparsity

    use m_sparsity_handling, only : Sp_union

    use m_region

    use ts_electrode_m
    use m_tbt_options, only : N_Elec, Elecs

    use m_tbt_hs, only: tTSHS, spin_idx, Volt, prep_next_HS
    
    use m_tbt_tri_init, only : tbt_tri_init, tbt_tri_print_opti
    use m_tbt_tri_init, only : DevTri, ElTri
    use m_tbt_trik, only: tbt_trik

    use m_tbt_contour

#ifdef NCDF_4
    use m_tbt_delta, only : read_delta_Sp
    use m_tbt_dH, only : use_dH, dH
    use m_tbt_dSE, only : use_dSE, dSE
    
    use m_tbt_kpoint, only : nkpnt, kpoint, kweight
    use m_tbt_options, only : save_DATA
    use m_tbt_options, only : cdf_fname, cdf_fname_sigma, cdf_fname_proj
    use m_tbt_regions, only : r_aDev, r_aBuf, sp_dev_sc, r_oDev
    use m_tbt_regions, only : r_aEl, r_oElpD

    use m_tbt_save
    use m_tbt_proj, only : N_mol, mols, init_proj_save
    use m_tbt_sigma_save, only : init_Sigma_save
#else
    use m_tbt_regions, only : r_oDev
    use m_tbt_kpoint, only : nkpnt
#endif
    use m_tbt_kregions, only : n_k, r_k, kregion_step, kregion_k

    use m_ts_gf, only : read_Green

    use dictionary

! ********************
! * INPUT variables  *
! ********************
    type(tTSHS), intent(inout) :: TSHS
    real(dp), intent(in) :: kT

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
! ************************************************************
#ifdef NCDF_4
    ! Temporary variables
    integer :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)
    real(dp) :: k(3)
#endif

! * local variables
    integer :: iEl, NEn, no_used, no_used2, ispin, ils, i
#ifdef NCDF_4
    type(Sparsity) :: sp1, sp_total
#endif

    ! Total number of energy-points...
    NEn = N_TBT_E()

#ifdef NCDF_4

    ! Ensure that the entire dH sparse elements exists in the
    ! final pattern
    if ( use_dH ) then
       call read_delta_Sp(dH, TSHS%no_u, sp_total)
    end if

    ! Ensure that the entire dSE sparse elements exists in the
    ! final pattern
    if ( use_dSE ) then
       call read_delta_Sp(dSE, TSHS%no_u, sp1)
       call Sp_union(TSHS%dit, sp_total, sp1, sp_total)
       call delete(sp1)
    end if

    ! Make union
    call Sp_union(TSHS%dit, sp_total, TSHS%sp, sp_total)

    ! Initialize the tri-diagonal matrices!
    if ( N_mol > 0 ) then
       call tbt_tri_init( TSHS%dit, sp_total, TSHS%cell, &
            TSHS%na_u, TSHS%xa, TSHS%lasto, &
            mols(:)%orb )
    else
       call tbt_tri_init( TSHS%dit, sp_total, TSHS%cell, &
            TSHS%na_u, TSHS%xa, TSHS%lasto)
    end if

    ! Clean-up
    call delete(sp_total)

#else
    call tbt_tri_init( TSHS%dit, TSHS%sp, TSHS%cell, &
         TSHS%na_u, TSHS%xa, TSHS%lasto)
#endif

    ! Suggest to the user an optimal device region for
    ! fastest calculation
    call tbt_tri_print_opti(TSHS%na_u,TSHS%lasto,r_oDev,N_Elec)

    ! Open GF files...
    ! Read-in header of Green's functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    ! We need to initialize tbtrans
    if ( IONode ) write(*,*) ! new-line

    call timer('TBT',1)

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

       ! If we using bulk electrodes, we need not the Hamiltonian, 
       ! nor the overlap...
       ! This is because we are downfolding the self-energies and thus the
       ! Gamma is constructed in the device sub-region.
       if ( .not. Elecs(iEl)%Bulk ) then
          call re_alloc(Elecs(iEl)%HA,1,no_used,1,no_used,1,nq(iEl),routine='tbtrans')
          call re_alloc(Elecs(iEl)%SA,1,no_used,1,no_used,1,nq(iEl),routine='tbtrans')
       end if

       ! We need to retain the maximum size of Gamma
       ! I.e. take into account that the down-projected region
       ! could be larger than the read in self-energy
       no_used = Elecs(iEl)%device_orbitals()
       no_used2 = no_used
       if ( Elecs(iEl)%pre_expand == 0 ) then
          no_used2 = Elecs(iEl)%no_used
       end if
       ! We need to assure that the entire down-folded Gamma
       ! can be saved
       if ( no_used * no_used2 < Elecs(iEl)%o_inD%n ** 2 ) then
          no_used2 = Elecs(iEl)%o_inD%n ** 2 / no_used + 1
       end if

       call re_alloc(Elecs(iEl)%Gamma,1,no_used*no_used2,routine='tbtrans')

       ! This seems stupid, however, we never use the expansion array and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       no_used2 = no_used
       if ( Elecs(iEl)%pre_expand == 0 ) no_used2 = Elecs(iEl)%no_used
       Elecs(iEl)%GA => Elecs(iEl)%Gamma(1:no_used*no_used2)

    end do

    call open_GF(N_Elec,Elecs,uGF,nkpnt,NEn,spin_idx)
    if ( n_k > 0 ) then
       ! We initialize the region that belongs to
       ! each electrode
       do iEl = 1 , N_Elec
          uGF(iEl) = -1
          do ils = 0 , n_k
             if ( in_rgn(r_k(ils)%atm,Elecs(iEl)%idx_a) ) then
                uGF(iEl) = ils
                exit
             end if
          end do
          if ( uGF(iEl) == -1 ) then
             call die('Error in setting up different k-regions, &
                  &all electrodes are not fully contained.')
          end if
       end do
    end if

    ! Print memory usage (after open_GF)
    call print_memory()
    
    do ils = 1 , TSHS%nspin

       ! Determine the actual spin-index
       if ( spin_idx == 0 ) then
          ispin = ils
       else
          ispin = spin_idx
       end if

#ifdef NCDF_4
       if ( use_dH ) then
          ! Force to re-read the following dH matrices
          dH%lvl = -1
          dH%ispin = ispin
       end if
       if ( use_dSE ) then
          ! Force to re-read the following dSE matrices
          dSE%lvl = -1
          dSE%ispin = ispin
       end if
#endif

       ! The initial spin has already been 
       ! setup for the first spin, hence we
       ! only need to re-read them for
       ! following spin calculations
       if ( ils > 1 ) then
          call prep_next_HS(ispin, Volt)
          do iEl = 1 , N_Elec
            ! Re-read in the electrode
            ! Hamiltonian and build the H/S
            ! for this spin.
            if ( Elecs(iEl)%out_of_core ) cycle
            call Elecs(iEl)%prepare_SE(IO=.false., spin=ispin)
          end do
       end if

#ifdef NCDF_4
       if ( n_k == 0 ) then
          nkpt =  nkpnt
          kpt  => kpoint(:,:)
          wkpt => kweight(:)
       else
          nullify(kpt,wkpt)
          nkpt = 1
          do i = 0 , n_k
             nkpt = nkpt * size(r_k(i)%wkpt)
          end do
          allocate(kpt(3,nkpt),wkpt(nkpt))
          r_k(:)%ik   = 1
          r_k(n_k)%ik = 0
          do i = 1 , nkpt
             call kregion_step( )
             call kregion_k(-1, k, w = wkpt(i) )
             call kpoint_convert(TSHS%cell,k,kpt(:,i),1)
          end do
       end if

       ! If the user has requested only to calculate
       ! the projections, we do not initialize the cdf_fname
       if ( ('proj-only'.nin.save_DATA).and.('Sigma-only'.nin.save_DATA) ) then

          ! Initialize data files
         call name_save( ispin, TSHS%nspin, cdf_fname, end = 'nc')
         call init_cdf_save(cdf_fname,TSHS,r_oDev,DevTri,ispin, &
             N_Elec, Elecs, r_aEl, r_oElpd, ElTri, &
             nkpt, kpt, wkpt, NEn, tbt_Eta, r_aDev, r_aBuf, sp_dev_sc, save_DATA )

       end if
       
       call name_save( ispin, TSHS%nspin, cdf_fname_sigma, end = 'SE.nc')
       call init_Sigma_save(cdf_fname_sigma,TSHS,r_oDev,DevTri,ispin, &
           N_Elec, Elecs, r_aEl, r_oElpd, ElTri, &
           nkpt, kpt, wkpt, NEn, tbt_Eta, r_aDev, r_aBuf )

       if ( ('Sigma-only'.nin.save_DATA) ) then
       
          call name_save( ispin, TSHS%nspin, cdf_fname_proj, end = 'Proj.nc' )
          call init_Proj_save( cdf_fname_proj,TSHS,r_oDev,DevTri,ispin, &
              N_Elec, Elecs, r_aEl, r_oElpd, ElTri, &
              nkpt, kpt, wkpt, NEn, tbt_Eta, r_aDev, r_aBuf, sp_dev_sc, save_DATA )
       end if

       if ( n_k /= 0 ) then
          deallocate(kpt,wkpt)
          nullify(kpt,wkpt)
       end if
#endif

       call tbt_trik(ispin, N_Elec, Elecs, TSHS, nq, uGF)

       ! the spin-index is zero for all,
       ! and one of the allowed spin indices if
       ! a specific one is requested...
       if ( ispin == spin_idx ) exit

    end do

    ! Close files
    do iEl = 1 , N_Elec
       if ( IONode .and. Elecs(iEl)%out_of_core ) then
          call io_close(uGF(iEl))
       end if
    end do
       
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
          call de_alloc(Elecs(iEl)%HA,routine='tbtrans')
          call de_alloc(Elecs(iEl)%SA,routine='tbtrans')
       end if
       call de_alloc(Elecs(iEl)%Gamma,routine='tbtrans')
       ! Clean pointer
       nullify(Elecs(iEl)%GA)
    end do

    deallocate(uGF,nq)

    call timer('TBT',2)

  contains

    subroutine print_memory()

      use class_dSpData2D, only: nnzs, spar_dim, size
      use precision, only: i8b
      use m_verbosity, only : verbosity
      use m_tbt_regions, only : sp_uc, sp_dev_sc
      use byte_count_m, only: byte_count_t

      integer :: nsize, nspin

      character(len=32) :: c_tmp
      type(byte_count_t) :: mem, spaux
      integer :: iEl

      if ( .not. IONode ) return
      if ( verbosity <= 4 ) return

      call mem%reset()

      ! Calculate size of electrodes
      nsize = 0
      do iEl = 1 , N_Elec
        ! These are doubles (not complex)
        ! hence, we need not count the overlap matrices as they are already
        ! doubled in the complex conversion
        if ( spar_dim(Elecs(iEl)%H00) == 1 ) then
          nspin = size(Elecs(iEl)%H00, 2)
        else
          nspin = size(Elecs(iEl)%H00, 1)
        end if
        nsize = nnzs(Elecs(iEl)%H00)
        call mem%add_type(0._dp, nsize, nspin + 1) ! H, S
        call mem%add_type(0, nsize) ! sparse pattern
        nsize = nnzs(Elecs(iEl)%H01)
        call mem%add_type(0._dp, nsize, nspin + 1) ! H, S
        call mem%add_type(0, nsize) ! sparse pattern

        if ( .not. Elecs(iEl)%Bulk ) then
          call mem%add_type(cmplx(0, 0, dp), size(Elecs(iEl)%HA), 2) ! H, S
        end if
        call mem%add_type(cmplx(0, 0, dp), size(Elecs(iEl)%Gamma))
      end do

      ! Electrode memory usage in KB
      call mem%get_string(c_tmp)
      write(*,'(2a)') 'tbt: [memory] electrode: ', trim(c_tmp)

      call spaux%reset()
      ! first the size of the real matrices (H and S)
      ! Since we immediately reduce to only one spin-component we never
      ! have spin-degeneracy as a memory requirement.
      nsize = nnzs(TSHS%sp)
      call spaux%add_type(0._dp, nsize, 2) ! H, S
      call spaux%add_type(0, nsize) ! sparse pattern
#ifdef NCDF_4
      if ( initialized(sp_dev_sc) ) then
        ! this is the sparse orbital currents/COOP/COHP
        nsize = nnzs(sp_dev_sc)
        call spaux%add_type(0._dp, nsize)
        call spaux%add_type(0, nsize) ! sparse pattern
      end if
#endif

      ! Now the complex sparse matrices H, S, these could essentially be removed
      nsize = nnzs(sp_uc)
      call spaux%add_type(cmplx(0, 0, dp), nsize) ! H, S
      call spaux%add_type(0, nsize) ! sparse pattern
      call spaux%get_string(c_tmp)
      write(*,'(2a)') 'tbt: [memory] H, S and auxiliary matrices: ', trim(c_tmp)

      ! Total memory
      mem = mem + spaux
      call mem%get_string(c_tmp)
      write(*,'(2a/)') 'tbt: [memory] sum of electrode and sparse: ', trim(c_tmp)

    end subroutine print_memory

    subroutine open_GF(N_Elec,Elecs,uGF,nkpnt,NEn,spin_idx)
      integer, intent(in) :: N_Elec
      type(electrode_t), intent(inout) :: Elecs(N_Elec)
      integer, intent(out) :: uGF(N_Elec)
      integer, intent(in) :: nkpnt, NEn, spin_idx

      ! Local variables
      integer :: iEl
      
      do iEl = 1 , N_Elec

         ! Initialize k-points
         Elecs(iEl)%bkpt_cur(:) = 2352345._dp

         if ( Elecs(iEl)%out_of_core ) then
            
            if ( IONode ) then
               call io_assign(uGF(iEl))
               open(file=Elecs(iEl)%GFfile,unit=uGF(iEl),form='unformatted')
            end if
            
            call read_Green(uGF(iEl), Elecs(iEl), nkpnt, NEn)

            ! There can only be two spin-components
            if ( IONode .and. spin_idx > 1 .and. .not. Elecs(iEl)%is_gamma ) then
              ! In case the user has requested to only use
              ! one of the spin-channels we step forward to that one
              ! Skip all H and S arrays
              do i = 1 , nkpnt
                read(uGF(iEl)) ! H
                read(uGF(iEl)) ! S
              end do
              ! Skip all header lines AND GS lines
              do i = 1 , nkpnt * NEn
                read(uGF(iEl)) ! Header line
                read(uGF(iEl)) ! GS
              end do
            end if

         else
            
           ! prepare the electrode to create the surface self-energy
           call Elecs(iEl)%prepare_SE(IO=.false., spin=max(1, spin_idx))
            
         end if

      end do
      
    end subroutine open_GF

  end subroutine tbt

end module m_tbtrans

