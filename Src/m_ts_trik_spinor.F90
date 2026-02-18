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
! * It has been heavily inspired by the original authors of the 
!   Transiesta code (hence the references here are still remaining) *

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.

module m_ts_trik_spinor

  use precision, only : dp
  use m_region

  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM, update_zDM
  use m_ts_dm_update, only : add_k_DM
  
  use m_ts_weight_spinor, only : weight_DM_spinor
  use m_ts_weight, only : TS_W_K_METHOD
  use m_ts_weight, only : TS_W_K_CORRELATED
  use m_ts_weight, only : TS_W_K_UNCORRELATED

  use m_ts_tri_init, only : c_Tri

  use m_ts_method, only : orb_offset, no_Buf, r_pvt
  use m_ts_method, only : ts_A_method, TS_BTD_A_COLUMN
  
  implicit none

  public :: ts_trik_spinor

  private
  
contains

  subroutine ts_trik_spinor(N_Elec,Elecs, &
       nq, uGF, ucell, nspin, na_u, lasto, &
       sp_dist, sparse_pattern, &
       no_u, n_nzs, &
       Hs, Ss, DM, EDM, Ef, &
       DE_NEGF)

    use units, only : Pi
    use parallel, only : Node, Nodes

#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D
    use class_dSpData2D
    use class_zSpData2D
    use class_dSpData3D
    use class_zSpData3D
    use class_zTriMat

    use m_spin, only : spin

    use ts_electrode_m
    ! Self-energy read
    use m_ts_gf
    ! Self-energy expansion
    use m_ts_elec_se

    use ts_kpoint_scf_m, only : ts_kpoint_scf

    use m_ts_options, only : Calc_Forces
    use m_ts_options, only : N_mu, mus

    use m_ts_options, only : IsVolt

    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc, ltsup_sc_pnt
    use m_ts_sparse, only : sc_off

    use m_ts_cctype
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq, c2energy
    use m_ts_contour_neq, only : nEq_E, has_cE_nEq
    use m_ts_contour_neq, only : N_nEq_ID, c2weight_neq
    use m_ts_contour_eq, only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use m_iterator
    use m_mat_invert

    use m_trimat_invert

    ! Gf calculation
    use m_ts_trimat_invert

    use m_ts_tri_common, only : GFGGF_needed_worksize, nnzs_tri

    ! Gf.Gamma.Gf
    use m_ts_tri_scat_spinor

    use ts_dq_m, only: ts_dq

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in) :: N_Elec
    type(electrode_t), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: nq(N_Elec), uGF(N_Elec)
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: nspin, na_u, lasto(0:na_u)
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    integer, intent(in)  :: no_u
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef
    real(dp), intent(inout) :: DE_NEGF

! ******************* Computational arrays *******************
    integer :: nzwork, n_s
    complex(dp), pointer :: zwork(:)
    type(zTriMat) :: zwork_tri, GF_tri
    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D) :: spS
    type(zSpData3D) :: spH
    ! local sparsity pattern in local SC pattern
    type(dSpData3D) :: spDM, spDMneq
    type(dSpData3D) :: spEDM ! only used if calc_forces
    ! The different sparse matrices that will surmount to the integral
    ! These two lines are in global update sparsity pattern (UC)
    type(zSpData3D) ::  spuDM
    type(zSpData3D) :: spuEDM ! only used if calc_forces
    ! To figure out which parts of the tri-diagonal blocks we need
    ! to calculate
    logical, pointer :: calc_parts(:) => null()
! ************************************************************

! ****************** Electrode variables *********************
    integer :: padding, GFGGF_size ! with IsVolt we need padding and work-array
    complex(dp), pointer :: GFGGF_work(:) => null()
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    integer :: NE
    integer :: index_dq !< Index for the current charge calculation @ E == mu
    real(dp) :: kw, kpt(3), bkpt(3), dq_mu
    complex(dp) :: W, ZW
    type(tRgn)  :: pvt
! ************************************************************

! ******************** Loop variables ************************
    type(itt1) :: Kp
    integer, pointer :: ispin, ikpt
    integer :: iEl, iID
    integer :: iE, imu, io, idx
    integer :: no
! ************************************************************
#ifdef TRANSIESTA_TIMING
    call timer('TS_pre_Econtours', 1)
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    ! Create the back-pivoting region
    call rgn_init(pvt,nrows_g(sparse_pattern),val=0)
    do io = 1 , r_pvt%n
       pvt%r(r_pvt%r(io)) = io
    end do

    ! Number of supercells
    n_s = size(sc_off,dim=2)

    ! We do need the full GF AND a single work array to handle the
    ! left-hand side of the inversion...
    ! We will provide all work arrays as single dimension arrays.
    ! This will make interfaces more stringent and allow for
    ! re-use in several other places.
    ! However, this comes at the cost of programmer book-keeping.
    ! It should be expected that the work arrays return GARBAGE
    ! on ALL routines, i.e. they are not used for anything other
    ! than, yes, work.

    ! The zwork is needed to construct the LHS for solving: G^{-1} G = I
    ! Hence, we will minimum require this...
    if ( ts_A_method == TS_BTD_A_COLUMN .and. IsVolt ) then
       call GFGGF_needed_worksize(c_Tri%n,c_Tri%r, &
            N_Elec, Elecs, padding, GFGGF_size)
      padding = 2 * padding
      GFGGF_size = 2 * GFGGF_size
    else
       padding = 0
       GFGGF_size = 0
    end if
    ! Now figure out the required worksize for SE expansion
    call UC_minimum_worksize(IsVolt, N_Elec, Elecs, idx)
    io = int( nnzs_tri(c_Tri%n, c_Tri%r) )
    padding = max(padding, 2*(idx - io))
    ! Note the doubling of the block sizes
    call newzTriMat(zwork_tri,c_Tri%n,2*c_Tri%r,'GFinv', &
         padding=padding)
    nzwork = elements(zwork_tri,all=.true.)

    ! Initialize the tri-diagonal inversion routine
    call init_TriMat_inversion(zwork_tri)

    call newzTriMat(GF_tri,c_Tri%n,2*c_Tri%r,'GF')

    ! initialize the matrix inversion tool
    call init_mat_inversion(maxval(c_Tri%r))

    ! Allocate the logical array to handle calculated
    ! entries in the block-matrix
    call re_alloc(calc_parts,1,c_Tri%n)
    ! initialize to calculate all blocks
    calc_parts(:) = .true.

    ! we use the GF as a placement for the self-energies
    no = 0
    zwork => val(GF_tri,all=.true.)
    do iEl = 1 , N_Elec

       ! This seems stupid, however, we never use the Sigma and
       ! GF at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called
       ! first the Sigma is assigned and then 
       ! it is required that prepare_GF_inv is called
       ! immediately (which it is)
       ! Hence the GF must NOT be used in between these two calls!
       io = (2*Elecs(iEl)%device_orbitals()) ** 2
       Elecs(iEl)%Sigma => zwork(no+1:no+io)
       no = no + io

       ! if we need the cross-terms we can not skip the blocks
       ! that are fully inside the electrode
       if ( Elecs(iEl)%DM_update /= 0 ) cycle

       io  = Elecs(iEl)%idx_o
       io  = io - orb_offset(io)
       idx = io + Elecs(iEl)%device_orbitals() - 1

    end do

    ! Save the work-space
    ! Now the programmer should "keep a straight tongue"
    ! The zwork points to the array in the zwork_tri
    ! tri-diagonal array. This means that there are two
    ! arrays that point to the same.
    ! Generally the zwork need only to retain the value in
    ! one call!
    zwork => val(zwork_tri,all=.true.)

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(no_u,MPI_Comm_Self,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1           ,fdist,name='TS-fake dist')
#endif
    
    ! The Hamiltonian and overlap matrices (in Gamma calculations
    ! we will not have any phases, hence, it makes no sense to
    ! have the arrays in complex)
    ! TODO move into a Data2D (could reduce overhead of COMM)
    call newzSpData3D(ts_sp_uc,2,2,fdist,spH,name='TS spH',sparsity_dim=3)
    call newzSpData1D(ts_sp_uc,    fdist,spS,name='TS spS')

    ! If we have a bias calculation we need additional arrays.
    ! If not bias we don't need the update arrays (we already have
    ! all information in tsup_sp_uc (spDMu))

    ! Allocate space for global sparsity arrays
    no = max(N_mu,N_nEq_id)
    call newzSpData3D(tsup_sp_uc, 4, no, fdist, spuDM, name='TS spuDM', &
                      sparsity_dim=2)
    if ( Calc_Forces ) then
       call newzSpData3D(tsup_sp_uc,4,N_mu,fdist, spuEDM, &
                         name='TS spuEDM', sparsity_dim=2)
    end if
    
    if ( IsVolt ) then
       ! Allocate space for update arrays, local sparsity arrays
       call newdSpData3D(ltsup_sp_sc,spin%DM ,N_mu,    sp_dist,spDM   , &
                         name='TS spDM', sparsity_dim=2)
       call newdSpData3D(ltsup_sp_sc,spin%DM,N_nEq_id,sp_dist,spDMneq, &
                         name='TS spDM-neq', sparsity_dim=2)
       if ( Calc_Forces ) then
          call newdSpData3D(ltsup_sp_sc,spin%EDM,N_mu, sp_dist,spEDM  , &
                            name='TS spEDM', sparsity_dim=2)
       end if
    end if

    if ( ts_A_method == TS_BTD_A_COLUMN .and. IsVolt ) then
       ! we need only allocate one work-array for
       ! Gf.G.Gf^\dagger
       call re_alloc(GFGGF_work,1,GFGGF_size,routine='transiesta')
    end if

#ifdef TRANSIESTA_WEIGHT_DEBUG
    do iel = 1 , n_mu
       print '(a20,tr1,i3)','k  '//trim(mus(iEl)%name),iel
    end do
#endif

    ! Initialize the charge correction scheme (will return if not used)
    call ts_dq%initialize_dq()

    ! Total number of energy points
    NE = N_Eq_E() + N_nEq_E()

    ! start the itterators
    call itt_init  (Kp,end=ts_kpoint_scf%N)
    ! point to the index iterators
    call itt_attach(Kp,cur=ikpt)

#ifdef TRANSIESTA_TIMING
    call timer('TS_pre_Econtours', 2)
    call timer('TS_Econtours', 1)
#endif

    call init_DM(sp_dist,sparse_pattern, &
        n_nzs, spin, DM, EDM, &
        tsup_sp_uc, Calc_Forces)
    do iEl = 1, N_Elec
      call reread_Gamma_Green(Elecs(iEl), uGF(iEl), NE, 1)
    end do
    
    do while ( .not. itt_step(Kp) )
        
       ! Include spin factor and 1/(2\pi)
       kpt(:) = ts_kpoint_scf%k(:,ikpt)
       ! create the k-point in reciprocal space
       call kpoint_convert(ucell,kpt,bkpt,1)

       ! Include spin factor and 1/\pi
       ! Since we are calculating G - G^\dagger in the equilibrium contour
       ! we need *half* the weight
       kw = ts_kpoint_scf%w(ikpt) / (Pi * 2._dp)
       
#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       call create_HS(sp_dist,sparse_pattern, &
            spin%H, Ef, &
            N_Elec, Elecs, no_u, n_s, & ! electrodes, SIESTA size
            n_nzs, Hs, Ss, sc_off, &
            spH, spS, kpt, &
            nzwork, zwork)

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif


#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',1)
#endif

       ! ***************
       ! * EQUILIBRIUM *
       ! ***************
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = Nodes - Node
       cE = Eq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(1, ikpt, bkpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .false. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, zwork_tri, r_pvt, pvt, &
               N_Elec, Elecs, &
               spH=spH , spS=spS )

          ! *******************
          ! * calc GF         *
          ! *******************
          if ( .not. cE%fake ) then
             call invert_TriMat(zwork_tri,GF_tri,calc_parts)
          end if
          
          ! ** At this point we have calculated the Green function

          ! ****************
          ! * save GF      *
          ! ****************
          do imu = 1 , N_mu
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution...
             call ID2idx(cE,mus(imu)%ID,idx)
             if ( idx < 1 ) cycle
             
             call c2weight_eq(cE,idx, kw, W ,ZW)

             ! Figure out if this point is a charge-correction energy
             index_dq = ts_dq%get_index(imu, iE)
             if ( index_dq > 0 ) then
               ! Accummulate charge at the electrodes chemical potentials
               ! Note that this dq_mu does NOT have the prefactor 1/Pi
               call add_DM( spuDM, W, spuEDM, ZW, &
                   GF_tri, r_pvt, pvt, &
                   N_Elec, Elecs, &
                   DMidx=mus(imu)%ID, &
                   spS=spS, q=dq_mu)
               ts_dq%mus(imu)%dq(index_dq) = ts_dq%mus(imu)%dq(index_dq) + dq_mu * kw
             else
               call add_DM( spuDM, W, spuEDM, ZW, &
                   GF_tri, r_pvt, pvt, &
                   N_Elec, Elecs, &
                   DMidx=mus(imu)%ID)
             end if
           end do

          ! step energy-point
          iE = iE + Nodes
          cE = Eq_E(iE,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',2)
#endif

#ifdef MPI
       ! We need to reduce all the arrays
       call MPI_Barrier(MPI_Comm_World,io)
       call timer('TS_comm',1)
       call AllReduce_SpData(spuDM,nzwork,zwork,4,N_mu)
       if ( Calc_Forces ) then
          call AllReduce_SpData(spuEDM,nzwork,zwork,4,N_mu)
       end if
       call timer('TS_comm',2)
#endif

       if ( .not. IsVolt ) then
          call update_zDM(sp_dist,sparse_pattern, n_nzs,&
               DM, spuDM, Ef, &
               EDM, spuEDM, kpt, n_s, sc_off)
          
          ! The remaining code segment only deals with 
          ! bias integration... So we skip instantly
          do iEl = 1, N_Elec
            call reread_Gamma_Green(Elecs(iEl), uGF(iEl), NE, 1)
          end do

          cycle

       end if

       ! *****************
       ! * only things with non-Equilibrium contour...
       ! *****************

       ! initialize to zero
       ! local sparsity update patterns
       if ( TS_W_K_METHOD == TS_W_K_UNCORRELATED ) then
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       else if ( itt_first(Kp) ) then
          ! we only need to initialize once per spin
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       end if

       ! transfer data to local sparsity arrays
       call add_k_DM(spDM, spuDM, spin%DM, N_mu, &
            spEDM, spuEDM, spin%EDM, N_mu, &
            n_s, sc_off, kpt, non_Eq = .false. )

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',1)
#endif

       ! *******************
       ! * NON-EQUILIBRIUM *
       ! *******************

       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = Nodes - Node
       cE = nEq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(1, ikpt, bkpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .true. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, zwork_tri, r_pvt, pvt, &
               N_Elec, Elecs, &
               spH =spH , spS =spS)
          
          ! *******************
          ! * prep GF         *
          ! *******************
          if ( .not. cE%fake ) then
             call invert_BiasTriMat_prep(zwork_tri,GF_tri, &
                  all_nn = .true. )
          end if

          ! ** At this point we have calculated the needed
          ! ** information to create the Green function column
          ! ** for all the electrodes

          ! ****************
          ! * save GF      *
          ! ****************
          do iEl = 1 , N_Elec
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution

             ! ******************
             ! * calc GF-column *
             ! ******************
             if ( ts_A_method == TS_BTD_A_COLUMN ) then
              call invert_BiasTriMat_rgn(GF_tri,zwork_tri, &
                   r_pvt, pvt, Elecs(iEl)%o_inD)

              call GF_Gamma_GF(zwork_tri, Elecs(iEl), Elecs(iEl)%o_inD%n, &
                   calc_parts, GFGGF_size, GFGGF_work)
#ifdef TRANSIESTA_WEIGHT_DEBUG
              print '(a7,tr1,i3,2(tr1,f10.5),tr5,2(tr1,f10.5))', &
                   trim(Elecs(iEl)%name),iE,zwork(index(zwork_tri,28,28)),cE%e
#endif

             else
              call dir_GF_Gamma_GF(Gf_tri, zwork_tri, r_pvt, pvt, &
                   Elecs(iEl), calc_parts)
             end if
             
             do iID = 1 , N_nEq_ID
                
                if ( .not. has_cE_nEq(cE,iEl,iID) ) cycle
                
                call c2weight_nEq(cE,iID,kw,W,imu,ZW)

#ifdef TRANSIESTA_WEIGHT_DEBUG
                print '(a20,2(tr1,i3),2(tr1,e12.5))', &
                     trim(Elecs(iEl)%name),iID,imu,W
#endif 

                call add_DM( spuDM, W, spuEDM, ZW, &
                     zwork_tri, r_pvt, pvt, &
                     N_Elec, Elecs, &
                     DMidx=iID, EDMidx=imu, is_eq = .false.)
             end do
          end do

          ! step energy-point
          iE = iE + Nodes
          cE = nEq_E(iE,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',2)
#endif

#ifdef MPI
       ! We need to reduce all the arrays
       call MPI_Barrier(MPI_Comm_World,io)
       call timer('TS_comm',1)
       call AllReduce_SpData(spuDM, nzwork, zwork, 4, N_nEq_id)
       if ( Calc_Forces ) then
          call AllReduce_SpData(spuEDM, nzwork, zwork, 4, N_mu)
       end if
       call timer('TS_comm',2)
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',1)
#endif

       ! 1. move from global UC to local SC
       ! 2. calculate the correct contribution by applying the weight
       ! 3. add the density to the real arrays
       call add_k_DM(spDMneq, spuDM, spin%DM, N_nEq_id, &
            spEDM, spuEDM, spin%EDM, N_mu, &
            n_s, sc_off, kpt, non_Eq = .true. )

       if ( TS_W_K_METHOD == TS_W_K_UNCORRELATED ) then
          call weight_DM_spinor( N_Elec, Elecs, N_mu, mus, na_u, lasto, &
               sp_dist, sparse_pattern, Ss, &
               spDM, spDMneq, spEDM, n_s, sc_off, DE_NEGF)
          
#ifdef TRANSIESTA_WEIGHT_DEBUG
          call die('')
#endif
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM, spDM, Ef=Ef, &
               EDM=EDM, spEDM=spEDM, ipnt=ltsup_sc_pnt)
       else if ( itt_last(Kp) ) then ! TS_W_K_METHOD == TS_W_K_CORRELATED
          call weight_DM_spinor( N_Elec, Elecs, N_mu, mus, na_u, lasto, &
               sp_dist, sparse_pattern, Ss, &
               spDM, spDMneq, spEDM, n_s, sc_off, DE_NEGF)
          
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM, spDM, Ef=Ef, &
               EDM=EDM, spEDM=spEDM, ipnt=ltsup_sc_pnt)          
       end if

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',2)
#endif

       do iEl = 1, N_Elec
         call reread_Gamma_Green(Elecs(iEl), uGF(iEl), NE, ispin)
       end do

    end do ! spin, k-point

#ifdef TRANSIESTA_TIMING
    call timer('TS_Econtours', 2)
    call timer('TS_post_Econtours', 1)
#endif

    call itt_destroy(Kp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA SCF'
#endif

!***********************
! CLEAN UP
!***********************

    call de_alloc(calc_parts)

    call delete(zwork_tri)
    call delete(GF_tri)

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    call delete(spuEDM)

    call delete(spDM)
    call delete(spDMneq)
    call delete(spEDM)

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    call clear_TriMat_inversion()
    if ( ts_A_method == TS_BTD_A_COLUMN .and. IsVolt ) then
       call de_alloc(GFGGF_work, routine='transiesta')
    end if
    call clear_mat_inversion()

    call rgn_delete(pvt)

    ! Nullify external pointers
    do iEl = 1, N_Elec
      nullify(Elecs(iEl)%Sigma)
    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif

#ifdef TRANSIESTA_TIMING
    call timer('TS_post_Econtours', 2)
#endif

  end subroutine ts_trik_spinor
                        

  subroutine add_DM(DM, DMfact, EDM, EDMfact, &
       GF_tri, r, pvt, &
       N_Elec, Elecs, &
       DMidx,EDMidx, &
       spS, q, &
       is_eq)

    use intrinsic_missing, only: SFIND

    use class_zSpData1D
    use class_zSpData3D
    use class_Sparsity
    use class_zTriMat

    use ts_electrode_m

    ! The DM and EDM equivalent matrices
    type(zSpData3D), intent(inout) :: DM
    complex(dp), intent(in) :: DMfact
    type(zSpData3D), intent(inout) :: EDM
    complex(dp), intent(in) :: EDMfact

    ! The Green function
    type(zTriMat), intent(inout) :: GF_tri
    type(tRgn), intent(in) :: r, pvt
    integer, intent(in) :: N_Elec
    type(electrode_t), intent(in) :: Elecs(N_Elec)
    ! the index of the partition
    integer, intent(in) :: DMidx
    integer, intent(in), optional :: EDMidx
    !< Overlap matrix setup for a k-point is needed for calculating q
    type(zSpData1D), intent(in), optional :: spS
    !< Charge calculated at this energy-point
    !!
    !! This does not contain the additional factor 1/Pi
    real(dp), intent(inout), optional :: q
    logical, intent(in), optional :: is_eq

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: s_ncol(:), s_ptr(:), s_col(:), sp_col(:)
    complex(dp), pointer :: D(:,:,:), E(:,:,:)
    complex(dp), pointer :: Gf(:)
    complex(dp) :: GG
    complex(dp), pointer :: Sk(:)
    integer :: sp, sind
    integer :: io, ind, iu, ju, ii, jj, idx(2,2), idxT(2,2), i1, i2
    logical :: calc_q
    logical :: hasEDM, lis_eq
#ifdef TRANSIESTA_TIMING
    call timer('TS_add_DM', 1)
#endif

    lis_eq = .true.
    if ( present(is_eq) ) lis_eq = is_eq

    calc_q = present(q) .and. present(spS)
    
    s => spar(DM)
    call attach(s, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    D => val(DM)
    hasEDM = initialized(EDM)
    if ( hasEDM ) E => val(EDM)
    
    i1 = DMidx
    i2 = i1
    if ( present(EDMidx) ) i2 = EDMidx

    ! Calculate new DM contribution from current GF.
    !
    ! (1) In equilibrium we use :
    !
    !         DM = \int GF Gamma GF^Dag dz = I/2 ( \int GF dz - (\int GF dz)^Dag )
    !
    !     and DM contribution of the current energy point is
    !
    !        Gf * DMfact - conjg(Gf * DMfact)
    !
    ! (2) In non-equilbrium we have to use
    !
    !         DM = \int GF Gamma GF^Dag dz
    !
    !     and GF already contains the triple-matrix-product: GF Gamma GF^Dag.
    !     The DM contribution is simply:
    !
    !         Gf * DMfact
    !
    ! (3) Structure of D(:,:)
    !
    !     We still need to multiply with complex the density-matrix with
    !     k-phase later, when transforming back from recirocal space.
    !     Since, the k-phase is complex, it is beneficial to keep the density
    !     matrix complex, and separate Real and Imaginary part, layer when
    !     copying to the real array (Dscf) in update_zDM at the end of the 
    !     energy contour.
    !    
    !     For now we store four complex numbers in the following order.
    !     TODO We could make D 4-dimensional?
    !    
    !            | D(1,:)    D(2,:) |     | DM^{uu}   DM^{ud} |
    !       DM = |                  |  =  |                   |
    !            | D(3,:)    D(4,:) |     | DM^{du}   DM^{dd} |

    Gf => val(Gf_tri)

    if ( lis_eq ) then
      ! *********** !
      ! Equilibrium !
      ! *********** !
      ! 
      !  This is the equilibrium case so we use:
      !
      !      GF Gamma GF^Dag = I/2 (GF - GF^Dag)

      if ( calc_q ) then
        q = 0._dp
        s => spar(spS)
        Sk => val(spS)
        call attach(s, n_col=s_ncol, list_ptr=s_ptr, list_col=s_col)
      end if

      if ( calc_q .and. hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ju,ii,jj,sp,sp_col,ind,idx,idxT,sind,GG) &
!$OMP&  reduction(-:q)
       do iu = 1 , r%n
          io = r%r(iu)
          if ( l_ncol(io) /= 0 ) then

          sp = s_ptr(io)
          sp_col => s_col(sp+1:sp+s_ncol(io))

          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

            ju = pvt%r(l_col(ind))
            do ii = 1, 2
              do jj = 1, 2
                idx (ii,jj) = index(Gf_tri,2*(iu-1)+ii,2*(ju-1)+jj)
                idxT(jj,ii) = index(Gf_tri,2*(ju-1)+jj,2*(iu-1)+ii)
              end do
            end do

            ! Search for overlap index
            sind = sp + SFIND(sp_col, l_col(ind))
            if ( sp < sind ) then
              ! For charge calculation it is Tr[(G-G^\dagger).S] i.e. the matrix
              ! multiplication
              q = q - aimag((Gf(idx(1,1)) - conjg(Gf(idxT(1,1)))) * Sk(sind))
              q = q - aimag((Gf(idx(2,2)) - conjg(Gf(idxT(2,2)))) * Sk(sind))
            end if
            D(1,ind,i1) = D(1,ind,i1) + Gf(idx(1,1)) * DMfact - conjg(Gf(idxT(1,1)) * DMfact)
            D(2,ind,i1) = D(2,ind,i1) + Gf(idx(1,2)) * DMfact - conjg(Gf(idxT(2,1)) * DMfact)
            D(3,ind,i1) = D(3,ind,i1) + Gf(idx(2,1)) * DMfact - conjg(Gf(idxT(1,2)) * DMfact)
            D(4,ind,i1) = D(4,ind,i1) + Gf(idx(2,2)) * DMfact - conjg(Gf(idxT(2,2)) * DMfact)

            E(1,ind,i2) = E(1,ind,i2) + Gf(idx(1,1)) * EDMfact - conjg(Gf(idxT(1,1)) * EDMfact)
            E(2,ind,i2) = E(2,ind,i2) + Gf(idx(1,2)) * EDMfact - conjg(Gf(idxT(2,1)) * EDMfact)
            E(3,ind,i2) = E(3,ind,i2) + Gf(idx(2,1)) * EDMfact - conjg(Gf(idxT(1,2)) * EDMfact)
            E(4,ind,i2) = E(4,ind,i2) + Gf(idx(2,2)) * EDMfact - conjg(Gf(idxT(2,2)) * EDMfact)
             
          end do
          end if
       end do
!$OMP end parallel do

     else if ( hasEDM ) then

!$OMP parallel do default(shared), private(io,iu,ind,idx,idxT)
       do iu = 1 , r%n
          io = r%r(iu)
          if ( l_ncol(io) /= 0 ) then
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

            ju = pvt%r(l_col(ind))
            do ii = 1, 2
              do jj = 1, 2
                idx (ii,jj) = index(Gf_tri,2*(iu-1)+ii,2*(ju-1)+jj)
                idxT(jj,ii) = index(Gf_tri,2*(ju-1)+jj,2*(iu-1)+ii)
              end do
            end do

             D(1,ind,i1) = D(1,ind,i1) + Gf(idx(1,1)) * DMfact - conjg(Gf(idxT(1,1)) * DMfact)
             D(2,ind,i1) = D(2,ind,i1) + Gf(idx(1,2)) * DMfact - conjg(Gf(idxT(2,1)) * DMfact)
             D(3,ind,i1) = D(3,ind,i1) + Gf(idx(2,1)) * DMfact - conjg(Gf(idxT(1,2)) * DMfact)
             D(4,ind,i1) = D(4,ind,i1) + Gf(idx(2,2)) * DMfact - conjg(Gf(idxT(2,2)) * DMfact)
 
             E(1,ind,i2) = E(1,ind,i2) + Gf(idx(1,1)) * EDMfact - conjg(Gf(idxT(1,1)) * EDMfact)
             E(2,ind,i2) = E(2,ind,i2) + Gf(idx(1,2)) * EDMfact - conjg(Gf(idxT(2,1)) * EDMfact)
             E(3,ind,i2) = E(3,ind,i2) + Gf(idx(2,1)) * EDMfact - conjg(Gf(idxT(1,2)) * EDMfact)
             E(4,ind,i2) = E(4,ind,i2) + Gf(idx(2,2)) * EDMfact - conjg(Gf(idxT(2,2)) * EDMfact)
          end do
          end if
       end do
!$OMP end parallel do

     else if ( calc_q ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,sp,sp_col,ind,idx,idxT,sind), &
!$OMP&  reduction(-:q)
       do iu = 1 , r%n
          io = r%r(iu)
          if ( l_ncol(io) /= 0 ) then

          sp = s_ptr(io)
          sp_col => s_col(sp+1:sp+s_ncol(io))

          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

            ju = pvt%r(l_col(ind))
            do ii = 1, 2
              do jj = 1, 2
                idx (ii,jj) = index(Gf_tri,2*(iu-1)+ii,2*(ju-1)+jj)
                idxT(jj,ii) = index(Gf_tri,2*(ju-1)+jj,2*(iu-1)+ii)
              end do
            end do

            sind = sp + SFIND(sp_col, l_col(ind))
            if ( sp < sind ) then
              q = q - aimag((Gf(idx(1,1)) - conjg(Gf(idxT(1,1)))) * Sk(sind))
              q = q - aimag((Gf(idx(2,2)) - conjg(Gf(idxT(2,2)))) * Sk(sind))
            end if
            D(1,ind,i1) = D(1,ind,i1) + Gf(idx(1,1)) * DMfact - conjg(Gf(idxT(1,1)) * DMfact)
            D(2,ind,i1) = D(2,ind,i1) + Gf(idx(1,2)) * DMfact - conjg(Gf(idxT(2,1)) * DMfact)
            D(3,ind,i1) = D(3,ind,i1) + Gf(idx(2,1)) * DMfact - conjg(Gf(idxT(1,2)) * DMfact)
            D(4,ind,i1) = D(4,ind,i1) + Gf(idx(2,2)) * DMfact - conjg(Gf(idxT(2,2)) * DMfact)
             
          end do
          end if
       end do
!$OMP end parallel do

     else

!$OMP parallel do default(shared), private(io,iu,ind,idx,idxT)
       do iu = 1 , r%n
          io = r%r(iu)
          if ( l_ncol(io) /= 0 ) then
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

            ju = pvt%r(l_col(ind))
            do ii = 1, 2
              do jj = 1, 2
                idx (ii,jj) = index(Gf_tri,2*(iu-1)+ii,2*(ju-1)+jj)
                idxT(jj,ii) = index(Gf_tri,2*(ju-1)+jj,2*(iu-1)+ii)
              end do
            end do

            D(1,ind,i1) = D(1,ind,i1) + Gf(idx(1,1)) * DMfact - conjg(Gf(idxT(1,1)) * DMfact)
            D(2,ind,i1) = D(2,ind,i1) + Gf(idx(1,2)) * DMfact - conjg(Gf(idxT(2,1)) * DMfact)
            D(3,ind,i1) = D(3,ind,i1) + Gf(idx(2,1)) * DMfact - conjg(Gf(idxT(1,2)) * DMfact)
            D(4,ind,i1) = D(4,ind,i1) + Gf(idx(2,2)) * DMfact - conjg(Gf(idxT(2,2)) * DMfact)
          end do
          end if
       end do
!$OMP end parallel do

     end if

    else
      ! *************** !
      ! Non-Equilibrium !
      ! *************** !
      !
      !  In the non-equilbrium case we have to use the triple matrix product
      !
      !      DM = \int GF Gamma GF^Dag dz
      !
      !  to calcualted the density matrix.
      !
      !  GF already contains the triple matrix product. So, we only
      !  need to add the elements of tri-diagonal matrix to the correct elements
      !  of the normal sparse matrix.
     
     if ( hasEDM ) then

!$OMP parallel do default(shared), private(io,iu,ind,idx)
       do iu = 1 , r%n
          io = r%r(iu)
          if ( l_ncol(io) /= 0 ) then
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

            do ii = 1, 2
              do jj = 1, 2
                idx (ii,jj) = index(Gf_tri,2*(iu-1)+ii,2*(ju-1)+jj)
              end do
            end do

            D(1,ind,i1) = D(1,ind,i1) + Gf(idx(1,1)) * DMfact
            D(2,ind,i1) = D(2,ind,i1) + Gf(idx(1,2)) * DMfact
            D(3,ind,i1) = D(3,ind,i1) + Gf(idx(2,1)) * DMfact
            D(4,ind,i1) = D(4,ind,i1) + Gf(idx(2,2)) * DMfact


            E(1,ind,i2) = E(1,ind,i2) + Gf(idx(1,1)) * EDMfact
            E(2,ind,i2) = E(2,ind,i2) + Gf(idx(1,2)) * EDMfact
            E(3,ind,i2) = E(3,ind,i2) + Gf(idx(2,1)) * EDMfact
            E(4,ind,i2) = E(4,ind,i2) + Gf(idx(2,2)) * EDMfact
          end do
          end if
       end do
!$OMP end parallel do

     else

!$OMP parallel do default(shared), private(io,iu,ind,idx)
       do iu = 1 , r%n
          io = r%r(iu)
          if ( l_ncol(io) /= 0 ) then
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

            do ii = 1, 2
              do jj = 1, 2
                idx (ii,jj) = index(Gf_tri,2*(iu-1)+ii,2*(ju-1)+jj)
              end do
            end do

            D(1,ind,i1) = D(1,ind,i1) + Gf(idx(1,1)) * DMfact
            D(2,ind,i1) = D(2,ind,i1) + Gf(idx(1,2)) * DMfact
            D(3,ind,i1) = D(3,ind,i1) + Gf(idx(2,1)) * DMfact
            D(4,ind,i1) = D(4,ind,i1) + Gf(idx(2,2)) * DMfact
          end do
          end if
       end do
!$OMP end parallel do

     end if

   end if

#ifdef TRANSIESTA_TIMING
    call timer('TS_add_DM', 2)
#endif

  end subroutine add_DM

  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, GFinv_tri, r, pvt, &
       N_Elec, Elecs, spH, spS)

    use class_Sparsity
    use class_zSpData1D
    use class_zSpData3D
    use class_zTriMat
    use ts_electrode_m
    use m_ts_tri_scat_spinor, only : insert_Self_Energies
    use m_ts_cctype, only : ts_c_idx

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    type(zTriMat), intent(inout) :: GFinv_tri
    type(tRgn), intent(in) :: r, pvt
    integer, intent(in) :: N_Elec
    type(electrode_t), intent(in) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D), intent(inout) :: spS
    type(zSpData3D), intent(inout) :: spH

    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: H(:,:,:), S(:)
    complex(dp), pointer :: Gfinv(:)
    integer :: io, iu, ju, ii, jj, ind, idx(2,2)

    if ( cE%fake ) return

#ifdef TRANSIESTA_TIMING
    call timer('TS-prep',1)
#endif

    Z = cE%e

    sp => spar(spH)
    H => val (spH)
    S => val (spS)

    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    Gfinv => val(Gfinv_tri)
    ! Initialize
    GFinv(:) = cmplx(0._dp,0._dp,dp)

!$OMP parallel default(shared), private(io,iu,ju,ii,jj,ind,idx)

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
!$OMP do
    do iu = 1, r%n
       io = r%r(iu)
       if ( l_ncol(io) /= 0 ) then

       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 

          ! Transpose to match phase convention in ts_sparse_helper (- phase) 
          ju = pvt%r(l_col(ind))
          do ii = 1, 2
            do jj = 1, 2
              idx(jj,ii) = index(GFinv_tri,2*(ju-1)+jj,2*(iu-1)+ii)
            end do
          end do

          GFinv(idx(1,1)) = Z * S(ind) - H(1,1,ind)
          GFinv(idx(2,1)) =            - H(2,1,ind)
          GFinv(idx(1,2)) =            - H(1,2,ind)
          GFinv(idx(2,2)) = Z * S(ind) - H(2,2,ind)
       end do
       end if
    end do
!$OMP end do

    do io = 1 , N_Elec
      call insert_Self_Energies(Gfinv_tri, Gfinv, pvt, Elecs(io))
    end do

!$OMP end parallel

#ifdef TRANSIESTA_TIMING
    call timer('TS-prep',2)
#endif

  end subroutine prepare_invGF
   
end module m_ts_trik_spinor
