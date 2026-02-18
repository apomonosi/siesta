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

module m_ts_fullk_spinor

  use precision, only : dp

  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM, update_zDM
  use m_ts_dm_update, only : add_k_DM
  
  use m_ts_weight_spinor, only : weight_DM_spinor
  use m_ts_weight, only : TS_W_K_METHOD
  use m_ts_weight, only : TS_W_K_CORRELATED
  use m_ts_weight, only : TS_W_K_UNCORRELATED

  use m_ts_method, only : orb_offset, no_Buf
  
  implicit none
  
  public :: ts_fullk_spinor
  
  private
  
contains
  
! ##################################################################
! ##                                                              ##       
! ##                       "TRANSIESTA"                           ##
! ##                                                              ##       
! ##          Non-equilibrium Density Matrix Subroutine           ##
! ##                   to be called from SIESTA                   ##
! ##                                                              ## 
! ## Originally:                By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##               Kurt Stokbro, ks@mic.dtu.dk                    ## 
! ##               Mikroelektronik Centret (MIC)                  ##
! ##           Technical University of Denmark (DTU)              ##
! ##                                                              ##
! ## Currently:                 By                                ##
! ##           Nick Papior Andersen, nickpapior@gmail.com         ##
! ##           Technical University of Denmark (DTU)              ##
! ##                                                              ##
! ## This code has been fully re-created to conform with the      ##
! ## sparsity patterns in SIESTA. Thus the memory requirements    ##
! ## has been greatly reduced.                                    ##
! ##                                                              ##
! ##################################################################
!
!
! Tight-binding density-matrix/transport program for the SIESTA
! package.
! Copyright by Mads Brandbyge, 1999, 2000, 2001, 2002.
! Copyright by Nick Papior Andersen, 2013.
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the authors.
!
  
  subroutine ts_fullk_spinor(N_Elec,Elecs, &
       nq,uGF, &
       ucell, nspin, na_u, lasto, &
       sp_dist, sparse_pattern, &
       no_u, n_nzs, &
       Hs, Ss, DM, EDM, Ef, DE_NEGF)

    use units, only : Pi, eV
    use parallel, only : Node, Nodes
    use m_spin, only : spin
#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_dSpData3D
    use class_zSpData1D
    use class_zSpData2D
    use class_zSpData3D

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
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq
    use m_ts_contour_neq, only : nEq_E, has_cE_nEq
    use m_ts_contour_neq, only : N_nEq_ID, c2weight_neq
    use m_ts_contour_eq, only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use m_iterator

    ! Gf calculation
    use m_ts_full_scat

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
    real(dp), intent(in) :: Hs(n_nzs,spin%H), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,spin%DM), EDM(n_nzs,spin%EDM)
    real(dp), intent(in) :: Ef
    real(dp), intent(inout) :: DE_NEGF

! ****************** Electrode variables *********************
    complex(dp), pointer :: GFGGF_work(:) => null()
! ************************************************************

! ******************* Computational arrays *******************
    integer :: nzwork, n_s
    complex(dp), allocatable, target :: zwork(:), GF(:)

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
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    integer :: NE
    integer :: index_dq !< Index for the current charge calculation @ E == mu
    real(dp) :: kw, kpt(3), bkpt(3), dq_mu
    complex(dp) :: W, ZW
! ************************************************************

! ******************** Loop variables ************************
    type(itt1) :: Kp
    integer, pointer :: ikpt
    integer :: iEl, iID, up_nzs
    integer :: iE, imu, io, idx
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr, no_u_TS, off, no, no_col, no_Els
! ************************************************************

#ifdef TS_SOC_DEBUG
    type(Sparsity), pointer :: s
    integer :: lnr, lio, lind, jo, col
    integer :: iu, iuT, ju, juT, i1, i2
    integer :: ind, nr
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: Hk(:,:,:), Sk(:)
    complex(dp), pointer :: D(:,:,:), E(:,:,:)
    character(len=1024) :: filename
#endif
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif
    if ( nspin /= 8 ) call die("ts_fullk_spinor: Must be called with spin=8")


    ! Number of supercells
    n_s = size(sc_off,dim=2)

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_Buf

    ! Number of elements that are transiesta updated
    up_nzs = nnzs(tsup_sp_uc)

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
    ! Hence, we will minimum require the full matrix...
    call UC_minimum_worksize(IsVolt, N_Elec, Elecs, no)
    nzwork = 4 * max(no_u_TS ** 2, no)
    allocate(zwork(nzwork),stat=ierr)
    if (ierr/=0) call die('Could not allocate space for zwork')
    call memory('A','Z',nzwork,'transiesta')

    ! We only need a partial size of the Green function
    io = 0
    no_Els = 0
    no_col = no_u_TS
    do iEl = 1 , N_Elec
      no = Elecs(iEl)%device_orbitals()
      if ( Elecs(iEl)%DM_update == 0 ) then
        no_col = no_col - no
      end if
      no_Els = no_Els + no
      io = io + no * no
    end do
    no = no_col
    ! when bias is needed we need the entire GF column
    ! for all the electrodes (at least some of the contour points needs this)
    if ( IsVolt ) then
      no = max(no, no_Els)
    end if
    no = no * no_u_TS
    no = 4 * max(no, io)
    allocate(GF(no),stat=ierr)
    if (ierr/=0) call die('Could not allocate space for GFpart')
    call memory('A','Z',no,'transiesta')

    no = 0
    do iEl = 1 , N_Elec

       ! This seems stupid, however, we never use the Sigma and
       ! GF at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called
       ! first the Sigma is assigned and then 
       ! it is required that prepare_GF_inv is called
       ! immediately (which it is)
       ! Hence the GF must NOT be used in between these two calls!
       io = 2 * Elecs(iEl)%device_orbitals()
       Elecs(iEl)%Sigma => GF(no+1:no+io**2)
       no = no + io ** 2

    end do

    if ( IsVolt ) then
       ! we need only allocate one work-array for
       ! Gf.G.Gf^\dagger
       call re_alloc(GFGGF_work,1,4*maxval(Elecs%device_orbitals())**2,routine='transiesta')
    end if

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(no_u,MPI_COMM_WORLD,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1            ,fdist,name='TS-fake dist')
#endif

    ! The Hamiltonian and overlap matrices (in Gamma calculations
    ! we will not have any phases, hence, it makes no sense to
    ! have the arrays in complex)
    call newzSpData3D(ts_sp_uc,2,2,fdist,spH,name='TS spH', sparsity_dim=3)
    call newzSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    ! If we have a bias calculation we need additional arrays.
    ! If not bias we don't need the update arrays (we already have
    ! all information in tsup_sp_uc (spDMu))

    ! Allocate space for global sparsity arrays
    ! These are complex arrays so we will need only 4 (not 8) components each
    ! [UP-UP, UP-DOWN, DOWN-UP, DOWN-DOWN]
    no = max(N_mu,N_nEq_id)
    call newzSpData3D(tsup_sp_uc,4,no,fdist, spuDM, name='TS spuDM', &
        sparsity_dim=2)
    if ( Calc_Forces ) then
       call newzSpData3D(tsup_sp_uc,4,N_mu,fdist, spuEDM, name='TS spuEDM', &
           sparsity_dim=2)
    end if
    
    if ( IsVolt ) then
       ! Allocate space for update arrays, local sparsity arrays
       call newdSpData3D(ltsup_sp_sc,spin%DM,N_mu,    sp_dist,spDM   ,name='TS spDM', &
           sparsity_dim=2)
       call newdSpData3D(ltsup_sp_sc,spin%DM,N_nEq_id,sp_dist,spDMneq,name='TS spDM-neq', &
           sparsity_dim=2)
       if ( Calc_Forces ) then
          call newdSpData3D(ltsup_sp_sc,spin%EDM,N_mu,sp_dist,spEDM  ,name='TS spEDM', &
              sparsity_dim=2)
       end if
    end if

    ! Initialize the charge correction scheme (will return if not used)
    call ts_dq%initialize_dq()

    ! Total number of energy points
    NE = N_Eq_E() + N_nEq_E()

     ! start the itterator
    call itt_init  (Kp,end=ts_kpoint_scf%N)
    ! point to the index iterator
    call itt_attach(Kp,cur=ikpt)

    call init_DM(sp_dist, sparse_pattern, &
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
            nspin, Ef, &
            N_Elec, Elecs, no_u, n_s, & ! electrodes, SIESTA size
            n_nzs, Hs, Ss, sc_off, &
            spH, spS, kpt, &
            nzwork, zwork)

#ifdef TS_SOC_DEBUG    
if (Node .eq. 0) then
write (filename, "('Hk-idx-',I0,'.dat')") ikpt
open(unit=340, file=filename)
write (filename, "('Hk-',I0,'.dat')") ikpt
open(unit=440, file=filename)
write (filename, "('Sk-',I0,'.dat')") ikpt
open(unit=540, file=filename)

s => spar(spH)
call attach(s,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, nrows=nr)
Hk => val(spH)
Sk => val(spS)
do io = 1 , nr
  if ( l_ncol(io) /= 0 ) then
    iu = io - orb_offset(io)
    iuT = iu - offset(N_Elec,Elecs,io)
    do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
      juT = l_col(ind) - orb_offset(l_col(ind))
      ju = juT - offset(N_Elec,Elecs,l_col(ind))
      write(340, "(4I10)") ind, iu, ju, ikpt
      write(440, "(8(spE30.15E4))") Hk(:,:,ind)
      write(540, "(2(spE30.15E4))") real (Sk(ind), dp), aimag(Sk(ind))
    end do
  end if
end do
    
close (340)
close (440)
close (540)
end if
#endif

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
          call prepare_invGF(cE, no_u_TS, zwork, &
               N_Elec, Elecs, &
               spH=spH , spS=spS)

#ifdef TS_SOC_DEBUG
write (filename, "('invGF-',I0,'-',I0'.dat')") iE, ikpt
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_u_TS
  write(340, "('('spE28.15E4 spE28.15E4'j)')") zwork(io) / eV
end do
close (340)
#endif
          ! *******************
          ! * calc GF         *
          ! *******************
          if ( all(Elecs(:)%DM_update /= 0) ) then
             call calc_GF(cE,2*no_u_TS, zwork, GF)

          else
             call calc_GF_Part_nc(cE, no_u, no_u_TS, no_col, &
                  N_Elec, Elecs, &
                  zwork, GF)

          end if

#ifdef TS_SOC_DEBUG
write (filename, "('GF-',I0,'-',I0,'.dat')") iE, ikpt
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_u_TS
write(340, "('('spE28.15E4 spE28.15E4'j)')") GF(io) * eV
end do
close (340)
#endif
          
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
                  no_u_TS, no_col, GF, &
                  N_Elec, Elecs, &
                  DMidx=mus(imu)%ID, spS=spS, q=dq_mu, &
                  iE=iE, ikpt=ikpt)
              ts_dq%mus(imu)%dq(index_dq) = ts_dq%mus(imu)%dq(index_dq) + dq_mu * kw
            else
              call add_DM( spuDM, W, spuEDM, ZW, &
                  no_u_TS, no_col, GF, &
                  N_Elec, Elecs, &
                  DMidx=mus(imu)%ID, &
                  iE=iE, ikpt=ikpt)
            end if
          end do

          ! step energy-point
          iE = iE + Nodes
          cE = Eq_E(iE,step=Nodes) ! we read them backwards
       end do ! eq. Energy contour

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

#ifdef TS_SOC_DEBUG    
if (Node .eq. 0) then
write (filename, "('zDM_eq-idx-',I0,'.dat')") ikpt
open(unit=340, file=filename)
do i1 = 1, N_mu
  write (filename, "('zDM_eq-',I0,'-',I0'.dat')") i1, ikpt
  open(unit=340+100*i1, file=filename)
end do

s => spar(spuDM)
call attach(s,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, nrows=nr)
D => val(spuDM)
do io = 1 , nr
  if ( l_ncol(io) /= 0 ) then
    iu = io - orb_offset(io)
    iuT = iu - offset(N_Elec,Elecs,io)
    do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
      juT = l_col(ind) - orb_offset(l_col(ind))
      ju = juT - offset(N_Elec,Elecs,l_col(ind))
      write(340, "(4I10)") ind, iu, ju, ikpt
      do i1 = 1, N_mu
        write(340+100*i1, "(8(spE30.15E4))") &
          real (D(1,ind,i1), dp), real (D(4,ind,i1), dp), &
          real (D(2,ind,i1), dp),-aimag(D(2,ind,i1)    ), &
          aimag(D(1,ind,i1)    ), aimag(D(4,ind,i1)    ), &
          real (D(3,ind,i1), dp), aimag(D(3,ind,i1)    )
      end do
    end do
  end if
end do
    
close (340)
do i1 = 1, N_mu
  close (340+100*i1)
end do
end if

#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,io)
#endif
#endif

       if ( .not. IsVolt ) then
          call update_zDM(sp_dist,sparse_pattern, n_nzs, &
               DM,  spuDM, Ef, &
               EDM, spuEDM, kpt, n_s, sc_off, &
               ikpt)

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
       ! if (tsweightmethod...)
       if ( TS_W_K_METHOD == TS_W_K_UNCORRELATED ) then
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       ! TODO NW if correlated this still have to be initiallized once!
       ! Is this done correctly?
       else if ( itt_first(Kp) ) then
          ! we only need to initialize once per spin
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       end if

#ifdef TS_SOC_DEBUG
write (filename, "('add_k_DMeq-idx-',I0,'-',I0,'.dat')") Node, ikpt
open(unit=1340+Node, file=filename)
write (filename, "('add_k_DMeq-',I0,'-',I0,'-',I0,'.dat')") Node, ikpt, 1
open(unit=1440+Node, file=filename)
write (filename, "('add_k_DMeq-',I0,'-',I0,'-',I0,'.dat')") Node, ikpt, 2
open(unit=1540+Node, file=filename)
#endif

       ! transfer equilibrium data to local sparsity arrays
       call add_k_DM(spDM, spuDM, spin%DM, N_mu, &
            spEDM, spuEDM, spin%EDM, N_mu, &
            n_s, sc_off, kpt, non_Eq = .false. )

#ifdef TS_SOC_DEBUG
    close (1340+Node)
    close (1440+Node)
    close (1540+Node)
#endif


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
          call prepare_invGF(cE, no_u_TS, zwork, &
               N_Elec, Elecs, &
               spH =spH , spS =spS)

#ifdef TS_SOC_DEBUG
write (filename, "('invGF_neq-',I0,'-',I0'.dat')") iE, ikpt
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_u_TS
  write(340, "('('spE28.15E4 spE28.15E4'j)')") zwork(io) / eV
end do
close (340)
#endif         

          ! *******************
          ! * calc GF         *
          ! *******************
          call calc_GF_Bias_nc(cE, no_u_TS, no_Els, &
               N_Elec, Elecs, &
               zwork, GF)

#ifdef TS_SOC_DEBUG
write (filename, "('GF_neq-',I0,'-',I0,'.dat')") iE, ikpt
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_Els
write(340, "('('spE28.15E4 spE28.15E4'j)')") GF(io) * eV
end do
close (340)
#endif

          ! ** At this point we have calculated the Green function

          ! ****************
          ! * save GF      *
          ! ****************
          off = 0
          do iEl = 1 , N_Elec
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution

             ! offset and number of orbitals
             no = Elecs(iEl)%device_orbitals() * 2

             call GF_Gamma_GF(Elecs(iEl), 2*no_u_TS, no, &
                  Gf(2*no_u_TS*off+1), zwork, size(GFGGF_work), GFGGF_work)

#ifdef TS_SOC_DEBUG
write (filename, "('Gamma-',I0,'-',I0,'-',I0,'.dat')") iEl, iE, ikpt
open(unit=340, file=filename)
do io = 1, no*no
  write(340, "('('spE28.15E4 spE28.15E4'j)')") Elecs(iEl)%Gamma(io) / eV
end do
close (340)
write (filename, "('GFGGF-',I0,'-',I0,'-',I0,'.dat')") iEl, iE, ikpt
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_u_TS
write(340, "('('spE28.15E4 spE28.15E4'j)')") zwork(io) * eV
end do
close (340)
#endif

             ! step to the next electrode position
             off = off + no
                
             do iID = 1 , N_nEq_ID
                
                if ( .not. has_cE_nEq(cE,iEl,iID) ) cycle
                
                call c2weight_neq(cE,iID,kw,W,imu,ZW)

                call add_DM( spuDM, W, spuEDM, ZW, &
                     no_u_TS, no_u_TS, zwork, &
                     N_Elec, Elecs, &
                     DMidx=iID, EDMidx=imu, &
                     is_eq = .false.)
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

#ifdef TS_SOC_DEBUG    
if (Node .eq. 0) then
  write (filename, "('zDM_neq-idx-',I0,'.dat')") ikpt
  open(unit=340, file=filename)

  do i1 = 1, N_mu
    write (filename, "('zDM_neq-',I0,'-',I0,'.dat')") i1, ikpt
    open(unit=340+100*i1, file=filename)
  end do

  s => spar(spuDM)
  call attach(s,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, nrows=nr)
  D => val(spuDM)
  do io = 1 , nr
    if ( l_ncol(io) /= 0 ) then
      iu = io - orb_offset(io)
      iuT = iu - offset(N_Elec,Elecs,io)
      do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
        juT = l_col(ind) - orb_offset(l_col(ind))
        ju = juT - offset(N_Elec,Elecs,l_col(ind))
        write(340, "(4I10)") ind, iu, ju, ikpt
        do i1 = 1, N_mu
          write(340+100*i1, "(8(spE30.15E4))"), &
             real(D(1,ind,i1), dp),  real(D(4,ind,i1), dp), &
             real(D(2,ind,i1), dp),-aimag(D(2,ind,i1)    ), &
            aimag(D(1,ind,i1)    ), aimag(D(4,ind,i1)    ), &
             real(D(3,ind,i1), dp), aimag(D(3,ind,i1)    )
        end do
      end do
    end if
  end do
      
  close (340)
  do i1 = 1, N_mu
    close (340+100*i1)
  end do
end if
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',1)
#endif

#ifdef TS_SOC_DEBUG
write (filename, "('add_k_DMneq-idx-',I0,'-',I0,'.dat')") Node, ikpt
open(unit=1340+Node, file=filename)
write (filename, "('add_k_DMneq-',I0,'-',I0,'-',I0,'.dat')") Node, ikpt, 1
open(unit=1440+Node, file=filename)
write (filename, "('add_k_DMneq-',I0,'-',I0,'-',I0,'.dat')") Node, ikpt, 2
open(unit=1540+Node, file=filename)
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
          
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM, spDM, Ef=Ef, &
               EDM=EDM, spEDM=spEDM, ipnt=ltsup_sc_pnt)
       else if ( itt_last(Kp) ) then ! TS_W_K_METHOD == TS_W_K_CORRELATED
          call weight_DM_spinor( N_Elec, Elecs, N_mu, mus, na_u, lasto, &
              sp_dist, sparse_pattern, Ss, &
              spDM, spDMneq, spEDM, n_s, sc_off, DE_NEGF)
          call mean_weight()
          
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM, spDM, Ef=Ef, &
               EDM=EDM, spEDM=spEDM, ipnt=ltsup_sc_pnt)          
       end if

#ifdef TS_SOC_DEBUG
    close (1340+Node)
    close (1440+Node)
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',2)
#endif

       do iEl = 1, N_Elec
         call reread_Gamma_Green(Elecs(iEl), uGF(iEl), NE, 1)
       end do

    end do ! kp

    call itt_destroy(Kp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA SCF'
#endif

!***********************
! CLEAN UP
!***********************

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    call delete(spuEDM)

    call delete(spDM)
    call delete(spDMneq)
    call delete(spEDM)

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    call memory('D','Z',size(zwork)+size(GF),'transiesta')
    deallocate(zwork,GF)

    ! In case of voltage calculations we need a work-array for
    ! handling the GF.Gamma.Gf^\dagger multiplication
    if ( IsVolt ) then
       call de_alloc(GFGGF_work, routine='transiesta')
    end if

    ! Nullify external pointers
    do iEl = 1, N_Elec
      nullify(Elecs(iEl)%Sigma)
    end do
   
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif
  contains

  ! Wrapper for simple weights
  subroutine mean_weight()
    use m_ts_contour_neq, only: id2mu

    real(dp), pointer :: eqDM(:,:,:), neqDM(:,:,:)
    real(dp), pointer :: qEDM(:,:,:)
    real(dp) :: w
    integer :: ID, imu

    eqDM => val(spDM)
    neqDM => val(spDMneq)
    qEDM => val(spEDM)

    w = 1._dp / real(size(eqDM, 3), dp)

    do ID = 1, size(neqDM, 3)
      imu = ID2mu(ID)
      eqDM(:,:,imu) = eqDM(:,:,imu) + neqDM(:,:,ID)
    end do
    eqDM(:,:,1) = sum(eqDM, 3) * w
    qEDM(:,:,1) = sum(qEDM, 3) * w

  end subroutine mean_weight

#ifdef TS_SOC_DEBUG
  pure function offset(N_Elec,Elecs,io)
  use ts_electrode_m

  integer, intent(in) :: N_Elec
  type(electrode_t), intent(in) :: Elecs(N_Elec)
  integer, intent(in) :: io
  integer :: offset
  offset = sum(Elecs(:)%device_orbitals(), &
      MASK=(Elecs(:)%DM_update==0) .and. Elecs(:)%idx_o <= io )
  end function offset
#endif

  end subroutine ts_fullk_spinor

  ! Update DM
  ! These routines are supplied for easy update of the update region
  ! sparsity patterns
  ! Note that these routines implement the usual rho(Z) \propto - GF
  ! TODO NW iE should be removed from arguments again
  subroutine add_DM(DM, DMfact,EDM, EDMfact, &
      no1,no2,GF, &
      N_Elec,Elecs, &
      DMidx, EDMidx, &
      spS, q, &
      is_eq, &
      iE, ikpt)

    use intrinsic_missing, only: SFIND

    use class_Sparsity
    use class_zSpData1D
    use class_zSpData3D
    use ts_electrode_m
#ifdef TS_SOC_DEBUG
    use parallel, only : Node
#endif

    ! The DM and EDM equivalent matrices
    type(zSpData3D), intent(inout) :: DM
    complex(dp), intent(in) :: DMfact
    type(zSpData3D), intent(inout) :: EDM
    complex(dp), intent(in) :: EDMfact
    ! The size of GF
    integer, intent(in) :: no1, no2
    ! The Green function
    complex(dp), intent(in) :: GF(2,no1,2,no2)
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
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: s_ncol(:), s_ptr(:), s_col(:), sp_col(:)
    complex(dp), pointer :: D(:,:,:), E(:,:,:), Sk(:)
    integer :: io, ind, nr
    integer :: sp, sind
    integer :: iu, iuT, ju, juT, i1, i2
    logical :: lis_eq, hasEDM, calc_q
    
    integer, intent(in), optional :: iE, ikpt

! #ifdef TS_SOC_DEBUG
! character(len=1024) :: filename
! write(*,*), iE, ikpt, 100+100*ikpt+iE, 500+100*ikpt+iE
! write (filename, "('DM-idx-',I2,'-',I2,'.dat')") iE, ikpt
! open(unit=100+100*ikpt+iE, file=filename)

! write (filename, "('zDM-',I2,'-',I2,'.dat')") iE, ikpt
! open(unit=500+100*ikpt+iE, file=filename)
! #endif
    lis_eq = .true.
    if ( present(is_eq) ) lis_eq = is_eq

    calc_q = present(q) .and. present(spS)

    ! Remember that this sparsity pattern HAS to be in Global UC
    s => spar(DM)
    call attach(s,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
        nrows=nr)
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

    if ( lis_eq ) then

      if ( calc_q ) then
        q = 0._dp
        s => spar(spS)
        Sk => val(spS)
        call attach(s, n_col=s_ncol, list_ptr=s_ptr, list_col=s_col)
      end if

      if ( calc_q .and. hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,sp,sp_col,iu,iuT,ind,ju,juT,sind), &
!$OMP&  reduction(-:q)
        do io = 1 , nr
          ! Quickly go past the buffer atoms...
          if ( l_ncol(io) /= 0 ) then

            sp = s_ptr(io)
            sp_col => s_col(sp+1:sp+s_ncol(io))

            ! The update region equivalent GF part
            iu = io - orb_offset(io)
            iuT = iu - offset(N_Elec,Elecs,io)

            do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
              
              juT = l_col(ind) - orb_offset(l_col(ind))
              ju = juT - offset(N_Elec,Elecs,l_col(ind))

              ! Search for overlap index, note it is transposed (- phase)
              ! So this is S(ju,iu)
              sind = sp + SFIND(sp_col, l_col(ind))
              if ( sp < sind ) then
                ! For charge calculation it is Tr[(G-G^\dagger).S] i.e. the matrix
                ! multiplication
                q = q - aimag((Gf(1,iu,1,ju) - conjg(Gf(1,juT,1,iuT))) * Sk(sind))
                q = q - aimag((Gf(2,iu,2,ju) - conjg(Gf(2,juT,2,iuT))) * Sk(sind))
              end if

              ! We still need to multiply with complex k-phase later, so 
              ! we do not split the imaginary and real part here.
              ! We will do this when copying to the real array in update_zDM
              ! at the end of the contour.
              !
              ! For now we store four complex numbers in the following order.
              ! TODO Make D 4-dimensional?
              !
              !        | D(1,:)    D(2,:) |     | GG11   GG12 |
              !   DM = |                  |  =  |             |
              !        | D(3,:)    D(4,:) |     | GG21   GG22 |
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)
              !
              !   where iu and ju are the orbital indices corresponding to ind
              !   and juT and iuT the transposed indices. 
              !
              !   The same for the energy density matrix.
              !

              D(1,ind,i1) = D(1,ind,i1) + Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
              D(2,ind,i1) = D(2,ind,i1) + Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
              D(3,ind,i1) = D(3,ind,i1) + Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
              D(4,ind,i1) = D(4,ind,i1) + Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)

              E(1,ind,i2) = E(1,ind,i2) + Gf(1,iu,1,ju) * EDMfact - conjg(Gf(1,juT,1,iuT) * EDMfact)
              E(2,ind,i2) = E(2,ind,i2) + Gf(1,iu,2,ju) * EDMfact - conjg(Gf(2,juT,1,iuT) * EDMfact)
              E(3,ind,i2) = E(3,ind,i2) + Gf(2,iu,1,ju) * EDMfact - conjg(Gf(1,juT,2,iuT) * EDMfact)
              E(4,ind,i2) = E(4,ind,i2) + Gf(2,iu,2,ju) * EDMfact - conjg(Gf(2,juT,2,iuT) * EDMfact)
              
            end do
            
          end if
        end do
!$OMP end parallel do

      else if ( hasEDM ) then

!$OMP parallel do default(shared), private(io,iu,iuT,ind,ju,juT)
        do io = 1 , nr
          if ( l_ncol(io) /= 0 ) then
            iu = io - orb_offset(io)
            iuT = iu - offset(N_Elec,Elecs,io)
            do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
              juT = l_col(ind) - orb_offset(l_col(ind))
              ju = juT - offset(N_Elec,Elecs,l_col(ind))

              ! We still need to multiply with complex k-phase later, so 
              ! we do not split the imaginary and real part here.
              ! We will do this when copying to the real array in update_zDM
              ! at the end of the contour.
              !
              ! For now we store four complex numbers in the following order.
              ! TODO Make D 4-dimensional?
              !
              !        | D(1,:)    D(2,:) |     | GG11   GG12 |
              !   DM = |                  |  =  |             |
              !        | D(3,:)    D(4,:) |     | GG21   GG22 |
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)
              !
              !   The same for the energy density matrix.
              !

              D(1,ind,i1) = D(1,ind,i1) + Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
              D(2,ind,i1) = D(2,ind,i1) + Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
              D(3,ind,i1) = D(3,ind,i1) + Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
              D(4,ind,i1) = D(4,ind,i1) + Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
#ifdef TS_SOC_DEBUG    
              ! write(100+100*ikpt+iE, "(4I10)") ind, iu, ju, Node
              ! write(500+100*ikpt+iE, "(10(spE30.15E4))") D(:,ind,i1), DMfact
#endif
              E(1,ind,i2) = E(1,ind,i2) + Gf(1,iu,1,ju) * EDMfact - conjg(Gf(1,juT,1,iuT) * EDMfact)
              E(2,ind,i2) = E(2,ind,i2) + Gf(1,iu,2,ju) * EDMfact - conjg(Gf(2,juT,1,iuT) * EDMfact)
              E(3,ind,i2) = E(3,ind,i2) + Gf(2,iu,1,ju) * EDMfact - conjg(Gf(1,juT,2,iuT) * EDMfact)
              E(4,ind,i2) = E(4,ind,i2) + Gf(2,iu,2,ju) * EDMfact - conjg(Gf(2,juT,2,iuT) * EDMfact)
            end do
          end if
        end do
!$OMP end parallel do

      else if ( calc_q ) then

!$OMP parallel do default(shared), &
!$OMP&   private(io,sp,sp_col,iu,iuT,ind,ju,juT,sind), &
!$OMP&   reduction(-:q)
        do io = 1 , nr
          if ( l_ncol(io) /= 0 ) then
            sp = s_ptr(io)
            sp_col => s_col(sp+1:sp+s_ncol(io))
            iu = io - orb_offset(io)
            iuT = iu - offset(N_Elec,Elecs,io)
            do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
              juT = l_col(ind) - orb_offset(l_col(ind))
              ju = juT - offset(N_Elec,Elecs,l_col(ind))
              sind = sp + SFIND(sp_col, l_col(ind))
              if ( sp < sind ) then
                q = q - aimag((Gf(1,iu,1,ju) - conjg(Gf(1,juT,1,iuT))) * Sk(sind))
                q = q - aimag((Gf(2,iu,2,ju) - conjg(Gf(2,juT,2,iuT))) * Sk(sind))
              end if

              ! We still need to multiply with complex k-phase later, so 
              ! we do not split the imaginary and real part here.
              ! We will do this when copying to the real array in update_zDM
              ! at the end of the contour.
              !
              ! For now we store four complex numbers in the following order.
              ! TODO Make D 4-dimensional?
              !
              !        | D(1,:)    D(2,:) |     | GG11   GG12 |
              !   DM = |                  |  =  |             |
              !        | D(3,:)    D(4,:) |     | GG21   GG22 |
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)

              D(1,ind,i1) = D(1,ind,i1) + Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
              D(2,ind,i1) = D(2,ind,i1) + Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
              D(3,ind,i1) = D(3,ind,i1) + Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
              D(4,ind,i1) = D(4,ind,i1) + Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
            end do
          end if
        end do
!$OMP end parallel do

      else

!$OMP parallel do default(shared), private(io,iu,iuT,ind,ju,juT)
        do io = 1 , nr
          if ( l_ncol(io) /= 0 ) then
            iu = io - orb_offset(io)
            iuT = iu - offset(N_Elec,Elecs,io)
            do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
              juT = l_col(ind) - orb_offset(l_col(ind))
              ju = juT - offset(N_Elec,Elecs,l_col(ind))

              ! We still need to multiply with complex k-phase later, so 
              ! we do not split the imaginary and real part here.
              ! We will do this when copying to the real array in update_zDM
              ! at the end of the contour.
              !
              ! For now we store four complex numbers in the following order.
              ! TODO Make D 4-dimensional?
              !
              !        | D(1,:)    D(2,:) |     | GG11   GG12 |
              !   DM = |                  |  =  |             |
              !        | D(3,:)    D(4,:) |     | GG21   GG22 |
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)

              D(1,ind,i1) = D(1,ind,i1) + Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
              D(2,ind,i1) = D(2,ind,i1) + Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
              D(3,ind,i1) = D(3,ind,i1) + Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
              D(4,ind,i1) = D(4,ind,i1) + Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
            end do
          end if
        end do
!$OMP end parallel do

      end if

    else ! lis_eq

#ifndef TS_NOCHECKS
      if ( no1 /= no2 ) call die("Error in matrix dimensions")
#endif
      ! *************** !
      ! Non-Equilibrium !
      ! *************** !

      if ( hasEDM ) then

!$OMP parallel do default(shared), private(io,iu,ind,ju)
        do io = 1 , nr
          if ( l_ncol(io) /= 0 ) then
            iu = io - orb_offset(io)
            do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
              ju = l_col(ind) - orb_offset(l_col(ind))

              ! We still need to multiply with complex k-phase later, so 
              ! we do not split the imaginary and real part here.
              ! We will do this when copying to the real array in update_zDM
              ! at the end of the contour.
              !
              ! For now we store four complex numbers in the following order.
              ! TODO Make D 4-dimensional?
              !
              !        | D(1,:)    D(2,:) |     | GG11   GG12 |
              !   DM = |                  |  =  |             |
              !        | D(3,:)    D(4,:) |     | GG21   GG22 |
              !
              !   GGij = Gf(i,iu ,j,ju ) * DMfact
              !
              !   The same for the energy density matrix.
              !

              D(1,ind,i1) = D(1,ind,i1) + Gf(1,iu,1,ju) * DMfact
              D(2,ind,i1) = D(2,ind,i1) + Gf(1,iu,2,ju) * DMfact
              D(3,ind,i1) = D(3,ind,i1) + Gf(2,iu,1,ju) * DMfact
              D(4,ind,i1) = D(4,ind,i1) + Gf(2,iu,2,ju) * DMfact


              E(1,ind,i2) = E(1,ind,i2) + Gf(1,iu,1,ju) * EDMfact
              E(2,ind,i2) = E(2,ind,i2) + Gf(1,iu,2,ju) * EDMfact
              E(3,ind,i2) = E(3,ind,i2) + Gf(2,iu,1,ju) * EDMfact
              E(4,ind,i2) = E(4,ind,i2) + Gf(2,iu,2,ju) * EDMfact
            end do
          end if
        end do
!$OMP end parallel do

      else

!$OMP parallel do default(shared), private(io,iu,ind,ju)
        do io = 1 , nr
          if ( l_ncol(io) /= 0 ) then
            iu = io - orb_offset(io)
            do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
              ju = l_col(ind) - orb_offset(l_col(ind))

              ! We still need to multiply with complex k-phase later, so 
              ! we do not split the imaginary and real part here.
              ! We will do this when copying to the real array in update_zDM
              ! at the end of the contour.
              !
              ! For now we store four complex numbers in the following order.
              ! TODO Make D 4-dimensional?
              !
              !        | D(1,:)    D(2,:) |     | GG11   GG12 |
              !   DM = |                  |  =  |             |
              !        | D(3,:)    D(4,:) |     | GG21   GG22 |
              !
              !   GGij = Gf(i,iu ,j,ju ) * DMfact
             !
              !   where iu and ju are the orbital indices corresponding to ind.

              D(1,ind,i1) = D(1,ind,i1) + Gf(1,iu,1,ju) * DMfact
              D(2,ind,i1) = D(2,ind,i1) + Gf(1,iu,2,ju) * DMfact
              D(3,ind,i1) = D(3,ind,i1) + Gf(2,iu,1,ju) * DMfact
              D(4,ind,i1) = D(4,ind,i1) + Gf(2,iu,2,ju) * DMfact
            end do
          end if
        end do
!$OMP end parallel do

      end if
    end if
! #ifdef TS_SOC_DEBUG    
! close (100+100*ikpt+iE)
! close (500+100*ikpt+iE)
! #endif
  contains
    
     function offset(N_Elec,Elecs,io)
      integer, intent(in) :: N_Elec
      type(electrode_t), intent(in) :: Elecs(N_Elec)
      integer, intent(in) :: io
      integer :: offset
      offset = sum(Elecs(:)%device_orbitals(), &
          MASK=(Elecs(:)%DM_update == 0) .and. Elecs(:)%idx_o <= io )
    end function offset

  end subroutine add_DM


  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, no_u,GFinv, &
       N_Elec, Elecs, spH, spS)

    use class_zSpData1D
    use class_zSpData3D
    use class_Sparsity
    use ts_electrode_m
    use m_ts_cctype, only : ts_c_idx
    use m_ts_full_scat, only : insert_Self_Energies_NC

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in) :: no_u
    complex(dp), intent(out) :: GFinv(2,no_u,2,no_u)
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
    integer :: io, iu, ju, ind, nr

    if ( cE%fake ) return

#ifdef TRANSIESTA_TIMING
    call timer('TS-prep',1)
#endif

    Z = cE%e
    
    sp => spar(spH)
    H  => val (spH)
    S  => val (spS)

    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows_g=nr)

    ! Initialize
    GFinv(:,:,:,:) = cmplx(0._dp,0._dp,dp)

!$OMP parallel default(shared), private(io,iu,ju,ind)

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
!$OMP do
    do io = 1, nr

       if ( l_ncol(io) /= 0 ) then

       iu = io - orb_offset(io)
       
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 

          ! Transpose to match phase convention in ts_sparse_helper (- phase)
          ju = l_col(ind) - orb_offset(l_col(ind))
          
          GFinv(1,ju,1,iu) = Z * S(ind) - H(1,1,ind)
          GFinv(2,ju,1,iu) =            - H(2,1,ind)
          GFinv(1,ju,2,iu) =            - H(1,2,ind)
          GFinv(2,ju,2,iu) = Z * S(ind) - H(2,2,ind)

       end do
       end if
    end do
!$OMP end do

    do io = 1 , N_Elec
      call insert_Self_Energies_NC(no_u, Gfinv, Elecs(io))
    end do

!$OMP end parallel

#ifdef TRANSIESTA_TIMING
    call timer('TS-prep',2)
#endif

  end subroutine prepare_invGF
   
end module m_ts_fullk_spinor
