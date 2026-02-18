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

module m_ts_fullg_spinor
  use precision, only : dp

  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM
  use m_ts_dm_update, only : add_Gamma_DM

  use m_ts_weight_spinor, only : weight_DM_spinor

  use m_ts_method, only: orb_offset, no_Buf

  implicit none

  public :: ts_fullg_spinor

  private

contains

  subroutine ts_fullg_spinor(N_Elec,Elecs, &
      nq, uGF, nspin, na_u, lasto, &
      sp_dist, sparse_pattern, &
      no_u, n_nzs, &
      Hs, Ss, DM, EDM, Ef, DE_NEGF)

    use units, only : eV, Pi
    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta
#endif
    use m_spin, only : spin

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D
    use class_dSpData3D

    use ts_electrode_m
    ! Self-energy read
    use m_ts_gf
    ! Self-energy expansion
    use m_ts_elec_se

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
    integer :: ndwork, nzwork, n_s
    real(dp), pointer :: dwork(:,:,:)
    complex(dp), allocatable, target :: zwork(:), GF(:)

    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData2D) :: spH
    type(dSpData1D) :: spS
    ! local sparsity pattern in local SC pattern
    type(dSpData3D) :: spDM, spDMneq
    type(dSpData3D) :: spEDM ! only used if calc_forces
    ! The different sparse matrices that will surmount to the integral
    ! These two lines are in global update sparsity pattern (UC)
    type(dSpData3D) ::  spuDM
    type(dSpData3D) :: spuEDM ! only used if calc_forces
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    integer :: index_dq !< Index for the current charge calculation @ E == mu
    real(dp) :: kw, dq_mu
    complex(dp) :: W, ZW
    logical :: eq_full_Gf
! ************************************************************

! ******************** Loop variables ************************
    integer :: iEl, iID
    integer :: iE, imu, io, idx
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr, no_u_TS, off, no, no_col, no_ELs
    real(dp), parameter :: bkpt(3) = (/0._dp,0._dp,0._dp/)
! ************************************************************

#ifdef TS_SOC_DEBUG
    type(Sparsity), pointer :: s
    integer :: lnr, lio, lind, jo, col
    integer :: iu, iuT, ju, juT, i1, i2
    integer :: ind, nr
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: D(:,:,:), E(:,:,:)
    character(len=1024) :: filename
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    if ( nspin /= 4 .and. nspin /=8 ) & 
      call die("ts_fullg_spinor: Must be called with spin=4 or 8")

    ! Number of supercells (even though its gamma we
    ! can have different schemes...)
    n_s = size(sc_off,dim=2)

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_Buf

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
      if ( Elecs(iEl)%DM_update == 0 ) then ! no elements in electrode are updated
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
       io = 4 * Elecs(iEl)%device_orbitals() ** 2
       Elecs(iEl)%Sigma => GF(no+1:no+io)
       no = no + io

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
    call newdSpData2D(ts_sp_uc,spin%H, fdist,spH,name='TS spH', sparsity_dim=2)
    call newdSpData1D(ts_sp_uc, fdist,spS,name='TS spS')

    ! If we have a bias calculation we need additional arrays.
    ! If not bias we don't need the update arrays (we already have
    ! all information in tsup_sp_uc (spDMu))

    ! Allocate space for global sparsity arrays
    no = max(N_mu,N_nEq_id)
    call newdSpData3D(tsup_sp_uc,spin%DM,no,fdist, spuDM, name='TS spuDM', &
        sparsity_dim=2)
    ! assign dwork, this will problably come at the expence of
    ! two full reductions, however, we save some memory...
    ndwork = nnzs(tsup_sp_uc)
    dwork => val(spuDM)
    if ( Calc_Forces ) then
      call newdSpData3D(tsup_sp_uc,spin%EDM,N_mu,fdist, spuEDM, name='TS spuEDM', &
          sparsity_dim=2)
    end if

    if ( IsVolt ) then
      ! Allocate space for update arrays, local sparsity arrays
      call newdSpData3D(ltsup_sp_sc,spin%DM,N_mu,sp_dist,spDM, name='TS spDM', &
          sparsity_dim=2)
      call newdSpData3D(ltsup_sp_sc,spin%DM,N_nEq_id,sp_dist,spDMneq,name='TS spDM-neq', &
          sparsity_dim=2)
      if ( nnzs(ltsup_sp_sc) > ndwork ) then
        ! only update if this array is larger (should only happen in
        ! few processor setups
        ndwork = nnzs(ltsup_sp_sc)
        dwork => val(spDMneq)
      end if
      if ( Calc_Forces ) then
        call newdSpData3D(ltsup_sp_sc,spin%EDM,N_mu,sp_dist,spEDM, name='TS spEDM', &
            sparsity_dim=2)
      end if
    end if

    ! Whether we should always calculate the full Green function
    eq_full_Gf = all(Elecs(:)%DM_update /= 0)

    ! TODO NW Initialize the charge correction scheme (see m_ts_fullg.F90)
    call ts_dq%initialize_dq()
    call init_DM(sp_dist, sparse_pattern, &
        n_nzs, spin, DM, EDM, &
        tsup_sp_uc, Calc_Forces)

    ! In NC/SO case we use DM = i\int((GF - GF^\dagger)/2
    ! We inlcude the factor of one-half here
    kw = 1._dp / (Pi * 2._dp)

#ifdef TRANSIESTA_TIMING
    call timer('TS_HS',1)
#endif

    ! Work-arrays are for MPI distribution...
    call create_HS(sp_dist,sparse_pattern, &
         nspin, Ef, &
         N_Elec, Elecs, no_u, & ! electrodes, SIESTA size
         n_nzs, Hs, Ss, &
         spH, spS, &
         ndwork, dwork(:,1,1)) ! annoyingly we can't pass the full array!!!!!

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
       call read_next_GS(1, 1, bkpt, &
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
write (filename, "('invGF-',I0,'.dat')") iE
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_u_TS
  write(340, "('('spE28.15E4 spE28.15E4'j)')") zwork(io) / eV
end do
close (340)
#endif

#ifdef TS_DEV
io=100+iE+node
open(io,form='unformatted')
write(io) cE%e / eV
write(io) no_u_TS,no_u_TS
write(io) zwork(1:no_u_TS**2) / eV
close(io)
#endif

       ! *******************
       ! * calc GF         *
       ! *******************
       if ( eq_full_Gf ) then
          call calc_GF(cE, 2*no_u_TS, zwork, GF)


#ifdef TS_DEV
io=300+iE+node
open(io,form='unformatted')
write(io) cE%e / eV
write(io) no_u_TS, no_u_TS
write(io) gf(1:no_u_TS**2) * eV
close(io)
#endif

       else
          call calc_GF_Part_nc(cE, no_u, no_u_TS, no_col,&
               N_Elec, Elecs, &
               zwork, GF)

       end if

#ifdef TS_SOC_DEBUG
write (filename, "('GF-',I0,'.dat')") iE
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
             DMidx=mus(imu)%ID, spS=spS, q=dq_mu)
         ts_dq%mus(imu)%dq(index_dq) = ts_dq%mus(imu)%dq(index_dq) + dq_mu * kw
        else
         call add_DM( spuDM, W, spuEDM, ZW, &
             no_u_TS, no_col, GF, &
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
    ! Here the work size is multiplied by 2 because zwork is declared
    ! as complex but the function accepts a real array
    call my_full_G_reduce_3D(spuDM,nzwork*2,zwork,spin%DM,N_mu)
    if ( Calc_Forces ) then
       call my_full_G_reduce_3D(spuEDM,nzwork*2,zwork,spin%EDM,N_mu)
    end if
    call timer('TS_comm',2)
#endif

#ifdef TS_SOC_DEBUG    
if (Node .eq. 0) then
  open(unit=340, file="DM_eq-idx.dat")
  
  do i1 = 1, N_mu
    write (filename, "('DM_eq-',I0,'.dat')") i1
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
        write(340, "(4I10)") ind, iu, ju, 1
        do i1 = 1, N_mu
          write(340+100*i1, "(8(spE30.15E4))") D(:,ind,i1)
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

    if ( .not. IsVolt ) then
       call update_DM(sp_dist,sparse_pattern, n_nzs, &
            DM(:,:), spuDM, Ef=Ef, &
            EDM=EDM(:,:), spEDM=spuEDM, &
            UpSpGlobal = .true.)
       ! The remaining code segment only deals with
       ! bias integration... So we skip instantly

    else
       ! *****************
       ! * only things with non-Equilibrium contour...
       ! *****************

       ! initialize to zero
       ! local sparsity update patterns
       call init_val(spDM)
       call init_val(spDMneq)
       if ( Calc_Forces ) call init_val(spEDM)

       ! transfer data to local sparsity arrays
       call add_Gamma_DM(spDM,   spuDM, D_dim2=spin%DM, D_dim3=N_mu, &
            spEDM=spEDM, spuEDM=spuEDM, E_dim2=spin%EDM, E_dim3=N_mu)

#ifdef TRANSIESTA_TIMING
        call timer('TS_NEQ',1)
#endif

       ! *******************
       ! * NON-EQUILIBRIUM *
       ! *******************

       ! We have the definition of: Gamma = i(\Sigma - \Sigma^\dagger)
       ! (not with one half)
       ! Hence we also use half the contribution for all contour points
       ! i.e. keep kw unchanged

       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = Nodes - Node
       cE = nEq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(1, 1, bkpt, &
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
               spH =spH , spS=spS)

#ifdef TS_SOC_DEBUG
write (filename, "('invGF_neq-',I0,'.dat')") iE
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
write (filename, "('GF_neq-',I0,'.dat')") iE
open(unit=340, file=filename)
do io = 1, 4*no_u_TS*no_Els
write(340, "('('spE28.15E4 spE28.15E4'j)')") GF(io) * eV
end do
close (340)
#endif

#ifdef TS_DEV
io = 500 + iE + Node
open(io,form='unformatted')
write(io) cE%e / eV
write(io) no_u_TS, TotUsedOrbs(Elecs(1))
write(io) GF(1:no_u_TS * totusedorbs(Elecs(1))) * eV
write(io) no_u_TS, TotUsedOrbs(Elecs(2))
write(io) GF(no_u_TS*totusedorbs(Elecs(1))+1:no_u_TS * sum(totusedorbs(Elecs))) * eV
close(io)
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
write (filename, "('Gamma-',I0,'-',I0,'.dat')") iEl, iE
open(unit=340, file=filename)
do io = 1, no*no
  write(340, "('('spE28.15E4 spE28.15E4'j)')") Elecs(iEl)%Gamma(io) / eV
end do
close (340)
write (filename, "('GFGGF-',I0,'-',I0,'.dat')") iEl, iE
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
      call my_full_G_reduce_3D(spuDM, nzwork*2, zwork, spin%DM, N_nEq_id)
      if ( Calc_Forces ) then
         call my_full_G_reduce_3D(spuEDM, nzwork*2, zwork, spin%EDM, N_mu)
      end if
      call timer('TS_comm',2)
#endif

#ifdef TS_SOC_DEBUG    
if (Node .eq. 0) then
  open(unit=340, file="DM_neq-idx.dat")

  do i1 = 1, N_mu
    write (filename, "('DM_neq-',I0,'.dat')") i1
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
        write(340, "(4I10)") ind, iu, ju, 1
        do i1 = 1, N_mu
          write(340+100*i1, "(8(spE30.15E4))") D(:,ind,i1)
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

      ! 1. move from global UC to local SC
      ! 2. calculate the correct contribution by applying the weight
      ! 3. add the density to the real arrays
      call add_Gamma_DM(spDMneq, spuDM, D_dim2=spin%DM, D_dim3=N_nEq_id, &
          spEDM=spEDM,  spuEDM=spuEDM, E_dim2=spin%EDM, E_dim3=N_mu)

      call weight_DM_spinor( N_Elec, Elecs, N_mu, mus, na_u, lasto, &
          sp_dist, sparse_pattern, Ss, &
          spDM, spDMneq, spEDM, n_s, sc_off, DE_NEGF)
      ! TODO NW : remove
      ! call mean_weight()

      call update_DM(sp_dist,sparse_pattern, n_nzs, &
           DM(:,:), spDM, Ef=Ef, &
           EDM=EDM(:,:), spEDM=spEDM, ipnt=ltsup_sc_pnt)

#ifdef TRANSIESTA_TIMING
      call timer('TS_weight',2)
#endif

      ! We don't need to do anything here..
  end if


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
    ! TODO NW: remove
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

  end subroutine ts_fullg_spinor

  subroutine add_DM(DM, DMfact,EDM, EDMfact, &
    no1,no2,GF, &
    N_Elec,Elecs, &
    DMidx, EDMidx, &
    spS, q, &
    is_eq)

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData3D
    use ts_electrode_m
    use m_spin, only: spin

    ! The DM and EDM equivalent matrices
    type(dSpData3D), intent(inout) :: DM
    complex(dp), intent(inout) :: DMfact
    type(dSpData3D), intent(inout) :: EDM
    complex(dp), intent(inout) :: EDMfact
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
     type(dSpData1D), intent(in), optional :: spS
     !< Charge calculated at this energy-point
     !!
     !! This does not contain the additional factor 1/Pi
     real(dp), intent(inout), optional :: q
    logical, intent(in), optional :: is_eq

    if (spin%Ncol) then
      call add_DM_ncol(DM, DMfact,EDM, EDMfact,no1,no2,GF,N_Elec,Elecs,&
          DMidx,EDMidx,spS,q,is_eq)
    else ! spin%SO
      call add_DM_so(DM, DMfact,EDM, EDMfact,no1,no2,GF,N_Elec,Elecs,&
          DMidx,EDMidx,spS,q,is_eq)
    end if
  end subroutine add_DM

  subroutine add_DM_ncol(DM, DMfact,EDM, EDMfact, &
      no1,no2,GF, &
      N_Elec,Elecs, &
      DMidx, EDMidx, &
      spS, q, &
      is_eq)

      use intrinsic_missing, only: SFIND

      use class_Sparsity
      use class_dSpData1D
      use class_dSpData3D
      use ts_electrode_m

#ifdef TS_SOC_DEBUG
      use units, only : eV
#endif

      ! The DM and EDM equivalent matrices
      type(dSpData3D), intent(inout) :: DM
      complex(dp), intent(inout) :: DMfact
      type(dSpData3D), intent(inout) :: EDM
      complex(dp), intent(inout) :: EDMfact
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
      type(dSpData1D), intent(in), optional :: spS
      !< Charge calculated at this energy-point
      !!
      !! This does not contain the additional factor 1/Pi
      real(dp), intent(inout), optional :: q
      logical, intent(in), optional :: is_eq

      ! Arrays needed for looping the sparsity
      type(Sparsity), pointer :: s
      integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
      integer, pointer :: s_ncol(:), s_ptr(:), s_col(:), sp_col(:)
      real(dp), pointer :: D(:,:,:), E(:,:,:), Sg(:)
      integer :: io, ind, nr
      integer :: sp, sind
      integer :: iu, iuT, ju, juT, i1, i2
      logical :: lis_eq, hasEDM, calc_q

      ! Conjugate of weighting factors
      complex(dp) :: GG

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
            Sg => val(spS)
            call attach(s, n_col=s_ncol, list_ptr=s_ptr, list_col=s_col)
         end if

         if ( calc_q .and. hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ind,ju,iuT,juT,GG,sp,sp_col,sind), &
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

                     ! Search for overlap index
                     sind = sp + SFIND(sp_col, l_col(ind))
                     if ( sp < sind ) then
                        ! Charge is given as Tr[(G-G^\dagger).S]
                        q = q - aimag(Gf(1,iu,1,ju) - conjg(Gf(1,juT,1,iuT))) * Sg(sind)
                        q = q - aimag(Gf(2,iu,2,ju) - conjg(Gf(2,juT,2,iuT))) * Sg(sind)
                     end if

                     !   In Gf the spin-boxes are expanded. We have to undo this here
                     !   when, copying into the DM array (D).
                     !
                     !        | D(1,:)              D(3,:) - iD(4,:) |
                     !   DM = |                                      |
                     !        | D(3,:) + iD(4,:)              D(2,:) |
                     !
                     !         | -Im(GG11)                -Im(GG12) + i Re(GG12) |
                     !       = |                                                 |
                     !         | -Im(GG12) - i Re(GG12)   -Im(GG22)              |
                     !
                     !   GGij =       Gf(i,iu ,j,ju ) * DMfact
                     !        - conjg(Gf(j,juT,i,iuT) * DMfact)
                     !
                     !   where iu and ju are the orbital indices corresponding to ind
                     !   and juT and iuT the transposed indices.
                     !
                     !   Currently we do not average the two off diagonal spin-box
                     !   elements G12 and G21.
                     !
                     !   The same for the energy density matrix.

                     GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
                     D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
                     GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
                     D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
                     GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
                     D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
                     D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)


                     GG = Gf(1,iu,1,ju) * EDMfact - conjg(Gf(1,juT,1,iuT) * EDMfact)
                     E(1,ind,i2) = E(1,ind,i2) - aimag(GG)
                     GG = Gf(2,iu,2,ju) * EDMfact - conjg(Gf(2,juT,2,iuT) * EDMfact)
                     E(2,ind,i2) = E(2,ind,i2) - aimag(GG)
                     GG = Gf(1,iu,2,ju) * EDMfact - conjg(Gf(2,juT,1,iuT) * EDMfact)
                     E(3,ind,i2) = E(3,ind,i2) - aimag(GG)
                     E(4,ind,i2) = E(4,ind,i2) - real(GG, dp)
                  end do
               end if
            end do
!$OMP end parallel do

         else if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ind,ju,iuT,juT,GG), &
            do io = 1 , nr
               ! Quickly go past the buffer atoms...
               if ( l_ncol(io) /= 0 ) then

                  ! The update region equivalent GF part
                  iu = io - orb_offset(io)
                  iuT = iu - offset(N_Elec,Elecs,io)

                  do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

                     juT = l_col(ind) - orb_offset(l_col(ind))
                     ju = juT - offset(N_Elec,Elecs,l_col(ind))

                     !   In Gf the spin-boxes are expanded. We have to undo this here
                     !   when, copying into the DM array (D).
                     !
                     !        | D(1,:)              D(3,:) - iD(4,:) |
                     !   DM = |                                      |
                     !        | D(3,:) + iD(4,:)              D(2,:) |
                     !
                     !        | -Im(GG11)                -Im(GG12) + i RE(GG12) |
                     !      = |                                                 |
                     !        | -Im(GG12) - i Re(GG12)   -Im(GG22)              |
                     !
                     !   GGij =       Gf(i,iu ,j,ju ) * DMfact
                     !        - conjg(Gf(j,juT,i,iuT) * DMfact)
                     !
                     !   where iu and ju are the orbital indices corresponding to ind
                     !   and juT and iuT the transposed indices.
                     !
                     !   Currently we do not average the two off diagonal spin-box
                     !   elements G12 and G21.
                     !
                     !   The same for the energy density matrix.

                     GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
                     D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
                     D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
                     GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
                     D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
                     D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
                     GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
                     D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
                     D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
                     GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
                     D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
                     D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)

                     GG = Gf(1,iu,1,ju) * EDMfact - conjg(Gf(1,juT,1,iuT) * EDMfact)
                     E(1,ind,i2) = E(1,ind,i2) - aimag(GG)
                     GG = Gf(2,iu,2,ju) * EDMfact - conjg(Gf(2,juT,2,iuT) * EDMfact)
                     E(2,ind,i2) = E(2,ind,i2) - aimag(GG)
                     GG = Gf(1,iu,2,ju) * EDMfact - conjg(Gf(2,juT,1,iuT) * EDMfact)
                     E(3,ind,i2) = E(3,ind,i2) - aimag(GG)
                     E(4,ind,i2) = E(4,ind,i2) - real(GG, dp)
                  end do
               end if
            end do
!$OMP end parallel do

         else if ( calc_q ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ind,ju,iuT,juT,GG,sp,sp_col,sind), &
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

                     ! Search for overlap index
                     sind = sp + SFIND(sp_col, l_col(ind))
                     if ( sp < sind ) then
                        ! Charge is given as Tr[(G-G^\dagger).S]
                        q = q - aimag(Gf(1,iu,1,ju) - conjg(Gf(1,juT,1,iuT))) * Sg(sind)
                        q = q - aimag(Gf(2,iu,2,ju) - conjg(Gf(2,juT,2,iuT))) * Sg(sind)
                     end if

                     !   In Gf the spin-boxes are expanded. We have to undo this here
                     !   when, copying into the DM array (D).
                     !
                     !        | D(1,:)              D(3,:) - iD(4,:) |
                     !   DM = |                                      |
                     !        | D(3,:) + iD(4,:)              D(2,:) |
                     !
                     !        | -Im(GG11)                -Im(GG12) + i RE(GG12) |
                     !      = |                                                 |
                     !        | -Im(GG12) - i Re(GG12)   -Im(GG22)              |
                     !
                     !   GGij =       Gf(i,iu ,j,ju ) * DMfact
                     !        - conjg(Gf(j,juT,i,iuT) * DMfact)
                     !
                     !   where iu and ju are the orbital indices corresponding to ind
                     !   and juT and iuT the transposed indices.
                     !
                     !   Currently we do not average the two off diagonal spin-box
                     !   elements G12 and G21.

                     GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
                     D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
                     D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
                     GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
                     D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
                     D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
                     GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
                     D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
                     D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
                     GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
                     D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
                     D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)
                  end do
               end if
            end do
!$OMP end parallel do

         else

!$OMP parallel do default(shared), &
!$OMP&private(io,iu,ind,ju,juT,GG)
            do io = 1 , nr
               ! Quickly go past the buffer atoms...
               if ( l_ncol(io) /= 0 ) then

                  ! The update region equivalent GF part
                  iu = io - orb_offset(io)
                  iuT = iu - offset(N_Elec,Elecs,io)

                  do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

                     juT = l_col(ind) - orb_offset(l_col(ind))
                     ju = juT - offset(N_Elec,Elecs,l_col(ind))
                     !   In Gf the spin-boxes are expanded. We have to undo this here
                     !   when, copying into the DM array (D).
                     !
                     !        | D(1,:)              D(3,:) - iD(4,:) |
                     !   DM = |                                      |
                     !        | D(3,:) + iD(4,:)              D(2,:) |
                     !
                     !        | -Im(GG11)                -Im(GG12) + i Re(GG12) |
                     !      = |                                                 |
                     !        | -Im(GG12) - i Re(GG12)   -Im(GG22)              |
                     !
                     !   GGij =       Gf(i,iu ,j,ju ) * DMfact
                     !        - conjg(Gf(j,juT,i,iuT) * DMfact)
                     !
                     !   where iu and ju are the orbital indices corresponding to ind
                     !   and juT and iuT the transposed indices.
                     !
                     !   Currently we do not average the two off diagonal spin-box
                     !   elements G12 and G21.

                     GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
                     D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
                     D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
                     GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
                     D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
                     D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
                     GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
                     D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
                     D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
                     GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
                     D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
                     D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)

                  end do
               end if
            end do
!$OMP end parallel do

         end if

      else ! lis_eq

         if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(io,iu,ind,ju)
            do io = 1 , nr
               if ( l_ncol(io) /= 0 ) then
                  iu = io - orb_offset(io)
                  do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                     ju = l_col(ind) - orb_offset(l_col(ind))

                     !   In non-equilibrium GF contains the product GF.Gamma.GF^Dag
                     !   In Gf the spin-boxes are expanded. We have to undo this here
                     !   when, copying into the DM array (D).
                     !
                     !        | D(1,:)              D(3,:) - iD(4,:) |
                     !   DM = |                                      |
                     !        | D(3,:) + iD(4,:)              D(2,:) |
                     !
                     !        | Re(GG11)                Re(GG12) + i Im(GG12) |
                     !      = |                                               |
                     !        | Re(GG12) - i Im(GG12)   Re(GG22)              |
                     !
                     !   GGij = Gf(i,iu ,j,ju ) * DMfact
                     !
                     !   where iu and ju are the orbital indices corresponding to ind
                     !   and juT and iuT the transposed indices.
                     !
                     !   Currently we do not average the two off diagonal spin-box
                     !   elements G12 and G21.
                     !
                     !   The same for the energy density matrix

                     D(1,ind,i1) = D(1,ind,i1) +  real( GF(1,iu,1,ju) * DMfact  ,dp)
                     D(2,ind,i1) = D(2,ind,i1) +  real( GF(2,iu,2,ju) * DMfact  ,dp)
                     D(3,ind,i1) = D(3,ind,i1) +  real( GF(1,iu,2,ju) * DMfact  ,dp)
                     D(4,ind,i1) = D(4,ind,i1) - aimag( GF(1,iu,2,ju) * DMfact)

                     E(1,ind,i2) = E(1,ind,i2) +  real( GF(1,iu,1,ju) * EDMfact ,dp)
                     E(2,ind,i2) = E(2,ind,i2) +  real( GF(2,iu,2,ju) * EDMfact ,dp)
                     E(3,ind,i2) = E(3,ind,i2) +  real( GF(1,iu,2,ju) * EDMfact ,dp)
                     E(4,ind,i2) = E(4,ind,i2) - aimag( GF(1,iu,2,ju) * EDMfact)
                  end do
               end if
            end do
!$OMP end parallel do

         else

!$OMP parallel do default(shared), &
!$OMP&private(io,iu,ind,ju)
            do io = 1 , nr
               if ( l_ncol(io) /= 0 ) then
                  iu = io - orb_offset(io)
                  do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                     ju = l_col(ind) - orb_offset(l_col(ind))

                     !   In non-equilibrium GF contains the product GF.Gamma.GF^Dag
                     !   In Gf the spin-boxes are expanded. We have to undo this here
                     !   when, copying into the DM array (D).
                     !
                     !        | D(1,:)              D(3,:) - iD(4,:) |
                     !   DM = |                                      |
                     !        | D(3,:) + iD(4,:)              D(2,:) |
                     !
                     !        | Re(GG11)                Re(GG12) + i Im(GG12) |
                     !      = |                                               |
                     !        | Re(GG12) - i Im(GG12)   Re(GG22)              |
                     !
                     !   GGij = Gf(i,iu ,j,ju )
                     !
                     !
                     !   where iu and ju are the orbital indices corresponding to ind
                     !   and juT and iuT the transposed indices.
                     !
                     !   Currently we do not average the two off diagonal spin-box
                     !   elements G12 and G21.

                     D(1,ind,i1) = D(1,ind,i1) +  real( GF(1,iu,1,ju) * DMfact  ,dp)
                     D(2,ind,i1) = D(2,ind,i1) +  real( GF(2,iu,2,ju) * DMfact  ,dp)
                     D(3,ind,i1) = D(3,ind,i1) +  real( GF(1,iu,2,ju) * DMfact  ,dp)
                     D(4,ind,i1) = D(4,ind,i1) - aimag( GF(1,iu,2,ju) * DMfact)
                  end do
               end if
            end do
!$OMP end parallel do

         end if
      end if

   contains

      pure function offset(N_Elec,Elecs,io)
         integer, intent(in) :: N_Elec
         type(electrode_t), intent(in) :: Elecs(N_Elec)
         integer, intent(in) :: io
         integer :: offset
         offset = sum(Elecs(:)%device_orbitals(), &
            MASK=(Elecs(:)%DM_update==0) .and. Elecs(:)%idx_o <= io )
      end function offset

   end subroutine add_DM_ncol

  subroutine add_DM_so(DM, DMfact,EDM, EDMfact, &
      no1,no2,GF, &
      N_Elec,Elecs, &
      DMidx, EDMidx, &
      spS, q, &
      is_eq)

      use intrinsic_missing, only: SFIND  

   use class_Sparsity
    use class_dSpData1D
   use class_dSpData3D
   use ts_electrode_m

#ifdef TS_SOC_DEBUG
   use units, only : eV
#endif

   ! The DM and EDM equivalent matrices
   type(dSpData3D), intent(inout) :: DM
   complex(dp), intent(inout) :: DMfact
   type(dSpData3D), intent(inout) :: EDM
   complex(dp), intent(inout) :: EDMfact
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
    type(dSpData1D), intent(in), optional :: spS
    !< Charge calculated at this energy-point
    !!
    !! This does not contain the additional factor 1/Pi
    real(dp), intent(inout), optional :: q
   logical, intent(in), optional :: is_eq

   ! Arrays needed for looping the sparsity
   type(Sparsity), pointer :: s
   integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: s_ncol(:), s_ptr(:), s_col(:), sp_col(:)
    real(dp), pointer :: D(:,:,:), E(:,:,:), Sg(:)
   integer :: io, ind, nr
    integer :: sp, sind
   integer :: iu, iuT, ju, juT, i1, i2
    logical :: lis_eq, hasEDM, calc_q

   ! Conjugate of weighting factors
   complex(dp) :: GG

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
        Sg => val(spS)
        call attach(s, n_col=s_ncol, list_ptr=s_ptr, list_col=s_col)
      end if

      if ( calc_q .and. hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ind,ju,iuT,juT,GG,sp,sp_col,sind), &
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

              ! Search for overlap index
              sind = sp + SFIND(sp_col, l_col(ind))
              if ( sp < sind ) then
                ! Charge is given as Tr[(G-G^\dagger).S]
                q = q - aimag(Gf(1,iu,1,ju) - conjg(Gf(1,juT,1,iuT))) * Sg(sind)
                q = q - aimag(Gf(2,iu,2,ju) - conjg(Gf(2,juT,2,iuT))) * Sg(sind)
              end if

              !   In Gf the spin-boxes are expanded. We have to undo this here
              !   when, copying into the DM array (D).
              !
              !        | D(1,:) + iD(5,:)    D(3,:) - iD(4,:) | 
              !   DM = |                                      |
              !        | D(7,:) + iD(8,:)    D(2,:) + iD(6,:) |
              !
              !         | -Im(GG11) + i Re(GG11)   -Im(GG12) + i Re(GG12) | 
              !       = |                                                 | 
              !         | -Im(GG21) + i Re(GG21)   -Im(GG22) + i Re(GG22) | 
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)
              !
              !   where iu and ju are the orbital indices corresponding to ind
              !   and juT and iuT the transposed indices. 
              !
              !   Similar for the energy density matrix, but we discard the
              !   imaginary parts on the spin-box diagonal and one of 
              !   the off-diagonal spin-box elements:
              !
              !         | E(1,:)              E(3,:) - iE(4,:) | 
              !   EDM = |                                      |
              !         | E(3,:) + iE(4,:)    E(2,:)           |
              !
              !         | -Im(GG11) + i Re(GG11)   -Im(GG12) + i Re(GG12) | 
              !       = |                                                 | 
              !         | -Im(GG12) - i Re(GG12)   -Im(GG22) + i Re(GG22) | 
              !
              !   GGij =       Gf(i,iu ,j,ju ) * EDMfact
              !        - conjg(Gf(j,juT,i,iuT) * EDMfact)


              GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
              D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
              D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
              GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
              D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
              D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
              GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
              D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
              D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
              GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
              D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
              D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)

              GG = Gf(1,iu,1,ju) * EDMfact - conjg(Gf(1,juT,1,iuT) * EDMfact)
              E(1,ind,i2) = E(1,ind,i2) - aimag(GG)
              GG = Gf(2,iu,2,ju) * EDMfact - conjg(Gf(2,juT,2,iuT) * EDMfact)
              E(2,ind,i2) = E(2,ind,i2) - aimag(GG)
              GG = Gf(1,iu,2,ju) * EDMfact - conjg(Gf(2,juT,1,iuT) * EDMfact)
              E(3,ind,i2) = E(3,ind,i2) - aimag(GG)
              E(4,ind,i2) = E(4,ind,i2) - real(GG, dp)      
            end do
          end if
        end do
!$OMP end parallel do

      else if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ind,ju,iuT,juT,GG), &
       do io = 1 , nr
         ! Quickly go past the buffer atoms...
         if ( l_ncol(io) /= 0 ) then

           ! The update region equivalent GF part
           iu = io - orb_offset(io)
           iuT = iu - offset(N_Elec,Elecs,io)

           do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

             juT = l_col(ind) - orb_offset(l_col(ind))
             ju = juT - offset(N_Elec,Elecs,l_col(ind))

              !   In Gf the spin-boxes are expanded. We have to undo this here
              !   when, copying into the DM array (D).
              !
              !        | D(1,:) + iD(5,:)    D(3,:) - iD(4,:) |
              !   DM = |                                      |
              !        | D(7,:) + iD(8,:)    D(2,:) + iD(6,:) |
              !
              !         | -Im(GG11) + i Re(GG11)   -Im(GG21) + i Re(GG21) | 
              !       = |                                                 | 
              !         | -Im(GG21) - i Re(GG21)   -Im(GG22) + i Re(GG22) |  
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)
              !
              !   where iu and ju are the orbital indices corresponding to ind
              !   and juT and iuT the transposed indices. 
              !
              !   Similar for the energy density matrix, but we discard the
              !   imaginary parts on the spin-box diagonal and one of 
              !   the off-diagonal spin-box elements:
              !
              !         | E(1,:)              E(3,:) - iE(4,:) | 
              !   EDM = |                                      |
              !         | E(3,:) + iE(4,:)    E(2,:)           |
              !
              !         | -Im(GG11) + i Re(GG11)   -Im(GG12) + i Re(GG12) | 
              !       = |                                                 | 
              !         | -Im(GG12) - i Re(GG12)   -Im(GG22) + i Re(GG22) | 
              !
              !   GGij =       Gf(i,iu ,j,ju ) * EDMfact
              !        - conjg(Gf(j,juT,i,iuT) * EDMfact)

             GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
             D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
             D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
             GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
             D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
             D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
             GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
             D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
             D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
             GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
             D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
             D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)
             
             GG = Gf(1,iu,1,ju) * EDMfact - conjg(Gf(1,juT,1,iuT) * EDMfact)
             E(1,ind,i2) = E(1,ind,i2) - aimag(GG)
             GG = Gf(2,iu,2,ju) * EDMfact - conjg(Gf(2,juT,2,iuT) * EDMfact)
             E(2,ind,i2) = E(2,ind,i2) - aimag(GG)
             GG = Gf(1,iu,2,ju) * EDMfact - conjg(Gf(2,juT,1,iuT) * EDMfact)
             E(3,ind,i2) = E(3,ind,i2) - aimag(GG)
             E(4,ind,i2) = E(4,ind,i2) - real(GG, dp)
           end do
          end if
        end do
!$OMP end parallel do
      
      else if ( calc_q ) then

!$OMP parallel do default(shared), &
!$OMP&  private(io,iu,ind,ju,iuT,juT,GG,sp,sp_col,sind), &
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

              ! Search for overlap index
              sind = sp + SFIND(sp_col, l_col(ind))
              if ( sp < sind ) then
                ! Charge is given as Tr[(G-G^\dagger).S]
                q = q - aimag(Gf(1,iu,1,ju) - conjg(Gf(1,juT,1,iuT))) * Sg(sind)
                q = q - aimag(Gf(2,iu,2,ju) - conjg(Gf(2,juT,2,iuT))) * Sg(sind)
              end if

              !   In Gf the spin-boxes are expanded. We have to undo this here
              !   when, copying into the DM array (D).
              !
              !        | D(1,:) + iD(5,:)    D(3,:) - iD(4,:) |
              !   DM = |                                      |
              !        | D(7,:) + iD(8,:)    D(2,:) + iD(6,:) |
              !
              !        | -Im(GG11) + i Re(GG11)   -Im(GG21) + i Re(GG21) | 
              !      = |                                                 | 
              !        | -Im(GG21) - i Re(GG21)   -Im(GG22) + i Re(GG22) |  
              !
              !   GGij =       Gf(i,iu ,j,ju ) * DMfact
              !        - conjg(Gf(j,juT,i,iuT) * DMfact)
              !
              !   where iu and ju are the orbital indices corresponding to ind
              !   and juT and iuT the transposed indices. 

              GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
              D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
              D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
              GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
              D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
              D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
              GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
              D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
              D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
              GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
              D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
              D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)
            end do
         end if
       end do
!$OMP end parallel do

     else

!$OMP parallel do default(shared), &
!$OMP&private(io,iu,ind,ju,juT,GG)
       do io = 1 , nr
         ! Quickly go past the buffer atoms...
         if ( l_ncol(io) /= 0 ) then

           ! The update region equivalent GF part
           iu = io - orb_offset(io)
           iuT = iu - offset(N_Elec,Elecs,io)

           do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

             juT = l_col(ind) - orb_offset(l_col(ind))
             ju = juT - offset(N_Elec,Elecs,l_col(ind))

             !   In Gf the spin-boxes are expanded. We have to undo this here
             !   when, copying into the DM array (D).
             !
             !        | D(1,:) + iD(5,:)    D(3,:) - iD(4,:) |
             !   DM = |                                      |
             !        | D(7,:) + iD(8,:)    D(2,:) + iD(6,:) |
             !
             !         | -Im(GG11) + i Re(GG11)   -Im(GG12) + i Re(GG12) | 
             !       = |                                                 | 
             !         | -Im(GG21) + i Re(GG21)   -Im(GG22) + i Re(GG22) |  
             !
             !   i
             !
             !   GGij = Gf(i,iu ,j,ju ) * DMfact
             !        - conjg(Gf(j,juT,i,iuT) * DMfact)
             !
             !   where iu and ju are the orbital indices corresponding to ind
             !   and juT and iuT the transposed indices. 

             GG = Gf(1,iu,1,ju) * DMfact - conjg(Gf(1,juT,1,iuT) * DMfact)
             D(1,ind,i1) = D(1,ind,i1) - aimag(GG)
             D(5,ind,i1) = D(5,ind,i1) + real(GG, dp)
             GG = Gf(2,iu,2,ju) * DMfact - conjg(Gf(2,juT,2,iuT) * DMfact)
             D(2,ind,i1) = D(2,ind,i1) - aimag(GG)
             D(6,ind,i1) = D(6,ind,i1) + real(GG, dp)
             GG = Gf(1,iu,2,ju) * DMfact - conjg(Gf(2,juT,1,iuT) * DMfact)
             D(3,ind,i1) = D(3,ind,i1) - aimag(GG)
             D(4,ind,i1) = D(4,ind,i1) - real(GG, dp)
             GG = Gf(2,iu,1,ju) * DMfact - conjg(Gf(1,juT,2,iuT) * DMfact)
             D(7,ind,i1) = D(7,ind,i1) - aimag(GG)
             D(8,ind,i1) = D(8,ind,i1) + real(GG, dp)

           end do
         end if
       end do
!$OMP end parallel do

     end if

   else ! lis_eq

     if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(io,iu,ind,ju)
       do io = 1 , nr
         if ( l_ncol(io) /= 0 ) then
           iu = io - orb_offset(io)
           do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             ju = l_col(ind) - orb_offset(l_col(ind))

             !   In non-equilibrium GF contains the product GF.Gamma.GF^Dag
             !   In Gf the spin-boxes are expanded. We have to undo this here
             !   when, copying into the DM array (D).
             !
             !        | D(1,:) + iD(5,:)    D(3,:) - iD(4,:) |
             !   DM = |                                      |
             !        | D(7,:) + iD(8,:)    D(2,:) + iD(6,:) |
             !
             !        | Re(GG11) + i Im(GG11)   Re(GG12) + i Im(GG12) | 
             !      = |                                               | 
             !        | Re(GG21) + i Im(GG21)   Re(GG22) + i Im(GG22) |  
             !
             !   GGij = Gf(i,iu ,j,ju ) * DMfact
             !
             !   where iu and ju are the orbital indices corresponding to ind
             !   and juT and iuT the transposed indices. 
             !
             !   Similar for the energy density matrix, but we discard the
             !   imaginary parts on the spin-box diagonal and one of 
             !   the off-diagonal spin-box elements:
             !
             !         | E(1,:)              E(3,:) - iE(4,:) | 
             !   EDM = |                                      |
             !         | E(3,:) + iE(4,:)    E(2,:)           |
             !
             !         | -Im(GG11) + i Re(GG11)   -Im(GG12) + i Re(GG12) | 
             !       = |                                                 | 
             !         | -Im(GG12) - i Re(GG12)   -Im(GG22) + i Re(GG22) | 
             !
             !   GGij =       Gf(i,iu ,j,ju ) * EDMfact
             !        - conjg(Gf(j,juT,i,iuT) * EDMfact)

             D(1,ind,i1) = D(1,ind,i1) +  real( GF(1,iu,1,ju) * DMfact  ,dp)
             D(2,ind,i1) = D(2,ind,i1) +  real( GF(2,iu,2,ju) * DMfact  ,dp)
             D(3,ind,i1) = D(3,ind,i1) +  real( GF(1,iu,2,ju) * DMfact  ,dp)
             D(4,ind,i1) = D(4,ind,i1) - aimag( GF(1,iu,2,ju) * DMfact)
             D(5,ind,i1) = D(5,ind,i1) + aimag( GF(1,iu,1,ju) * DMfact)
             D(6,ind,i1) = D(6,ind,i1) + aimag( GF(2,iu,2,ju) * DMfact)
             D(7,ind,i1) = D(7,ind,i1) +  real( GF(2,iu,1,ju) * DMfact  ,dp)
             D(8,ind,i1) = D(8,ind,i1) + aimag( GF(2,iu,1,ju) * DMfact)

             E(1,ind,i2) = E(1,ind,i2) +  real( GF(1,iu,1,ju) * EDMfact ,dp)
             E(2,ind,i2) = E(2,ind,i2) +  real( GF(2,iu,2,ju) * EDMfact ,dp)
             E(3,ind,i2) = E(3,ind,i2) +  real( GF(1,iu,2,ju) * EDMfact ,dp)
             E(4,ind,i2) = E(4,ind,i2) - aimag( GF(1,iu,2,ju) * EDMfact)
           end do
         end if
       end do
!$OMP end parallel do

     else

!$OMP parallel do default(shared), &
!$OMP&private(io,iu,ind,ju)
       do io = 1 , nr
         if ( l_ncol(io) /= 0 ) then
           iu = io - orb_offset(io)
           do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             ju = l_col(ind) - orb_offset(l_col(ind))

             !   In non-equilibrium GF contains the product GF.Gamma.GF^Dag
             !   In Gf the spin-boxes are expanded. We have to undo this here
             !   when, copying into the DM array (D).
             !
             !        | D(1,:) + iD(5,:)    D(3,:) - iD(4,:) |
             !   DM = |                                      |
             !        | D(7,:) + iD(8,:)    D(2,:) + iD(6,:) |
             !
             !        | Re(GG11) + i Im(GG11)   Re(GG12) + i Im(GG12) | 
             !      = |                                               | 
             !        | Re(GG21) + i Im(GG21)   Re(GG22) + i Im(GG22) |  
             !
             !   GGij = Gf(i,iu ,j,ju ) * DMfact
             !
             !   where iu and ju are the orbital indices corresponding to ind
             !   and juT and iuT the transposed indices. 

             D(1,ind,i1) = D(1,ind,i1) +  real( GF(1,iu,1,ju) * DMfact  ,dp)
             D(2,ind,i1) = D(2,ind,i1) +  real( GF(2,iu,2,ju) * DMfact  ,dp)
             D(3,ind,i1) = D(3,ind,i1) +  real( GF(1,iu,2,ju) * DMfact  ,dp)
             D(4,ind,i1) = D(4,ind,i1) - aimag( GF(1,iu,2,ju) * DMfact)
             D(5,ind,i1) = D(5,ind,i1) + aimag( GF(1,iu,1,ju) * DMfact)
             D(6,ind,i1) = D(6,ind,i1) + aimag( GF(2,iu,2,ju) * DMfact)
             D(7,ind,i1) = D(7,ind,i1) +  real( GF(2,iu,1,ju) * DMfact  ,dp)
             D(8,ind,i1) = D(8,ind,i1) + aimag( GF(2,iu,1,ju) * DMfact)
           end do
         end if
       end do
!$OMP end parallel do

     end if
   end if

 contains

   pure function offset(N_Elec,Elecs,io)
     integer, intent(in) :: N_Elec
     type(electrode_t), intent(in) :: Elecs(N_Elec)
     integer, intent(in) :: io
     integer :: offset
     offset = sum(Elecs(:)%device_orbitals(), &
         MASK=(Elecs(:)%DM_update==0) .and. Elecs(:)%idx_o <= io )
      end function offset

   end subroutine add_DM_so

   ! creation of the GF^{-1}.
   ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF
   subroutine prepare_invGF(cE, no_u, GFinv, &
      N_Elec, Elecs, spH, spS)

   use class_dSpData1D
      use class_dSpData2D
      use class_Sparsity
      use m_spin, only: spin
      use ts_electrode_m
      use m_ts_cctype, only : ts_c_idx
#ifdef TS_DEV
      use parallel,only:ionode
#endif
      use m_ts_full_scat, only : insert_Self_Energies_NC

      ! the current energy point
      type(ts_c_idx), intent(in) :: cE
      ! Remark that we need the left buffer orbitals
      ! to calculate the actual orbital of the sparse matrices...
      integer, intent(in) :: no_u
      complex(dp), intent(out) :: GFinv(2,no_u,2,no_u)
      integer, intent(in) :: N_Elec
      type(electrode_t), intent(in) :: Elecs(N_Elec)
      ! The Hamiltonian and overlap sparse matrices
      type(dSpData1D), intent(inout) :: spS
      type(dSpData2D), intent(inout) :: spH

      ! Local variables
      complex(dp) :: Z
      type(Sparsity), pointer :: sp
      integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
      real(dp), pointer :: H(:,:), S(:)
      integer :: io, iu, ju, ind, nr

#ifdef TS_DEV
      logical, save :: hasSaved = .false.
      integer :: i
#endif

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

      GFinv(:,:,:,:) = cmplx(0._dp,0._dp, dp)

!$OMP parallel default(shared), private(io,iu,ju,ind)

      ! We will only loop in the central region
      ! We have constructed the sparse array to only contain
      ! values in this part...
   if (spin%Ncol) then
!$OMP do
   do io = 1, nr
     if ( l_ncol(io) == 0 ) cycle

     iu = io - orb_offset(io)

     do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

       ju = l_col(ind) - orb_offset(l_col(ind))

       GFinv(1,ju,1,iu) = Z * S(ind) - cmplx(H(1,ind),    0._dp, dp)
       GFinv(2,ju,1,iu) =            - cmplx(H(3,ind), H(4,ind), dp)
       GFinv(1,ju,2,iu) =            - cmplx(H(3,ind),-H(4,ind), dp)
       GFinv(2,ju,2,iu) = Z * S(ind) - cmplx(H(2,ind),    0._dp, dp)
     end do
   end do
!$OMP end do
    else ! spin%SO
!$OMP do
   do io = 1, nr
     if ( l_ncol(io) == 0 ) cycle

     iu = io - orb_offset(io)

     do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

       ju = l_col(ind) - orb_offset(l_col(ind))

       GFinv(1,ju,1,iu) = Z * S(ind) - cmplx(H(1,ind), H(5,ind), dp)
       GFinv(2,ju,1,iu) =            - cmplx(H(7,ind), H(8,ind), dp)
       GFinv(1,ju,2,iu) =            - cmplx(H(3,ind),-H(4,ind), dp)
       GFinv(2,ju,2,iu) = Z * S(ind) - cmplx(H(2,ind), H(6,ind), dp)
     end do
   end do
!$OMP end do
   end if

   do io = 1 , N_Elec
      call insert_Self_Energies_NC(no_u, Gfinv, Elecs(io))
   end do

!$OMP end parallel

#ifdef TRANSIESTA_TIMING
   call timer('TS-prep',2)
#endif

 end subroutine prepare_invGF

end module m_ts_fullg_spinor
