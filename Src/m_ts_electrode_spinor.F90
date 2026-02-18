! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_ts_electrode_spinor
!
! Routines that are used for Electrodes GFs calculations
! Heavily updated by Nick Papior Andersen, 2012
!

  use precision, only : dp
  use m_ts_electrode, only: SSR_sGreen_DOS, SSR_sGreen_NoDOS, print_Elec_Green

  implicit none

  public :: create_Green_spinor
  public :: calc_next_GS_Elec_spinor

  private

  ! BLAS parameters
  complex(dp), parameter :: z_0 = cmplx(0._dp,0._dp,dp)
  complex(dp), parameter :: z_1 = cmplx(1._dp,0._dp,dp)
  complex(dp), parameter :: z_m1 = cmplx(-1._dp,0._dp,dp)

  interface set_HS_transfer
     module procedure set_HS_Transfer_1d
     module procedure set_HS_Transfer_4d
  end interface set_HS_transfer

contains

! ##################################################################
! ## Driver subroutine for calculating the (ideal)                ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                 Updated by : Nick Papior Andersen            ##
! ## It has now been parallelized to speed up electrode           ##
! ## surface Green function generation.                           ##
! ## It generates the surface Green function by handling          ##
! ## repetition as well.                                          ##
! ##################################################################
  subroutine create_Green_spinor(El, &
       ucell,nkpnt,kpoint,kweight, &
       NEn,ce, &
       DOS,T)

    use precision,  only : dp
    use parallel  , only : Node, Nodes, IONode
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Sum, MPI_Max, MPI_integer
    use mpi_siesta, only : MPI_Wait,MPI_Status_Size
    use mpi_siesta, only : MPI_double_complex
    use mpi_siesta, only : MPI_double_precision
#endif
    use ts_electrode_m
    use m_mat_invert

    use m_ts_elec_se, only : update_UC_expansion_A1D

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use m_iterator

! ***********************
! * INPUT variables     *
! ***********************
    type(electrode_t), intent(inout) :: El  ! The electrode
    integer,  intent(in)      :: nkpnt ! Number of k-points
    real(dp), intent(in)      :: kpoint(3,nkpnt) ! k-points
    real(dp), intent(in)      :: kweight(nkpnt) ! weights of kpoints
    real(dp), dimension(3,3)  :: ucell ! The unit cell of the CONTACT
    integer, intent(in)       :: NEn ! Number of energy points
    complex(dp), intent(in)   :: ce(NEn) ! the energy points

! ***********************
! * OUTPUT variables    *
! ***********************
    real(dp), intent(inout), optional :: DOS(El%no_u*2,NEn)
    real(dp), intent(inout), optional :: T(NEn)

! ***********************
! * LOCAL variables     *
! ***********************
    ! Array for holding converted k-points
    real(dp) :: bkpt(3), kpt(3), kq(3), wq, rcell(3,3)
    real(dp), allocatable :: lDOS(:)
    
    ! Dimensions
    integer :: nq, nspin, n_s
    integer :: nuo_E, nS, nuou_E, nuS, no_X, n_X

    ! Electrode transfer and hamiltonian matrix
    complex(dp), pointer :: H00(:) => null()
    complex(dp), pointer :: S00(:) => null()
    complex(dp), pointer :: H01(:) => null()
    complex(dp), pointer :: S01(:) => null()
    complex(dp), pointer :: zwork(:) => null()
    complex(dp), pointer :: zHS(:) => null()
    real(dp), allocatable :: sc_off(:,:)

    ! Expanded arrays
    complex(dp), pointer :: X(:) => null()

    ! Green function variables
    complex(dp), pointer :: GS(:)
    complex(dp), pointer :: Hq(:), Sq(:), Gq(:)
    complex(dp) :: ZEnergy

    ! In order to print information about the recursize algorithm
    integer, allocatable :: iters(:,:,:,:)
    real(dp) :: i_mean, i_std

    integer :: uGF
    ! Big loop counters
    type(itt1) :: it1
    integer, pointer :: ikpt
    integer :: iEn, iqpt
    ! Counters
    integer :: i, j, io, jo, off

    logical :: CalcDOS, CalcT, pre_expand
    logical :: is_left, Gq_allocated, reduce_size

#ifdef MPI
    integer :: MPIerror, curNode
    integer :: req, status(MPI_Status_Size)
    integer, allocatable :: reqs(:)
#endif

#ifdef TS_SOC_DEBUG
   character(len=1024) :: filename
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE create_Green_spinor' )
#endif

    call timer('TS_SE',1)

    CalcDOS = present(DOS)
    CalcT = present(T)

    ! Check input for what to do
    if( El%inf_dir == INF_NEGATIVE ) then
       is_left = .true.
    else if( El%inf_dir == INF_POSITIVE ) then
       is_left = .false.
    else
       call die("init electrode has received wrong job ID [L,R].")
    endif

    ! Initialize TSGF-file
    call init_TSGF()

    ! capture information from the electrode
    nspin  = El%nspin
    nuo_E  = 2 * El%no_u
    nS     = nuo_E ** 2
    nuou_E = 2 * El%no_used
    nuS    = nuou_E ** 2
    ! create expansion k-points (weight of q-points)
    nq     = El%Bloch%size()
    wq     = 1._dp / real(nq,dp)
    ! We also need to invert to get the contribution in the
    reduce_size = nuo_E /= nuou_E
    no_X = nuou_E * nq
    n_X  = no_X ** 2
    pre_expand = El%pre_expand > 0 .and. nq > 1

    ! Calculate offsets
    n_s = size(El%isc_off,dim=2)
    allocate(sc_off(3,n_s))
    do i = 1, n_s
      sc_off(1,i) = El%cell(1,1) * El%isc_off(1,i) + &
          El%cell(1,2) * El%isc_off(2,i) + &
          El%cell(1,3) * El%isc_off(3,i)
      sc_off(2,i) = El%cell(2,1) * El%isc_off(1,i) + &
          El%cell(2,2) * El%isc_off(2,i) + &
          El%cell(2,3) * El%isc_off(3,i)
      sc_off(3,i) = El%cell(3,1) * El%isc_off(1,i) + &
          El%cell(3,2) * El%isc_off(2,i) + &
          El%cell(3,3) * El%isc_off(3,i)
    end do

    ! Print information on file-size and electrode type.
    call print_Elec_Green(El, NEn, nkpnt)

    ! Initialize Green function and Hamiltonian arrays
    nullify(GS)
    if ( nS /= nuS ) then
       allocate(GS(nS))
       call memory('A','Z',nS,'create_green')
    !else
    !  the regions are of same size, so we can just point
    !  to the correct memory segment
    end if

    ! Allocate work array
    i = max(nS*9,nuS*nq*2)
    if ( pre_expand ) then
       i = max(i,n_X)
    end if
    allocate(zwork(i))
    call memory('A','Z',i,'create_green')

    ! Point the hamiltonian and the overlap to the work array
    ! The work-array is only used for calculation the surface
    ! Green function and
    Hq => zwork(1:nuS*nq)
    Sq => zwork(nuS*nq+1:nuS*nq*2)
    if ( size(zwork) >= nS * 9 + nuS*nq ) then
       Gq => zwork(nS*9+1:nS*9+nuS*nq)
       Gq_allocated = .false.
    else
       nullify(Gq)
       allocate(Gq(nuS*nq))
       call memory('A','Z',size(Gq),'create_green')
       Gq_allocated = .true.
    end if

    if ( pre_expand ) then
       ! We allocate space for pre-expansion of the arrays
       allocate(X(n_X))
    end if

    ! all the Hamiltonian and overlaps
    allocate(zHS(nS * nq * 4))
    call memory('A','Z',nS * nq * 4,'create_green')

    ! Prepare for the inversion
    i = max(no_X,nuo_E)
    call init_mat_inversion(i)

    ! Reset bulk DOS
    if ( CalcDOS ) then
       allocate(lDOS(nuo_E))
       DOS(:,:) = 0._dp
    end if
    if ( CalcT ) then
       T = 0._dp
    end if

!******************************************************************
!           Start Green function calculation
!******************************************************************
    
#ifdef MPI
    if ( IONode ) then
       allocate(reqs(Nodes-1))
       call memory('A','I',Nodes-1,'create_green')
       ! Create request handles for communication
       ! This is a rather new feature which enhances communication times.
       ! However, this is perhaps overkill as we never have VERY many 
       ! contour points. Say NEn > 1000
       ! Look in the loop for MPI_Start(...) for where this is used
       do i = 1 , Nodes - 1
          if ( pre_expand ) then
             call MPI_Recv_Init(X(1),n_X,MPI_double_complex, &
                  i,i,MPI_Comm_World,reqs(i),MPIerror)
          else
             call MPI_Recv_Init(Gq(1),nuS*nq,MPI_double_complex, &
                  i,i,MPI_Comm_World,reqs(i),MPIerror)
          end if
       end do
    else
       ! Create request handles for communication
       if ( pre_expand ) then
          call MPI_Send_Init(X(1),n_X,MPI_double_complex, &
               0,Node,MPI_Comm_World,req,MPIerror)
       else
          call MPI_Send_Init(Gq(1),nuS*nq,MPI_double_complex, &
               0,Node,MPI_Comm_World,req,MPIerror)
       end if
    end if
#endif

    ! prepare the iteration counter
    allocate(iters(nq,NEn,nkpnt,2))
    if ( IONode ) then
       ! TODO when adding new surface-Green functions schemes, please update here
       write(*,'(1x,a)') 'Lopez Sancho, Lopez Sancho & Rubio recursive &
            &surface self-energy calculation...'
    end if

    ! start up the iterators
    call itt_init  (it1,end=nkpnt)
    call itt_attach(it1,cur=ikpt)

    call reclat(El%cell,rcell,1)

    iters(:,:,:,:) = 0

    ! do k-point loop
    do while ( .not. itt_step(it1) )

       ! Init kpoint, in reciprocal vector units ( from CONTACT ucell)
       call El%kpoint_convert(ucell,kpoint(:,ikpt),bkpt, opt = 2)
       ! We need to save the k-point for the "expanded" super-cell
       El%bkpt_cur = bkpt
       
       ! loop over the repeated cell...
       HSq_loop: do iqpt = 1 , nq
             
          ! point to the correct segment of memory
          H00 => zHS((     iqpt-1)*nS+1:      iqpt *nS)
          S00 => zHS((  nq+iqpt-1)*nS+1:(  nq+iqpt)*nS)
          H01 => zHS((2*nq+iqpt-1)*nS+1:(2*nq+iqpt)*nS)
          S01 => zHS((3*nq+iqpt-1)*nS+1:(3*nq+iqpt)*nS)

          ! init qpoint in reciprocal lattice vectors
          kpt(:) = bkpt(:) + El%unfold_k(iqpt)
          ! Convert to 1/Bohr
          call kpoint_convert(rcell,kpt(1),kq(1),-2)

          ! Setup the transfer matrix and the intra cell at the k-point and q-point
          ! Calculate transfer matrices @Ef (including the chemical potential)
          call set_HS_transfer(El, n_s,sc_off, kq, &
              El%no_u, H00,S00,H01,S01)

          i = (iqpt-1)*nuS
          if ( reduce_size ) then
             if( is_left ) then
                ! Left, we use the last orbitals
                off = nuo_E - nuou_E + 1
                do jo = off - 1 , nuo_E - 1
                   do io = off , nuo_E
                      i = i + 1
                      Hq(i) = H00(nuo_E*jo+io)
                      Sq(i) = S00(nuo_E*jo+io)
                   end do
                end do
             else
                ! Right, the first orbitals
                do jo = 0 , nuou_E - 1
                   do io = 1 , nuou_E
                      i = i + 1
                      Hq(i) = H00(nuo_E*jo+io)
                      Sq(i) = S00(nuo_E*jo+io)
                   end do   ! io
                end do      ! jo
             end if
          end if

       end do HSq_loop
#ifdef TS_SOC_DEBUG
if (is_left) then
   write (filename, "('S00_L-',I0'.dat')") ikpt
   open(unit=340, file=filename)
   write (filename, "('H00_L-',I0'.dat')") ikpt
   open(unit=341, file=filename)
   write (filename, "('S01_L-',I0'.dat')") ikpt
   open(unit=342, file=filename)
   write (filename, "('H01_L-',I0'.dat')") ikpt
   open(unit=343, file=filename)
else
   write (filename, "('S00_R-',I0'.dat')") ikpt
   open(unit=340, file=filename)
   write (filename, "('H00_R-',I0'.dat')") ikpt
   open(unit=341, file=filename)
   write (filename, "('S01_R-',I0'.dat')") ikpt
   open(unit=342, file=filename)
   write (filename, "('H01_R-',I0'.dat')") ikpt
   open(unit=343, file=filename)
end if

do io = 1, nuou_E*nuou_E
   write(340, "('('spE28.15E4 spE28.15E4'j)')") S00(io)
   write(341, "('('spE28.15E4 spE28.15E4'j)')") H00(io) / eV
   write(342, "('('spE28.15E4 spE28.15E4'j)')") S01(io)
   write(343, "('('spE28.15E4 spE28.15E4'j)')") H01(io) / eV
end do
close (340)
close (341)
close (342)
close (343)
#endif
       
       ! Save Hamiltonian and overlap
       call store_HS()
  	    
       Econtour_loop: do iEn = 1, NEn

#ifdef MPI
          ! Every node takes one energy point
          ! This asserts that IONode = Node == 0 will have iEn == 1
          ! Important !
          curNode = MOD(iEn-1,Nodes)
          E_Nodes: if ( curNode == Node ) then
#endif
             ! as we already have shifted H,S to Ef + mu, and ZEnergy is
             ! wrt. mu, we don't need to subtract mu again
             ZEnergy = ce(iEn)
             i_mean = 0._dp
             
             ! loop over the repeated cell...
             q_loop: do iqpt = 1 , nq

                H00 => zHS((     iqpt-1)*nS+1:      iqpt *nS)
                S00 => zHS((  nq+iqpt-1)*nS+1:(  nq+iqpt)*nS)
                H01 => zHS((2*nq+iqpt-1)*nS+1:(2*nq+iqpt)*nS)
                S01 => zHS((3*nq+iqpt-1)*nS+1:(3*nq+iqpt)*nS)
                if ( nS == nuS ) then
                   ! instead of doing a copy afterward, we can
                   ! put it the correct place immediately
                   GS => Gq((    iqpt-1)*nS+1:      iqpt *nS)
                end if

                ! Calculate the surface Green function
                ! Zenergy is wrt. to the system Fermi-level
                if ( CalcDOS ) then
                   lDOS = 0._dp
                   call SSR_sGreen_DOS(nuo_E,ZEnergy,H00,S00,H01,S01, &
                        El%accu, GS, &
                        lDOS,i_mean,9*nS,zwork, &
                        iterations=iters(iqpt,iEn,ikpt,1), final_invert = reduce_size)
                   
                   ! We also average the k-points.
                   DOS(:,iEn) = DOS(:,iEn) + lDOS * wq * kweight(ikpt)
                   if ( CalcT ) T(iEn) = T(iEn) + i_mean

                else
                   call SSR_sGreen_NoDos(nuo_E,ZEnergy,H00,S00,H01,S01, &
                        El%accu, GS, &
                        8*nS,zwork, &
                        iterations=iters(iqpt,iEn,ikpt,1), final_invert = reduce_size)
                   
                end if
                  
#ifdef TS_SOC_DEBUG
if (is_left) then
   write (filename, "('SE_L-',I0,'-',I0,'.dat')") iEn, ikpt
   open(unit=315, file=filename)
else
   write (filename, "('SE_R-',I0,'-',I0,'.dat')") iEn, ikpt
   open(unit=315, file=filename)
end if
do io = 1, nuou_E*nuou_E
write(315, "('('spE28.15E4 spE28.15E4'j)')") GS(io) / eV
end do
close (315)
#endif

                ! Copy over surface Green function
                i = (iqpt-1)*nuS
                if ( reduce_size ) then
                   if ( is_left ) then
                      ! Left, we use the last orbitals
                      off = nuo_E - nuou_E + 1
                      do jo = off - 1 , nuo_E - 1
                         do io = off , nuo_E
                            i = i + 1
                            Gq(i) = GS(nuo_E*jo+io)
                         end do           ! io
                      end do              ! jo
                   else
                      ! Right, the first orbitals
                      do jo = 0 , nuou_E-1
                         do io = 1 , nuou_E
                            i = i + 1
                            Gq(i) = GS(nuo_E*jo+io)
                         end do           ! io
                      end do              ! jo
                   end if

                   if ( nq == 1 ) then
                      ! We invert back here, instead of in
                      ! the SCF (this is important as the
                      ! decreased size of the surface-Green function
                      ! would otherwise yield a different result)
                      call mat_invert(Gq(1:nuS),zwork(1:nuS),&
                           nuou_E, &
                           MI_IN_PLACE_LAPACK)

                   end if

                end if
                
             end do q_loop

             if ( pre_expand ) then
                ! Expand this energy-point
                call update_UC_expansion_A1D(nuou_E,no_X,El,nq, &
                     El%na_used,El%lasto_used,Gq,X)
                if ( reduce_size ) then
                   call mat_invert(X(1:n_X),zwork(1:n_X),&
                        no_X, &
                        MI_IN_PLACE_LAPACK)
                end if
             end if
                
#ifdef MPI
             ! If not IONode we should send message
             ! This message parsing is directly connected to 
             ! a predefined size of the message, see right before
             ! spin loop.
             ! It communicates the Gq array to the Gq array
             if ( .not. IONode ) then
                call MPI_Start(req,MPIerror)
                call MPI_Wait(req,status,MPIerror)
             end if

          end if E_Nodes

#endif

          ! Save the surface Green function file
          call store_GS()

       end do Econtour_loop

    end do
!*******************************************************************
!          Green function calculation is done
!*******************************************************************

#ifdef MPI
    call MPI_Reduce(iters(1,1,1,1), iters(1,1,1,2), nq*NEn*nkpnt, &
        MPI_Integer, MPI_Sum, 0, MPI_Comm_World, MPIerror)
!$OMP parallel default(shared), private(j,i,iqpt)
#else
!$OMP parallel default(shared), private(j,i,iqpt)

    iters(:,:,:,2) = iters(:,:,:,1)
#endif
    if ( IONode ) then
      i_mean = sum(iters(:,:,:,2)) / real(nq*NEn*nkpnt,dp)

!$OMP single
      i_std = 0._dp
!$OMP end single ! keep barrier

!$OMP do reduction(+:i_std)
      do j = 1 , nkpnt
        do i = 1 , NEn
          do iqpt = 1 , nq
            i_std = i_std + ( iters(iqpt,i,j,2) - i_mean ) ** 2
          end do
        end do
      end do
!$OMP end do

!$OMP master
      i_std = sqrt(i_std/real(NEn*nq*nkpnt,dp))
      ! TODO if new surface-Green function scheme is implemented, fix here
      write(*,'(1x,a,f10.4,'' / '',f10.4)') 'Lopez Sancho, Lopez Sancho & Rubio: &
          &Mean/std iterations: ', i_mean             , i_std
      write(*,'(1x,a,i10,'' / '',i10)')     'Lopez Sancho, Lopez Sancho & Rubio: &
          &Min/Max iterations : ', minval(iters(:,:,:,2)) , maxval(iters(:,:,:,2))
!$OMP end master
             
    end if
!$OMP end parallel

    deallocate(iters)

#ifdef MPI
    ! Free requests made for the communications
    if ( IONode ) then
       do i = 1 , Nodes - 1 
          call MPI_Request_Free(reqs(i),MPIerror)
       end do
       call memory('D','I',Nodes-1,'create_green')
       deallocate(reqs)
    else
       call MPI_Request_Free(req,MPIerror)
    end if
#endif

    ! Close and finish
    call finish_TSGF()
    
    ! Clean up computational arrays
    if ( nS /= nuS ) then
       call memory('D','Z',size(GS),'create_green')
       deallocate(GS)
    end if

    if ( Gq_allocated ) then
       call memory('D','Z',size(Gq),'create_green')
       deallocate(Gq)
    end if

    ! Work-arrays
    call memory('D','Z',size(zwork),'create_green')
    deallocate(zwork)
    call memory('D','Z',size(zHS),'create_green')
    deallocate(zHS)

    deallocate(sc_off)

    if ( CalcDOS ) deallocate(lDOS)

    if ( pre_expand ) deallocate(X)

    call itt_destroy(it1)

    call clear_mat_inversion()

#ifdef MPI
    if ( CalcDOS ) then
       ! Sum the bulkdensity of states
       ! Here we can safely use the array as temporary (Gq)
       ! TODO NW this should be changed to 4 DOS components eventually
       allocate(lDOS(nuo_E*NEn))
    else if ( CalcT ) then
       allocate(lDOS(NEn))
    end if
    if ( allocated(lDOS) ) then
       call memory('A','D',size(lDOS),'create_green')
    end if
    
    if ( CalcDOS ) then
       call MPI_AllReduce(DOS(1,1),lDOS(1),nuo_E*NEn, &
            MPI_double_precision, &
            MPI_Sum,MPI_Comm_World,MPIerror)
       i = 0
        do io = 1 , NEn
            DOS(1:nuo_E,io) = lDOS(i+1:i+nuo_E)
            i = i + nuo_E
        end do

    end if
    if ( CalcT ) then
       call MPI_AllReduce(T(1),lDOS(1),NEn, &
            MPI_double_precision, &
            MPI_Sum,MPI_Comm_World,MPIerror)
       T(1:NEn) = lDOS(1:NEn)

    end if

    if ( allocated(lDOS) ) then
       call memory('D','D',size(lDOS),'create_green')
       deallocate(lDOS)
    end if
#endif

    call timer('TS_SE',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS create_Green' )
#endif

  contains


    subroutine init_TSGF()
      
      real(dp), allocatable :: kE(:,:)

      if ( .not. IONode ) return
      
      call io_assign(uGF)
      open(FILE=El%GFfile,UNIT=uGF,FORM='UNFORMATTED')
      
      ! Electrode information
      write(uGF) El%nspin, El%cell
      write(uGF) El%na_u, El%no_u
      write(uGF) El%na_used, El%no_used
      write(uGF) El%xa_used, El%lasto_used
      write(uGF) El%repeat, El%Bloch%B(:), El%pre_expand
      write(uGF) El%mu%mu
      
      ! Write out explicit information about this content
      write(uGF) nkpnt
      ! Notice that we write the k-points for the ELECTRODE
      ! They will be stored in units of the reciprocal lattice vector
      !   1/b 
      allocate(kE(3,nkpnt))
      do i = 1 , nkpnt
         ! Store the k-points in units of reciprocal lattice
         call El%kpoint_convert(ucell,kpoint(:,i),kE(:,i), opt = 2)
      end do
      write(uGF) kE, kweight
      deallocate(kE)
      
      ! write out the contour information
      write(uGF) NEn
      write(uGF) ce ! energy points
      
    end subroutine init_TSGF

    subroutine store_HS()
      
      if ( .not. IONode ) return
      ! k-point and energy-point is in front of Hamiltonian

      write(uGF) ikpt, 1, ce(1) ! k-point and energy point
      if ( reduce_size ) then
         if ( pre_expand .and. El%pre_expand > 1 ) then
            call update_UC_expansion_A1D(nuou_E,no_X,El,nq,&
                 El%na_used,El%lasto_used,Hq,X)
            write(uGF) X
            call update_UC_expansion_A1D(nuou_E,no_X,El,nq,&
                 El%na_used,El%lasto_used,Sq,X)
            write(uGF) X
         else
            write(uGF) Hq
            write(uGF) Sq
         end if
      else
         H00 => zHS(      1:nq*nS  )
         S00 => zHS(nq*nS+1:nq*nS*2)
         if ( pre_expand .and. El%pre_expand > 1 ) then
            call update_UC_expansion_A1D(nuo_E,no_X,El,nq, &
                 El%na_used,El%lasto_used,H00,X)
            write(uGF) X
            call update_UC_expansion_A1D(nuo_E,no_X,El,nq,&
                 El%na_used,El%lasto_used,S00,X)
            write(uGF) X
         else
            write(uGF) H00
            write(uGF) S00
         end if
      end if

    end subroutine store_HS

    subroutine store_GS()
      
      ! Save Surface Green function
      if ( .not. IONode ) return

#ifdef MPI
      if ( curNode /= Node ) then
         call MPI_Start(reqs(curNode),MPIerror)
         call MPI_Wait(reqs(curNode),status,MPIerror)
      end if
#endif
         
      ! Write out calculated information at E point
      if ( iEn /= 1 ) write(uGF) ikpt, iEn, ce(iEn)
      if ( pre_expand ) then
         write(uGF) X
      else
         write(uGF) Gq
      end if
      
    end subroutine store_GS

    subroutine finish_TSGF()
      
      ! Close file
      if ( .not. IONode ) return
      
      call io_close(uGF)
      write(*,'(a,/)') "Done creating '"//trim(El%GFfile)//"'."  

    end subroutine finish_TSGF
    
  end subroutine create_Green_spinor

!**********
! Create the Hamiltonian for the electrode as well
! as creating the transfer matrix.
!**********
  subroutine set_HS_Transfer_1d(El,n_s,sc_off,kq, &
      no,Hk,Sk,Hk_T,Sk_T)

    use iso_c_binding

    use sys, only : die
    use precision, only : dp
    use ts_electrode_m
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in) :: no
    type(electrode_t), intent(inout) :: El
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    real(dp), intent(in) :: kq(3)   ! k + q-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(:), target, intent(inout) :: Hk,Sk,Hk_T,Sk_T

    complex(dp), dimension(:,:,:,:), pointer :: Hk4D,Sk4D,Hk_T4D,Sk_T4D

    call c_f_pointer(c_loc(Hk), Hk4D, [2, no, 2, no])
    call c_f_pointer(c_loc(Sk), Sk4D, [2, no, 2, no])
    call c_f_pointer(c_loc(Hk_T), Hk_T4D, [2, no, 2, no])
    call c_f_pointer(c_loc(Sk_T), Sk_T4D, [2, no, 2, no])

    call set_HS_Transfer_4d(El,n_s,sc_off,kq, &
        no,Hk4D,Sk4D,Hk_T4D,Sk_T4D)
    
  end subroutine set_HS_Transfer_1d
  
  subroutine set_HS_Transfer_4d(El,n_s,sc_off,kq, &
       no,Hk,Sk,Hk_T,Sk_T)
    use sys, only : die
    use precision, only : dp
    use ts_electrode_m
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in) :: no
    type(electrode_t), intent(inout) :: El
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    real(dp), intent(in) :: kq(3) ! k + q-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(:,:,:,:) :: Hk,Sk,Hk_T,Sk_T

    select case ( El%nspin )
    case ( 4 )
      call set_HS_Transfer_4d_NC(El,n_s,sc_off,kq, &
          no,Hk,Sk,Hk_T,Sk_T)
    case ( 8 )
      call set_HS_Transfer_4d_SO(El,n_s,sc_off,kq, &
          no,Hk,Sk,Hk_T,Sk_T)
    case default
      write(*,*) "Electrode spin:", El%nspin
      call die("set_hs_transfer_4d: Electrode calculation should non-collinear or spin-orbit")
    end select

  end subroutine set_HS_Transfer_4d

  subroutine set_HS_Transfer_4d_nc(El,n_s,sc_off,kq, &
      no,Hk,Sk,Hk_T,Sk_T)
    use sys, only : die
    use precision, only : dp
    use ts_electrode_m
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    type(electrode_t), intent(inout) :: El
    integer, intent(in) :: n_s, no
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    real(dp), intent(in) :: kq(3)   ! k + q-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(:,:,:,:) :: Hk,Sk,Hk_T,Sk_T

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: Ef
    complex(dp) :: ph(0:n_s-1)
    integer :: i, j, io, jo, ind, is
    integer, pointer :: ncol00(:), l_ptr00(:), l_col00(:)
    integer, pointer :: ncol01(:), l_ptr01(:), l_col01(:)
    real(dp), pointer :: H00(:,:) , S00(:), H01(:,:), S01(:)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE elec_HS_Transfer' )
#endif

    ! we need to subtract as the below code shifts to Ef
    Ef = El%Ef - El%mu%mu

    if ( El%no_u /= no ) call die('Wrong size of the electrode array')

    ! retrieve values
    call attach(El%sp00,n_col=ncol00,list_ptr=l_ptr00,list_col=l_col00)
    call attach(El%sp01,n_col=ncol01,list_ptr=l_ptr01,list_col=l_col01)
    ! point to the data-segments...
    H00 => val(El%H00)
    H01 => val(El%H01)
    S00 => val(El%S00)
    S01 => val(El%S01)

    ! The algorithm outside should take care of the
    ! nullification of the k-point in the semi-infinite direction
    ! While the Siesta convention matches xijg = R(j) - R(i)
    ! and this means ij -> (i,j) we cannot simply take the negative phase here.
    ! This is because the transfer matrix (H01) is not symmetric
    ! and therefore we *must* setup the Hamiltonian in the
    ! correct order (not relying on symmetries)
    do is = 0 , n_s - 1
       ph(is) = exp(cmplx(0._dp, -dot_product(kq,sc_off(:,is)),kind=dp))
    end do

    ! Initialize arrays
    Hk(:,:,:,:) = z_0
    Sk(:,:,:,:) = z_0
    Hk_T(:,:,:,:) = z_0
    Sk_T(:,:,:,:) = z_0

!$OMP parallel default(shared), private(io,jo,ind,is)

    ! We will not have any data-race condition here
!$OMP do 
    do io = 1 , no

       ! Create 00
       do ind = l_ptr00(io) + 1 , l_ptr00(io) + ncol00(io)
          jo = ucorb(l_col00(ind),no)
          is = (l_col00(ind)-1) / no
          
          Hk(1,jo,1,io) = Hk(1,jo,1,io) + (H00(ind,1) - Ef * S00(ind)) * ph(is)
          Hk(2,jo,1,io) = Hk(2,jo,1,io) + cmplx(H00(ind,3), H00(ind,4), dp) * ph(is)
          Hk(1,jo,2,io) = Hk(1,jo,2,io) + cmplx(H00(ind,3),-H00(ind,4), dp) * ph(is)
          Hk(2,jo,2,io) = Hk(2,jo,2,io) + (H00(ind,2) - Ef * S00(ind)) * ph(is)
          Sk(1,jo,1,io) = Sk(1,jo,1,io) + S00(ind) * ph(is)
          Sk(2,jo,2,io) = Sk(2,jo,2,io) + S00(ind) * ph(is)
       end do

       ! Create 01
       do ind = l_ptr01(io) + 1 , l_ptr01(io) + ncol01(io)
          jo = ucorb(l_col01(ind),no)
          is = (l_col01(ind)-1) / no

          Hk_T(1,jo,1,io) = Hk_T(1,jo,1,io) + (H01(ind,1) - Ef * S01(ind)) * ph(is)
          Hk_T(2,jo,1,io) = Hk_T(2,jo,1,io) + cmplx(H01(ind,3), H01(ind,4), dp) * ph(is)
          Hk_T(1,jo,2,io) = Hk_T(1,jo,2,io) + cmplx(H01(ind,3),-H01(ind,4), dp) * ph(is)
          Hk_T(2,jo,2,io) = Hk_T(2,jo,2,io) + (H01(ind,2) - Ef * S01(ind)) * ph(is)
          Sk_T(1,jo,1,io) = Sk_T(1,jo,1,io) + S01(ind) * ph(is)
          Sk_T(2,jo,2,io) = Sk_T(2,jo,2,io) + S01(ind) * ph(is)
          
       end do

    end do
!$OMP end do nowait

!$OMP end parallel

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS elec_HS_Transfer' )
#endif

  end subroutine set_HS_Transfer_4d_nc
  
  subroutine set_HS_Transfer_4d_so(El,n_s,sc_off,kq, &
       no,Hk,Sk,Hk_T,Sk_T)
    use sys, only : die
    use precision, only : dp
    use ts_electrode_m
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in) :: no
    type(electrode_t), intent(inout) :: El
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    real(dp), intent(in) :: kq(3) ! k + q-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(:,:,:,:) :: Hk, Sk, Hk_T, Sk_T

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: Ef
    complex(dp) :: ph(0:n_s-1)
    integer :: i, j, io, jo, ind, is
    integer, pointer :: ncol00(:), l_ptr00(:), l_col00(:)
    integer, pointer :: ncol01(:), l_ptr01(:), l_col01(:)
    real(dp), pointer :: H00(:,:) , S00(:), H01(:,:), S01(:)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE elec_HS_Transfer' )
#endif

    ! we need to subtract as the below code shifts to Ef
    Ef = El%Ef - El%mu%mu

    if ( El%no_u /= no ) call die('Wrong size of the electrode array')

    ! retrieve values
    call attach(El%sp00,n_col=ncol00,list_ptr=l_ptr00,list_col=l_col00)
    call attach(El%sp01,n_col=ncol01,list_ptr=l_ptr01,list_col=l_col01)
    ! point to the data-segments...
    H00 => val(El%H00)
    H01 => val(El%H01)
    S00 => val(El%S00)
    S01 => val(El%S01)

    ! The algorithm outside should take care of the
    ! nullification of the k-point in the semi-infinite direction
    ! While the Siesta convention matches xijg = R(j) - R(i)
    ! and this means ij -> (i,j) we cannot simply take the negative phase here.
    ! This is because the transfer matrix (H01) is not symmetric
    ! and therefore we *must* setup the Hamiltonian in the
    ! correct order (not relying on symmetries)
    do is = 0 , n_s - 1
       ph(is) = exp(cmplx(0._dp, -dot_product(kq,sc_off(:,is)),kind=dp))
    end do

    ! Initialize arrays
    Hk(:,:,:,:) = z_0
    Sk(:,:,:,:) = z_0
    Hk_T(:,:,:,:) = z_0
    Sk_T(:,:,:,:) = z_0


!$OMP parallel default(shared), private(io,jo,ind,is)

    ! We will not have any data-race condition here
    ! We create the transfer matrices in 2x2 blocks for each combination of
    ! orbitals ((up   up, up   down),
    !           (down up, down down))
!$OMP do 
    do io = 1 , no

       ! Create 00
       do ind = l_ptr00(io) + 1 , l_ptr00(io) + ncol00(io)
          jo = ucorb(l_col00(ind),no)
          is = (l_col00(ind)-1) / no
          
          Hk(1,jo,1,io) = Hk(1,jo,1,io) + cmplx(H00(ind,1)-Ef*S00(ind), H00(ind,5), dp) * ph(is)
          Hk(2,jo,1,io) = Hk(2,jo,1,io) + cmplx(H00(ind,7), H00(ind,8), dp) * ph(is)
          Hk(1,jo,2,io) = Hk(1,jo,2,io) + cmplx(H00(ind,3),-H00(ind,4), dp) * ph(is)
          Hk(2,jo,2,io) = Hk(2,jo,2,io) + cmplx(H00(ind,2)-Ef*S00(ind), H00(ind,6), dp) * ph(is)
          Sk(1,jo,1,io) = Sk(1,jo,1,io) + S00(ind) * ph(is)
          Sk(2,jo,2,io) = Sk(2,jo,2,io) + S00(ind) * ph(is)
       end do

       ! Create 01
       do ind = l_ptr01(io) + 1 , l_ptr01(io) + ncol01(io)
          jo = ucorb(l_col01(ind),no)
          is = (l_col01(ind)-1) / no

          Hk_T(1,jo,1,io) = Hk_T(1,jo,1,io) + cmplx(H01(ind,1)-Ef*S01(ind), H01(ind,5), dp) * ph(is)
          Hk_T(2,jo,1,io) = Hk_T(2,jo,1,io) + cmplx(H01(ind,7), H01(ind,8), dp) * ph(is)
          Hk_T(1,jo,2,io) = Hk_T(1,jo,2,io) + cmplx(H01(ind,3),-H01(ind,4), dp) * ph(is)
          Hk_T(2,jo,2,io) = Hk_T(2,jo,2,io) + cmplx(H01(ind,2)-Ef*S01(ind), H01(ind,6), dp) * ph(is)
          Sk_T(1,jo,1,io) = Sk_T(1,jo,1,io) + S01(ind) * ph(is)
          Sk_T(2,jo,2,io) = Sk_T(2,jo,2,io) + S01(ind) * ph(is)

       end do

    end do
!$OMP end do nowait

!$OMP end parallel

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS elec_HS_Transfer' )
#endif

  end subroutine set_HS_Transfer_4d_so

  subroutine calc_next_GS_Elec_spinor(El,bkpt,Z,nzwork,in_zwork,DOS,T)
    use iso_c_binding
    use precision,  only : dp

    use ts_electrode_m
    use m_mat_invert

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use alloc, only : re_alloc, de_alloc

    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    type(electrode_t), intent(inout) :: El
    ! the k-point in reciprocal units of the electrode
    ! also with / Rep
    real(dp), intent(in) :: bkpt(3)
    complex(dp), intent(in) :: Z
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: in_zwork(nzwork)
    ! Possibly the bulk density of states from the electrode
    ! If the DOS, also BULK transmission
    real(dp), intent(inout), optional :: DOS(:), T

    ! ***********************
    ! * LOCAL variables     *
    ! ***********************
    integer  :: iq
    real(dp) :: kpt(3), kq(3), rcell(3,3)
    
    ! Dimensions
    integer :: nq, nw
    integer :: nuo_E, nS, nuou_E, nuS, nuouT_E

    ! Electrode transfer and hamiltonian matrix
    complex(dp), pointer :: H00(:), H01(:), S00(:), S01(:)

    complex(dp), pointer :: zwork(:)
    complex(dp), pointer :: zHS(:) => null()
    real(dp), allocatable :: sc_off(:,:)

    ! Green function variables
    complex(dp), pointer :: GS(:)

    ! size requirement
    integer :: size_req(2)
    ! Counters
    integer :: i, ios, ioe, off, n_s
    logical :: is_left, reduce_size
    logical :: zHS_allocated
    logical :: same_k, calc_DOS

    ! Check input for what to do
    is_left = El%inf_dir == INF_NEGATIVE
    calc_DOS = present(DOS)
    if ( calc_DOS .and. .not. present(T) ) then
       call die('Need both DOS and T')
    end if

    zHS_allocated = .false.

    ! constants for this electrode
    nuo_E  = 2 * El%no_u
    nS     = nuo_E ** 2
    nuou_E = 2 * El%used_orbitals()
    nuS    = nuou_E ** 2
    ! create expansion k-points
    nq     = El%Bloch%size()
    ! We also need to invert to get the contribution in the
    ! reduced region
    reduce_size = nuo_E /= nuou_E
    nuouT_E = El%device_orbitals()

    if ( calc_DOS ) then
       if ( El%no_u > size(DOS) ) &
            call die('Error in DOS size for calculation bulk DOS')
       ! Initialize density of states
       DOS(1:El%no_u) = 0._dp
       T = 0._dp
    end if

    n_s = size(El%isc_off,dim=2)
    allocate(sc_off(3,n_s))
    do i = 1, n_s
      sc_off(1,i) = El%cell(1,1) * El%isc_off(1,i) + &
          El%cell(1,2) * El%isc_off(2,i) + &
          El%cell(1,3) * El%isc_off(3,i)
      sc_off(2,i) = El%cell(2,1) * El%isc_off(1,i) + &
          El%cell(2,2) * El%isc_off(2,i) + &
          El%cell(2,3) * El%isc_off(3,i)
      sc_off(3,i) = El%cell(3,1) * El%isc_off(1,i) + &
          El%cell(3,2) * El%isc_off(2,i) + &
          El%cell(3,3) * El%isc_off(3,i)
    end do
    
    ! whether we already have the H and S set correctly, 
    ! update accordingly, it will save a bit of time, but not much
    same_k = abs( bkpt(1) - El%bkpt_cur(1) ) < 1.e-8_dp
    same_k = same_k .and. abs( bkpt(2) - El%bkpt_cur(2) ) < 1.e-8_dp
    same_k = same_k .and. abs( bkpt(3) - El%bkpt_cur(3) ) < 1.e-8_dp
    if ( .not. same_k ) then
      El%bkpt_cur(:) = bkpt

      ! In case we do not need the hamiltonian
      ! This will be the case for non-bias points and when using bulk electrode
      same_k = .not. associated(El%HA)
    end if

    ! determine whether there is room enough
    if ( reduce_size ) then
       size_req(1) = 5 * nS
    else
       size_req(1) = 4 * nS
    end if
    if ( calc_DOS ) then
       size_req(2) = 9 * nS
    else
       size_req(2) = 8 * nS
    end if
    if ( size_req(1) + size_req(2) <= nzwork ) then

       ! we have enough room in the regular work-array for everything
       i = 0
       H00 => in_zwork(i+1:i+nS)
       i = i + nS
       S00 => in_zwork(i+1:i+nS)
       i = i + nS
       H01 => in_zwork(i+1:i+nS)
       i = i + nS
       S01 => in_zwork(i+1:i+nS)
       i = i + nS
       if ( reduce_size ) then
          GS => in_zwork(i+1:i+nS)
          i = i + nS
       end if
       zwork => in_zwork(i+1:nzwork)

    else if ( size_req(2) <= nzwork ) then

       ! we will allocate H00,H01,S00,S01,GS arrays
       call re_alloc(zHS,1,size_req(1),routine='next_GS')
       zHS_allocated = .true.

       i = 0
       H00 => zHS(i+1:i+nS)
       i = i + nS
       S00 => zHS(i+1:i+nS)
       i = i + nS
       H01 => zHS(i+1:i+nS)
       i = i + nS
       S01 => zHS(i+1:i+nS)
       i = i + nS
       if ( reduce_size ) then
          GS => zHS(i+1:i+nS)
       end if
       
       ! the work-array fits the input work-array
       zwork => in_zwork(1:nzwork)

    else if ( size_req(1) <= nzwork ) then
       ! we will allocate 8*nS work array

       i = 0
       H00 => in_zwork(i+1:i+nS)
       i = i + nS
       S00 => in_zwork(i+1:i+nS)
       i = i + nS
       H01 => in_zwork(i+1:i+nS)
       i = i + nS
       S01 => in_zwork(i+1:i+nS)
       i = i + nS
       if ( reduce_size ) then
          GS => in_zwork(i+1:i+nS)
       end if

       call re_alloc(zHS,1,size_req(2),routine='next_GS')
       zHS_allocated = .true.
       zwork => zHS(:)

    else

       call die('Your electrode is too large compared &
            &to your system in order to utilize the in-core &
            &calculation of the self-energies.')

    end if
    ! Get actual size of work-array
    nw = size(zwork)

    call init_mat_inversion(nuo_E)

    ! prepare the indices for the Gamma array
    ios = 1
    ioe = nuS 

    ! create the offset to be used for copying over elements
    off = nuo_E - nuou_E + 1

    call reclat(El%cell,rcell,1)

    ! loop over the repeated cell...
    q_loop: do iq = 1 , nq

       ! init qpoint in reciprocal lattice vectors
       kpt(:) = bkpt(:) + El%unfold_k(iq)
      
       ! Convert to 1/Bohr
       call kpoint_convert(rcell,kpt(1),kq(1),-2)

       ! Calculate transfer matrices @Ef (including the chemical potential)
       call set_HS_Transfer(El, n_s,sc_off, kq, &
           El%no_u, H00,S00,H01,S01)
       
       if ( .not. same_k ) then
         if ( reduce_size ) then
           ! we only need to copy over the data if we don't already have it calculated
!$OMP parallel default(shared)
           call copy_over(is_left,nuo_E,H00,nuou_E,El%HA(:,:,iq),off)
           call copy_over(is_left,nuo_E,S00,nuou_E,El%SA(:,:,iq),off)
!$OMP end parallel
         else
           ! means that it is the full thing
           call zcopy(nuS, H00(1), 1, El%HA(1,1,iq), 1)
           call zcopy(nuS, S00(1), 1, El%SA(1,1,iq), 1)
         end if
       end if

       if ( .not. reduce_size ) then
          ! Instead of doing a copy, we store it directly
          GS => El%GA(ios:ioe)
       end if

       ! calculate the contribution for this q-point
       if ( calc_DOS ) then
          call SSR_sGreen_DOS(nuo_E,Z,H00,S00,H01,S01, &
               El%accu, GS, &
               DOS(1:El%no_u), T, &
               nw, zwork, &
               final_invert = reduce_size, non_col=.true.)
       else
          call SSR_sGreen_NoDOS(nuo_E,Z,H00,S00,H01,S01, &
               El%accu, GS, &
               nw, zwork, &
               final_invert = reduce_size)
       end if

       if ( reduce_size ) then
          ! Copy over surface Green function
          ! first we need to determine the correct placement
!$OMP parallel default(shared)
          call copy_over(is_left,nuo_E,GS,nuou_E,El%GA(ios:ioe),off)
!$OMP end parallel

          ! we need to invert back as we don't need to
          ! expand. And the algorithm expects it to be in correct format
          if ( nq == 1 ) then
             call mat_invert(El%GA(ios:ioe),zwork(1:nuS),&
                  nuou_E, &
                  MI_IN_PLACE_LAPACK)
          end if
       end if

       ! correct indices of Gamma-array
       ios = ios + nuS
       ioe = ioe + nuS

    end do q_loop

    ! We normalize DOS as this will be comparable to a bulk
    ! calculation.
    if ( calc_DOS .and. nq > 1 ) then
       DOS(1:El%no_u) = DOS(1:El%no_u) / nq
    end if
       
    if ( zHS_allocated ) then
       call de_alloc(zHS, routine='next_GS')
    end if

    deallocate(sc_off)
    call clear_mat_inversion()

  contains
    
    subroutine copy_over(is_left,fS,from,tS,to,off)
      logical, intent(in) :: is_left
      integer, intent(in) :: fS, tS, off
      complex(dp), intent(in) :: from(fS,fS)
      complex(dp), intent(inout) :: to(tS,tS)

      integer :: i, j, ioff

      if ( is_left ) then
         ! Left, we use the last orbitals
         ioff = 1 - off ! ioff is private in OMP orphaned routines
!$OMP do private(j,i)
         do j = off , fS
            do i = off , fS
               to(ioff+i,ioff+j) = from(i,j)
            end do
         end do
!$OMP end do nowait
      else
         ! Right, the first orbitals
!$OMP do private(j,i)
         do j = 1 , tS
            do i = 1 , tS
               to(i,j) = from(i,j)
            end do
         end do
!$OMP end do nowait
      end if

    end subroutine copy_over

  end subroutine calc_next_GS_Elec_spinor


end module m_ts_electrode_spinor
