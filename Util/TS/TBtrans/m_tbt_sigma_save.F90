! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---


! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication

! Creation of the hssigma files.
#include "mpi_macros.f"

module m_tbt_sigma_save

  use precision, only : dp

  use m_tbt_hs, only : tTSHS
  use m_tbt_save, only : tNodeE
#ifdef NCDF_4
  use m_tbt_save, only : tbt_cdf_precision
#endif
  
  implicit none

  private 

  public :: init_Sigma_options, print_Sigma_options

#ifdef NCDF_4
  logical, save :: sigma_save      = .false.
  logical, save :: sigma_mean_save = .false.
  logical, save :: sigma_parallel  = .false.
  integer, save :: cmp_lvl    = 0

  public :: open_cdf_Sigma
  public :: init_Sigma_save
  public :: state_Sigma_save
  public :: state_Sigma2mean
#endif

contains


  subroutine init_Sigma_options(save_DATA)

    use dictionary
    use fdf

    type(dictionary_t), intent(inout) :: save_DATA

#ifdef NCDF_4

    sigma_save = fdf_get('TBT.CDF.SelfEnergy.Save',.false.)
    sigma_save = fdf_get('TBT.SelfEnergy.Save',sigma_save)
    if ( sigma_save ) then
      sigma_mean_save = fdf_get('TBT.CDF.SelfEnergy.Save.Mean',.false.)
      sigma_mean_save = fdf_get('TBT.SelfEnergy.Save.Mean',sigma_mean_save)
    end if
    cmp_lvl = fdf_get('CDF.Compress',0)
    cmp_lvl = fdf_get('TBT.CDF.Compress',cmp_lvl)
    cmp_lvl = fdf_get('TBT.CDF.SelfEnergy.Compress',cmp_lvl)
    if ( cmp_lvl < 0 ) cmp_lvl = 0
    if ( cmp_lvl > 9 ) cmp_lvl = 9
#ifdef NCDF_PARALLEL
    sigma_parallel = fdf_get('TBT.CDF.MPI',.false.)
    sigma_parallel = fdf_get('TBT.CDF.SelfEnergy.MPI',sigma_parallel)
    if ( sigma_parallel ) then
       cmp_lvl = 0
    end if
#endif

    if ( sigma_save .and. fdf_get('TBT.SelfEnergy.Only',.false.) ) then
       save_DATA = save_DATA // ('Sigma-only'.kv.1)
    end if
#endif
    
  end subroutine init_Sigma_options

  subroutine print_Sigma_options( save_DATA )

    use parallel, only: IONode
    use dictionary

    type(dictionary_t), intent(inout) :: save_DATA

    character(len=*), parameter :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'
    character(len=*), parameter :: f12='(''tbt: '',a,t53,''='',tr2,i0)'
    character(len=*), parameter :: f11='(''tbt: '',a)'

    if ( .not. IONode ) return
    
#ifdef NCDF_4
    write(*,f1)'Saving downfolded self-energies',sigma_save
    if ( .not. sigma_save ) return
#ifdef NCDF_PARALLEL
    write(*,f1)'Use parallel MPI-IO for self-energy file', sigma_parallel
#endif

    write(*,f1)'Only calc downfolded self-energies', &
         ('Sigma-only'.in.save_DATA)
    if ( cmp_lvl > 0 ) then
       write(*,f12)'Compression level of TBT.SE.nc files',cmp_lvl
    else
       write(*,f11)'No compression level of TBT.SE.nc files'
    end if
    write(*,f1)'k-average downfolded self-energies',sigma_mean_save
#else
    write(*,f11)'Saving downfolded self-energies not enabled (NetCDF4)'
#endif

  end subroutine print_Sigma_options

#ifdef NCDF_4

  subroutine open_cdf_Sigma(fname, ncdf, N_Elec, Elecs)
    use netcdf_ncdf, ncdf_parallel => parallel
    use ts_electrode_m, only: electrode_t

#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
#endif

    character(len=*), intent(in) :: fname
    type(hNCDF), intent(inout) :: ncdf

    integer, intent(in) :: N_Elec
    type(electrode_t), intent(in) :: Elecs(:)

#ifdef NCDF_PARALLEL
    integer :: iEl
    type(hNCDF) :: grp
#endif

    if ( .not. sigma_save ) return

#ifdef NCDF_PARALLEL
    if ( sigma_parallel ) then
      call ncdf_open(ncdf,fname, mode=ior(NF90_WRITE,NF90_MPIIO), &
          comm = MPI_COMM_ID(MPI_COMM_WORLD) )

      ! Assign all writes to be collective
      ! Collective is faster since we don't need syncronization
      call ncdf_par_access(ncdf, access=NF90_COLLECTIVE)
      do iEl = 1 , N_Elec
        call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)
        call ncdf_par_access(grp, access=NF90_COLLECTIVE)
      end do

     else
#endif
       call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#ifdef NCDF_PARALLEL
    end if
#endif
    
  end subroutine open_cdf_Sigma

  ! Save the self-energies of the electrodes and
  subroutine init_Sigma_save(fname, TSHS, r, btd, ispin, &
      N_Elec, Elecs, raEl, roElpd, btd_El, &
      nkpt, kpt, wkpt, NE, Eta, &
      a_Dev, a_Buf)

    use parallel, only : IONode
    use m_os, only : file_exist
    use byte_count_m, only: byte_count_t

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast, MPI_Logical, MPI_Barrier
#endif
    use ts_electrode_m
    use m_region
    use dictionary
    use m_tbt_save, only: add_cdf_common

    ! The file name that we save in
    character(len=*), intent(in) :: fname
    ! The full Hamiltonian and system at present investigation.
    ! Note the H have been shifted to zero energy
    type(tTSHS), intent(in) :: TSHS
    ! The device region that we are checking
    ! This is the device regions pivot-table!
    type(tRgn), intent(in) :: r, btd
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(electrode_t), intent(in) :: Elecs(:)
    type(tRgn), intent(in) :: raEl(:), roElpd(:), btd_El(:)
    integer, intent(in) :: nkpt
    real(dp), intent(in) :: kpt(:,:), wkpt(:)
    integer, intent(in) :: NE
    real(dp), intent(in) :: Eta

    ! Device atoms
    type(tRgn), intent(in) :: a_Dev
    ! Buffer atoms
    type(tRgn), intent(in) :: a_Buf

    type(hNCDF) :: ncdf, grp
    type(dictionary_t) :: dic
    type(tRgn) :: r_tmp

    logical :: prec_Sigma
    logical :: exist, same_dims
    character(len=256) :: c_tmp
    type(byte_count_t) :: mem
    integer :: i, iEl, no_e
    real(dp), allocatable :: r2(:,:)
#ifdef MPI
    integer :: MPIerror
#endif

    if ( .not. sigma_save ) return

    exist = file_exist(fname, Bcast = .true. )

    call tbt_cdf_precision('SelfEnergy','single',prec_Sigma)

    ! in case it already exists...
    if ( exist ) then

      ! Create a dictionary to check that the sigma file is the
      ! same
      call ncdf_open(ncdf,fname)

      dic = ('no_u'.kv.TSHS%no_u) // ('na_u'.kv.TSHS%na_u) // &
          ('nkpt'.kv.nkpt ) // ('no_d'.kv.r%n) // &
          ('ne'.kv. NE ) // ('n_btd'.kv.btd%n)
      dic = dic // ('na_d'.kv. a_Dev%n)
      if ( a_Buf%n > 0 ) then
        dic = dic // ('na_b'.kv.a_Buf%n)
      end if
      call ncdf_assert(ncdf,exist,dims=dic)
      call delete(dic)
#ifdef MPI
      call MPI_Bcast(same_dims,1,MPI_Logical,0, &
          MPI_Comm_World,MPIerror)
#endif
      if ( .not. same_dims ) then
        call die('Dimensions in the TBT.SE.nc file does not conform &
            &to the current simulation.')
      end if

      do iEl = 1 , N_Elec
        call ncdf_open_grp(ncdf,Elecs(iEl)%name,grp)
        dic = dic // ('no_e'.kv.Elecs(iEl)%o_inD%n)
        call ncdf_assert(grp,exist,dims=dic)
        if ( .not. exist ) then
          write(*,*) 'Assertion of dimensions in file: '//trim(fname)//' failed.'

          call die('We could not assert the dimensions TBT.SE.nc file.')
        end if
      end do
      call delete(dic)

      ! Check the variables
      ! Check the variables
      dic = ('lasto'.kvp. TSHS%lasto(1:TSHS%na_u) ) // &
          ('pivot'.kvp. r%r ) // ('btd'.kvp.btd%r)
      call rgn_copy(a_Dev, r_tmp)
      call rgn_sort(r_tmp)
      dic = dic // ('a_dev'.kvp.r_tmp%r )
      dic = dic // ('xa'.kvp. TSHS%xa)
      if ( a_Buf%n > 0 )then
        dic = dic // ('a_buf'.kvp.a_Buf%r )
      end if
      call ncdf_assert(ncdf,same_dims,vars=dic, d_EPS = 1.e-4_dp )
      call delete(dic,dealloc=.false.) ! we have them pointing...
#ifdef MPI
      call MPI_Bcast(same_dims,1,MPI_Logical,0, &
          MPI_Comm_World,MPIerror)
#endif
      if ( .not. same_dims ) then
        call die('pivot, lasto, xa or a_buf in the TBT.nc file does &
            &not conform to the current simulation.')
      end if
      call rgn_delete(r_tmp)

      ! Check the k-points
      allocate(r2(3,nkpt))
      do i = 1 , nkpt
        call kpoint_convert(TSHS%cell,kpt(1,i),r2(1,i),1)
      end do
      dic = ('kpt'.kvp. r2) // ('wkpt'.kvp. wkpt)
      call ncdf_assert(ncdf,same_dims,vars=dic, d_EPS = 1.e-7_dp )
      if ( .not. same_dims ) then
        call die('k-points or k-weights are not the same')
      end if
      call delete(dic,dealloc = .false. )
      deallocate(r2)

      call die('Currently the TBT.SE.nc file exists, &
          &we do not currently implement a continuation scheme.')

      ! We currently overwrite the Sigma-file
      if ( IONode ) then
        write(*,'(2a)')'tbt: Overwriting self-energy file: ',trim(fname)
      end if

    else

      if ( IONode ) then
        write(*,'(2a)')'tbt: Initializing self-energy file: ',trim(fname)
      end if

    end if

    ! We need to create the file
#ifdef NCDF_PARALLEL
    if ( sigma_parallel ) then
      call ncdf_create(ncdf,fname, mode=ior(NF90_NETCDF4,NF90_MPIIO), overwrite=.true., &
          comm = MPI_COMM_ID(MPI_COMM_WORLD), &
          parallel = .true. )
    else
      call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)
    end if
#else
    call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)
#endif

    ! Reset memory counter
    call mem%reset()

    ! Define all default variables
    call add_cdf_common(ncdf, TSHS, ispin, r, btd, &
        N_Elec, Elecs, raEl, roElpd, btd_El, &
        nkpt, kpt, wkpt, NE, Eta, &
        a_Dev, a_Buf, mem)
    
    dic = ('info'.kv.'Downfolded self-energy')
#ifdef TBT_PHONON
    dic = dic//('unit'.kv.'Ry**2')
#else
    dic = dic//('unit'.kv.'Ry')
#endif
    do iEl = 1 , N_Elec

      call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

      no_e = Elecs(iEl)%o_inD%n

      ! Chunking greatly reduces IO cost
      call ncdf_def_var(grp,'SelfEnergy', prec_Sigma, &
          (/'no_e','no_e','ne  ','nkpt'/), compress_lvl = cmp_lvl, &
          atts = dic , chunks = (/no_e,no_e,1,1/))
      call delete(dic)

      call mem%add_cdf(prec_Sigma, no_e, no_e, NE, nkpt)

      if ( sigma_mean_save ) then
        ! Add mean sigma
        call mem%add_cdf(prec_Sigma, no_e, no_e, NE)
      end if

    end do

    call delete(dic)

    call ncdf_close(ncdf)

#ifdef MPI
    ! Ensure that the processors are aligned
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

    if ( IONode ) then
      call mem%get_string(c_tmp)
      write(*,'(4a/)') 'tbt: Estimated file size of ', trim(fname), ': ', trim(c_tmp)
    end if

  end subroutine init_Sigma_save

  subroutine state_Sigma_save(ncdf, ikpt, nE, N_Elec, Elecs, nzwork,zwork)

    use iso_c_binding
    use parallel, only : Node, Nodes

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer
    use mpi_siesta, only : Mpi_double_precision
    USE_MPI_ONLY_STATUS
#endif
    use ts_electrode_m

    ! The file name we save too
    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_Elec
    type(electrode_t), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: zwork(:)

    type(hNCDF) :: grp
    integer :: iEl, iN, no_e, n_e
    type(c_ptr) :: sigma_ptr
    complex(dp), pointer :: Sigma2D(:,:)
#ifdef MPI
    integer :: MPIerror
    MPI_STATUS_TYPE :: status
#endif

    if ( .not. sigma_save ) return

    ! Save the energy-point
    ikpt_if: if ( ikpt == 1 ) then
    if ( parallel_io(ncdf) ) then
      if ( nE%iE(Node) > 0 ) then
        call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
      else
        call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/1/), count=(/0/))
      end if
    else
      do iN = 0 , Nodes - 1
        if ( nE%iE(iN) <= 0 ) cycle
        call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
      end do
    end if
    end if ikpt_if

#ifdef MPI
    if ( .not. sigma_parallel .and. Nodes > 1 ) then
      no_e = 0
      do iEl = 1 , N_Elec
        no_e = max(no_e,Elecs(iEl)%o_inD%n)
      end do
      n_e = no_e ** 2
      if ( n_e > nzwork ) then
        call die('Could not re-use the work array for Sigma &
            &communication.')
      end if
    end if
#endif

    do iEl = 1 , N_Elec
       
      call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

      no_e = Elecs(iEl)%o_inD%n

      sigma_ptr = c_loc(Elecs(iEl)%Sigma)
      call c_f_pointer(sigma_ptr, Sigma2D, [no_e, no_e])

      if ( nE%iE(Node) > 0 ) then
        call ncdf_put_var(grp,'SelfEnergy', Sigma2D, &
            start = (/1,1,nE%iE(Node),ikpt/) )
      else if ( sigma_parallel ) then
        ! Collective!!!
        call ncdf_put_var(grp,'SelfEnergy', Sigma2D, &
            start = (/1,1,1,ikpt/), count=(/0,0,0,0/))
      end if

#ifdef MPI
      if ( .not. sigma_parallel .and. Nodes > 1 ) then
        n_e = no_e ** 2
        if ( Node == 0 ) then
          sigma_ptr = c_loc(zwork)
          ! Because we are using a work-array to retrieve data
          call c_f_pointer(sigma_ptr, Sigma2D, [no_e, no_e])
          do iN = 1 , Nodes - 1
            if ( nE%iE(iN) <= 0 ) cycle
            call MPI_Recv(Sigma2D(1,1),n_e,MPI_Double_Complex,iN,iN, &
                Mpi_comm_world,status,MPIerror)
            call ncdf_put_var(grp,'SelfEnergy',Sigma2D, &
                start = (/1,1,nE%iE(iN),ikpt/) )
          end do
        else if ( nE%iE(Node) > 0 ) then
          call MPI_Send(Sigma2D(1,1),n_e,MPI_Double_Complex,0,Node, &
              Mpi_comm_world,MPIerror)
        end if
      end if
#endif

    end do

  end subroutine state_Sigma_save

  subroutine state_Sigma2mean(fname, N_Elec, Elecs)

    use iso_c_binding

    use parallel, only : IONode

    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
    use mpi_siesta, only : MPI_Barrier
#endif
    use ts_electrode_m

    ! The file name we save too
    character(len=*), intent(in) :: fname
    integer, intent(in) :: N_Elec
    type(electrode_t), intent(in) :: Elecs(N_Elec)

    type(dictionary_t) :: dic
    type(hNCDF) :: ncdf, grp
    integer :: iEl, iE, ikpt
    integer :: NE, nkpt, no_e
    real(dp), allocatable :: rwkpt(:)
    complex(dp), allocatable :: c2(:,:)
    type(c_ptr) :: sigma_ptr
    complex(dp), pointer :: Sigma(:,:)

#ifdef MPI
    integer :: MPIerror
#endif

    ! If we should not save the mean, we return immediately.
    if ( .not. sigma_save ) return
    if ( .not. sigma_mean_save ) return

    if ( .not. IONode ) then
#ifdef MPI
       call MPI_Barrier(Mpi_comm_world,MPIerror)
#endif
       return
    end if

    call timer('SE-mean', 1)

    ! We do this on one processor
    call ncdf_open(ncdf,fname, mode=NF90_WRITE)

    ! We read in the dimensions
    call ncdf_inq_dim(ncdf,'ne',len=NE)
    call ncdf_inq_dim(ncdf,'nkpt',len=nkpt)

    ! Allocate space
    allocate(rwkpt(nkpt))
    call ncdf_get_var(ncdf,'wkpt',rwkpt)

    ! When taking the mean of self-energies
    ! we need the transpose, hence we need half the
    ! contribution from Sigma and Sigma^T
    rwkpt(:) = 0.5_dp * rwkpt(:)

    ! Loop over all electrodes
    do iEl = 1 , N_Elec

      call delete(dic)

       ! We need to extend the netcdf file with the SigmaMean
       ! variable

       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       ! Get size of Sigma
       call ncdf_inq_dim(grp,'no_e',len=no_e)

       dic = ('info'.kv.'Downfolded self-energy, k-averaged')
       dic = dic//('unit'.kv.'Ry')
       ! Chunking greatly reduces IO cost
       call ncdf_def_var(grp,'SelfEnergyMean',NF90_DOUBLE_COMPLEX, &
            (/'no_e','no_e','ne  '/), chunks = (/no_e,no_e,1/) , &
            atts = dic ,compress_lvl = cmp_lvl )
       call delete(dic)

       ! Allocate space for the self-energy mean
       allocate(c2(no_e,no_e))

       ! Point the sigma
       ! This is a hack to ease the processing
       sigma_ptr = c_loc(Elecs(iEl)%Sigma)
       call c_f_pointer(sigma_ptr, Sigma, [no_e, no_e])

       ! loop over all energy points
       do iE = 1 , NE

         ! Loop over k-points to average
         call ncdf_get_var(grp,'SelfEnergy',Sigma, &
             start=(/1,1,iE,1/) )
         
         c2(:,:) = rwkpt(1) * ( Sigma + transpose(Sigma) )

         do ikpt = 2 , nkpt
           
           ! Loop over k-points to average
           call ncdf_get_var(grp,'SelfEnergy',Sigma, &
               start=(/1,1,iE,ikpt/) )
           
           c2(:,:) = c2(:,:) + rwkpt(ikpt) * ( Sigma + transpose(Sigma) )
           
         end do

         call ncdf_put_var(grp,'SelfEnergyMean',c2, start=(/1,1,iE/) )

       end do

       deallocate(c2)

    end do

    deallocate(rwkpt)
    
    call ncdf_close(ncdf)

    call timer('SE-mean', 2)
        
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

  end subroutine state_Sigma2mean

#endif

end module m_tbt_sigma_save
