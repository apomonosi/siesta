!
! Copyright (C) 1996-2016 The SIESTA group
! This file is distributed under the terms of the
! GNU General Public License: see COPYING in the top directory
! or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
!> Re-designed module to allow a pre-computation of the needed matrix elements
!> in parallel, followed by a globalization of the data among all the MPI processes.
!> Once the interpolation tables are setup, further calls to the matrix-element
!> evaluator (here renamed 'get_matel') are cheap. This has a dramatic effect in
!> some routines (such as nlefsm) that had to perform the table-building operations
!> under conditions that did not scale in parallel.

!> Concept: Rogeli Grima (BSC) and Alberto Garcia (ICMAB)
!> Initial implementation: Rogeli Grima (BSC)

module matel_table_m

  use matel_params_m
  use matel_ylm_m, only: spher_harm_t
  
  use precision, only : dp
  use alloc,     only : re_alloc, de_alloc, alloc_default, allocDefaults

  use parallel,  only : Node, Nodes
  use m_radfft,  only : radfft
  use m_matel_registry, only : RCUT
  use spher_harm, only : RLYLM, YLMEXP, ylmylm, lofilm
#ifdef MPI
  use mpi_siesta,  only : MPI_Allgatherv, MPI_Allgather
  use mpi_siesta,  only : MPI_INTEGER, MPI_COMM_WORLD
  use mpi_siesta,  only : MPI_DOUBLE_PRECISION, MPI_IN_PLACE
#endif      
  private

  integer, parameter, public :: MODE_S = 1, MODE_T = 2, MODE_XYZ = 3

  !> The main type that holds the interpolation tables and all the indexes needed to find
  !> the right radial tables. Note that it contains pointers to "spherical-harmonic decompositions",
  !> with two versions (1 and 2), but '2' is not always needed.
  type, public :: matel_t
    integer :: MODE           ! One of MODE_S, MODE_T or MODE_XYZ

    integer :: min_ig1        ! Registry index of first "F_1" function
    integer :: max_ig1        ! Registry index of last "F_1" function
    integer :: min_ig2        ! Registry index of first "F_2" function
    integer :: max_ig2        ! Registry index of last "F_2" function

    integer :: MAT_ROWS       ! Leading size of matrix (in num of harmonics)
    integer :: MAT_COLS       ! Columns of matrix (in num of harmonics)

    type(spher_harm_t), pointer :: SpHa1 => null()  ! Orbitals & projectors
    type(spher_harm_t), pointer :: SpHa2 => null()
    integer :: BASE_SPHA1     ! Base index for block of interest in SpHa1
    integer :: BASE_SPHA2     ! Base index for block of interest in SpHa2

    real(dp), pointer :: FFR(:,:,:) => null()   ! Radial functions. To save space we only store
    integer           :: MFFR, NFFR             ! linear independent functions.
    integer,  pointer :: INDFFR(:) => null()    !

    real(dp), pointer :: FFY(:,:) => null()    ! Expansion of the product of two spherical harmonics
    integer,  pointer :: ILMFF(:) => null()    ! Spherical harmonics indexes
    integer           :: MFFY, NFFY            !
    integer,  pointer :: INDFFY(:) => null()   ! Number of radial functions (after COMPUTE_RADEXP)
                                               ! Index of radial functions (after REDUCE_RADEXP)

    integer           :: MILM
    !> Work space for computation of rlylm
    real(dp), pointer :: Y(:) => null(), DYDR(:,:) => null()

    contains
      procedure :: init
      procedure :: get_matel
      procedure :: delete   
      procedure :: write_info_matel_table
      ! called only internally
      procedure :: compute_radexp
      procedure :: reduce_radexp
      procedure :: remove_redundant_data

  end type matel_t

  real(dp),         parameter :: EXPAND    =  1.20_dp
  integer,          parameter :: MINEXPAND =  32
  real(dp),         parameter :: FFTOL     =  1.e-8_dp
  character(len=*), parameter :: MYNAME    =  'MATEL_TABLE'

  external :: die

contains
  
  subroutine init(this, IOPER, SpHa1, lo1, hi1, SpHa2, lo2, hi2)
    !     Initialize a MATEL table structure using the previously calculated expansions
    !     in spherical harmonics of the fourier transforms of the functions.

    !      M(R) =  < F_1(r) | OP | F_2(r-R) >
    !      M(R) =  < F_1(r) | G_2(r-R) >

    !           G_2(r) is F_2(r) when OP is the identity (simple overlap)
    !           G_2(r) can be  {x,y,z}*F_2  when OP = {x,y,z}
    !           When OP is the laplacian, this routine will just insert a factor of Q^2
    !           in the product of F_1(Q) and F_2(Q).

    implicit none

    class(matel_t)        :: this

    !> One of MODE_S, MODE_T or MODE_XYZ
    integer, intent(in) :: IOPER

    !> Structures with the expansions of F_1 and F_2 in spherical harmonics
    type(spher_harm_t), intent(in), target :: SpHa1, SpHa2

    !> First and last entries to consider in SpHa1
    integer, intent(in) :: lo1, hi1
    !> First and last entries to consider in SpHa2
    integer, intent(in) :: lo2, hi2

    integer :: n_harms_hi, n_harms_lo
    
    this%MODE  = IOPER

    this%min_ig1 = lo1 + spha1%base
    this%max_ig1 = hi1 + spha1%base
    this%min_ig2 = lo2 + spha2%base
    this%max_ig2 = hi2 + spha2%base

    ! Note simple pointing, without allocation
    ! The SpHa objects should be already allocated and configured
    this%SpHa1 => SpHa1
    this%SpHa2 => SpHa2

    ! Find total number of harmonics in the block of SpHa1 we are
    ! interested in

    n_harms_hi = this%SpHa1%get_Dim(hi1)
    n_harms_lo = 0
    if (lo1 > 1) then
       n_harms_lo = this%SpHa1%get_Dim(lo1 - 1)
    endif
    this%MAT_ROWS = n_harms_hi - n_harms_lo
    this%BASE_SPHA1 = n_harms_lo

    ! Find total number of harmonics in the block of SpHa2 we are
    ! interested in

    n_harms_hi = this%SpHa2%get_Dim(hi2)
    n_harms_lo = 0
    if (lo2 > 1) then
       n_harms_lo = this%SpHa2%get_Dim(lo2 - 1)
    endif
    this%MAT_COLS = n_harms_hi - n_harms_lo
    this%BASE_SPHA2 = n_harms_lo

    ! Compute and reduce Radial Expansion
    call this%compute_RadExp()
    call this%reduce_RadExp()

    ! If we have used several MPI processes to build
    ! the tables, we need to remove the redundant data
    ! that was not handled initially by the same process.

    if (Nodes > 1) then
       call this%remove_redundant_data()
    endif
    
#ifdef SIESTA__MATEL_INIT_DEBUG    
    call this%write_info_matel_table(6,4)
#endif
    
  end subroutine init

  subroutine compute_RadExp(this)
    ! Fill a matrix with the radial expansion functions
    ! Compute the radial expansion between all harmonics in parallel.
    ! For very harmonic permutation, we compute its Radial functions and
    ! the expansion of the product of spherical harmonics
    
    use interpolation, only: spline
    implicit none
    
    class(matel_t) :: this
    integer :: dim1, dim2, nlocal, first, last, i, j, &
               IR, JR, JG, ih1, ih2, ig1, ig2, IQ, &
               L1, L2, L3, L1L2, JLM, NILM
    real(dp) :: C, Q, R, DFFR0, DFFRMX, CPROP
    integer, pointer :: IFFR(:) => null()
    real(dp), pointer :: FQ1(:) => null(), FQ2(:) => null()
    real(dp), pointer :: FFQ(:) => null(), FFL(:) => null(), CFFR(:) => null()
    logical, external :: propor

    ! Initialize temporal arrays
    nullify(this%Y,this%DYDR)
    this%MILM = (MINEXPAND*MINEXPAND+1)**2
    call RE_ALLOC(this%Y, 1, this%MILM, 'Y', MYNAME)
    call RE_ALLOC(this%DYDR, 1, 3, 1, this%MILM, 'DYDR', MYNAME)

    dim1 = this%MAT_ROWS
    dim2 = this%MAT_COLS
    nullify( this%INDFFY )
    call RE_ALLOC(this%INDFFY, 1, dim1*dim2+1, 'INDFFY', MYNAME)

    nlocal = GET_LOOP_LIMITS(dim1*dim2, Node, first, last)

    nullify(this%FFR,this%INDFFR)
    this%MFFR = MAX(INT(dim1*EXPAND), dim1+MINEXPAND)
    call RE_ALLOC(this%FFR, 0, NRTAB, 1, 2, 1, this%MFFR, 'FFR', MYNAME)
    this%NFFR = 0

    nullify(this%FFY,this%ILMFF,this%INDFFR)
    this%MFFY = MAX(INT(nlocal*EXPAND), nlocal+MINEXPAND)
    call RE_ALLOC(this%FFY, 1, 1, 1, this%MFFY, 'FFY', MYNAME)
    call RE_ALLOC(this%INDFFR, 1, this%MFFY, 'INDFFR', MYNAME)
    call RE_ALLOC(this%ILMFF, 1, this%MFFY, 'ILMFF', MYNAME)
    this%NFFY = 0

    nullify(IFFR,CFFR,FFQ,FFL)
    call RE_ALLOC(FFQ, 0, NQ, 'FFQ', MYNAME)
    call RE_ALLOC(FFL, 0, NQ, 'FFL', MYNAME)
    L1L2 = MINEXPAND
    call RE_ALLOC(IFFR, 0, L1L2, 'IFFR', MYNAME)
    call RE_ALLOC(CFFR, 0, L1L2, 'CFFR', MYNAME)

    ! Iterate all combinations of harmonics
    do i = first, last
      ! Get the harmonic terms IDS
      ! Note accounting for block of interest
      ih1 = MOD(i-1, dim1) + 1 + this%BASE_SPHA1
      ih2 = (i-1)/dim1 + 1 + this%BASE_SPHA2
      ! Get the matel-registry IDS of the functions
      ! from which the harmonic terms originate
      ig1 = this%SpHa1%Harm2Orb(ih1)
      ig2 = this%SpHa2%Harm2Orb(ih2)

      ! Check interaction range
      if (RCUT(IG1)+RCUT(IG2) > RMAX) then
        call die('MATEL: NQ too small for required cutoff.')
      endif
      FQ1(0:) => this%SpHa1%F(:,ih1)
      FQ2(0:) => this%SpHa2%F(:,ih2)

      ! Find orbitals convolution by multiplication in k-space
      C = (2.0_dp * PI)**1.5_dp
      do IQ = 0, NQ
        FFQ(IQ) = C * FQ1(IQ) * FQ2(IQ)
        if (this%MODE==MODE_T) then
          Q = IQ * DQ
          FFQ(IQ) = FFQ(IQ) * Q*Q
        endif
      enddo
      
      ! Loop on possible values of l quantum number of product
      L1 = LOFILM(this%SpHa1%ILM(IH1))
      L2 = LOFILM(this%SpHa2%ILM(IH2))
      if (L1+L2 > L1L2) then
        ! Reallocate more memory if necessary
        L1L2 = L1+L2
        call RE_ALLOC(IFFR, 0, L1L2, 'IFFR', MYNAME)
        call RE_ALLOC(CFFR, 0, L1L2, 'CFFR', MYNAME)
      endif
      
      do L3 = ABS(L1-L2), L1+L2, 2
        ! Return to real space
        call RADFFT(L3, NQ, NQ*PI/RMAX, FFQ, FFL)
        ! FFL(NQ) = 0._dp
        if (MOD(ABS(L1-L2-L3)/2, 2) /= 0) then
          do IR = 0, NR
            FFL(IR) = -FFL(IR)
          enddo
        endif
        ! Divide by R**L
        if (L3 /= 0) then
          do IR = 1, NR
            R = IR * DR
            FFL(IR) = FFL(IR) / R**L3
          enddo
          ! Parabolic extrapolation to R=0
          FFL(0) = (4.0_dp * FFL(1) - FFL(2)) / 3.0_dp
        endif
        ! Select NRTAB out of NR points
        if (MOD(NR, NRTAB) /= 0) then
          call DIE('matel ERROR: NQ must be multiple of NRTAB')
        endif
        do IR = 0, NRTAB
          JR = IR * NR / NRTAB
          FFL(IR) = FFL(JR)
        enddo

        ! To save space, we check if the new radial function has
        ! been seen before (save a global factor of proportionality)
        IFFR(L3) = 0
        do JG = 1, this%NFFR
          if (PROPOR(NRTAB, FFL(1), this%FFR(1,1,JG), FFTOL, CPROP)) then
            ! If found, save the index and the coefficient of proportionality
            IFFR(L3) = JG
            CFFR(L3) = CPROP
            exit
          endif
        enddo
        ! If not found, store new radial function
        if (IFFR(L3)==0) then
          this%NFFR = this%NFFR + 1
          if (this%NFFR > this%MFFR) then
            this%MFFR = EXPAND * this%NFFR
            call RE_ALLOC(this%FFR, 0, NRTAB, 1, 2, 1, this%MFFR, 'FFR', MYNAME)
          endif
          IFFR(L3) = this%NFFR
          CFFR(L3) = 1._dp
          do IR = 0, NRTAB
            this%FFR(IR,1,this%NFFR) = FFL(IR)
          enddo
          ! Setup spline interpolation
          ! Force derivative, rather than second derivative, to zero
          ! DFFR0 = HUGE(1.0_dp)
          DFFR0 = 0.0_dp
          DFFRMX = 0.0_dp
          call SPLINE(RMAX/NRTAB, this%FFR(0:NRTAB,1,this%NFFR), &
                      NRTAB+1, DFFR0, DFFRMX, &
                      this%FFR(0:NRTAB,2,this%NFFR))
        endif
      enddo

      ! Reallocate some arrays for the angular expansion
      NILM = (L1+L2+1)**2
      if (NILM > this%MILM) then
        call RE_ALLOC(this%Y, 1, NILM, 'Y', MYNAME, .FALSE.)
        call RE_ALLOC(this%DYDR, 1, 3, 1, NILM, 'DYDR', MYNAME, .FALSE.)
        this%MILM = NILM
      endif
      if (this%NFFY+NILM > this%MFFY) then
        this%MFFY = EXPAND * (this%NFFY+NILM)
        call RE_ALLOC(this%FFY, 1, 1, 1, this%MFFY, 'FFY', MYNAME)
        call RE_ALLOC(this%ILMFF, 1, this%MFFY, 'ILMFF', MYNAME)
        call RE_ALLOC(this%INDFFR, 1, this%MFFY, 'INDFFR', MYNAME)
      endif
      
      ! Expand the product of two spherical harmonics (SH) also in SH
      call YLMEXP(L1+L2, RLYLM, YLMYLM, this%SpHa1%ILM(IH1), &
                 this%SpHa2%ILM(IH2), 1, 1, 1.0_dp, NILM, &
                 this%ILMFF(this%NFFY+1:), this%FFY(:,this%NFFY+1:))
                 
      ! Loop on possible lm values of orbital product
      do J = 1, NILM
        this%NFFY = this%NFFY + 1
        JLM = this%ILMFF(this%NFFY)
        L3 = LOFILM(JLM)
        this%INDFFR(this%NFFY) = IFFR(L3)
        this%FFY(1,this%NFFY) = this%FFY(1,this%NFFY) * CFFR(L3)
      enddo
      this%INDFFY(I) = NILM
    enddo

    call DE_ALLOC(CFFR, 'CFFR', MYNAME)
    call DE_ALLOC(IFFR, 'IFFR', MYNAME)
    call DE_ALLOC(FFQ, 'FFQ', MYNAME)
    call DE_ALLOC(FFL, 'FFL', MYNAME)
  end subroutine compute_RadExp

  subroutine reduce_RadExp(this)
    ! Reduce Radial functions and expansion of spherical harmonics
    ! Reduce the values of the different arrays of MATEL
    ! At input, INDFFY contains the number of harmonics of every
    ! interaction. At exit, it can be used as and index.
    ! NOTE: TO DO!!! The reduction of FFR can be improved. We are saving
    ! functions that are not linearly independent (since we did the
    ! work on different processors). We should choose between saving time
    ! or saving space.
    
    implicit none
    class(matel_t) :: this
    
    ! Local variables
    integer :: I, PREV, NACUM, DIM1, DIM2
#ifdef MPI
    integer :: J, FIRST, LAST, MPIERR, IG1, IG2
    integer,  pointer :: COUNT(:) => null(), DISPL(:) => null()
    integer,  pointer :: COUN2(:) => null(), DISP2(:) => null()
    integer,  pointer :: G_ILMFF(:) => null(), G_INDFFR(:) => null()
    real(dp), pointer :: G_FFR(:,:,:) => null(), G_FFY(:,:) => null()
#endif

    dim1 = this%MAT_ROWS
    dim2 = this%MAT_COLS
    
#ifdef MPI
    if (NODES > 1) then
      ! Allocate space for temporal arrays
       nullify(count,displ,coun2,disp2)
      call RE_ALLOC(COUNT, 0, NODES-1, 'COUNT,', MYNAME)
      call RE_ALLOC(DISPL, 0, NODES, 'DISPL,', MYNAME)
      call RE_ALLOC( COUN2, 0, NODES-1, 'COUN2,', MYNAME )
      call RE_ALLOC( DISP2, 0, NODES, 'DISP2,', MYNAME )
      DISPL(1) = 0
      do I = 0, NODES-1
        COUNT(I) = GET_LOOP_LIMITS(dim1*dim2, I, first, last)
        DISPL(I+1) = DISPL(I) + COUNT(I)
      enddo
      ! Reduce INDFFY
      call MPI_Allgatherv(MPI_IN_PLACE, COUNT(NODE), MPI_INTEGER, &
                         this%INDFFY, COUNT, DISPL, MPI_INTEGER, &
                         MPI_COMM_WORLD, MPIERR)
    endif
#endif

    ! Create the acumulated version of INDFFY
    PREV = this%INDFFY(1)
    this%INDFFY(1) = 1
    do I = 1, dim1*dim2
      NACUM = PREV + this%INDFFY(I)
      PREV = this%INDFFY(I+1)
      this%INDFFY(I+1) = NACUM
    enddo
    
#ifdef MPI
    if (NODES > 1) then
      DISP2(0) = 0
      do I = 0, NODES-1
        COUN2(I) = this%INDFFY(DISPL(I+1)+1)-this%INDFFY(DISPL(I)+1)
        DISP2(I+1) = DISP2(I) + COUN2(I)
      enddo

      ! Get the accumulated NFFR
      call MPI_AllGather(this%NFFR, 1, MPI_INTEGER, COUNT, &
                        1, MPI_INTEGER, MPI_COMM_WORLD, MPIERR)
      DISPL(0) = 0
      do I = 0, NODES-1
        DISPL(I+1) = DISPL(I) + COUNT(I)
      enddo

      ! Update INDFFR to local indices
      this%INDFFR(1:this%NFFY) = this%INDFFR(1:this%NFFY)+DISPL(NODE)
      ! Get the global NFFY
      this%NFFY = DISP2(NODES)

      ! Reduce FFY
      nullify(G_FFY)
      call RE_ALLOC(G_FFY, 1, 1, 1, this%NFFY, 'FFY', MYNAME)
      call MPI_Allgatherv(this%FFY, COUN2(NODE), MPI_DOUBLE_PRECISION, &
        G_FFY, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, &
        MPIERR)
      call DE_ALLOC(this%FFY, 'FFY', MYNAME)
      this%FFY => G_FFY

      ! Reduce ILMFF
      nullify(G_ILMFF)
      call RE_ALLOC(G_ILMFF, 1, this%NFFY, 'ILMFF', MYNAME)
      call MPI_Allgatherv(this%ILMFF, COUN2(NODE), MPI_INTEGER, &
        G_ILMFF, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD, &
        MPIERR)
      call DE_ALLOC(this%ILMFF, 'ILMFF', MYNAME)
      this%ILMFF => G_ILMFF

      ! Reduce INDFFR
      nullify(G_INDFFR)
      call RE_ALLOC(G_INDFFR, 1, this%NFFY, 'INDFFR', MYNAME)
      call MPI_Allgatherv(this%INDFFR, COUN2(NODE), MPI_INTEGER, &
        G_INDFFR, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD, &
        MPIERR)
      call DE_ALLOC(this%INDFFR, 'INDFFR', MYNAME)
      this%INDFFR => G_INDFFR

      this%NFFR = DISPL(NODES)
      COUN2 = COUNT*2*(NRTAB+1)
      DISP2 = DISPL*2*(NRTAB+1)
      nullify(G_FFR)
      call RE_ALLOC(G_FFR, 0, NRTAB, 1, 2, 1, this%NFFR, 'FFR', MYNAME)
      call MPI_Allgatherv(this%FFR, COUN2(NODE), MPI_DOUBLE_PRECISION, &
        G_FFR, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, &
        MPIERR)
      call DE_ALLOC(this%FFR, 'FFR', MYNAME)
      this%FFR => G_FFR

      call DE_ALLOC(COUNT, 'COUNT', MYNAME)
      call DE_ALLOC(DISPL, 'DISPL', MYNAME)
      call DE_ALLOC(COUN2, 'COUN2', MYNAME)
      call DE_ALLOC(DISP2, 'DISP2', MYNAME)
    endif
#endif
  end subroutine reduce_RadExp

  !> Utility function      
  function GET_LOOP_LIMITS(N, ID, FIRST, LAST) RESULT(SIZE)
    ! Get the loop limits adequate to split N elements among Nodes
    ! INPUT:
    !   integer N  : Number of elements to distribute
    !   integer ID : Id of the Node
    ! OUTPUT:
    !   integer FIRST : first element of the loop
    !   integer LAST  : last element of the loop
    !   integer SIZE  : Number of elements of the loop
    
    implicit none
    integer, intent(in)  :: N, ID
    integer, intent(out) :: FIRST, LAST
    integer              :: size, di, mo
    
    if (NODES == 1) then
      FIRST = 1
      LAST  = N
    else
      di    = N / NODES
      mo    = MOD(N, NODES)
      FIRST = 1 + ID*di + MIN(ID, mo)
      LAST  = (ID+1)*di + MIN(ID+1, mo)
    endif
    SIZE = LAST - FIRST + 1
  end function GET_LOOP_LIMITS

  !>  Get the matrix element and its derivative
  subroutine get_matel(THIS, IG1, IG2, R12, S12, DSDR)
    use interpolation, only: splint

    implicit none

    !> Matel_table object
    class(matel_t) :: THIS

    !> Global index of 1st function
    integer, intent(in)      :: IG1
    !> Global index of 2nd function
    integer, intent(in)      :: IG2

    !> Vector from first to second atom
    real(dp), intent(in)    :: R12(3)

    !> Matrix element
    real(dp), intent(out)    :: S12
    !> Derivative (gradient) of S12 with respect to R12.
    real(dp), intent(out)    :: DSDR(3)

    real(dp), parameter :: TINY = 1.e-12_dp

    integer  :: IX, IH, LMAX, IFFY, JLM, JFFR, IFLM1, IFLM2
    real(dp) :: X12(3), R, SR, DSRDR

    if ((ig1 < this%min_ig1) .or. (ig1 > this%max_ig1)) then
       call die("IG1 out of range in get_matel")
    endif
    if ((ig2 < this%min_ig2) .or. (ig2 > this%max_ig2)) then
       call die("IG2 out of range in get_matel")
    endif

    ! Initialize output
    S12 = 0.0_dp
    DSDR(1) = 0.0_dp
    DSDR(2) = 0.0_dp
    DSDR(3) = 0.0_dp

    ! Avoid R12=0
    X12(1) = R12(1)
    X12(2) = R12(2)
    X12(3) = R12(3)
    R = SQRT(X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3))
    if (R < TINY) then
      X12(3) = TINY
      R = SQRT(X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3))
    endif

    ! Find if orbitals are far (out of range)
    if (R <= RCUT(IG1)+RCUT(IG2)) then
      ! Find spherical harmonics times R**L
      do IFLM1 = this%SpHa1%INDF(IG1), this%SpHa1%INDF(IG1+1)-1
        do IFLM2 = this%SpHa2%INDF(IG2), this%SpHa2%INDF(IG2+1)-1
          ! Note
          IH = (IFLM1-this%BASE_SPHA1) + &
               (IFLM2-this%BASE_SPHA2-1)*this%MAT_ROWS
          LMAX = 0
          do IFFY = this%INDFFY(IH), this%INDFFY(IH+1)-1
            JLM = this%ILMFF(IFFY)
            LMAX = MAX(LMAX, LOFILM(JLM))
          enddo
          call RLYLM(LMAX, X12, this%Y, this%DYDR)
          
          ! Interpolate radial functions and obtain SH expansion
          do IFFY = this%INDFFY(IH), this%INDFFY(IH+1)-1
            JFFR = this%INDFFR(IFFY)
            call SPLINT(RMAX/NRTAB, this%FFR(0:NRTAB,1,JFFR), &
                       this%FFR(0:NRTAB,2,JFFR), NRTAB+1, R, SR, DSRDR)
            JLM = this%ILMFF(IFFY)
            S12 = S12 + SR * this%FFY(1,IFFY) * this%Y(JLM)
            do IX = 1, 3
              DSDR(IX) = DSDR(IX) + &
                DSRDR * this%FFY(1,IFFY) * this%Y(JLM) * X12(IX) / R + &
                SR * this%FFY(1,IFFY) * this%DYDR(IX,JLM)
            enddo
          enddo
        enddo
      enddo
    endif
  end subroutine get_matel

  subroutine delete(this)
    ! Deallocates all allocatable arrays in the MATEL type
    ! Note: SpHa1 and SpHa2 are just pointed to, not owned
    
    implicit none
    class(matel_t) :: this
    
    ! Local variables
    character(len=*), parameter :: MYNAME = 'MATEL_DELETE'

    ! Deallocate FFR, FFY arrays and associated index arrays
    if (associated(this%FFR)) then
      call de_alloc(this%FFR, 'FFR', MYNAME)
    endif
    nullify(this%FFR)
    
    if (associated(this%FFY)) then
      call de_alloc(this%FFY, 'FFY', MYNAME)
    endif
    nullify(this%FFY)
    
    if (associated(this%ILMFF)) then
      call de_alloc(this%ILMFF, 'ILMFF', MYNAME)
    endif
    nullify(this%ILMFF)
    
    if (associated(this%INDFFR)) then
      call de_alloc(this%INDFFR, 'INDFFR', MYNAME)
    endif
    nullify(this%INDFFR)
    
    if (associated(this%INDFFY)) then
      call de_alloc(this%INDFFY, 'INDFFY', MYNAME)
    endif
    nullify(this%INDFFY)

    ! Clean up workspace arrays
    if (associated(this%Y)) then
      call de_alloc(this%Y, 'Y', MYNAME)
    endif
    nullify(this%Y)
    
    if (associated(this%DYDR)) then
      call de_alloc(this%DYDR, 'DYDR', MYNAME)
    endif
    nullify(this%DYDR)

    ! Just nullify pointers to spherical harmonic objects
    nullify(this%SpHa1)
    nullify(this%SpHa2)

    ! Reset scalars
    this%MODE = 0
    this%min_ig1 = 0 
    this%max_ig1 = 0
    this%min_ig2 = 0
    this%max_ig2 = 0
    this%MAT_ROWS = 0
    this%MAT_COLS = 0
    this%BASE_SPHA1 = 0
    this%BASE_SPHA2 = 0
    this%MFFR = 0
    this%NFFR = 0
    this%MFFY = 0
    this%NFFY = 0
    this%MILM = 0

  end subroutine delete

  SUBROUTINE WRITE_INFO_MATEL_TABLE(THIS, LUN, INDENT)
    ! Dumps information about a MATEL table in YAML-like format
    use precision, only: dp
    use parallel,  only: IONode
    implicit none
    ! Arguments
    CLASS(matel_t) :: THIS

    integer, intent(in) :: LUN      ! Logical unit number for output
    integer, intent(in), optional :: INDENT  ! Indentation level

    ! Local variables 
    integer :: IND
    character(len=4)  :: INDSTR
    character(len=80) :: SPACES

    if (.not. IONode) return

    ! Set indentation
    IND = 0
    if (present(INDENT)) IND = INDENT
    write(INDSTR,'(I4)') IND
    INDSTR = adjustl(INDSTR)
    SPACES = ' '

    ! Write header info
    write(LUN,'(A,A)') SPACES(1:IND), 'matel_table:'

    IND = IND + 2

    ! Write base info using a simple format
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'mode: ', THIS%MODE
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'min_ig1: ', THIS%MIN_IG1  
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'max_ig1: ', THIS%MAX_IG1
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'min_ig2: ', THIS%MIN_IG2
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'max_ig2: ', THIS%MAX_IG2
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'mat_rows: ', THIS%MAT_ROWS
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'mat_cols: ', THIS%MAT_COLS
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'mffr: ', THIS%MFFR
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'nffr: ', THIS%NFFR
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'mffy: ', THIS%MFFY
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'nffy: ', THIS%NFFY

    ! Dump SpHa1 info
    write(LUN,'(A,A)') SPACES(1:IND), 'spha1:'
    IND = IND + 2

    write(LUN,'(A,A,I0)') SPACES(1:IND), 'n: ', THIS%SpHa1%N
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'base: ', THIS%SpHa1%BASE
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'size of ilm array: ',     &
                                         size(THIS%SpHa1%ILM)

    ! Dump SpHa2 info  
    IND = IND - 2
    write(LUN,'(A,A)') SPACES(1:IND), 'spha2:'
    IND = IND + 2

    write(LUN,'(A,A,I0)') SPACES(1:IND), 'n: ', THIS%SpHa2%N
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'base: ', THIS%SpHa2%BASE
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'size of ilm array: ',   &
                          size(THIS%SpHa2%ILM)

    ! Record some stats about the FFR array
    IND = IND - 2
    write(LUN,'(A,A)') SPACES(1:IND), 'ffr_stats:'
    IND = IND + 2

    write(LUN,'(A,A,I0)') SPACES(1:IND), 'total elements: ', size(THIS%FFR)
    write(LUN,'(A,A,I0)') SPACES(1:IND), 'non-zero elements: ', &
                          count(abs(THIS%FFR) > 1.e-10_dp)

    RETURN
  end subroutine WRITE_INFO_MATEL_TABLE

  subroutine REMOVE_REDUNDANT_DATA(this)
    use precision,     only: dp
    use parallel,      only: IONode
    use alloc,         only: re_alloc, de_alloc

    implicit none
    class(MATEL_t)       :: this

    integer              :: i, j, nkeep, orig_size
    real(dp)             :: cprop
    integer,  pointer    :: redundant_map(:) => null()
    real(dp), pointer    :: scale_factors(:) => null()
    real(dp), pointer    :: tmp_FFR(:,:,:) => null()
    logical              :: is_redundant
    logical, external    :: propor

    character(len=132)   :: message
    real(dp)            :: reduction_pct

    ! Store original size for statistics
    orig_size = this%NFFR

    ! Allocate temporary arrays
    nullify(redundant_map, scale_factors)
    call re_alloc( redundant_map, 1, this%NFFR, &
         'redundant_map', 'REMOVE_REDUNDANT_DATA' )
    call re_alloc( scale_factors, 1, this%NFFR, &
         'scale_factors', 'REMOVE_REDUNDANT_DATA' )
    redundant_map = 0
    scale_factors = 1.0_dp

    ! Allocate temporary array for unique functions
    nullify(tmp_FFR)
    call re_alloc( tmp_FFR, 0, NRTAB, 1, 2, 1, this%NFFR, &
         'tmp_FFR', 'REMOVE_REDUNDANT_DATA' )

    ! Identify redundant functions
    nkeep = 0
    do i = 1, this%NFFR
       is_redundant = .false.
       ! Compare with previously kept functions
       do j = 1, nkeep
          if ( propor(NRTAB, this%FFR(1,1,i), tmp_FFR(1,1,j), FFTOL, cprop) ) then
             is_redundant = .true.
             redundant_map(i) = j
             scale_factors(i) = cprop
             EXIT
          endif
       enddo

       if (.not. is_redundant) then
          nkeep = nkeep + 1
          ! Copy this function to the kept set
          tmp_FFR(:,:,nkeep) = this%FFR(:,:,i)
          redundant_map(i) = nkeep
          scale_factors(i) = 1.0_dp
       endif
    enddo

    ! Update array sizes
    this%NFFR = nkeep
    this%MFFR = nkeep

    ! Reallocate and update FFR with kept functions
    call de_alloc( this%FFR, 'FFR', 'REMOVE_REDUNDANT_DATA' )
    call re_alloc( this%FFR, 0, NRTAB, 1, 2, 1, this%MFFR,    &
         'FFR', 'REMOVE_REDUNDANT_DATA' )
    this%FFR(:,:,1:nkeep) = tmp_FFR(:,:,1:nkeep)

    ! Update indices and coefficients in FFY using the map and scale factors
    do i = 1, size(this%INDFFR)
       if (this%INDFFR(i) .gt. 0) then
          j = this%INDFFR(i)
          this%INDFFR(i) = redundant_map(j)
          this%FFY(1,i) = this%FFY(1,i) * scale_factors(j)
       endif
    enddo

#ifdef SIESTA__MATEL_INIT_DEBUG    
    ! Print statistics from root node only

    if (IONode) then
       if (orig_size .gt. 0) then
          reduction_pct = 100.0_dp * (orig_size - nkeep) / orig_size
       else
          reduction_pct = 0.0_dp
       endif

       write(6,'(/,a)') 'Matrix element table cleanup statistics:'
       write(6,'(a,i8)')  '  Original number of functions: ', orig_size
       write(6,'(a,i8)')  '  Final number of functions:    ', nkeep
       write(6,'(a,i10)') '  Array elements saved:         ', &
            (orig_size - nkeep) * (NRTAB + 1) * 2
       write(6,'(a,f8.2,a)') '  Size reduction:             ', reduction_pct, ' %'
       write(6,'(a)') ' '
    endif
#endif
    
    ! Cleanup
    call de_alloc( redundant_map, 'redundant_map', 'REMOVE_REDUNDANT_DATA' )
    call de_alloc( scale_factors, 'scale_factors', 'REMOVE_REDUNDANT_DATA' )
    call de_alloc( tmp_FFR, 'tmp_FFR', 'REMOVE_REDUNDANT_DATA' )

  end subroutine REMOVE_REDUNDANT_DATA

end module matel_table_m
