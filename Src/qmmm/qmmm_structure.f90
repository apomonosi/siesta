module QMMM_structure
  !! This module contains the data structure for MM atoms, or, more
  !! correctly, point charges, including also the forces acting over
  !! said charges.
  use precision, only : dp

  implicit none

  private
  type mm_struct_t
    integer :: n = 0
      !! Total number of classical point charges.
    real(dp), pointer :: r(:,:) => null()
      !! XYZ positions of said charges.
    real(dp), pointer :: pc(:)  => null()
      !! The value of the partial charge in e-.
    real(dp), pointer :: f(:,:) => null()
      !! The forces acting over said charges.
    real(dp), pointer :: stress(:,:) => null()
      !! Components added to the MM stress.
    real(dp), pointer :: Vpc(:) => null()
      !! The electrostatic potential on grid generated
      !! by the MM partial charges.
  contains
    procedure :: update      => update_charges
    procedure :: receive_rpc => read_rpc
    procedure :: send_forces => write_forces
    procedure :: clean       => clean_structure
    procedure :: add_potential
    procedure :: subtract_potential
  end type mm_struct_t

  type(mm_struct_t), public :: mm_charges
    !! The main structure containing all partial charges information.

contains

  subroutine update_charges( this, na_qm, rqm, ucell )
    use alloc         , only : re_alloc, de_alloc
    use precision     , only : dp
    use QMMM_core     , only : QMMM_cutoff, doing_QMMM
    use QMMM_neighbour, only : update_QMMM_mneighb

    !! Reads input QMMM charges and updates neighbour lists.
    implicit none
    class(mm_struct_t), intent(inout) :: this
      !! The structure containing MM information.
    integer           , intent(in)    :: na_qm
      !! Number of QM atoms.
    real(dp)          , intent(in)    :: rqm(3,na_qm)
      !! Positions of QM atoms.
    real(dp)          , intent(inout) :: ucell(3,3)
      !! Periodic unit cell.

    integer  :: ntot
    real(dp) :: QMMMcut
    real(dp), pointer :: rtot(:,:)

    if ( .not. doing_QMMM() ) return

    call this%receive_rpc( )
    if ( this%n < 1 ) return

    ntot = na_qm + this%n

    nullify( rtot )
    call re_alloc( rtot, 1, 3, 1, ntot, 'rtot', 'update_charges' )

    rtot(1:3,1:na_qm)      = rqm(1:3,1:na_qm)
    rtot(1:3,na_qm+1:ntot) = this%r(1:3,1:this%n)

    QMMMcut = QMMM_cutoff( )
    call update_QMMM_mneighb( na_qm, this%n, ntot, rtot, ucell, &
                              this%pc, QMMMcut )

    call de_alloc( rtot, 'rtot', 'update_charges' )
  end subroutine update_charges

  subroutine read_rpc( this )
    !! Reads the total number of partial charges, their positions
    !! and their charge values from the fdf file. Also allocates
    !! arrays.
    !!
    !! Input coordinates must be in Bohr and charges in e-.
    use alloc         , only : de_alloc, re_alloc
    use fdf           , only : block_fdf, parsed_line, fdf_bline , fdf_block, &
                               fdf_block_linecount   , fdf_bmatch, fdf_breals
    use m_io          , only : io_getout
    use parallel      , only : IOnode
    use precision     , only : dp
    use QMMM_core     , only : QMMM_uses_external_driver
    use QMMM_interface, only : QMMM_get_partialcharges
    use sys           , only : die
    use units         , only : Ang

    implicit none
    class(mm_struct_t), intent(inout) :: this
      !! The structure containing MM information.

    integer                    :: iCrg, stdOut
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    if ( .not. associated(this%stress) ) then
      nullify( this%stress )
      call re_alloc( this%stress, 1, 3, 1, 3, 'stress' , 'read_rpc' )
    endif


    if ( associated(this%r)  ) call de_alloc( this%r , 'r' , 'read_rpc' )
    if ( associated(this%f)  ) call de_alloc( this%f , 'f' , 'read_rpc' )
    if ( associated(this%pc) ) call de_alloc( this%pc, 'pc', 'read_rpc' )
    nullify( this%r, this%f, this%pc )

    if ( QMMM_uses_external_driver() ) then
      ! Early return when using an external driver.
      call QMMM_get_partialcharges( this%n, this%r, this%pc )

    else
      ! Counts the number of partial charges found.
      this%n = 0
      if ( fdf_block( 'PartialCharges', bfdf ) ) then
        this%n = fdf_block_linecount( 'PartialCharges', 'rrrr' )
      endif

      call io_getout( stdOut )
      if ( this%n < 1 ) then
        if ( IOnode ) &
          write( stdOut,'(A)' ) 'WARNING: QM/MM calculation requested '//&
                                'but no partial charges found.'
        return
      endif
      if ( IOnode ) &
        write( stdOut,'(A31,I7,A16)' ) 'Running QM/MM calculation with ',&
                                       this%n, ' partial charges.'

      ! Allocates and initalises arrays to zero.
      call re_alloc( this%r , 1, 3, 1, this%n, 'r' , 'read_rpc' )
      call re_alloc( this%pc,       1, this%n, 'pc', 'read_rpc' )
      this%r(:,:) = 0.0_dp
      this%pc(:)  = 0.0_dp

      iCrg = 0
      do while( fdf_bline( bfdf, pline ) )
        if ( .not. fdf_bmatch( pline, 'rrrr' ) ) cycle

        iCrg = iCrg +1
        if ( iCrg > this%n ) &
          call die( 'Error while reading PartialCharges block.' )

        this%r(1,iCrg) = fdf_breals( pline, 1 ) * Ang
        this%r(2,iCrg) = fdf_breals( pline, 2 ) * Ang
        this%r(3,iCrg) = fdf_breals( pline, 3 ) * Ang
        this%pc(iCrg)  = fdf_breals( pline, 4 )
      enddo

    endif

    ! Allocates forces array.
    call re_alloc( this%f , 1, 3, 1, this%n, 'f' , 'read_rpc' )
    this%f(:,:) = 0.0_dp

  end subroutine read_rpc

  subroutine write_forces( this )
    !! Writes forces acting over partial charges to an output file.
    use files, only : slabel
    use m_io , only : io_assign, io_close
    use units, only : Ang, eV

    implicit none
    class(mm_struct_t), intent(in) :: this
      !! The structure containing MM information.

    integer :: funit, iCrg

    if ( this%n < 1 ) return

    call io_assign( funit )
    open( unit = funit, file = trim(slabel)//'.FAPC' )

    do iCrg = 1, this%n
      write( funit, '(F14.7,F14.7,F14.7)') &
        this%f(1,iCrg) * Ang / eV, this%f(2,iCrg) * Ang / eV, &
        this%f(3,iCrg) * Ang / eV
    enddo

    call io_close( funit )
  end subroutine write_forces

  subroutine add_potential( this, ntpl, nspin, V_scf )
    !! Adds the partial charges grid potential to a potential
    !! depending on both the grid and spin.
    use precision, only : dp

    implicit none
    class(mm_struct_t), intent(in)    :: this
      !! The structure containing MM information.
    integer           , intent(in)    :: ntpl
      !! Total number of points in grid.
    integer           , intent(in)    :: nspin
      !! Number of spin coordinates.
    real(dp)          , intent(inout) :: V_scf(ntpl, nspin)
      !! The grid potential over which the partial charges potential
      !! is added.

    integer :: ip, ispin

    if ( this%n < 1 ) return
!$OMP parallel default(shared), private(ip,ispin)
    do ispin = 1, nspin
!$OMP do
      do ip = 1, ntpl
        V_scf(ip,ispin) = V_scf(ip,ispin) + this%Vpc(ip)
      enddo
!$OMP end do
    enddo
!$OMP end parallel
  end subroutine add_potential

  subroutine subtract_potential( this, ntpl, nspin, V_scf )
    !! Subtract the partial charges grid potential from a potential
    !! depending on both the grid and spin.
    use precision, only : dp

    implicit none
    class(mm_struct_t), intent(in)    :: this
      !! The structure containing MM information.
    integer           , intent(in)    :: ntpl
      !! Total number of points in grid.
    integer           , intent(in)    :: nspin
      !! Number of spin coordinates.
    real(dp)          , intent(inout) :: V_scf(ntpl, nspin)
      !! The grid potential from which the partial charges potential
      !! is subtracted.

    integer :: ip, ispin

    if ( this%n < 1 ) return
!$OMP parallel default(shared), private(ip,ispin)
    do ispin = 1, nspin
!$OMP do
      do ip = 1, ntpl
        V_scf(ip,ispin) = V_scf(ip,ispin) - this%Vpc(ip)
      enddo
!$OMP end do
    enddo
!$OMP end parallel
  end subroutine subtract_potential

  subroutine clean_structure( this )
    !! Clears all pointers within the structure.
    use alloc, only : de_alloc

    implicit none
    class(mm_struct_t), intent(inout) :: this
      !! The structure containing MM information.

    if ( associated(this%r)   ) &
      call de_alloc( this%r  , 'r'  , 'clean_structure' )
    if ( associated(this%f)   ) &
      call de_alloc( this%f  , 'f'  , 'clean_structure' )
    if ( associated(this%pc)  ) &
      call de_alloc( this%pc , 'pc' , 'clean_structure' )
    if ( associated(this%Vpc) ) &
      call de_alloc( this%Vpc, 'Vpc', 'clean_structure' )

  end subroutine clean_structure
end module QMMM_structure
