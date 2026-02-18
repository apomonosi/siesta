! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module hsx_m

!
! Experimental module to process HSX files
!
implicit none

integer, parameter, private :: dp = selected_real_kind(14,100)
integer, parameter, private :: sp = selected_real_kind(6,30)

public  :: read_hsx_file, write_hsx_file
public  :: read_hs_file, write_hs_file

! Set derived type hsx_t to hold info of HSX file, containing:
!   nspecies                : number of chemical species
!   na_u                    : number of atoms in unit cell
!   no_u                    : number of orbitals in unit cell
!   no_s                    : number of orbitals in supercell
!   nspin                   : number of spin components
!   nh                      : dimension of arrays hamilt, Sover, and xij
!   gamma                   : was this a gamma-only calculation?
!   has_xij                 : does the file contain xij vectors?
!   no(nspecies)            : number of atomic orbitals of each species
!   nquant(nspecies,naoatx) : principal quantum number of each atomic orbital
!                             with naoatx=max(no)
!   lquant(nspecies,naoatx) : ang. momentum number of each atomic orbital
!   zeta(nspecies,naoatx)   : zeta-index of each atomic orbital
!   iaorb(no_u)             : atom to which each orbital belongs
!   iphorb(no_u)            : index of each orbital within its atom
!   label(nspecies)         : atomic label (symbol) of each species
!   numh(no_u)              : num of nonzero elements in each row of hamiltonian
!   listhptr(no_u)          : row-start index in sparse-matrix arrays
!   listh(nh)               : orbital index of nonzero matrix elements
!   indxuo(no_s)            : index of equivalent orbital in first unit cell
!   hamilt(nh,nspin)        : hamiltonian matrix elements in sparse format
!   Sover(nh)               : overlap matrix elements in sparse format
!   xij(3,nh)               : vector between each pair of connected orbitals
!   isa(na_u)               : species index of each atom
!   zval(nspecies)          : atomic number of each species
! To transform from sparse to full format, for a given k point:
!   S(1:no_u,1:no_u) = 0                       ! full complex overlap matrix
!   H(1:no_u,1:no_u,1:nspin) = 0               ! full complex hamiltonian
!   do io = 1,no_u                             ! loop on unit cell orbitals
!     do j = 1,numh(io)                        ! loop on connected orbitals
!        ij = listhptr(io)+j                   ! sparse-matrix array index
!        jos = listh(ij)                       ! index of connected orbital
!        jo = indxuo(jos)                      ! equiv. orbital in unit cell
!        phase = exp(-i*sum(k(:)*xij(:,ij)))   ! phase factor between orbs.
!        H(jo,io,:) += phase*hamilt(ij,:)      ! hamiltonian matrix element
!                           
!        S(jo,io) += phase*Sover(ij)           ! overlap matrix element
!     enddo
!   enddo
! Notice that io,jo are within unit cell, and jos is within supercell

type, public :: hsx_t
  integer :: nspecies
  integer :: na_u
  integer :: no_u
  integer :: nsc(3)
  integer :: no_s
  integer :: nspin
  integer :: nh
  logical :: gamma
  logical :: has_xij = .false.
  real(dp) :: ucell(3,3) = 0._dp
  integer, pointer :: no(:) => null()
  integer, pointer :: nquant(:,:) => null()
  integer, pointer :: lquant(:,:) => null()
  integer, pointer :: zeta(:,:) => null()
  integer, pointer :: iaorb(:) => null()
  integer, pointer :: iphorb(:) => null()
  character(len=20), pointer :: label(:) => null()
  integer, pointer  :: numh(:) => null() 
  integer, pointer  :: listhptr(:) => null()
  integer, pointer  :: listh(:) => null()  
  integer, pointer  :: indxuo(:) => null()
  real(dp), pointer :: xa(:,:) => null()
  real(dp), pointer :: hamilt(:,:) => null()
  real(dp), pointer :: Sover(:) => null()
  real(dp), pointer :: xij(:,:) => null()
  integer, pointer :: isc_off(:,:) => null()
  integer, pointer  :: isa(:) => null()
  real(dp), pointer :: zval(:) => null()
  real(dp)          :: Ef=0._dp, qtot=0._dp, temp=0._dp

  ! K-point information for the run, only present in version >=2
  real(dp) :: kdispl(3)
  integer :: kcell(3,3) = 0
  integer :: version = 0
end type

private

CONTAINS

  subroutine read_hsx_file(hsx, fname)
    type(hsx_t), intent(out)  :: hsx
    character(len=*), intent(in) :: fname

    integer :: version

    version = HSX_version(fname)
    if ( version == 0 ) then
      call read_hsx_file_version0(hsx, fname)
    else if ( version == 1 .or. version == 2 ) then
      call read_hsx_file_version1_2(hsx, fname)
    else
      STOP "unknown HSX file version [0, 1, 2]"
    end if
  end subroutine read_hsx_file

  function HSX_version(fname) result(version)
    character(len=*), intent(in) :: fname
    integer :: version
    integer :: iu
    integer :: na_u, no_u, no_s, nspin, n_nzs, err

    ! Initialize
    version = 0

    call get_unit_number(iu)

    ! Open file
    open( iu, file=fname, form='unformatted', status='unknown' )

    read(iu,iostat=err) na_u, no_s, nspin, n_nzs
    if ( err == 0 ) then
       ! we can successfully read 4 integers
       version = 0
    else
       rewind(iu)
       read(iu,iostat=err) version
    end if

    close(iu)

  end function HSX_version


  subroutine read_hsx_file_version1_2(hsx, fname)
    type(hsx_t), intent(out)  :: hsx
    character(len=*), intent(in) :: fname

    integer :: hs_u, iostat
    integer :: io, jo, is, ind, ia, ja
    integer :: max_atom_orbs
    logical :: is_dp
    integer, allocatable :: lasto(:)
    real(sp), allocatable :: rbuf(:)

    call get_unit_number(hs_u)

    write(6,"(1x,a)",advance='no') trim(fname)
    open(hs_u,file=trim(fname),status='old',form='unformatted')

    ! skip version
    read(hs_u, iostat=iostat) hsx%version
    if (iostat /= 0) STOP "version"

    if ( .not. any(hsx%version == [1, 2]) ) then
      STOP "incorrect call [version/=1|2]"
    end if

    ! now read whether this is double precision or not
    read(hs_u, iostat=iostat) is_dp
    if ( iostat /= 0 ) STOP "is_dp"

    read(hs_u,iostat=iostat) hsx%na_u, hsx%no_u, hsx%nspin, hsx%nspecies, hsx%nsc
    if ( iostat /= 0 ) STOP "geometry dimensions"
    hsx%no_s = hsx%no_u * product(hsx%nsc)
    hsx%gamma = hsx%no_s == hsx%no_u

    read(hs_u,iostat=iostat) hsx%ucell, hsx%Ef, hsx%qtot, hsx%temp

    ! Allocate the arrays for the atomic information
    ! First local arrays that won't be shared
    allocate(hsx%xa(3,hsx%na_u), lasto(0:hsx%na_u), hsx%isc_off(3, product(hsx%nsc)))
    lasto(0) = 0
    ! allocate shared variables
    allocate(hsx%isa(hsx%na_u), hsx%iaorb(hsx%no_u), hsx%iphorb(hsx%no_u))
    allocate(hsx%label(hsx%nspecies), hsx%zval(hsx%nspecies), hsx%no(hsx%nspecies))

    ! Now read data
    read(hs_u,iostat=iostat) hsx%isc_off, hsx%xa, hsx%isa, lasto(1:hsx%na_u)
    if ( iostat /= 0 ) STOP "geometry coordinates"
    read(hs_u,iostat=iostat) (hsx%label(is), hsx%zval(is), hsx%no(is), is=1,hsx%nspecies)
    if ( iostat /= 0 ) STOP "species information"

    ! Allocate remaning orbital information
    max_atom_orbs = maxval(hsx%no)
    allocate(hsx%nquant(hsx%nspecies,max_atom_orbs), hsx%lquant(hsx%nspecies,max_atom_orbs))
    allocate(hsx%zeta(hsx%nspecies,max_atom_orbs))
    do is = 1, hsx%nspecies
      read(hs_u,iostat=iostat) (hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io), io=1,hsx%no(is))
      if ( iostat /= 0 ) STOP "specific specie information"
    end do

    if ( hsx%version == 2 ) then
       read(hs_u,iostat=iostat) hsx%kcell, hsx%kdispl
       if ( iostat /= 0 ) STOP "k-information"
    end if

    ! Populate iaorb and iphorb
    ind = 0
    do ia = 1, hsx%na_u
      do io = 1, lasto(ia) - lasto(ia-1)
        ind = ind + 1
        hsx%iaorb(ind) = ia
        hsx%iphorb(ind) = io
      end do
    end do

    ! Now create indxuo (siesta always produces the same order)
    allocate(hsx%indxuo(hsx%no_s))
    ind = 0
    do is = 1 , product(hsx%nsc)
      do io = 1, hsx%no_u
        ind = ind + 1
        hsx%indxuo(ind) = io
      end do
    end do

    allocate(hsx%numh(hsx%no_u))
    read(hs_u,iostat=iostat) hsx%numh
    if (iostat /= 0) STOP "numh(io)"

    ! Create pointer
    allocate(hsx%listhptr(hsx%no_u))
    hsx%listhptr(1) = 0
    do io=2,hsx%no_u
      hsx%listhptr(io) = hsx%listhptr(io-1) + hsx%numh(io-1)
    end do

    ! Now read sparse data
    hsx%nh = hsx%listhptr(hsx%no_u) + hsx%numh(hsx%no_u)
    allocate(hsx%listh(hsx%nh))

    do io=1,hsx%no_u
      read(hs_u,iostat=iostat) hsx%listh(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io))
      if (iostat /= 0) STOP "listh"
    end do


    ! Now we need to re-create the xij array
    allocate(hsx%xij(3,hsx%nh))

    ! Create xij and dij
    do ia = 1 , hsx%na_u
      do io = lasto(ia-1) + 1, lasto(ia)
        do ind = hsx%listhptr(io) + 1, hsx%listhptr(io) + hsx%numh(io)
          is = (hsx%listh(ind) - 1) / hsx%no_u + 1
          ja = hsx%iaorb(ucorb(hsx%listh(ind), hsx%no_u))

          hsx%xij(:, ind) = hsx%xa(:,ja) - hsx%xa(:,ia) + &
              hsx%isc_off(1,is) * hsx%ucell(:,1) + &
              hsx%isc_off(2,is) * hsx%ucell(:,2) + &
              hsx%isc_off(3,is) * hsx%ucell(:,3)

        end do
      end do
    end do

    ! Clean-up
    deallocate(lasto)

    ! Read H and S
    allocate(hsx%hamilt(hsx%nh,hsx%nspin))
    allocate(hsx%Sover(hsx%nh))

    if ( is_dp ) then
      do is = 1, hsx%nspin
        do io = 1, hsx%no_u
          read(hs_u,iostat=iostat) hsx%hamilt(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io),is)
          if (iostat /= 0) STOP "H(dp)"
        end do
      end do

      do io = 1, hsx%no_u
        read(hs_u,iostat=iostat) hsx%Sover(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io))
        if (iostat /= 0) STOP "S(dp)"
      end do

    else
      allocate(rbuf(maxval(hsx%numh)))

      do is = 1, hsx%nspin
        do io = 1, hsx%no_u
          read(hs_u,iostat=iostat) rbuf(1:hsx%numh(io))
          if (iostat /= 0) STOP "H(sp)"
          hsx%hamilt(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io),is) = rbuf(1:hsx%numh(io))
        end do
      end do

      do io = 1, hsx%no_u
        read(hs_u,iostat=iostat) rbuf(1:hsx%numh(io))
        if (iostat /= 0) STOP "S(sp)"
        hsx%Sover(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io)) = rbuf(1:hsx%numh(io))
      end do

      deallocate(rbuf)
    end if

    close(hs_u)

  contains

    elemental function UCORB(a,p)
      integer, intent(in) :: a,p
      integer :: UCORB
      UCORB = MOD(a-1,p) + 1
    end function

  end subroutine read_hsx_file_version1_2


  subroutine read_hsx_file_version0(hsx,fname)
    type(hsx_t), intent(out)  :: hsx
    character(len=*), intent(in) :: fname

    !
    ! Reads HSX file "fname" and stores the info in the hsx data structure
    ! (Real arrays are stored in double precision)

    integer, allocatable  :: ibuff(:)
    real(sp), allocatable  :: hbuff(:)
    real(sp), allocatable  :: buff3(:,:)

    integer numx, ind, no_u, nnz, na_u, nspecies, nspin, nh, i
    integer :: im, is, hsx_u, ia, io, iostat, k, naoatx, no_s
    logical  :: debug = .false.

    call get_unit_number(hsx_u)
    print *, "Using unit: ", hsx_u

    hsx%version = 0

    open(hsx_u,file=trim(fname),status='old',form='unformatted')

    read(hsx_u,iostat=iostat) hsx%no_u, hsx%no_s, hsx%nspin, hsx%nh
    if (iostat /= 0) STOP "nnao, no_s..."

    no_u = hsx%no_u
    no_s = hsx%no_s

    read(hsx_u,iostat=iostat) hsx%gamma
    if (iostat /= 0) STOP "gamma"
    IF (DEBUG) PRINT *, "GAMMA=", hsx%gamma
    if (.not. hsx%gamma) then
      allocate(hsx%indxuo(no_s))
      read(hsx_u) (hsx%indxuo(i),i=1,hsx%no_s)
    else
      allocate(hsx%indxuo(hsx%no_u))
      do i=1,hsx%no_u
        hsx%indxuo(i) = i
      enddo
    endif

    nh  = hsx%nh
    nspin = hsx%nspin
    print *, "nh: ", nh
    allocate (hsx%numh(no_u), hsx%listhptr(no_u), hsx%listh(nh))

    allocate (hsx%xij(3,nh),stat=iostat)
    allocate (hsx%hamilt(nh,nspin),stat=iostat)
    allocate (hsx%Sover(nh),stat=iostat)

    read(hsx_u,iostat=iostat) (hsx%numh(io), io=1,no_u)      
    if (iostat /= 0) STOP "numh"

    numx = maxval(hsx%numh(1:no_u))
    allocate(ibuff(numx), hbuff(numx), buff3(3,numx))

    nnz = sum(hsx%numh(1:hsx%no_u))
    if (nnz > nh) STOP "nh overflow in HS"
    ! Create listhptr 
    hsx%listhptr(1)=0
    do io=2,hsx%no_u
      hsx%listhptr(io)=hsx%listhptr(io-1)+hsx%numh(io-1)
    enddo

    do io=1,hsx%no_u
      read(hsx_u,iostat=iostat) (ibuff(im), im=1,hsx%numh(io))
      if (iostat /= 0) STOP "listh"
      do im=1,hsx%numh(io)
        hsx%listh(hsx%listhptr(io)+im) = ibuff(im)
      enddo
    enddo

    do is=1,hsx%nspin
      do io=1,hsx%no_u
        read(hsx_u,iostat=iostat) (hbuff(im), im=1,hsx%numh(io))
        if (iostat /= 0) STOP "Hamilt"
        do im=1,hsx%numh(io)
          hsx%hamilt(hsx%listhptr(io)+im,is) = hbuff(im)
          if (debug) print *, "Hamilt ", io, im, hbuff(im)
        enddo
      enddo
    enddo
    !
    !       Read overlap matrix
    !
    do io=1,hsx%no_u
      read(hsx_u,iostat=iostat) (hbuff(im), im=1,hsx%numh(io))
      if (iostat /= 0) STOP "Overlap matrix read error"
      do im=1,hsx%numh(io)
        hsx%Sover(hsx%listhptr(io)+im) = hbuff(im)
        if (debug) print *, "S ", io, im, hbuff(im)
      enddo
    enddo

    read(hsx_u,iostat=iostat) hsx%qtot, hsx%temp           ! fossils
    if (iostat /= 0) STOP "Qtot, temp, read error"

    !
    !        Always read xijk
    !
    do io=1,hsx%no_u
      read(hsx_u,iostat=iostat) ((buff3(k,im), k=1,3), im=1,hsx%numh(io))
      if (iostat /= 0) STOP "xij(k) read error"
      do im=1,hsx%numh(io)
        ind = hsx%listhptr(io)+im
        if (debug) print *, "xijk ", buff3(:,im)
        hsx%xij(1:3,ind) = buff3(1:3,im)
      enddo
    enddo
    hsx%has_xij = .true.

    !
    !        Read auxiliary info
    !
    read(hsx_u) hsx%nspecies
    nspecies = hsx%nspecies
    print *, "nspecies: ", nspecies
    allocate(hsx%label(nspecies), hsx%zval(nspecies), hsx%no(nspecies))
    read(hsx_u) (hsx%label(is),hsx%zval(is),hsx%no(is), is=1,nspecies)
    naoatx = maxval(hsx%no(1:nspecies))
    allocate (hsx%nquant(nspecies,naoatx), hsx%lquant(nspecies,naoatx), &
        hsx%zeta(nspecies,naoatx))
    do is=1, nspecies
      do io=1, hsx%no(is)
        read(hsx_u) hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io)
      enddo
    enddo
    read(hsx_u) hsx%na_u
    na_u = hsx%na_u
    allocate(hsx%isa(na_u))
    allocate(hsx%iaorb(no_u), hsx%iphorb(no_u))
    read(hsx_u) (hsx%isa(ia), ia=1,na_u)
    read(hsx_u) (hsx%iaorb(io), hsx%iphorb(io), io=1,no_u)

    close(hsx_u)
    deallocate(ibuff, hbuff, buff3)

  end subroutine read_hsx_file_version0


!--------------------------------------------------------------

  subroutine write_hsx_file(hsx,fname, version)
    type(hsx_t), intent(in)  :: hsx
    character(len=*), intent(in) :: fname
    integer, intent(in), optional :: version

    if ( present(version) ) then
      if ( version == 0 ) then
        call write_hsx_file_version0(hsx,fname)
      else if ( version == 1 ) then
        call write_hsx_file_version1(hsx,fname)
      else if ( version == 2 ) then
        call write_hsx_file_version2(hsx,fname)
      else
        print *, 'Unknown version specifier for HSX file [0, 1, 2]: ', version
        stop 'Unknown version specifier [0, 1, 2]?'
      end if
    else if ( any(hsx%kcell > 0) ) then
      ! Check for variables only in 2
      call write_hsx_file_version2(hsx,fname)
    else if ( associated(hsx%isc_off) ) then
      ! Check for variables only in 1
      call write_hsx_file_version1(hsx,fname)
    else
      call write_hsx_file_version0(hsx,fname)
    end if

  end subroutine write_hsx_file

  subroutine write_hsx_file_version2(hsx, fname)
    type(hsx_t), intent(in)  :: hsx
    character(len=*), intent(in) :: fname

    integer :: iu, ia, io, ind, is
    integer, allocatable :: lasto(:)

    call get_unit_number(iu)

    ! Open file
    open( iu, file=trim(fname), form='unformatted', status='unknown' )

    ! Write version specification (to easily distinguish between different versions)
    write(iu) 1
    ! And what precision
    write(iu) .true.

    ! Write overall data
    write(iu) hsx%na_u, hsx%no_u, hsx%nspin, hsx%nspecies, hsx%nsc
    write(iu) hsx%ucell, hsx%Ef, hsx%qtot, hsx%temp

    ! Recreate lasto for easier storage
    allocate(lasto(hsx%na_u))
    lasto(1) = hsx%no(hsx%isa(1))
    do ia = 2, hsx%na_u
      lasto(ia) = lasto(ia-1) + hsx%no(hsx%isa(ia))
    end do
    write(iu) hsx%isc_off, hsx%xa, hsx%isa, lasto
    deallocate(lasto)

    ! Write other useful info
    write(iu) (hsx%label(is),hsx%zval(is),hsx%no(is),is=1,hsx%nspecies)
    do is = 1, hsx%nspecies
      write(iu) (hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io),io=1,hsx%no(is))
    end do

    ! This is the only difference from version == 1
    write(iu) hsx%kcell, hsx%kdispl

    write(iu) hsx%numh

    do io = 1, hsx%no_u
      write(iu) hsx%listh(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io))
    end do

    do is = 1, hsx%nspin
      do io = 1, hsx%no_u
        write(iu) hsx%hamilt(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io),is)
      end do
    end do

    do io = 1, hsx%no_u
      write(iu) hsx%Sover(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io))
    end do

    close( iu )

  end subroutine write_hsx_file_version2


  subroutine write_hsx_file_version1(hsx, fname)
    type(hsx_t), intent(in)  :: hsx
    character(len=*), intent(in) :: fname

    integer :: iu, ia, io, ind, is
    integer, allocatable :: lasto(:)

    call get_unit_number(iu)

    ! Open file
    open( iu, file=trim(fname), form='unformatted', status='unknown' )

    ! Write version specification (to easily distinguish between different versions)
    write(iu) 1
    ! And what precision
    write(iu) .true.

    ! Write overall data
    write(iu) hsx%na_u, hsx%no_u, hsx%nspin, hsx%nspecies, hsx%nsc
    write(iu) hsx%ucell, hsx%Ef, hsx%qtot, hsx%temp

    ! Recreate lasto for easier storage
    allocate(lasto(hsx%na_u))
    lasto(1) = hsx%no(hsx%isa(1))
    do ia = 2, hsx%na_u
      lasto(ia) = lasto(ia-1) + hsx%no(hsx%isa(ia))
    end do
    write(iu) hsx%isc_off, hsx%xa, hsx%isa, lasto
    deallocate(lasto)

    ! Write other useful info
    write(iu) (hsx%label(is),hsx%zval(is),hsx%no(is),is=1,hsx%nspecies)
    do is = 1, hsx%nspecies
      write(iu) (hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io),io=1,hsx%no(is))
    end do

    write(iu) hsx%numh

    do io = 1, hsx%no_u
      write(iu) hsx%listh(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io))
    end do

    do is = 1, hsx%nspin
      do io = 1, hsx%no_u
        write(iu) hsx%hamilt(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io),is)
      end do
    end do

    do io = 1, hsx%no_u
      write(iu) hsx%Sover(hsx%listhptr(io)+1:hsx%listhptr(io)+hsx%numh(io))
    end do

    close( iu )

  end subroutine write_hsx_file_version1

  subroutine write_hsx_file_version0(hsx,fname)
    type(hsx_t), intent(in)  :: hsx
    character(len=*), intent(in) :: fname

    !
    ! Writes HSX file "fname" 

    integer :: no_u, nnz, na_u, nspecies, nspin, nh, i
    integer :: im, is, hsx_u, ia, io, iostat, k, no_s

    if ( .not. hsx%has_xij) then
      print *, "Cannot generate an HSX file without Xij information"
      STOP
    endif
    call get_unit_number(hsx_u)
    print *, "Using unit: ", hsx_u

    open(hsx_u,file=trim(fname),status='unknown',form='unformatted')

    write(hsx_u,iostat=iostat) hsx%no_u, hsx%nspin, hsx%nsc
    if (iostat /= 0) STOP "nnao, no_s..."

    no_u = hsx%no_u
    no_s = hsx%no_s

    write(hsx_u,iostat=iostat) hsx%gamma
    if (.not. hsx%gamma) then
      write(hsx_u) (hsx%indxuo(i),i=1,hsx%no_s)
    endif

    nh  = hsx%nh
    nspin = hsx%nspin
    print *, "nh: ", nh

    write(hsx_u,iostat=iostat) (hsx%numh(io), io=1,no_u)      
    if (iostat /= 0) STOP "numh"

    nnz = sum(hsx%numh(1:hsx%no_u))
    if (nnz /= nh) STOP "nnz /= nh"

    do io=1,hsx%no_u
      write(hsx_u,iostat=iostat)  &
          (hsx%listh(hsx%listhptr(io)+im), im=1,hsx%numh(io))
    enddo

    do is=1,hsx%nspin
      do io=1,hsx%no_u
        write(hsx_u,iostat=iostat)   &
            (real(hsx%hamilt(hsx%listhptr(io)+im,is),kind=sp),   &
            im=1,hsx%numh(io))
      enddo
    enddo
    !
    !   overlap matrix
    !
    do io=1,hsx%no_u
      write(hsx_u,iostat=iostat)    &
          (real(hsx%Sover(hsx%listhptr(io)+im),kind=sp),im=1,hsx%numh(io))
    enddo

    write(hsx_u,iostat=iostat) hsx%qtot, hsx%temp           ! fossils
    !
    !        Always write xijk
    !
    do io=1,hsx%no_u
      write(hsx_u,iostat=iostat)    &
          ((real(hsx%xij(k,hsx%listhptr(io)+im),kind=sp), &
          k=1,3), im=1,hsx%numh(io))
    enddo
    !
    !        Write auxiliary info
    !
    write(hsx_u) hsx%nspecies
    nspecies = hsx%nspecies
    print *, "nspecies: ", nspecies
    write(hsx_u) (hsx%label(is),hsx%zval(is),hsx%no(is), is=1,nspecies)

    do is=1, nspecies
      do io=1, hsx%no(is)
        write(hsx_u) hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io)
      enddo
    enddo
    write(hsx_u) hsx%na_u

    na_u = hsx%na_u

    write(hsx_u) (hsx%isa(ia), ia=1,na_u)
    write(hsx_u) (hsx%iaorb(io), hsx%iphorb(io), io=1,no_u)

    close(hsx_u)

  end subroutine write_hsx_file_version0

!----------------------------------------------------------
subroutine get_unit_number(lun)
integer, intent(out) :: lun

logical :: used
integer :: iostat

do lun= 10, 99
   inquire(unit=lun, opened=used, iostat=iostat)
   if (iostat .ne. 0) used = .true.
   if (.not. used) return
enddo
STOP "Cannot get unit"
end subroutine get_unit_number


subroutine read_hs_file(hsx,fname)
type(hsx_t), intent(out)  :: hsx
character(len=*), intent(in) :: fname

!
! Reads HS file "fname" and stores the info in the hsx data structure

  integer, allocatable  :: ibuff(:)
  real(dp), allocatable  :: hbuff(:)
  real(dp), allocatable  :: buff3(:,:)

  integer numx, ind, no_u, nnz, na_u, nspecies, nspin, nh, i
  integer :: im, is, hs_u, ia, io, iostat, k, naoatx, no_s
  logical  :: debug = .false.

  call get_unit_number(hs_u)
  print *, "Reading ", trim(fname), " using unit: ", hs_u

  open(hs_u,file=trim(fname),status='old',form='unformatted')

  read(hs_u,iostat=iostat) hsx%no_u, hsx%no_s, hsx%nspin, hsx%nh
  if (iostat /= 0) STOP "nnao, no_s..."

  no_u = hsx%no_u
  no_s = hsx%no_s

  read(hs_u,iostat=iostat) hsx%gamma
  if (iostat /= 0) STOP "gamma"
  IF (DEBUG) PRINT *, "GAMMA=", hsx%gamma
  if (.not. hsx%gamma) then
     allocate(hsx%indxuo(no_s))
     read(hs_u) (hsx%indxuo(i),i=1,hsx%no_s)
  else
     allocate(hsx%indxuo(hsx%no_u))
     do i=1,hsx%no_u
        hsx%indxuo(i) = i
     enddo
  endif

  nh  = hsx%nh
  nspin = hsx%nspin
  print *, "nh: ", nh
  allocate (hsx%numh(no_u), hsx%listhptr(no_u), hsx%listh(nh))

       allocate (hsx%hamilt(nh,nspin),stat=iostat)
       allocate (hsx%Sover(nh),stat=iostat)

       do io = 1,no_u
          read(hs_u,iostat=iostat) hsx%numh(io)
          if (iostat /= 0) STOP "numh"
       enddo

  numx = maxval(hsx%numh(1:no_u))
  allocate(ibuff(numx), hbuff(numx), buff3(3,numx))

  nnz = sum(hsx%numh(1:hsx%no_u))
   if (nnz > nh) STOP "nh overflow in HS"
  ! Create listhptr 
  hsx%listhptr(1)=0
  do io=2,hsx%no_u
     hsx%listhptr(io)=hsx%listhptr(io-1)+hsx%numh(io-1)
  enddo

  do io=1,hsx%no_u
     do im = 1, hsx%numh(io)
        read(hs_u,iostat=iostat) ibuff(im)
        if (iostat /= 0) STOP "listh"
     enddo
     do im=1,hsx%numh(io)
        hsx%listh(hsx%listhptr(io)+im) = ibuff(im)
     enddo
  enddo

  do is=1,hsx%nspin
     do io=1,hsx%no_u
        do im = 1, hsx%numh(io)
           read(hs_u,iostat=iostat) hbuff(im)
           if (iostat /= 0) STOP "Hamilt"
        enddo
        do im=1,hsx%numh(io)
           hsx%hamilt(hsx%listhptr(io)+im,is) = hbuff(im)
           if (debug) print *, "Hamilt ", io, im, hbuff(im)
        enddo
     enddo
  enddo
  !
  !       Read overlap matrix
  !
  do io=1,hsx%no_u
     do im = 1, hsx%numh(io)
        read(hs_u,iostat=iostat) hbuff(im)
        if (iostat /= 0) STOP "Overlap matrix read error"
     enddo
     do im=1,hsx%numh(io)
        hsx%Sover(hsx%listhptr(io)+im) = hbuff(im)
        if (debug) print *, "S ", io, im, hbuff(im)
     enddo
  enddo

  read(hs_u,iostat=iostat) hsx%qtot, hsx%temp           ! fossils
  if (iostat /= 0) STOP "Qtot, temp, read error"

  !
  !   read xijk if not gamma
  !
  if (.not. hsx%gamma) then
    allocate (hsx%xij(3,nh),stat=iostat)
    do io=1,hsx%no_u
     do im = 1, hsx%numh(io)
        read(hs_u,iostat=iostat) (buff3(k,im), k=1,3)
        if (iostat /= 0) STOP "xij(k) read error"
     enddo
     do im=1,hsx%numh(io)
        ind = hsx%listhptr(io)+im
        if (debug) print *, "xijk ", buff3(:,im)
        hsx%xij(1:3,ind) = buff3(1:3,im)
     enddo
    enddo
    hsx%has_xij = .true.
  else
    print *, "Cannot read Xij from old HS file if gamma-only ..."
    print *, "... will mark the data structure as missing Xij"
    print *, "... in the future you might be able to reconstruct it"
    hsx%has_xij = .false.
  endif
  !
  !        Read auxiliary info
  !
  read(hs_u) hsx%nspecies
  nspecies = hsx%nspecies
  print *, "nspecies: ", nspecies
  allocate(hsx%label(nspecies), hsx%zval(nspecies), hsx%no(nspecies))
  read(hs_u) (hsx%label(is),hsx%zval(is),hsx%no(is), is=1,nspecies)
  naoatx = maxval(hsx%no(1:nspecies))
  allocate (hsx%nquant(nspecies,naoatx), hsx%lquant(nspecies,naoatx), &
       hsx%zeta(nspecies,naoatx))
  do is=1, nspecies
     do io=1, hsx%no(is)
        read(hs_u) hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io)
     enddo
  enddo
  read(hs_u) hsx%na_u
  na_u = hsx%na_u
  allocate(hsx%isa(na_u))
  allocate(hsx%iaorb(no_u), hsx%iphorb(no_u))
  read(hs_u) (hsx%isa(ia), ia=1,na_u)
  read(hs_u) (hsx%iaorb(io), hsx%iphorb(io), io=1,no_u)

  close(hs_u)
  deallocate(ibuff, hbuff, buff3)

end subroutine read_hs_file

subroutine write_hs_old( filename, gamma, no_u, no_s, nspin, indxuo,     &
                   maxnh, numh, listhptr, listh, H, S, qtot, temp, &
                   xij, nspecies, iphorb, iaorb, na_u, isa, no,    &
                   label, zval, nquant, lquant, zeta)
!
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible

! This routine is deprecated. New version iohsx is more compact.
! *************************** INPUT **********************************
! logical       gamma         : Is only gamma point used?
! ******************** INPUT 
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer nspin               : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
!
  implicit          none

  integer, parameter ::  dp = selected_real_kind(14,100)

!
  character(len=*), intent(in) :: filename

  logical, intent(in) ::  gamma
  integer, intent(in) ::  maxnh, no_u, no_s, nspin
  integer, intent(in) ::  indxuo(no_s), listh(maxnh),  &
                          numh(*), listhptr(*)
  real(dp), intent(in)  ::  H(maxnh,nspin), S(maxnh),  &
                            qtot, temp, xij(3,maxnh)
!
  integer, intent(in) ::   nspecies
  integer, intent(in) ::   iphorb(:), iaorb(:)
  integer, intent(in) ::   na_u, isa(:)
  integer, intent(in) ::   no(:)
  character(len=20), intent(in)  :: label(:)
  real(dp), intent(in)  :: zval(:)

  integer, intent(in) :: nquant(:,:) 
  integer, intent(in) :: lquant(:,:) 
  integer, intent(in) :: zeta(:,:) 

! Internal variables and arrays
  integer    im, is, iu, ju, k, ns, ia, io
  integer    ih,hl
!
  logical  write_xijk
!
!---------------------------------------
  iu = 3
  open( iu, file=filename, form='unformatted', status='unknown' )      
  write(iu) no_u, no_s, nspin, maxnh
  write(iu) gamma
!
! Write out indxuo
  if (.not.gamma) then
     write(iu) (indxuo(ih),ih=1,no_s)
  endif

  do ih = 1,no_u
     write(iu) numh(ih)
  enddo
!
! Write listh
      do ih = 1,no_u
         hl = ih
         do im = 1,numh(hl)
            write(iu) listh(listhptr(hl)+im)
         enddo
      enddo
!
! Write Hamiltonian
      do is=1,nspin
         do ih=1,no_u
            hl = ih
            do im=1,numh(hl)
               write(iu) H(listhptr(hl)+im,is)
            enddo
         enddo
      enddo
!
! Write Overlap matrix
      do ih = 1,no_u
         hl = ih
         do im = 1,numh(hl)
            write(iu) S(listhptr(hl)+im)
         enddo
      enddo
!
      write(iu) qtot,temp
!
      ! Write xij if requested
      write_xijk = .not. gamma
!
      if (write_xijk) then
         do ih = 1,no_u
            hl = ih
            do im = 1,numh(hl)
               write(iu) (xij(k,listhptr(hl)+im),k=1,3)
            enddo
         enddo
      endif   ! write_xijk
!

!       Write other useful info
!
      write(iu) nspecies
      write(iu) (label(is), zval(is),no(is),is=1,nspecies)
      do is = 1, nspecies
         do io=1,no(is)
            write(iu) nquant(is,io), lquant(is,io), zeta(is,io)
         enddo
      enddo
      write(iu) na_u
      write(iu) (isa(ia),ia=1,na_u)
      write(iu) (iaorb(io), iphorb(io), io=1,no_u)
!
! Close file
      close(iu)

end subroutine write_hs_old

!-----------------------------------------------------------

subroutine write_hs_file( h, filename)
type(hsx_t), intent(in) :: h
character(len=*), intent(in) :: filename


! Internal variables and arrays
  integer    im, is, iu, ju, k, ns, ia, io, no_u
  integer    ih,hl, nspecies
!
  logical  write_xijk
!
!---------------------------------------
  iu = 3
  open( iu, file=filename, form='unformatted', status='unknown' )      
  write(iu) h%no_u, h%no_s, h%nspin, h%nh
  write(iu) h%gamma
!
! Write out indxuo
  if (.not. h%gamma) then
     write(iu) (h%indxuo(ih),ih=1,h%no_s)
  endif

  no_u = h%no_u
  do ih = 1,no_u
     write(iu) h%numh(ih)
  enddo
!
! Write listh
      do ih = 1,no_u
         hl = ih
         do im = 1,h%numh(hl)
            write(iu) h%listh(h%listhptr(hl)+im)
         enddo
      enddo
!
! Write Hamiltonian

      do is=1,h%nspin
         do ih=1,no_u
            hl = ih
            do im=1,h%numh(hl)
               write(iu) h%Hamilt(h%listhptr(hl)+im,is)
            enddo
         enddo
      enddo
!
! Write Overlap matrix
      do ih = 1,no_u
         hl = ih
         do im = 1,h%numh(hl)
            write(iu) h%Sover(h%listhptr(hl)+im)
         enddo
      enddo
!
      write(iu) h%qtot, h%temp
!
      ! Write xij if not gamma
      ! (By convention, HS files only have X if not gamma...)

      write_xijk = .not. h%gamma
!
      if (write_xijk) then
         do ih = 1,no_u
            hl = ih
            do im = 1,h%numh(hl)
               write(iu) (h%xij(k,h%listhptr(hl)+im),k=1,3)
            enddo
         enddo
      else
         if (h%has_xij) then
            print *, "Your hsx record contains Xij info..."
            print *, "... even if it came from a Gamma-only calc..."
            print *, "... you will lose it by writing an HS file."
         endif
      endif   ! write_xijk
!

!       Write other useful info
!
      nspecies = h%nspecies
      write(iu) h%nspecies
      write(iu) (h%label(is), h%zval(is), h%no(is),is=1,nspecies)
      do is = 1, nspecies
         do io=1,h%no(is)
            write(iu) h%nquant(is,io), h%lquant(is,io), h%zeta(is,io)
         enddo
      enddo
      write(iu) h%na_u
      write(iu) (h%isa(ia),ia=1,h%na_u)
      write(iu) (h%iaorb(io), h%iphorb(io), io=1,no_u)
!
! Close file
      close(iu)

end subroutine write_hs_file

end module hsx_m
