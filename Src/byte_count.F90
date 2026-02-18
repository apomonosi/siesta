!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior, 2020, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for accummulating byte counts from array dimensions.

! This module is primarily intended for determining file
! sizes by passing sizes of arrays stored.
!
! Its basic usage may be illustrated by this:
! \code
! character(len=32) :: mem_unit
! type(byte_count_t) :: byte
!
! call byte%add(4, 10, 20, 30) ! add array(10, 20, 30) with each element 4 bytes
!
! call byte%get_string(mem_unit) ! retrieve a string representing the total size (with unit)
!
! \endcode

module byte_count_m

  use precision, only: sp, dp, i4b, i8b

  implicit none

  private

  type :: byte_count_t

    !< Mega bytes accummulated
    real(dp) :: MB = 0._dp

  contains

    procedure, pass :: reset => reset_
    procedure, pass :: get_string => get_string_
#ifdef NCDF_4
    procedure, pass, private :: add_cdf_basic_, add_cdf_complex_
    generic :: add_cdf => add_cdf_basic_, add_cdf_complex_
#endif
    procedure, pass, private :: add_i4_, add_i8_, add_i8i8_
    procedure, pass, private :: add_ai4_, add_ai8_, add_ai8i8_
    procedure, pass, private :: add_bc_, add_r8i4_, add_r8i8_
    generic :: add => add_bc_, add_i4_, add_i8_, add_i8i8_, &
        add_ai4_, add_ai8_, add_ai8i8_, add_r8i4_, add_r8i8_

    ! Type based (using storage_size)
    procedure, pass, private :: add_type_l_i4_, add_type_l_i8_
    procedure, pass, private :: add_type_i4_i4_, add_type_i4_i8_
    procedure, pass, private :: add_type_i8_i4_, add_type_i8_i8_
    procedure, pass, private :: add_type_r4_i4_, add_type_r4_i8_
    procedure, pass, private :: add_type_r8_i4_, add_type_r8_i8_
    procedure, pass, private :: add_type_c4_i4_, add_type_c4_i8_
    procedure, pass, private :: add_type_c8_i4_, add_type_c8_i8_

    generic :: add_type => &
        add_type_l_i4_, add_type_l_i8_, &
        add_type_i4_i4_, add_type_i4_i8_, &
        add_type_i8_i4_, add_type_i8_i8_, &
        add_type_r4_i4_, add_type_r4_i8_, &
        add_type_r8_i4_, add_type_r8_i8_, &
        add_type_c4_i4_, add_type_c4_i8_, &
        add_type_c8_i4_, add_type_c8_i8_

    ! Mathematical stuff
    procedure, pass, private :: assign_, add_

    generic :: assignment(=) => assign_
    generic :: operator(+) => add_

  end type byte_count_t

  public :: byte_count_t

  !< Conversion from bytes to mega-bytes
  real(dp), private, parameter :: B2MB = 1._dp / 1024._dp ** 2
  real(dp), private, parameter :: Bi2B = 1._dp / 8._dp
  real(dp), private, parameter :: Bi2MB = Bi2B * B2MB

contains

  subroutine assign_(this, other)
    class(byte_count_t), intent(inout) :: this
    class(byte_count_t), intent(in) :: other
    this%MB = other%MB
  end subroutine assign_

  function add_(this, other) result(res)
    class(byte_count_t), intent(in) :: this
    class(byte_count_t), intent(in) :: other
    type(byte_count_t) :: res
    res%MB = this%MB + other%MB
  end function add_

  subroutine reset_(this)
    class(byte_count_t), intent(inout) :: this
    this%MB = 0._dp
  end subroutine reset_

  subroutine add_bc_(this, o1, o2, o3, o4, o5, o6)
    class(byte_count_t), intent(inout) :: this
    class(byte_count_t), intent(in) :: o1
    class(byte_count_t), intent(in), optional :: o2, o3, o4, o5, o6

    ! Number of elements
    this%MB = this%MB + o1%MB
    if ( present(o2) ) this%MB = this%MB + o2%MB
    if ( present(o3) ) this%MB = this%MB + o3%MB
    if ( present(o4) ) this%MB = this%MB + o4%MB
    if ( present(o5) ) this%MB = this%MB + o5%MB
    if ( present(o6) ) this%MB = this%MB + o6%MB

  end subroutine add_bc_

  subroutine add_i4_(this, bytes, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    this%MB = this%MB + MB

  end subroutine add_i4_

  subroutine add_ai4_(this, bytes, n)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i4b), intent(in) :: n(:)

    integer :: i
    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n(1), dp) * B2MB
    do i = 2, size(n)
      MB = MB * real(n(i), dp)
    end do
    this%MB = this%MB + MB

  end subroutine add_ai4_

  subroutine add_ai8_(this, bytes, n)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i8b), intent(in) :: n(:)

    integer :: i
    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n(1), dp) * B2MB
    do i = 2, size(n)
      MB = MB * real(n(i), dp)
    end do
    this%MB = this%MB + MB

  end subroutine add_ai8_

  subroutine add_ai8i8_(this, bytes, n)
    class(byte_count_t), intent(inout) :: this
    integer(i8b), intent(in) :: bytes
    integer(i8b), intent(in) :: n(:)

    integer :: i
    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n(1), dp) * B2MB
    do i = 2, size(n)
      MB = MB * real(n(i), dp)
    end do
    this%MB = this%MB + MB

  end subroutine add_ai8i8_

  subroutine add_i8_(this, bytes, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    this%MB = this%MB + MB

  end subroutine add_i8_

  subroutine add_i8i8_(this, bytes, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i8b), intent(in) :: bytes
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    this%MB = this%MB + MB

  end subroutine add_i8i8_

  subroutine add_r8i4_(this, bytes, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    real(dp), intent(in) :: bytes
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    this%MB = this%MB + MB

  end subroutine add_r8i4_

  subroutine add_r8i8_(this, bytes, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    real(dp), intent(in) :: bytes
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = bytes * real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    this%MB = this%MB + MB

  end subroutine add_r8i8_

#ifdef NCDF_4
  subroutine add_cdf_basic_(this, nf_var, n1, n2, n3, n4, n5, n6)
    use netcdf, only: NF90_BYTE, NF90_UBYTE
    use netcdf, only: NF90_SHORT, NF90_USHORT
    use netcdf, only: NF90_INT, NF90_UINT
    use netcdf, only: NF90_INT64, NF90_UINT64
    use netcdf, only: NF90_FLOAT, NF90_DOUBLE

    class(byte_count_t), intent(inout) :: this
    integer, intent(in) :: nf_var
    integer, intent(in) :: n1
    integer, intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    select case ( nf_var )
    case ( NF90_BYTE, NF90_UBYTE )
      ! do nothing
    case ( NF90_SHORT, NF90_USHORT )
      MB = MB * 2
    case ( NF90_INT, NF90_UINT, NF90_FLOAT )
      MB = MB * 4
    case ( NF90_INT64, NF90_UINT64, NF90_DOUBLE )
      MB = MB * 8
    case default
      !error
    end select

    this%MB = this%MB + MB

  end subroutine add_cdf_basic_

  subroutine add_cdf_complex_(this, nf_var, n1, n2, n3, n4, n5, n6)
    use netcdf_ncdf, only: NF90_FLOAT_COMPLEX, NF90_DOUBLE_COMPLEX

    class(byte_count_t), intent(inout) :: this
    logical, intent(in) :: nf_var
    integer, intent(in) :: n1
    integer, intent(in), optional :: n2, n3, n4, n5, n6

    real(dp) :: MB

    ! Number of elements
    MB = real(n1, dp) * B2MB
    if ( present(n2) ) MB = MB * real(n2, dp)
    if ( present(n3) ) MB = MB * real(n3, dp)
    if ( present(n4) ) MB = MB * real(n4, dp)
    if ( present(n5) ) MB = MB * real(n5, dp)
    if ( present(n6) ) MB = MB * real(n6, dp)

    if ( nf_var .eqv. NF90_FLOAT_COMPLEX ) then
      MB = MB * 8
    else if ( nf_var .eqv. NF90_DOUBLE_COMPLEX ) then
      MB = MB * 16
    end if

    this%MB = this%MB + MB

  end subroutine add_cdf_complex_
#endif

  ! Type additions

  subroutine add_type_l_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    logical, intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_l_i4_

  subroutine add_type_l_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    logical, intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_l_i8_


  subroutine add_type_i4_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_i4_i4_

  subroutine add_type_i4_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_i4_i8_

  subroutine add_type_i8_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i8b), intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_i8_i4_

  subroutine add_type_i8_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    integer(i8b), intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_i8_i8_

  subroutine add_type_r4_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    real(sp), intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_r4_i4_

  subroutine add_type_r4_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    real(sp), intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_r4_i8_

  subroutine add_type_r8_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    real(dp), intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_r8_i4_

  subroutine add_type_r8_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    real(dp), intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_r8_i8_

  subroutine add_type_c4_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    complex(sp), intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_c4_i4_

  subroutine add_type_c4_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    complex(sp), intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_c4_i8_

  subroutine add_type_c8_i4_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    complex(dp), intent(in) :: typ
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_c8_i4_

  subroutine add_type_c8_i8_(this, typ, n1, n2, n3, n4, n5, n6)
    class(byte_count_t), intent(inout) :: this
    complex(dp), intent(in) :: typ
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5, n6

    call this%add(storage_size(typ) * Bi2B, n1, n2, n3, n4, n5, n6)

  end subroutine add_type_c8_i8_


  subroutine get_string_(this, mem_str, unit)
    use m_char, only: lcase 
    !< Memory class containing the amount of memory
    class(byte_count_t), intent(in) :: this
    !< Output string for the memory with an appropriate unit.
    !!
    !! If `unit` is passed the unit will be in the requested unit.
    !! Otherwise the unit will be automatically determined [kB, MB, GB, TB, PB, EB]
    !! The format will be [es*.2] depending on the length of `mem_str`
    character(len=*), intent(inout) :: mem_str
    character(len=2), intent(in), optional :: unit

    character(len=2) :: lunit
    character(len=24) :: fmt

    real(dp) :: mem
    integer :: fw, fd, fe

    ! Determine length available
    fw = len(mem_str)

    ! Determine default unit
    lunit = 'kB'
    mem = this%MB * 1024._dp
    if ( mem >= 1000._dp ) then
      lunit = 'MB'
      mem = mem / 1024._dp
      if ( mem >= 1000._dp ) then
        lunit = 'GB'
        mem = mem / 1024._dp
        if ( mem >= 1000._dp ) then
          lunit = 'TB'
          mem = mem / 1024._dp
          if ( mem >= 1000._dp ) then
            lunit = 'PB'
            mem = mem / 1024._dp
            if ( mem >= 1000._dp ) then
              lunit = 'EB'
              mem = mem / 1024._dp
            end if
          end if
        end if
      end if
    end if

    if ( present(unit) ) then
      lunit = lcase(unit)
      select case ( lunit )
      case ( 'b' )
        lunit = 'B'
        mem = this%MB * 1024._dp ** 2
      case ( 'kb' )
        lunit = 'kB'
        mem = this%MB * 1024._dp
      case ( 'mb' )
        lunit = 'MB'
        mem = this%MB
      case ( 'gb' )
        lunit = 'GB'
        mem = this%MB / 1024._dp
      case ( 'tb' )
        lunit = 'TB'
        mem = this%MB / 1024._dp ** 2
      case ( 'pb' )
        lunit = 'PB'
        mem = this%MB / 1024._dp ** 3
      case ( 'eb' )
        lunit = 'EB'
        mem = this%MB / 1024._dp ** 4
      case default
        lunit = 'GB'
        mem = this%MB / 1024._dp
      end select
    end if

    ! Now determine format
    !  ENw.dEe
    ! w = 3 + [.] + d + [E] + [+-] + e + [ kB]

    if ( mem < 1000._dp ) then

      write(mem_str, '(f9.3,tr1,a)') mem, lunit
      ! Move to the left
      mem_str = adjustl(mem_str)

      return

    end if

    ! Subtract ' *B'
    fw = fw - 3

    fe = 1
    if ( abs(exponent(mem)) > 12 ) fe = 2

    ! Maximally 3 digits after '.'
    fd = min(fw - 3 - 1 - 1 - 1 - fe, 3)
    write(fmt,'(''('',a,i0,a,i0,a,i0,'',tr1,a)'')') 'EN',fw,'.',fd,'E',fe
    write(mem_str, trim(fmt)) mem, lunit
    ! Move to the left
    mem_str = adjustl(mem_str)

  end subroutine get_string_

end module byte_count_m
