module units_common_m

   implicit none
   public

   integer, parameter, private :: dp = selected_real_kind(14,100)

contains

  ! Returns information about a unit in the units table
  ! 
  ! Unit specifications might include an optional 'physical dimension'
  ! qualifier (e.g. 'bfield:g')
  ! In this case, 'phys_dim' returns the physical dimension, and the
  ! qualifier is used to match the unit.
  ! This version is case-insensitive (e.g. 'g' and 'G' could stand for 'Gauss').
  ! As the above example indicates, in the absence of a physical dimension qualifier,
  ! 'g' might be ambiguous ('bfield' or 'mass'?). The routine will return 'stat=-1'
  ! in this case.
  ! Units might be ambiguous in a more serious way: 'meV' and 'MeV' could both be
  ! present in the table. In this case, it might be advisable to use a case-sensitive
  ! version of this routine (replacing 'leqi' by 'leqi_strict' below).
  ! If the unit is not found in the table, the routine returns 'stat=-2'.

  subroutine inquire_unit_table(unit_str, stat, phys_dim, unit_name, unit_value, &
     nu, dimm, name, unit)
    
    character(len=*), intent(in)   :: unit_str   ! unit specification
    character(len=*), intent(out)  :: phys_dim   ! physical dimension (e.g. 'mass')
    character(len=*), intent(out)  :: unit_name  ! unit name (e.g. 'g')
    real(dp), intent(out)          :: unit_value ! actual value (e.g. 1.e-3)
    integer, intent(out)           :: stat       ! status code
    integer, intent(in) :: nu
    character(*), intent(in) :: dimm(nu)
    character(*), intent(in) :: name(nu)
    real(dp), intent(in) :: unit(nu)

    integer           :: idx_colon, iu, idx
    logical           :: phys_dim_specified, match
    
    idx_colon = index(unit_str,":")
    if (idx_colon /= 0) then
       ! spec includes dimension prefix
       phys_dim = unit_str(1:idx_colon-1)
       unit_name = unit_str(idx_colon+1:)
       phys_dim_specified = .true.
    else
       phys_dim = ""
       unit_name = unit_str
       phys_dim_specified = .false.
    endif

    stat = 0
    idx = 0

    do iu= 1, nu
         match = .false.
         if (leqi(name(iu), unit_name)) then
            if (phys_dim_specified) then
               if (leqi(dimm(iu), phys_dim)) then
                  match = .true.
               endif
            else
               match = .true.
            endif
         endif
         if (match) then
            if (idx /= 0) then  ! ambiguous
               stat = 1
               RETURN
            endif
            idx = iu
         endif
      enddo
      
      if (idx == 0) then
         stat = -2    ! not found
      else
         phys_dim = trim(dimm(idx))
         unit_value = unit(idx)
      endif
      



  end subroutine

!
!   Case-insensitive lexical equal-to comparison
!
    FUNCTION leqi(string1, string2)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: string1, string2

!-------------------------------------------------------------- Output Variables
      logical          :: leqi

!--------------------------------------------------------------- Local Variables
      logical          :: completed
      character        :: char1, char2
      integer          :: i, len1, len2, lenc

!------------------------------------------------------------------------- BEGIN
      len1 = len(string1)
      len2 = len(string2)
      lenc = min(len1, len2)

      i = 1
      leqi      = .TRUE.
      completed = .FALSE.
      do while((.not. completed) .and. (i .le. lenc))
        char1 = string1(i:i)
        char2 = string2(i:i)
        call chrcap(char1, 1)
        call chrcap(char2, 1)
        if (char1 .ne. char2) then
          leqi      = .FALSE.
          completed = .TRUE.
        endif

        i = i + 1
      enddo

      if (leqi) then
        if ((len1 .gt. lenc) .and. (string1(lenc+1:len1) .ne. ' '))     &
          leqi = .FALSE.
        if ((len2 .gt. lenc) .and. (string2(lenc+1:len2) .ne. ' '))     &
          leqi = .FALSE.
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION leqi

!
!   Examples of eq_func's for search function (Case sensitive)
!
    FUNCTION leqi_strict(str1, str2)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: str1, str2

!-------------------------------------------------------------- Output Variables
      logical          :: leqi_strict

!------------------------------------------------------------------------- BEGIN
      leqi_strict = (str1 .eq. str2)
      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION leqi_strict

!
!   CHRCAP accepts a STRING of NCHAR characters and replaces
!   any lowercase letters by uppercase ones.
!
    SUBROUTINE chrcap(string, nchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer  :: nchar

!-------------------------------------------------------------- Output Variables
      character(*) :: string

!--------------------------------------------------------------- Local Variables
      integer  :: i, itemp, ncopy

!------------------------------------------------------------------------- BEGIN
      if (nchar .le. 0) then
        ncopy = LEN(string)
      else
        ncopy = nchar
      endif

      do i= 1, ncopy
        if (LGE(string(i:i),'a') .and. LLE(string(i:i),'z')) then
          itemp = ICHAR(string(i:i)) + ICHAR('A') - ICHAR('a')
          string(i:i) = CHAR(itemp)
        endif
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE chrcap

!
!   CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
!   the length of the string up to the last NONBLANK, NONNULL.
!     
    SUBROUTINE chrlen(string, nchar, lchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*) :: string
      integer  :: nchar

!-------------------------------------------------------------- Output Variables
      integer  :: lchar

!------------------------------------------------------------------------- BEGIN
      lchar = nchar
      if (lchar .le. 0) lchar = LEN(string)

      do while(((string(lchar:lchar) .eq. ' ') .or. (string(lchar:lchar) &
               .eq. CHAR(0))) .and. (lchar .gt. 0))
        lchar = lchar - 1
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE chrlen

end module
