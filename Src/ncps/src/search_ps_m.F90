! SPECIFICATION

! We are going to allow an optional "ps file spec" in the ChemicalSpecies Block, and
! at the same time provide a mechanism for searching for ps files in a "search path" as specified
! in PSEUDO_PATH or SIESTA_PS_PATH. 

! The behavior is as follows:

! The spec can be:

! - A single label, .e.g: Si

!   In this case, the routines search for Si.{ext} for ext in {"vps", "psf", "psml"}

! - A ps file name with an extension. e.g. Fe.psml

!   In this case, the search is only performed for "Fe.psml"

! - A full path, or a path ending in a label, e.g.   /home/user/pseudos/Si

!   In this case, the search is done only in the full path specified, completing the extension
!   if needed. 

! - A relative path, such as "sub/Si" or "../Mo.psml". This case is similar to the one above, except
!   that the search continues with the spec added to each section of PSEUDO_PATH. For example, if

!   PSEUDO_PATH=/home/pseudos:/data/pseudos

!   and spec=pbe/Si

!   the search might succeeded with /home/pseudos/pbe/Si.psf

!   Note that in this case the final piece of the path is NOT searched in the current directory.


! A further extension would be to allow environmental variable specifications. For example:

! %block Chemical-Species
! 1 14  Si_surf  $ENV{SURF_DATA_FILES}/pseudos/Si
! 2 201  H_star   $ENV{AUX_DATA_FILES}/terminators/H-0.9.psf
! %endblock Chemical-Species


! =====================================
! IMPLEMENTATION NOTES

! For convenience, we use "allocatable character variables".

! Note that the setting of the "automatic reallocation behavior" flag
! in the compiler CAN affect also allocatable strings.

! For example, if we compile with:

!   ifort -nostandard-realloc-lhs -check bounds -o tt env_utils_m.o search_path.F90

! standard numeric arrays cannot be dynamically allocated by assignment,
! but character variables can.

! For gfortran, with the option:

!   gfortran -fno-realloc-lhs -fcheck=all -o tt env_utils_m.o search_path.F90

! both numeric arrays AND allocatable strings fail to allocate dynamically.

! So, to be completely safe, either the option must be avoided (and
! F2003 and above standard-conformance required), at least for
! string-handling modules, or string-handling modules should be
! re-written in f95 style (i.e., with explicit allocations).
! We have taken the second route, allocating explicitly.

! See:
! https://stackoverflow.com/questions/42140832/automatic-array-allocation-upon-assignment-in-fortran


module search_ps_m

  public :: search_ps
  
CONTAINS
  
  subroutine search_ps(label, env_var, path, stat, extensions, debug_in)
    use env_utils, only: get_env_var
  
    implicit none
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: env_var
    character(len=:), allocatable, intent(out) :: path
    integer, intent(out) :: stat
    character(len=*), intent(in) :: extensions(:)
    logical, intent(in), optional :: debug_in

    character(len=:), allocatable :: pseudo_path
    character(len=:), allocatable :: try_path
    integer :: i, j, init, status_env, endpos, colon_idx, max_len
    integer :: max_len_ext
    logical :: exists, debug
    character(len=256) :: file

    stat = -1
    debug = .false.
    if (present(debug_in)) then
       debug = debug_in
    endif

    ! Extensions could be simply [""]
    max_len_ext = 0
    do j = 1, size(extensions)
       max_len_ext = max(max_len_ext,len_trim(extensions(j)))
    enddo
    
    ! Try raw entry, which could be a single file name (
    !   and hence will be tried in the current directory)
    !   or a more complex path (../File.psf, /Some/path/to/Si.psf)
    !
    do j = 1, size(extensions)
       file = trim(label)//trim(extensions(j))
       if (debug) print *, "Trying: ", trim(file)
       inquire(file=trim(file),exist=exists)
       if (exists) then
          allocate(character(len_trim(file)) :: path)
          path(:) = trim(file)
          stat = 0
          return
       endif
    enddo

    ! If we used an absolute path, do not try the environment search path
    if (label(1:1) == "/") RETURN
    
    ! If not found as is, try environment variable
    call get_env_var(env_var, pseudo_path, status_env)
    if (status_env /= 0) RETURN
    
    ! Unpack sections of PSEUDO_PATH
    init = 1
    ! Allocate work string
    max_len = len(pseudo_path) + 1 + len_trim(label) + max_len_ext
    allocate(character(len=max_len) :: try_path)
    
    do
      if (debug) print *, "Remaining: ", pseudo_path(init:)
      colon_idx = index(pseudo_path(init:), ':')
      
       if (colon_idx == 0) then
          endpos = len(pseudo_path)
       else
          if (colon_idx == 1) then
             init = init + 1
             cycle  !  Skip over initial ":" or "::" sections
          endif
          endpos = init + colon_idx - 2
       endif

       do j=1,size(extensions)
          ! note that using string section in lhs inhibits reallocation
          try_path(:) = pseudo_path(init:endpos) // &
                        "/" // trim(label) //       &
                        trim(extensions(j))
          if (debug) print *, "Trying: ", try_path
          inquire(file=try_path,exist=exists)
          if (exists) then
             allocate(character(len_trim(try_path)) :: path)
             path(:) = trim(try_path)
             stat = 0
             return
          endif
       enddo
       
       if (colon_idx == 0) exit
       init = init + colon_idx 
    enddo
    
  end subroutine search_ps
end module search_ps_m

