! Module for mingling with the OS

! Fully implemented by Nick Papior Andersen
#include "mpi_macros.f"

module m_os
  implicit none
  
  private
  
  public :: file_exist
  public :: dir_exist
  public :: file_delete


contains

  function file_exist(file, Bcast, Comm, all) result(exist)
#ifdef MPI
    use mpi_siesta
#endif
    implicit none
    character(len=*), intent(in) :: file
    logical, intent(in), optional :: Bcast ! whether it should be b-casted
    MPI_COMM_TYPE, intent(in), optional  :: Comm ! the communicator (default World)
    logical, intent(in), optional :: all ! whether all Nodes sees the file
    logical :: exist
    
#ifdef MPI
    logical :: all_exist
    integer :: MPIerror
    MPI_COMM_TYPE :: Com
#endif

    ! Now we have all required information
    inquire(file=file, exist=exist)
    
#ifdef MPI
    Com = MPI_Comm_World
    if ( present(Comm) ) Com = Comm
    
    if ( present(Bcast) ) then

     if ( Bcast ) then
       
        call MPI_Bcast(exist,1,MPI_Logical, 0, Com, MPIerror)
       
     end if

    else if ( present(all) ) then

     if ( all ) then

        call MPI_AllReduce(exist, all_exist, 1, MPI_Logical, &
             MPI_LAnd, Com, MPIerror )
        
        exist = all_exist
        
     end if

    end if
#endif

  end function file_exist

  
  function dir_exist(dir, Bcast, Comm, all) result(exist)
#ifdef MPI
    use mpi_siesta
#endif
    implicit none
    character(len=*), intent(in) :: dir
    logical, intent(in), optional :: Bcast ! whether it should be b-casted
    MPI_COMM_TYPE, intent(in), optional  :: Comm ! the communicator (default World)
    logical, intent(in), optional :: all ! whether all Nodes sees the file
    logical :: exist
    
    integer :: ldir
#ifdef MPI
    logical :: all_exist
    integer :: MPIerror
    MPI_COMM_TYPE :: Com
#endif

    ! A directory of length 0 is the "top" directory,
    ! of course it exists
    ldir = len_trim(dir)
    if ( ldir == 0 ) then
       exist = .true.
       return
    else if ( ldir == 1 .and. dir(1:1) == '.' ) then
       exist = .true.
       return
    end if
    
    ! Now we have all required information
#ifdef __INTEL_COMPILER
    inquire(directory=dir, exist=exist)
#else
    if ( ldir == 1 ) then
       if ( dir(ldir:ldir) == '/' ) then
          inquire(file=trim(dir)//'.' , exist=exist)
       else
          inquire(file=trim(dir)//'/.', exist=exist)
       end if
    else
       if ( dir(ldir-1:ldir) == '/.' ) then
          inquire(file=trim(dir)      , exist=exist)
       else if ( dir(ldir:ldir) == '/' ) then
          inquire(file=trim(dir)//'.' , exist=exist)
       else
          inquire(file=trim(dir)//'/.', exist=exist)
       end if
    end if
#endif

#ifdef MPI
    Com = MPI_Comm_World
    if ( present(Comm) ) Com = Comm

    if ( present(Bcast) ) then

     if ( Bcast ) then
       
        call MPI_Bcast(exist,1,MPI_Logical, 0, Com, MPIerror)
       
     end if

    else if ( present(all) ) then

     if ( all ) then

        call MPI_AllReduce(exist, all_exist, 1, MPI_Logical, &
             MPI_LAnd, Com, MPIerror )
        
        exist = all_exist

     end if

    end if
#endif

  end function dir_exist


  function file_delete(file) result(existed)

    character(len=*), intent(in) :: file
    logical :: existed
    
    integer :: lfile, iu, iostat
    logical :: is_open

    ! A directory of length 0 is the "top" directory,
    ! of course it exists
    lfile = len_trim(file)
    if ( lfile == 0 ) then
      existed = .false.
      return
    else if ( lfile == 1 .and. file(1:1) == '.' ) then
      existed = .true.
      return
    end if

    existed = file_exist(file)
    if ( .not. existed ) return
    
    do iu = 1000, 10000
      inquire(unit=iu, opened=is_open)
      if ( is_open ) cycle
    end do

    ! Open file and close (delete it)
    open(unit=iu, file=file, iostat=iostat, status='old')
    if ( iostat == 0 ) then
      close(iu, status='delete')
    end if

  end function file_delete

end module m_os
