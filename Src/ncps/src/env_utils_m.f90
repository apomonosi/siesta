module env_utils

  public :: get_env_var
  
CONTAINS
  
  subroutine get_env_var(name,value,status)
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: value
    integer, intent(out)         :: status
    
    integer len
    
    call get_environment_variable (name, length=len, status=status, trim_name=.true.)

    if (status .eq. 1) then
       !write (*,*) 'env var does not exist'
       return
    end if
    if (len .eq. 0) then
       !write (*,*) 'env var exists  but has no value'
       status = -2
       return
    end if

    allocate(character(len) :: value)
    call get_environment_variable (name=name, value=value, trim_name=.true.)
    
  end subroutine get_env_var
  
end module env_utils
