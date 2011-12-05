!----------------------------------------------------------------------------- best with 100 columns

!***********
! utilities
!***********
module moduleUtility
contains
  
  !---------------------
  ! show error and stop
  !---------------------
  subroutine showError(string)
    character(*),intent(in)::string
    
    write(*,'(/,a,a)'),'ERROR: ',string
    stop
  end subroutine
end module
