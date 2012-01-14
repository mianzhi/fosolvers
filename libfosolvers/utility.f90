!----------------------------------------------------------------------------- best with 100 columns

!***********
! utilities
!***********
module moduleUtility
  private
  
  ! procedures
  public showError
  public showWarning
  public showNoAdv
  
contains
  
  !---------------------
  ! show error and stop
  !---------------------
  subroutine showError(string)
    character(*),intent(in)::string
    
    write(*,'(/,a,a)'),'ERROR: ',string
    stop
  end subroutine
  
  !--------------
  ! show warning
  !--------------
  subroutine showWarning(string)
    character(*),intent(in)::string
    
    write(*,'(/,a,a)'),'WARNING: ',string
  end subroutine
  
  !----------------------------------------
  ! show string without advancing the line
  !----------------------------------------
  subroutine showNoAdv(string)
    character(*),intent(in)::string
    
    ! clean line
    write(*,'(80a,80a,$)'),(char(8),i=1,80),(char(32),i=1,80)
    ! write line
    write(*,'(80a,a,$)'),(char(8),i=1,80),string
  end subroutine
end module
