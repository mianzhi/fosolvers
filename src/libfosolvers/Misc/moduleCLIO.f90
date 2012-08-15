!----------------------------------------------------------------------------- best with 100 columns

!> command line input/output
module moduleCLIO
  private
  
  ! public procedures
  public::showNoAdv
  public::showProg
  
contains
  
  !> show string without advancing the line
  subroutine showNoAdv(str)
    character(*),intent(in)::str !< the string to be written
    
    write(*,'(80a,80a,$)'),(char(8),i=1,80),(char(32),i=1,80)
    write(*,'(80a,a,$)'),(char(8),i=1,80),str
  end subroutine
  
  !> show progress
  subroutine showProg(prog)
    double precision,intent(in)::prog !< the progress
    character(80) str
    
    write(str,'(a,f5.1,a)'),'progress:',prog*100d0,'%'
    call showNoAdv(trim(str))
  end subroutine
  
end module
