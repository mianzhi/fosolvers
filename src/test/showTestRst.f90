!----------------------------------------------------------------------------- best with 100 columns

module utilityTest
  
  double precision,parameter::TOLERANCE=1d-6
  
contains
  
  subroutine showPass(string)
    character(*),intent(in)::string
    write(*,'(a,a)'),'Pass: ',string
  end subroutine
  
  subroutine showFail(string)
    character(*),intent(in)::string
    write(*,'(a,a)'),'Fail: ',string
  end subroutine
  
end module
