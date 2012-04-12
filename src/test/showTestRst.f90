!----------------------------------------------------------------------------- best with 100 columns

module showTestRst
contains
  
  subroutine showPass(string)
    character(*),intent(in)::string
    write(*,'(/,a,a)'),'Pass: ',string
  end subroutine
  
  subroutine showFail(string)
    character(*),intent(in)::string
    write(*,'(/,a,a)'),'Fail: ',string
  end subroutine
  
end module
