!----------------------------------------------------------------------------- best with 100 columns

program test
  integer,external::otGrid1
  integer,external::otGrid2
  
  call try('otGrid1',otGrid1)
  call try('otGrid2',otGrid2)
  
contains
  
  subroutine try(str,item)
    character(*),intent(in)::str
    integer,external::item
    
    i=item()
    if(i==0)then
      write(*,'(a,a,a)'),'[P]: "',trim(adjustl(str)),'"'
    else
      write(*,'(a,a,a,i6)'),'[F]: "',trim(adjustl(str)),'", ERR=',i
    end if
  end subroutine
  
end program
