!----------------------------------------------------------------------------- best with 100 columns

program test
  integer,external::otGrid1,otGrid2,&
  &                 polyX1,polyMesh1,polyGrid1,&
  &                 readGTS1,readVTK1
  
  call try('otGrid1',otGrid1)
  call try('otGrid2',otGrid2)
  call try('polyX1',polyX1)
  call try('polyMesh1',polyMesh1)
  call try('polyGrid1',polyGrid1)
  call try('readGTS1',readGTS1)
  call try('readVTK1',readVTK1)
  
contains
  
  subroutine try(str,item)
    character(*),intent(in)::str
    integer,external::item
    
    i=item()
    if(i==0)then
      write(*,'(a,a,a)'),'[P]: "',trim(adjustl(str)),'"'
    else
      write(*,'(a,a,a,i6)'),'[FAIL]: "',trim(adjustl(str)),'", ERR=',i
    end if
  end subroutine
  
end program
