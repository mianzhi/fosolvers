!----------------------------------------------------------------------------- best with 100 columns

program test
  integer,external::otGrid1,otGrid2,&
  &                 polyX1,polyMesh1,polyGrid1,polyFvGrid1,&
  &                 readGTS1,readVTK1,writeVTK1,&
  &                 readCondTab1,&
  &                 gradient1,advection1,advection2,euler,eulerJac,&
  &                 diffusion1,&
  &                 UDF1,&
  &                 fixPt1
  
  call try('otGrid1',otGrid1)
  call try('otGrid2',otGrid2)
  call try('polyX1',polyX1)
  call try('polyMesh1',polyMesh1)
  call try('polyGrid1',polyGrid1)
  call try('polyFvGrid1',polyFvGrid1)
  call try('readGTS1',readGTS1)
  call try('readVTK1',readVTK1)
  call try('writeVTK1',writeVTK1)
  call try('readCondTab',readCondTab1)
  call try('gradient1',gradient1)
  call try('advection1',advection1)
  call try('advection2',advection2)
  call try('diffusion1',diffusion1)
  call try('euler',euler)
  call try('eulerJac',eulerJac)
  call try('UDF1',UDF1)
  call try('fixPt1',fixPt1)
  
  write(*,'(a)')'check md5 of output files:'
  call system('md5sum -c sums.md5')
  
contains
  
  subroutine try(str,item)
    character(*),intent(in)::str
    integer,external::item
    
    i=item()
    if(i<0)then
      write(*,'(a,a,a,i3)')'[P]: "',trim(adjustl(str)),'", epsilon < 1E',i
    else if(i==0)then
      write(*,'(a,a,a)')'[P]: "',trim(adjustl(str)),'"'
    else
      write(*,'(a,a,a,i6)')'[FAIL]: "',trim(adjustl(str)),'", ERR=',i
    end if
  end subroutine
  
end program
