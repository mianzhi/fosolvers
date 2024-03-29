!----------------------------------------------------------------------------- best with 100 columns

program test
  integer,external::otGrid1,otGrid2,&
  &                 polyX1,polyMesh1,polyGrid1,polyFeGrid1,&
  &                 readGTS1,readVTK1,writeVTK1,&
  &                 readCondTab1,UDF1,&
  &                 fixPt1,NewtonKrylov1,BDFNewtonKrylov1,ILU1,multiFront1,matVecMul,&
  &                 basicFEM1,edgeFEM1
  
  call try('otGrid1',otGrid1)
  call try('otGrid2',otGrid2)
  call try('polyX1',polyX1)
  call try('polyMesh1',polyMesh1)
  call try('polyGrid1',polyGrid1)
  call try('polyFeGrid1',polyFeGrid1)
  call try('readGTS1',readGTS1)
  call try('readVTK1',readVTK1)
  call try('writeVTK1',writeVTK1)
  call try('readCondTab',readCondTab1)
  call try('UDF1',UDF1)
  call try('fixPt1',fixPt1)
  call try('NewtonKrylov1',NewtonKrylov1)
  call try('BDFNewtonKrylov1',BDFNewtonKrylov1)
  call try('ILU1',ILU1)
  call try('multiFront1',multiFront1)
  call try('matVecMul',matVecMul)
  call try('basicFEM1',basicFEM1)
  call try('edgeFEM1',edgeFEM1)
  
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
