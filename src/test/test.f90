!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridInspection
  use moduleGridOperation
  use moduleGrid1D
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleCondition
  use moduleNonlinearSolve
  use moduleMPIComm
  use moduleSpatialHashing
  type(typeSHT)::sht
  integer,allocatable::rst(:)
  double precision pos(3,1000)
  
  call initMPI()
  if(pidMPI==0)then
    l=0
    do i=1,10
      do j=1,10
        do k=1,10
          l=l+1
          pos(:,l)=dble([i,j,k])
        end do
      end do
    end do
    call sht%fill(pos,3000)
    rst=sht%findNeib([4d0,7d0,8d0],10)
    do i=1,size(rst)
      write(*,*),pos(:,rst(i))
    end do
  else
  end if
  call finalMPI()
end program
