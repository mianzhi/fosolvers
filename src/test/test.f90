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
  integer rst(10)
  double precision pos(3,1000),dist(10)
  
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
    call findNearestSHT(pos,sht,[3d0,6d0,7d0],10,rst,dist)
    do i=1,size(rst)
      write(*,*),rst(i),dist(i)
    end do
  else
  end if
  call finalMPI()
end program
