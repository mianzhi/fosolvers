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
  
  call initMPI()
  if(pidMPI==0)then
    allocate(sht%list(3,3,3))
    sht%SideLength=1d0
    sht%BoundBox(:,1)=1.3d0
    sht%BoundBox(:,2)=3.8d0
    do i=1,3
      do j=1,3
        do k=1,3
          call pushArr(sht%list(i,j,k)%dat,i*100+j*10+k)
        end do
      end do
    end do
    rst=sht%lookup(sht%hash([1.5d0,2.5d0,3.5d0]))
    write(*,*),rst
  else
  end if
  call finalMPI()
end program
