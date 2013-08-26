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
  double precision src(3,1000),tar(3,27),dats(1000)
  double precision,allocatable::datt(:)
  
  call initMPI()
  if(pidMPI==0)then
    l=0
    do i=1,10
      do j=1,10
        do k=1,10
          l=l+1
          src(:,l)=dble([i,j,k])
          dats(l)=norm2(src(:,l))
        end do
      end do
    end do
    call sht%fill(src,3000)
    l=0
    do i=1,3
      do j=1,3
        do k=1,3
          l=l+1
          tar(:,l)=dble([i,j,k])*1.5d0+0.3d0
        end do
      end do
    end do
    datt=itplScat2Scat(dats,src,sht,tar)
    do i=1,size(datt)
      write(*,*),i,datt(i)-norm2(tar(:,i))
    end do
  else
  end if
  call finalMPI()
end program
