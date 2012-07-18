!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleMPIComm
  type(typeGrid)::grid
  double precision,allocatable::v(:),vi(:)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH4.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    allocate(v(grid%nNode))
    forall(i=1:grid%nNode)
      v(i)=grid%NodePos(1,i)
    end forall
    vi=itplNode2Intf(v,grid)
    do i=1,grid%nIntf
      write(*,*),i,':',vi(i)
    end do
  else
  end if
  call finalMPI()
end program
