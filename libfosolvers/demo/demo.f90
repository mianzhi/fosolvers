!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  use moduleCond
  use moduleWrite
  use moduleMPIvar
  
  character(100) fnameGrid,fnameCond,fnameData,fnameRst
  double precision,target,allocatable::v(:),vv(:,:)
  double precision g(3)
  integer,parameter::FGRID_ID=10
  integer,parameter::FCOND_ID=11
  integer,parameter::FDATA_ID=12
  integer,parameter::FRST_ID=15
  
  call initMPI()
  call genPrtDataTab(0,0,0,1,1,0)
  
  if(pidMPI==ROOT_PID)then
    call getfNames(fnameGrid,fnameCond,fnameData,fnameRst)
    call readmsh(fnameGrid,FGRID_ID)
    call updateFacetPara()
    call updateElePara()
    
    ! prepare data space
    allocate(v(nEle))
    v(:)=0d0
    transEleScal(1)%ptr=>v
    allocate(vv(nEle,DIMS))
    vv(:,:)=0d0
    transEleVect(1)%ptr=>vv
    
    ! define conditions
    call readcod(fnameCond,FCOND_ID)
    
    ! read data tables
    call readdata(fnameData,FDATA_ID)
    
    ! broadcast static information
    call bcastStatic()
    
    ! assign task
    call distriPrt(2,1)
    
    ! gather data
    call gathData(i,j)
    
    ! write result
    call initWriteEnv(0,0,0,1,1,0)
    rstEleScal(1)%ptr=>v
    rstEleVect(1)%ptr=>vv
    call writerst(fnameRst,FRST_ID,.false.)
    
  else
  
    ! receive static informaion
    call recvStatic()
    
    ! receive task
    call recvPrt()
    allocate(v(nEle))
    v=transEleScal(1)%ptr
    deallocate(transEleScal(1)%ptr)
    allocate(vv(nEle,DIMS))
    vv=transEleVect(1)%ptr
    deallocate(transEleVect(1)%ptr)
    
    ! do work
    do i=1,nNode
      if(allocated(CondNode(i)%Cond))then
        write(*,*),'Node',CondNode(i)%Cond(1)%what,CondNode(i)%Cond(1)%Val
      end if
    end do
    do i=1,nFacet
      if(allocated(CondFacet(i)%Cond))then
        write(*,*),'Facet',CondFacet(i)%Cond(1)%what,CondFacet(i)%Cond(1)%Tab
      end if
    end do
    do i=1,nEle
      if(allocated(CondEle(i)%Cond))then
        write(*,*),'Ele',CondEle(i)%Cond(1)%what,CondEle(i)%Cond(1)%Val,CondEle(i)%Cond(1)%Tab
      end if
    end do
    do i=1,nEle
      v(i)=sin(t+10d0*Ele(i)%PC(1))+cos(5d0*Ele(i)%PC(2))-sin(5d0*Ele(i)%PC(3))
    end do
    do i=1,nEle
      call findEleGradScal(i,v,g)
      vv(i,:)=g(:)
    end do
    
    ! return data
    transEleScal(1)%ptr=>v
    transEleVect(1)%ptr=>vv
    call retnData()
  end if
  
  call finalMPI()
end program
