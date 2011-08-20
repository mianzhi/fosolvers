program libtest
  use moduleGrid
  use moduleWrite
  
  double precision,allocatable::v(:),vv(:,:),vvv(:,:,:)
  double precision g(3),gg(3,3)
  
  call readmsh('../heat3dFVM/grid.msh',10)
  call sortEle()
  call updateFacetPara()
  call updateElePara()
  call initWriteEnv(1,1,1,0,0,0)
  
  allocate(v(nNode))
  allocate(vv(nNode,3))
  allocate(vvv(nNode,3,3))
  
  tFinal=5d0
  do l=0,50
    t=dble(l)/10d0
    do i=1,nNode
      v(i)=(sin(t+10d0*Node(i)%Pos(1))+cos(t+5d0*Node(i)%Pos(2))-sin(t+5d0*Node(i)%Pos(3)))
    end do
    do i=1,nNode
      call findNodeGradScal(i,v,g)
      vv(i,:)=g(:)
    end do
    do i=1,nNode
      call findNodeGradScal(i,vv,gg)
      vvv(i,:,:)=gg(:,:)
    end do
    rstNodeScal(1,:)=v(:)
    rstNodeVect(1,:,:)=vv(:,:)
    rstNodeTens(1,:,1:3)=vvv(:,1,:)
    rstNodeTens(1,:,4:6)=vvv(:,2,:)
    rstNodeTens(1,:,7:9)=vvv(:,3,:)
    call writerstSpan('rst.msh',11)
  end do
end program
