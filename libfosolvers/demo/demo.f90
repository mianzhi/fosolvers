program libtest
  use moduleGrid
  use moduleWrite
  
  double precision,allocatable::v(:),vv(:,:),vvv(:,:,:)
  double precision g(3),gg(3,3)
  
  call readmsh('./grid.msh',10)
  call sortEle()
  call updateFacetPara()
  call updateElePara()
  call initWriteEnv(0,0,0,1,1,1)
  
  allocate(v(nEle))
  allocate(vv(nEle,3))
  allocate(vvv(nEle,3,3))
  
  tFinal=5d0
  do l=0,50
    t=dble(l)/10d0
    do i=1,nEle
      v(i)=sin(t+10d0*Ele(i)%PC(1))+cos(5d0*Ele(i)%PC(2))-sin(5d0*Ele(i)%PC(3))
    end do
    do i=1,nEle
      call findEleGradScal(i,v,g)
      vv(i,:)=g(:)
    end do
    do i=1,nEle
      call findEleGradVect(i,vv,gg)
      vvv(i,:,:)=gg(:,:)
    end do
    rstEleScal(1,:)=v(:)
    rstEleVect(1,:,:)=vv(:,:)
    rstEleTens(1,:,1:3)=vvv(:,1,:)
    rstEleTens(1,:,4:6)=vvv(:,2,:)
    rstEleTens(1,:,7:9)=vvv(:,3,:)
    call writerst('rst.msh',11,.true.)
  end do
  call writerst('snap.msh',12,.false.)
end program
