!----------------------------------------------------------------------------- best with 100 columns

!> fons main program
program fons
  use moduleBasicDataStruct
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  use moduleInterpolation
  use moduleFVMGrad
  use moduleFVMConvect
  use moduleNonlinearSolve
  use moduleCLIO
  use miscNS
  double precision,allocatable::uBlock(:,:) !< velocity at block
  double precision,allocatable::uIntf(:,:) !< velocity at interface
  double precision,allocatable::pIntf(:) !< pressure at interface
  double precision,allocatable::gradRho(:,:) !< gradient of rho
  double precision,allocatable::gradRhou(:,:,:) !< gradient of rhou
  double precision,allocatable::gradRhoE(:,:) !< gradient of rhoE
  double precision pWork !< pressure work done on block surface
  external::resMom
  double precision,allocatable::u1d(:) !< unwrapped velocity
  
  ! read simulation control file
  open(11,file='bin/sim',status='old')
  call readSim(11)
  close(11)
  ! read grid
  open(11,file='bin/gridGMSH5.msh',status='old')
  call readGMSH(11,grid)
  call grid%updateIntf()
  close(11)
  ! read conditions
  open(11,file='bin/cond',status='old')
  call readCondition(11,condition)
  close(11)
  ! allocate storage
  allocate(rho(grid%nBlock))
  allocate(u(DIMS,grid%nNode))
  allocate(p(grid%nBlock))
  allocate(E(grid%nBlock))
  allocate(rhou(DIMS,grid%nNode))
  allocate(rhoE(grid%nBlock))
  allocate(rhoNode(grid%nNode))
  allocate(uBlock(DIMS,grid%nBlock))
  allocate(uIntf(DIMS,grid%nIntf))
  allocate(pIntf(grid%nIntf))
  allocate(Mass(grid%nBlock))
  allocate(Mom(DIMS,grid%nNode))
  allocate(Energy(grid%nBlock))
  allocate(gradRho(DIMS,grid%nBlock))
  allocate(gradRhou(DIMS,DIMS,grid%nNode))
  allocate(gradRhoE(DIMS,grid%nBlock))
  allocate(u1d(DIMS*grid%nNode))
  ! simulation control
  dt=1d-5
  ! initial value of variables
  call grid%updateBlockPos()
  gamm=1.4d0
  u(:,:)=0d0
  rhou(:,:)=0d0
  forall(i=1:grid%nBlock)
    rho(i)=merge(1d0,0.125d0,grid%BlockPos(1,i)<0.5d0)
    p(i)=merge(1d5,1d4,grid%BlockPos(1,i)<0.5d0)
    E(i)=p(i)/rho(i)/(gamm-1d0)
    rhoE(i)=rho(i)*E(i)
  end forall
  t=0d0
  iWrite=0
  ! write initial states
  open(13,file='rstNS.msh',status='replace')
  call writeGMSH(13,grid)
  call writeGMSH(13,rho,grid,BIND_BLOCK,'rho',iWrite,t)
  call writeGMSH(13,u,grid,BIND_NODE,'u',iWrite,t)
  call writeGMSH(13,p,grid,BIND_BLOCK,'p',iWrite,t)
  ! advance in time
  do while(t<tFinal)
    call grid%updateDualBlock()
    call grid%updateBlockVol()
    call grid%updateFacetNorm()
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    ! Lagrangian step
    Mass(:)=rho(:)*grid%BlockVol(:)
    forall(i=1:grid%nNode)
      Mom(:,i)=rhou(:,i)*grid%NodeVol(i)
    end forall
    Energy(:)=rhoE(:)*grid%BlockVol(:)
    rhoNode=itplBlock2Node(rho,grid)
    u1d(:)=reshape(u,[DIMS*grid%nNode])
    if(maxval(abs(u1d))<=tiny(1d0))then
      u1d(:)=1d0
    end if
    ProblemFunc=>resMom
    call solveNonlinear(u1d)
    u=reshape(u1d,[DIMS,grid%nNode])
    forall(i=1:grid%nNode)
      rhou(:,i)=u(:,i)*rhoNode(i)
      Mom(:,i)=rhou(:,i)*grid%NodeVol(i)
    end forall
    uIntf=itplNode2Intf(u,grid)
    pIntf=itplBlock2Intf(p,grid)
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      pWork=dt*pIntf(i)*grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),uIntf(:,i))
      Energy(m)=Energy(m)-pWork
      Energy(n)=Energy(n)+pWork
    end do
    ! apply BC
    do i=1,grid%nFacet
      do j=1,grid%Facet(i)%nNode
        k=grid%Facet(i)%iNode(j)
        Mom(:,k)=Mom(:,k)-grid%FacetNorm(:,i)*dot_product(grid%FacetNorm(:,i),Mom(:,k))
      end do
    end do
    ! recover intensive state
    rhoNode=itplBlock2Node(rho,grid)
    forall(i=1:grid%nNode)
      rhou(:,i)=Mom(:,i)/grid%NodeVol(i)
      u(:,i)=rhou(:,i)/rhoNode(i)
    end forall
    uBlock=itplNode2Block(u,grid) !FIXME: seems can be removed
    forall(i=1:grid%nBlock)
      rhoE(i)=Energy(i)/grid%BlockVol(i)
    end forall
    rho(:)=Mass(:)/grid%BlockVol(:)
    call mvGrid(grid,dt*u)
    ! Euler rezoning
    gradRho=findGrad(rho,grid,BIND_BLOCK)
    gradRhou=findGrad(rhou,grid,BIND_NODE)
    gradRhoE=findGrad(rhoE,grid,BIND_BLOCK)
    Mass=Mass+findDispConvect(rho,BIND_BLOCK,-dt*u,grid,gradRho,limiter=vanLeer)
    Mom=Mom+findDispConvect(rhou,BIND_NODE,-dt*u,grid,gradRhou,limiter=vanLeer)
    Energy=Energy+findDispConvect(rhoE,BIND_BLOCK,-dt*u,grid,gradRhoE,limiter=vanLeer)
    call mvGrid(grid,-dt*u)
    ! recover state
    call grid%updateDualBlock()
    call grid%updateBlockVol()
    rhoNode=itplBlock2Node(rho,grid)
    forall(i=1:grid%nNode)
      rhou(:,i)=Mom(:,i)/grid%NodeVol(i)
      u(:,i)=rhou(:,i)/rhoNode(i)
    end forall
    uBlock=itplNode2Block(u,grid)
    forall(i=1:grid%nBlock)
      rhoE(i)=Energy(i)/grid%BlockVol(i)
      E(i)=rhoE(i)/rho(i)
      p(i)=(E(i)-dot_product(uBlock(:,i),uBlock(:,i))/2d0)*rho(i)*(gamm-1d0)
    end forall
    rho(:)=Mass(:)/grid%BlockVol(:)
    t=t+dt
    ! write results
    if(t/tWrite>=iWrite)then
      iWrite=iWrite+1
      call writeGMSH(13,rho,grid,BIND_BLOCK,'rho',iWrite,t)
      call writeGMSH(13,u,grid,BIND_NODE,'u',iWrite,t)
      call writeGMSH(13,p,grid,BIND_BLOCK,'p',iWrite,t)
    end if
    call showProg(t/tFinal)
  end do
  write(*,*),''
  deallocate(rho)
  deallocate(u)
  deallocate(p)
  deallocate(E)
  deallocate(rhou)
  deallocate(rhoE)
  deallocate(rhoNode)
  deallocate(uBlock)
  deallocate(uIntf)
  deallocate(pIntf)
  deallocate(Mass)
  deallocate(Mom)
  deallocate(Energy)
  deallocate(gradRho)
  deallocate(gradRhou)
  deallocate(gradRhoE)
  deallocate(u1d)
  close(13)
end program
