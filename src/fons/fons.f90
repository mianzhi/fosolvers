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
  double precision pWork !< pressure work done on block surface
  external::resMom
  external::resEnergy
  external::resPressure
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
  call setEnv()
  allocate(u1d(DIMS*grid%nNode))
  ! simulation control
  dt=1d-5
  ! initial value of variables
  call grid%updateBlockPos()
  gamm=1.4d0
  u(:,:)=0d0
  rhou(:,:)=0d0
  forall(i=1:grid%nBlock)
    p(i)=merge(1d5,1d4,grid%BlockPos(1,i)<0.5d0)
    Temp(i)=500d0
    rho(i)=p(i)/200d0/Temp(i) !TODO:rho=rho(p,T), Ru=200
    IE(i)=200d0*Temp(i)/(gamm-1d0) !TODO:IE=IE(p,T), Ru=200
    E(i)=IE(i) !zero velocity
    rhoE(i)=rho(i)*E(i)
  end forall
  call grid%updateDualBlock()
  call grid%updateBlockVol()
  Mass(:)=rho(:)*grid%BlockVol(:)
  forall(i=1:grid%nNode)
    Mom(:,i)=rhou(:,i)*grid%NodeVol(i)
  end forall
  IEnergy(:)=IE(:)*Mass(:)
  Energy(:)=E(:)*Mass(:)
  visc(:)=1d-3 !TODO: mu, k are functions of T,p
  viscRate(:)=-2d0/3d0
  thermK(:)=1d-4
  t=0d0
  iWrite=0
  ! write initial states
  open(13,file='rstNS.msh',status='replace')
  call writeGMSH(13,grid)
  call writeGMSH(13,rho,grid,BIND_BLOCK,'rho',iWrite,t)
  call writeGMSH(13,u,grid,BIND_NODE,'u',iWrite,t)
  call writeGMSH(13,p,grid,BIND_BLOCK,'p',iWrite,t)
  call writeGMSH(13,Temp,grid,BIND_BLOCK,'T',iWrite,t)
  ! advance in time
  do while(t<tFinal)
    call grid%updateDualBlock()
    call grid%updateBlockVol()
    call grid%updateFacetNorm()
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    ! Lagrangian step
    preU(:,:)=u(:,:)
    oldU(:,:)=u(:,:)
    oldP(:)=p(:)
    oldMom(:,:)=Mom(:,:)
    oldEnergy(:)=Energy(:)
    do l=1,10
      ! solve momentum equation for velocity using assumed pressure
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
      ! solve energy equation for temperature, exclude pressure work
      tao=findTao(u)
      uBlock=itplNode2Block(u,grid)
      uIntf=itplNode2Intf(u,grid)
      taoIntf=itplBlock2Intf(tao,grid)
      ProblemFunc=>resEnergy
      call solveNonlinear(Temp)
      forall(i=1:grid%nBlock)
        IE(i)=Temp(i)*200d0/(gamm-1) !TODO:IE=IE(p,T)
        E(i)=IE(i)+dot_product(uBlock(:,i),uBlock(:,i))/2d0
        rhoE(i)=rho(i)*E(i)
        IEnergy(i)=IE(i)*Mass(i)
        Energy(i)=E(i)*Mass(i)
      end forall
      ! remove from the momentum equation the effect of assumed pressure
      do i=1,grid%nNode
        do j=1,size(grid%NodeNeibBlock(i)%dat)
          Mom(:,i)=Mom(:,i)+dt*grid%NBAreaVect(i)%dat(:,j)*p(grid%NodeNeibBlock(i)%dat(j))
        end do
      end do
      ! couple pressure with fluid displacement, add pressure effects on momentum and energy
      ProblemFunc=>resPressure
      call solveNonlinear(p)
      do i=1,grid%nNode
        do j=1,size(grid%NodeNeibBlock(i)%dat)
          Mom(:,i)=Mom(:,i)-dt*grid%NBAreaVect(i)%dat(:,j)*p(grid%NodeNeibBlock(i)%dat(j))
          rhou(:,i)=Mom(:,i)/grid%NodeVol(i)
          u(:,i)=rhou(:,i)/rhoNode(i)
        end do
      end do
      pIntf=itplBlock2Intf(p,grid)
      do i=1,grid%nIntf
        m=grid%IntfNeibBlock(1,i)
        n=grid%IntfNeibBlock(2,i)
        pWork=dt*pIntf(i)*grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),uIntf(:,i))
        Energy(m)=Energy(m)-pWork
        Energy(n)=Energy(n)+pWork
      end do
      rhoE(:)=Energy(:)/grid%BlockVol(:)
      E(:)=rhoE(:)/rho(:)
      uBlock=itplNode2Block(u,grid)
      forall(i=1:grid%nBlock)
        IE(i)=E(i)-dot_product(uBlock(:,i),uBlock(:,i))/2d0
        IEnergy(i)=IE(i)*Mass(i)
      end forall
      write(*,*),l,maxval(norm2(u(:,:)-preU(:,:),1))
      if(maxval(norm2(u(:,:)-preU(:,:),1))<1d0)then
        exit
      else
        preU(:,:)=u(:,:)
        u(:,:)=oldU(:,:)
        Mom(:,:)=oldMom(:,:)
        Energy(:)=oldEnergy(:)
      end if
    end do
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
      IE(i)=E(i)-dot_product(uBlock(:,i),uBlock(:,i))/2d0
      IEnergy(i)=IE(i)*Mass(i)
      p(i)=IE(i)*rho(i)*(gamm-1d0) !TODO:p=p(rho,T)
      Temp(i)=IE(i)*(gamm-1)/200d0 !TODO:T=T(IE,p)
    end forall
    rho(:)=Mass(:)/grid%BlockVol(:)
    t=t+dt
    ! write results
    if(t/tWrite>=iWrite)then
      iWrite=iWrite+1
      call writeGMSH(13,rho,grid,BIND_BLOCK,'rho',iWrite,t)
      call writeGMSH(13,u,grid,BIND_NODE,'u',iWrite,t)
      call writeGMSH(13,p,grid,BIND_BLOCK,'p',iWrite,t)
      call writeGMSH(13,Temp,grid,BIND_BLOCK,'T',iWrite,t)
    end if
    call showProg(t/tFinal)
  end do
  write(*,*),''
  call clearEnv()
  deallocate(u1d)
  close(13)
end program
