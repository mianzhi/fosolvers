!----------------------------------------------------------------------------- best with 100 columns

!> foburn1d main program
program foburn1d
  use miscBurn1D
  use moduleGrid1D
  use moduleFVMConvect
  use moduleFVMDiffus
  use moduleFVMGrad
  use moduleInterpolation
  use moduleCLIO
  
  call grid%genUniform(0d0,1d0,100)
  call setEnv()
  t=0d0
  dt=1d-6
  tFinal=2d-2!8d-4
  ! initial state
  gamm=1.4d0
  Dm=1d-4
  nu=1d-4
  alpha=1d-4
  Y(:)=1d0
  u(:)=0d0
  p(:)=1d5 !< inital pressure
  !p(1:2)=1d6
  Temp(:)=300d0 !< initial temperature
  Temp(1:4)=2000d0
  rho(:)=p(:)/300d0/Temp(:) !TODO:rho=rho(p,T), R=300
  IE(:)=300d0*Temp(:)/(gamm-1d0) !TODO:IE=IE(p,T), R=300
  uCell=itplNode2Cell(u,grid)
  E(:)=IE(:)+uCell(:)**2d0/2d0
  rhoNode=itplCell2Node(rho,grid)
  rhou(:)=rhoNode(:)*u(:)
  rhoE(:)=rho(:)*E(:)
  Mass(:)=rho(:)*grid%CellWidth(:)
  Mom(:)=rhou(:)*grid%NodeWidth(:)
  Energy(:)=rhoE(:)*grid%CellWidth(:)
  IEnergy(:)=IE(:)*Mass(:)
  ! advance in time
  do while(t<tFinal)
    ! boundary sources
    forall(i=2:grid%nNode-1)
      Mom(i)=Mom(i)+dt*(p(NlC(i))-p(NrC(i)))
    end forall
    Mom(:)=Mom(:)+dt/3d0*findDiffus(nu*[(1,i=1,grid%nCell)],BIND_NODE,u,grid)
    pNode=itplCell2Node(p,grid)
    forall(i=1:grid%nCell)
      Energy(i)=Energy(i)+dt*&
      &  ((pNode(ClN(i))+nu*(uCell(ClC(i))-uCell(i))/grid%NodeWidth(ClN(i))/3d0)*u(ClN(i))&
      &  -(pNode(CrN(i))+nu*(uCell(CrC(i))-uCell(i))/grid%NodeWidth(CrN(i))/3d0)*u(CrN(i)))
    end forall
    burnR(:)=1d9*Y(:)*exp(-1.8d4/Temp(:))
    burnR(:)=min(burnR(:),Y(:)/dt)
    Energy(:)=Energy(:)+dt*rho(:)*burnR(:)*grid%CellWidth*1.2d6
    !write(*,*),rho*grid%CellWidth
    !write(*,*),Energy
    !write(*,*),burnR
    !stop
    Energy(:)=Energy(:)+dt*findDiffus(alpha*[(1,i=1,grid%nCell)],BIND_CELL,Temp,grid)
    Y(:)=(Y(:)*Mass(:)+dt*findDiffus(Dm*rho(:),BIND_CELL,Y,grid))/Mass(:)
    Y(:)=Y(:)-dt*burnR(:)
    Y(:)=max(Y(:),[(0d0,i=1,grid%nCell)])
    ! closed boundary
    u(1)=0d0
    u(grid%nNode)=0d0
    Mom(1)=0d0
    Mom(grid%nNode)=0d0
    ! recover states
    rhoNode=itplCell2Node(rho,grid)
    rhou(:)=Mom(:)/grid%NodeWidth(:)
    u(:)=rhou(:)/rhoNode(:)
    uCell(:)=itplNode2Cell(u,grid)
    rhoE(:)=Energy(:)/grid%CellWidth(:)
    ! convection
    rhox=findGrad(rho,BIND_CELL,grid)
    rhoux=findGrad(rhou,BIND_NODE,grid)
    rhoEx=findGrad(rhoE,BIND_CELL,grid)
    Mass(:)=Mass(:)+dt*findConvect(rho,u,BIND_NODE,grid,rhox)
    do i=2,grid%nNode-1 !TODO: needs a node convection scheme
      if(rhou(i)>0d0)then
        Mom(i)=Mom(i)+dt*(rhou(i-1)-rhou(i))
      else
        Mom(i)=Mom(i)+dt*(rhou(i)-rhou(i+1))
      end if
    end do
    Energy(:)=Energy(:)+dt*findConvect(rhoE,u,BIND_NODE,grid,rhoEx)
    Y(:)=Y(:)+dt*findConvect(Y,u,BIND_NODE,grid,Yx)
    ! update the new intensive state and internal energy
    rhou(:)=Mom(:)/grid%NodeWidth(:)
    rhoE(:)=Energy(:)/grid%CellWidth(:)
    rhoNode=itplCell2Node(rho,grid)
    u(:)=rhou(:)/rhoNode(:)
    E(:)=rhoE(:)/rho(:)
    uCell=itplNode2Cell(u,grid)
    IEnergy(:)=Energy(:)-Mass(:)*uCell(:)**2d0/2d0
    IE(:)=IEnergy/Mass(:)
    Temp(:)=IE(:)/300d0*(gamm-1d0)
    p(:)=rho(:)*300d0*Temp(:)
    rho(:)=Mass(:)/grid%CellWidth(:)
    !call showProg(t/tFinal)
    t=t+dt
  end do
  
  do i=1,grid%nCell
    write(*,'(10g11.3)'),grid%CellPos(i),p(i),rho(i),u(i),Temp(i),Y(i)
  end do
  ! clean up
  call clearEnv()
  
end program
