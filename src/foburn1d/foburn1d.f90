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
  dt=1d-7
  tFinal=1d0!8d-4
  ! initial state
  gamm=1.4d0
  Dm=1d-4
  nu=1d-4
  alpha=1d-4
  Y(:)=1d0
  u(:)=0d0
  p(:)=1d5 !< inital pressure
  p(1:50)=6d5
  Temp(:)=300d0 !< initial temperature
  !Temp(1:2)=600d0
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
    ! find new conserved extensive variable and Y of reactant
    rhox=findGrad(rho,BIND_CELL,grid)
    Mass(:)=Mass(:)+dt*findConvect(rho,u,BIND_NODE,grid,rhox)
    rhoux=findGrad(rhou,BIND_NODE,grid)
    !TODO: needs node convection scheme
    !Mom(:)=itplCell2Node(itplNode2Cell(Mom(:),grid)+dt*findConvect(rho,u,BIND_NODE,grid,rhoux),grid)
    do i=2,grid%nCell-1
      if(rhou(i)>0d0)then
        Mom(i)=Mom(i)+dt*(rhou(i-1)-rhou(i))
      else
        Mom(i)=Mom(i)+dt*(rhou(i)-rhou(i+1))
      end if
    end do
    forall(i=1:grid%nNode)
      Mom(i)=Mom(i)+dt*(p(NlC(i))-p(NrC(i)))
    end forall
    u(1)=0d0
    u(grid%nNode)=0d0
    Mom(1)=0d0
    Mom(grid%nNode)=0d0
    pNode=itplCell2Node(p,grid)
    rhoEx=findGrad(rhoE,BIND_CELL,grid)
    Energy(:)=Energy(:)+dt*findConvect(rhoE,u,BIND_NODE,grid,rhoEx)
    forall(i=1:grid%nCell)
      Energy(i)=Energy(i)+dt*(pNode(ClN(i))*u(ClN(i))-pNode(CrN(i))*u(CrN(i)))
    end forall
    !Energy(:)=Energy(:)+dt*findDiffus(alpha*[(1,i=1,grid%nCell)],BIND_CELL,Temp,grid)
    Y(:)=(Y(:)*Mass(:)+dt*findDiffus(Dm*rho(:),BIND_CELL,Y,grid))/Mass(:)
    ! update the new intensive state and internal energy
    rho(:)=Mass(:)/grid%CellWidth(:)
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
    !call showProg(t/tFinal)
    t=t+dt
  end do
  
  do i=1,grid%nCell
    write(*,'(10g11.3)'),grid%CellPos(i),p(i),rho(i),u(i)
  end do
  ! clean up
  call clearEnv()
  
end program
