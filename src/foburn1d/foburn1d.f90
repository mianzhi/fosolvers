!----------------------------------------------------------------------------- best with 100 columns

!> foburn1d main program
program foburn1d
  use miscBurn1D
  use moduleGrid1D
  use moduleFVMConvect
  use moduleFVMDiffus
  use moduleFVMGrad
  
  call grid%genUniform(0d0,1d0,100)
  call setEnv()
  t=0d0
  dt=1d-5
  tFinal=1d0
  ! initial state
  gamm=1.4d0
  Y(:)=1d0
  u(:)=0d0
  p(:)=1d6 !< inital pressure
  Temp(:)=500d0 !< initial temperature
  rho(:)=p(:)/300d0/Temp(:) !TODO:rho=rho(p,T), R=300
  IE(:)=300d0*Temp(:)/(gamm-1d0) !TODO:IE=IE(p,T), R=300
  E(:)=IE(:)+u(:)**2d0/2d0
  rhou(:)=rho(:)*u(:)
  rhoE(:)=rho(:)*E(:)
  Mass(:)=rho(:)*grid%CellWidth(:)
  Mom(:)=rhou(:)*grid%CellWidth(:)
  Energy(:)=rhoE(:)*grid%CellWidth(:)
  IEnergy(:)=IE(:)*Mass(:)
  ! advance in time
  do while(t<tFinal)
    ! find new conserved extensive variable and Y of reactant
    rhox=findGrad(rho,BIND_CELL,grid)
    Mass(:)=Mass(:)+dt*findConvect(rho,u,BIND_CELL,grid,rhox)
    Mom(:)=Mom(:)+dt*0d0
    Energy(:)=Energy(:)+dt*0d0
    Y(:)=(Y(:)*Mass(:)+dt*0d0)/Mass(:)
    ! update the new intensive state and internal energy
    rho(:)=Mass(:)/grid%CellWidth(:)
    rhou(:)=Mom(:)/grid%CellWidth(:)
    rhoE(:)=Energy(:)/grid%CellWidth(:)
    u(:)=rhou(:)/rho(:)
    E(:)=rhoE(:)/rho(:)
    IEnergy(:)=Energy(:)-Mass(:)*u(:)**2d0/2d0
    IE(:)=IEnergy/Mass(:)
    Temp(:)=IE(:)/300d0*(gamm-1d0)
    p(:)=rho(:)*300d0*Temp(:)
    t=t+dt
  end do
  ! clean up
  call clearEnv()
  
end program
