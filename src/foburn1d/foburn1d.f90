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
  
  width=0.5d0
  call grid%genUniform(0d0,width,100)
  call setEnv()
  t=0d0
  dt=1d-3
  tFinal=50d0
  ! initial state
  gamm=1.4d0
  Dm=5d-4 ![m^2/s]
  alpha=10d-4 ![m^2/s]
  mw=30d-3 ![kg/mol]
  R=RU/mw ![J/kg/K]
  Cv=R/(gamm-1d0) ![J/kg/K]
  Cp=Cv*gamm ![J/kg/K]
  Q=50d3 ![J/mol reactant]
  Y(:)=1d0
  u(:)=0d0
  p=1d5 !< inital pressure
  !p(1:2)=1d6
  Temp(:)=300d0 !< initial temperature
  !Temp(1:4)=2000d0
  rho(:)=p/R/Temp(:)
  Mass(:)=rho(:)*grid%CellWidth(:)
  Temp(1)=2d3
  ! advance in time
  do while(t<tFinal)
    ! diffusion and burn
    burnR(:)=1d9*(Y(:)*rho(:)/mw)*exp(-1.5d4/Temp(:))
    burnR(:)=min(burnR(:),Y(:)*rho(:)/mw/dt)
    Temp(:)=Temp(:)+dt*findDiffus(alpha*[(1d0,i=1,grid%nCell)],BIND_CELL,Temp,grid)
    Temp(:)=Temp(:)+dt*Q/Cp*burnR(:)/rho(:)
    Y(:)=Y(:)+dt*findDiffus(Dm*[(1d0,i=1,grid%nCell)],BIND_CELL,Y,grid)
    Y(:)=Y(:)-dt*mw*burnR(:)/rho(:)
    ! move gas
    rho(:)=p/R/Temp(:)
    grid%CellWidth(:)=Mass(:)/rho(:)
    do i=2,grid%nNode
      grid%NodePos(i)=grid%NodePos(NlN(i))+grid%CellWidth(NlC(i))
    end do
    p=p*(grid%NodePos(grid%nNode)/width)**gamm
    Temp(:)=Temp(:)*(grid%NodePos(grid%nNode)/width)**(gamm-1d0)
    grid%NodePos(:)=grid%NodePos(:)*width/grid%NodePos(grid%nNode)
    call grid%update()
    rho(:)=Mass(:)/grid%CellWidth(:)
    t=t+dt
  end do
  
  do i=1,grid%nCell
    write(*,'(10g11.3)'),grid%CellPos(i),Temp(i),Y(i),rho(i)
  end do
  ! clean up
  call clearEnv()
  
end program
