!----------------------------------------------------------------------------- best with 100 columns

!> the pressure-based coupled solver
program fopbc
  use modPbc
  character(20)::tmpStr
  double precision::dt !< time step size [s]
  integer::ier !< error signal
  
  call init()
  write(tmpStr,*),iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    call deriveState(rho,rhou,rhoE,Y,p,u,temp)
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    !call preSolve(dt,xscale,rscale)
    dt=1d-4
    call fillVar(p,u,temp,x)
    !call fkinsol(x,1,nocale,nocale,ier)
    call extractVar(x,p,u,temp)
    ! TODO solve Y transport here
    call recoverState(p,u,temp,Y,rho,rhou,rhoE)
    t=t+dt
    if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      !write(tmpStr,*),iOut
      !call syncState(y)
      !call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
    stop
  end do
  call clear()
end program
