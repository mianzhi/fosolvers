!----------------------------------------------------------------------------- best with 100 columns

!> the pressure-based coupled solver
program fopbc
  use modPbc
  character(20)::tmpStr
  integer::ier !< error signal
  
  call init()
  write(tmpStr,*),iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  !do while(t<tFinal)
    call deriveState(rho,rhou,rhoE,Y,p,u,temp)
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    call preSolve(dt,xscale,rscale)
    call fillVar(p,u,temp,x)
    write(*,*),'start'
    call fkinsol(x,1,noscale,noscale,ier)
    write(*,*),'done solve'
    call extractVar(x,p,u,temp)
    ! TODO solve Y transport here
    call recoverState(p,u,temp,Y,rho,rhou,rhoE)
    t=t+dt
    !if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      write(tmpStr,*),iOut
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    !end if
  !end do
  call clear()
end program
