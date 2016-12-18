!----------------------------------------------------------------------------- best with 100 columns

!> the pressure-based coupled solver
program fopbc
  use modPbc
  use modFileIO
  character(20)::tmpStr
  integer::ier !< error signal
  double precision,allocatable::r(:),rRho(:),rRhou(:,:),rRhoE(:)
  
  call init()
  allocate(r(nEq))
  allocate(rRho(grid%nE))
  allocate(rRhou(DIMS,grid%nE))
  allocate(rRhoE(grid%nE))
  rRho=0d0
  rRhou=0d0
  rRhoE=0d0
  write(tmpStr,*)iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    call setBC()
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    call preSolve()
    write(*,*)'start',t
    call fkinsol(x,1,xscale,rscale,ier)
    write(*,*)'done'
    call extractVar(x,p,u,temp)
    ! TODO solve Y transport here
    call recoverState(p,u,temp,Y,rho,rhou,rhoE)
    call setBC()
    t=t+dt
    if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      write(*,*)'writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      call fkfun(x,r,ier)
      do i=1,grid%nC
        j=(i-1)*5
        rRho(i)=r(j+1)
        rRhou(:,i)=r(j+2:j+4)
        rRhoE(i)=r(j+5)
      end do
      open(10,file='r_'//trim(adjustl(tmpStr))//'.vtk',action='write')
        call writeVTK(10,grid)
        call writeVTK(10,grid,E_DATA)
        call writeVTK(10,'rRho',rRho)
        call writeVTK(10,'rRhoU',rRhou)
        call writeVTK(10,'rRhoE',rRhoE)
      close(10)
      tNext=tNext+tInt
    end if
  end do
  call clear()
end program
