!----------------------------------------------------------------------------- best with 100 columns

!> the compressible PISO solver
program fopiso
  use modPiso
  use modFileIO
  use modDiffusion
  character(20)::tmpStr
  
  call init()
  write(tmpStr,*)iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do l=1,40!while(t<tFinal)
    write(*,*)'start',t
    call setBC()
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    call preSolve()
    call predictMomentum()
    call findDiff(grid,p,[(1d0,i=1,grid%nC)],laP)
    temp1(:)=temp(:)
    
    !iOut=iOut+1
    !write(tmpStr,*)iOut
    !write(*,*)'writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
    !call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
    
    do k=1,5!ITMAX_PISO
      laP1(:)=laP(:)
      presF1(:,:)=presF(:,:)
      p1(:)=p(:)
      rho1(:)=rho(:)
      call solvePressure()
      call correctMomentum()
      !call predictMomentum()
      !rho(1:grid%nC)=p(1:grid%nC)/287.058d0/temp(1:grid%nC) ! TODO gas property
      call solveDensity()
      !p(1:grid%nC)=rho(1:grid%nC)*287.058d0*temp(1:grid%nC)
      forall(i=1:grid%nC)
        u(:,i)=rhou(:,i)/rho(i)
      end forall
      temp1(:)=temp(:)
      call solveEnergy()
      
      !iOut=iOut+1
      !write(tmpStr,*)iOut
      !write(*,*)'writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      !call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      
    end do
    t=t+dt
    write(*,*)'done'
    if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      write(*,*)'writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
  end do
  call clear()
end program
