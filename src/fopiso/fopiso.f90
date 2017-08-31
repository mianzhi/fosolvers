!----------------------------------------------------------------------------- best with 100 columns

!> the compressible PISO solver
program fopiso
  use modPiso
  use modFileIO
  use modDiffusion
  character(20)::tmpStr
  
  call init()
  write(tmpStr,*)iOut
  write(*,'(a)')'[I] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    if(needRetry)then
      nRetry=nRetry+1
      needRetry=.false.
    end if
    if(nRetry>0)then
      call loadState0()
    end if
    call setBC()
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    call preSolve()
    call predictMomentum()
    if(needRetry) cycle
    call findDiff(grid,p,[(1d0,i=1,grid%nC)],laP)
    temp1(:)=temp(:)
    do nItPISO=1,MAXIT_PISO
      laP1(:)=laP(:)
      presF1(:,:)=presF(:,:)
      p1(:)=p(:)
      rho1(:)=rho(:)
      call solvePressure()
      if(needRetry) exit
      call correctMomentum()
      call solveDensity()
      if(needRetry) exit
      forall(i=1:grid%nC)
        u(:,i)=rhou(:,i)/rho(i)
      end forall
      temp1(:)=temp(:)
      call solveEnergy()
      if(needRetry) exit
      write(*,*)'PISO',maxval(abs(rho(1:grid%nC)-rho1(1:grid%nC)))/rhoScale
      if(maxval(abs(rho(1:grid%nC)-rho1(1:grid%nC)))/rhoScale<=RTOL_DENSITY)then
        exit
      elseif(nItPISO==MAXIT_PISO)then
        needRetry=.true.
      end if
    end do
    if(needRetry)then
      write(*,'(a,i2,a,g12.6)')'[W] start retry No. ',nRetry,' at t: ',dt
      cycle
    else
      nRetry=0
    end if
    t=t+dt
    call postSolve()
    if(t*(1d0+tiny(1d0))>=tNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      write(*,'(a)')'[I] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
  end do
  call clear()
end program
