!----------------------------------------------------------------------------- best with 100 columns

!> the compressible PISO solver
program fopiso
  use modPiso
  use modFileIO
  use modDiffusion
  character(20)::tmpStr
  
  call init()
  write(tmpStr,*)iOut
  write(*,'(a)')'[i] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    if(needRetry)then
      nRetry=nRetry+1
      needRetry=.false.
      call loadState0()
    else
      nRetry=0
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
    do nItPISO=1,MAXIT_PISO
      rho1(:)=rho(:)
      laP1(:)=laP(:)
      presF1(:,:)=presF(:,:)
      p1(:)=p(:)
      call solvePressure()
      if(needRetry) exit
      call correctMomentum()
      write(*,*)'rho and p error: ',maxval(abs(rho(1:grid%nC)-rho1(1:grid%nC)))/rhoScale,&
      &          maxval(abs(p(1:grid%nC)-p1(1:grid%nC)))/pScale
      temp1(:)=temp(:)
      call solveEnergy()
      if(needRetry) exit
      if(maxval(abs(rho(1:grid%nC)-rho1(1:grid%nC)))/rhoScale<=RTOL_DENSITY)then
        exit
      else if(nItPISO==MAXIT_PISO)then
        needRetry=.true.
        exit
      end if
    end do
    if(needRetry) cycle
    p(:)=rho(:)*Rgas*temp(:)
    t=t+dt
    if(t*(1d0+tiny(1d0))>=tNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      write(*,'(a)')'[i] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
    call postSolve()
  end do
  call clear()
end program
