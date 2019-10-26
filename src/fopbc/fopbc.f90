!----------------------------------------------------------------------------- best with 100 columns

!> the pressure-based coupled solver
program fopbc
  use modPbc
  use modFileIO
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
    visc(:)=2d-4! TODO update transport properties according to state0
    cond(:)=3d-1! TODO update transport properties according to state0
    call preSolve()
    do nItOuter=1,MAXIT_OUTER
      rho1(:)=rho(:)
      p1(:)=p(:)
      call solvePBC()
      if(needRetry) exit
      call solveEnergy()
      if(needRetry) exit
      if(maxval(abs(rho(1:grid%nC)-rho1(1:grid%nC)))/rhoScale<=1d-4)then
        exit
      else if(nItOuter==MAXIT_OUTER)then
        needRetry=.true.
        exit
      end if
    end do
    if(needRetry)then
      write(*,'(a,i2,a,i2,a,i2,a,i2,a)')'[W] failed step, nIt[PBC,Energy,Outer]: [',&
      &                                 nItPBC,',',nItEnergy,',',nItOuter,']'
      cycle
    end if
    t=t+dt
    call postSolve()
    if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      write(*,'(a)')'[i] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
  end do
  call clear()
end program
