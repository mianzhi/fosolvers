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
    visc(:)=20d-6! TODO update transport properties according to state0
    cond(:)=30d-3! TODO update transport properties according to state0
    call preSolve()
    do nItOuter=1,3!MAXIT_OUTER
      call solvePBC()
      if(needRetry) exit
      call solveEnergy()
      if(needRetry) exit
      ! TODO exit check
    end do
    if(needRetry) cycle
    t=t+dt
    if(t+tiny(1d0)>=tNext)then
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
