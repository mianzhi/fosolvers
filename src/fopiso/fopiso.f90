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
  do l=1,20!while(t<tFinal)
    write(*,*)'start',t
    call setBC()
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    call preSolve()
    call predictMomentum()
    call findDiff(grid,p,[(1d0,i=1,grid%nC)],laP)
    do k=1,5!ITMAX_PISO
      laP1(:)=laP(:)
      presF1(:,:)=presF(:,:)
      p1(:)=p(:)
      rho1(:)=rho(:)
      call solvePressure()
      call correctMomentum()
      rho(1:grid%nC)=p(1:grid%nC)/286.9d0/temp(1:grid%nC) ! TODO gas property
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
