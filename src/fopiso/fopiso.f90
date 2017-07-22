!----------------------------------------------------------------------------- best with 100 columns

!> the compressible PISO solver
program fopiso
  use modPiso
  use modFileIO
  character(20)::tmpStr
  integer::ier !< error signal
  
  call init()
  write(tmpStr,*)iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    write(*,*)'start',t
    call setBC()
    call recordState0()
    ! TODO update transport properties according to state0
    visc(:)=20d-6
    cond(:)=30d-3
    call preSolve()
    call predictMomentum()
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
