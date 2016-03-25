!----------------------------------------------------------------------------- best with 100 columns

!> the pressure-based coupled solver
program fopbc
  use modPbc
  
  !call init()
  !write(tmpStr,*),iOut
  !call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    !call fcvode(tNext,t,y,1,ier)
    if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      !write(tmpStr,*),iOut
      !call syncState(y)
      !call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
  end do
  !call clear()
end program
