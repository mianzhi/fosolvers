!----------------------------------------------------------------------------- best with 100 columns

function otGrid2()
  use modOtGrid
  integer otGrid2
  type(otGrid)::grid
  
  otGrid2=0
  call grid%init([0d0,0d0,0d0],1d0,lvl=2)
  if(norm2(grid%p(6)-[0.375d0,0.125d0,0.375d0])>tiny(1d0))then
    otGrid2=otGrid2+1
  end if
  do i=1,grid%nC
    if(abs(grid%h(i)-0.25d0)>tiny(1d0))then
      otGrid2=otGrid2+10
      exit
    end if
  end do
  do i=1,grid%nC
    if(abs(grid%a(i)-0.25d0**2)>tiny(1d0))then
      otGrid2=otGrid2+100
      exit
    end if
  end do
  do i=1,grid%nC
    if(abs(grid%v(i)-0.25d0**3)>tiny(1d0))then
      otGrid2=otGrid2+1000
      exit
    end if
  end do
  call grid%clear()
end function
