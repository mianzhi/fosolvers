!----------------------------------------------------------------------------- best with 100 columns

function otGrid1()
  use modOtGrid
  integer otGrid1
  type(otGrid)::grid
  
  otGrid1=0
  call grid%init([0d0,0d0,0d0],1d0,lvl=5)
  if(grid%nC/=32768)then
    otGrid1=otGrid1+1
  end if
  do i=1,grid%nC
    if(grid%lvl(i)/=5)then
      otGrid1=otGrid1+10
      exit
    end if
  end do
  do i=1,grid%nC
    if(grid%oid(i)/=int(i+o'400000'-1,kind=8))then
      otGrid1=otGrid1+100
      exit
    end if
  end do
  if(count(grid%neib(:,:)==0)/=6144)then
    otGrid1=otGrid1+1000
  end if
  if(maxval(grid%neib(:,:))/=32768)then
    otGrid1=otGrid1+2000
  end if
  if(minval(grid%neib(:,:))<0)then
    otGrid1=otGrid1+4000
  end if
  call grid%clear()
end function
