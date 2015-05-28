!----------------------------------------------------------------------------- best with 100 columns

function readCondTab1()
  use modFileIO
  use modCondition
  integer readCondTab1
  type(condTab)::cond
  
  readCondTab1=0
  open(10,file='data/condition_simple',action='read')
  call readCondTab(10,cond)
  if(any(cond%gid/=[10,20]))then
    readCondTab1=readCondTab1+1
  end if
  if(any(cond%t/=[1,2]))then
    readCondTab1=readCondTab1+10
  end if
  if(any(abs(cond%p-reshape([300d0,1d5,20d0,0d0,0d0,0.8d5,0d0,0d0,0d0,0d0],&
  &                         shape(cond%p)))>1d-13))then
    readCondTab1=readCondTab1+100
  end if
  call cond%clear()
  close(10)
end function
