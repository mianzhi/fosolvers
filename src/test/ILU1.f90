!----------------------------------------------------------------------------- best with 100 columns

function ILU1()
  use modSparse
  integer::ILU1
  type(ILU)::p
  integer::iA(7),jA(7)
  double precision::A(7)
  double precision::rhs(4),x(4)
  
  ILU1=0
  call p%init(4,7,4)
  iA=[1,2,2,3,3,4,4]
  jA=[1,2,1,3,4,4,2]
  A=[1d0,2d0,-1d0,5d0,-4d0,2d0,-1d0]
  rhs=[1d0,2d0,3d0,4d0]
  call p%setCOO(iA,jA,A,7)
  call p%fact(0d0)
  call p%solve(rhs,x)
  if(norm2(x-[1d0,3d0/2d0,14d0/5d0,11d0/4d0])>1d-9)then
    ILU1=1
  end if
end function
