!----------------------------------------------------------------------------- best with 100 columns

function matVecMul()
  use modSparse
  integer::matVecMul
  type(linEq)::p
  integer::iA(7),jA(7)
  double precision::A(7)
  double precision::y(4),x(4)
  
  matVecMul=0
  call p%initLinEq(4,7)
  iA=[1,2,2,3,3,4,4]
  jA=[1,2,1,3,4,4,2]
  A=[1d0,2d0,-1d0,5d0,-4d0,2d0,-1d0]
  x=[1d0,3d0/2d0,14d0/5d0,11d0/4d0]
  call p%setCOO(iA,jA,A,7)
  call p%mulVec(x,y)
  if(norm2(y-[1d0,2d0,3d0,4d0])>1d-9)then
    matVecMul=1
  end if
end function
