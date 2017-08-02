!----------------------------------------------------------------------------- best with 100 columns

function fixPt()
  use modNumerics
  integer fixPt
  type(fixPtEq)::p
  
  fixPt=0
  call p%init(3,func,maa=2)
  call p%clear()
  
contains
  
  subroutine func(x,y)
    double precision,intent(in)::x(*)
    double precision,intent(out)::y(*)
    
    y(1:3)=2d0*x(1:3)-[1d0,2d0,3d0]
  end subroutine
  
end function
