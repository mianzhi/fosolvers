!----------------------------------------------------------------------------- best with 100 columns

function fixPt()
  use modNumerics
  integer fixPt
  type(fixPtEq)::p
  
  fixPt=0
  call p%init(3,func,maa=2)
  call p%clear()
  
contains
  
  function func(cx,cy,dat)
    use iso_c_binding
    use modNumerics
    type(C_PTR),value::cx,cy,dat
    integer(C_INT)::func
    double precision::x(3)
    double precision::y(3)
    
    call extractVector(cx,x)
    y(1:3)=2d0*x(1:3)-[1d0,2d0,3d0]
    call exportVector(y,cy)
    func=0
  end function
  
end function
