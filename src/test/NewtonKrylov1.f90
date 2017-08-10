!----------------------------------------------------------------------------- best with 100 columns

function NewtonKrylov1()
  use modNumerics
  integer::NewtonKrylov1
  type(NewtonKrylov)::p
  double precision::s(3)
  
  NewtonKrylov1=0
  call p%init(3,func)
  s(:)=[0d0,0d0,0d0]
  call p%solve(s)
  if(norm2(s-[-0.5d0,-0.5d0,3d0/8d0])>1d-9)then
    NewtonKrylov1=1
  end if
  call p%clear()
  
contains
  
  function func(cx,cy,dat)
    use iso_c_binding
    use modNumerics
    type(C_PTR),value::cx,cy,dat
    integer(C_INT)::func
    double precision,pointer::x(:)
    double precision,pointer::y(:)
    type(C_PTR)::foo
    
    call associateVector(cx,x)
    call associateVector(cy,y)
    y(1:3)=[2d0,4d0,-8d0]*x(1:3)+[1d0,2d0,3d0]
    func=0
    foo=dat
  end function
  
end function
