!----------------------------------------------------------------------------- best with 100 columns

function BDFNewtonKrylov1()
  use modNumerics
  use iso_c_binding
  integer::BDFNewtonKrylov1
  type(BDFNewtonKrylov)::p
  double precision::t,s(2)
  
  BDFNewtonKrylov1=0
  call p%init(2,func)
  call p%setTol(1d-8,1d-12)
  call p%setMaxSteps(-1)
  t=0d0
  s(:)=[1d0,0d0]
  call p%setIV(t,s)
  do while(t<2d0)
    call p%step(2d0,t,s)
  end do
  if(norm2(s-[cos(t),-sin(t)])>1d-6)then
    BDFNewtonKrylov1=1
  end if
  call p%clear()
  
contains
  
  function func(ct,cx,cdx,dat)
    use iso_c_binding
    use modNumerics
    real(C_DOUBLE),value::ct
    type(C_PTR),value::cx,cdx,dat
    integer(C_INT)::func
    double precision::t
    double precision,pointer::x(:)
    double precision,pointer::dx(:)
    type(C_PTR)::foo
    
    t=ct
    call associateVector(cx,x)
    call associateVector(cdx,dx)
    dx(1:2)=[x(2),-x(1)]
    func=0
    foo=dat
  end function
  
end function
