!----------------------------------------------------------------------------- best with 100 columns

!> UDF module
module modUDF
  use iso_c_binding
  private
  
  integer,parameter::DIMS=3 !< dimensions
  
  !> table of UDFs
  type,public::UDFTab
    type(C_PTR),allocatable::ptr(:) !< UDF pointers
  contains
    procedure,public::init=>initUDFTab
    procedure,public::clear=>clearUDFTab
    procedure,public::eval=>evalUDFTab
    final::purgeUDFTab
  end type
  
  !> interface to libmatheval C functions
  interface
  
    !> create an evaluator
    function evaluator_create(string) bind(c,name='evaluator_create')
      use iso_c_binding
      character(kind=C_CHAR)::string(*)
      type(C_PTR)::evaluator_create
    end function
    
    !> destroy an evaluator
    subroutine evaluator_destroy(evaluator) bind(c,name='evaluator_destroy')
      use iso_c_binding
      type(C_PTR),value::evaluator
    end subroutine
    
    !> evaluate an evaluator
    function evaluator_evaluate(evaluator,count,names,values) bind(c,name='evaluator_evaluate')
      use iso_c_binding
      type(C_PTR),value::evaluator
      integer(C_INT),value::count
      type(C_PTR)::names(*)
      real(C_DOUBLE)::values(*)
      real(C_DOUBLE)::evaluator_evaluate
    end function
    
  end interface
  
  public::evaluator_create
  
contains
  
  !> initialize this UDFTab
  subroutine initUDFTab(this,n)
    class(UDFTab),intent(inout)::this !< this UDFTab
    integer,intent(in)::n !< number UDFs
    
    call this%clear()
    allocate(this%ptr(n))
    this%ptr(:)=C_NULL_PTR
  end subroutine
  
  !> clear this UDFTab
  subroutine clearUDFTab(this)
    class(UDFTab),intent(inout)::this !< this UDFTab
    
    if(allocated(this%ptr))then
      do i=1,size(this%ptr)
        call evaluator_destroy(this%ptr(i))
      end do
      deallocate(this%ptr)
    end if
  end subroutine
  
  !> evaluate the k th UDF in this UDFTab
  function evalUDFTab(this,k,x,t)
    class(UDFTab),intent(in)::this !< this UDFTab
    integer,intent(in)::k !< index of the UDF
    double precision,intent(in)::x(DIMS) !< position
    double precision,intent(in)::t !< time
    double precision::evalUDFTab !< result
    character(len=2,kind=C_CHAR),target::xName,yName,zName,tName
    type(C_PTR)::inputNames(4)
    real(C_DOUBLE)::inputVals(4)
    
    xName='x'//C_NULL_CHAR
    yName='y'//C_NULL_CHAR
    zName='z'//C_NULL_CHAR
    tName='t'//C_NULL_CHAR
    inputNames=[c_loc(xName),c_loc(yName),c_loc(zName),c_loc(tName)]
    inputVals=[x(1),x(2),x(3),t]
    evalUDFTab=dble(evaluator_evaluate(this%ptr(k),4,inputNames,inputVals))
  end function
  
  !> destructor of UDFTab
  subroutine purgeUDFTab(this)
    type(UDFTab),intent(inout)::this !< this UDFTab
    
    call this%clear()
  end subroutine
  
end module
