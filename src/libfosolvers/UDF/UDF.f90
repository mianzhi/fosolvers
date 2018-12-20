!----------------------------------------------------------------------------- best with 100 columns

!> UDF module
module modUDF
  use iso_c_binding
  private
  
  integer,parameter::DIMS=3 !< dimensions
  
  !> table of UDFs
  type,public::UDFTab
    integer(C_LONG),allocatable::ptr(:) !< UDF pointers
  contains
    procedure,public::init=>initUDFTab
    procedure,public::clear=>clearUDFTab
    procedure,public::eval=>evalUDFTab
    final::purgeUDFTab
  end type
  
contains
  
  !> initialize this UDFTab
  subroutine initUDFTab(this,n)
    class(UDFTab),intent(inout)::this !< this UDFTab
    integer,intent(in)::n !< number UDFs
    
    call this%clear()
    allocate(this%ptr(n))
    this%ptr(:)=0
  end subroutine
  
  !> clear this UDFTab
  subroutine clearUDFTab(this)
    class(UDFTab),intent(inout)::this !< this UDFTab
    external::evaluator_destroy_
    
    if(allocated(this%ptr))then
      do i=1,size(this%ptr)
        call evaluator_destroy_(this%ptr(i))
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
    real(C_DOUBLE)::evaluator_evaluate_
    
    evalUDFTab=evaluator_evaluate_(this%ptr(k),4,'x y z t',[x(1),x(2),x(3),t])
  end function
  
  !> destructor of UDFTab
  subroutine purgeUDFTab(this)
    type(UDFTab),intent(inout)::this !< this UDFTab
    
    call this%clear()
  end subroutine
  
end module
