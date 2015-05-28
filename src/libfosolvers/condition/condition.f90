!----------------------------------------------------------------------------- best with 100 columns

!> condition module
module modCondition
  private
  
  !> table of conditions
  type,public::condTab
    integer,allocatable::gid(:) !< geometric group identifier
    integer,allocatable::t(:) !< type of condition
    double precision,allocatable::p(:,:) !< parameters of condition
  contains
    procedure,public::init=>initCondTab
    procedure,public::clear=>clearCondTab
    !FIXME:final::purgeCondTab
  end type
  
contains
  
  !> initialize this condTab
  elemental subroutine initCondTab(this,n,m)
    class(condTab),intent(inout)::this !< this condTab
    integer,intent(in)::n !< number conditions
    integer,intent(in)::m !< maximum number of parameters per condition
    
    call this%clear()
    allocate(this%gid(n))
    allocate(this%t(n))
    allocate(this%p(m,n))
    this%gid(:)=0
    this%t(:)=0
    this%p(:,:)=0d0
  end subroutine
  
  !> clear this condTab
  elemental subroutine clearCondTab(this)
    class(condTab),intent(inout)::this !< this condTab
    
    if(allocated(this%gid)) deallocate(this%gid)
    if(allocated(this%t)) deallocate(this%t)
    if(allocated(this%p)) deallocate(this%p)
  end subroutine
  
  !> destructor of condTab
  elemental subroutine purgeCondTab(this)
    type(condTab),intent(inout)::this !< this condTab
    
    call this%clear()
  end subroutine
  
end module
