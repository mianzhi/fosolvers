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
    final::purgeCondTab
  end type
  
  !> map table of conditions to generic grid
  interface mapCondTab
    module procedure::mapCondTabPolyFvGrid
  end interface
  public::mapCondTab
  
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
  
  !> map table of conditions to polyFvGrid
  pure subroutine mapCondTabPolyFvGrid(grid,cTab,iCond)
    use modPolyFvGrid
    class(polyFvGrid),intent(in)::grid !< grid
    class(condTab),intent(in)::cTab !< condition table
    integer,allocatable,intent(inout)::iCond(:) !< index of condition at each element in grid
    
    if(allocated(iCond)) deallocate(iCond)
    allocate(iCond(grid%nE))
    iCond(:)=0
    do i=1,grid%nE
      do j=1,size(cTab%gid)
        if(grid%gid(i)==cTab%gid(j))then
          iCond(i)=j
          exit
        end if
      end do
    end do
  end subroutine
  
end module
