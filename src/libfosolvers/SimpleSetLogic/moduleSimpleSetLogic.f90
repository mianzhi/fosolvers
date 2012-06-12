!----------------------------------------------------------------------------- best with 100 columns

!> simple (mostly integer) set operations
module moduleSimpleSetLogic
  private
  
  ! public procedures
  public findIntersection
  public findUnion
  public applIntersection
  public applUnion
  
contains
  
  !> find the intersection of array a and b
  pure function findIntersection(a,b)
    use moduleBasicDataStruct
    integer,intent(in)::a(:) !< array a
    integer,intent(in)::b(:) !< array b
    integer,allocatable::findIntersection(:) !< intersection of a and b
    
    do i=1,size(b)
      if(any(a(:)==b(i)))then
        call pushArr(findIntersection,b(i))
      end if
    end do
    if(.not.allocated(findIntersection))then
      allocate(findIntersection(0))
    end if
  end function
  
  !> find the union of array a and b
  pure function findUnion(a,b)
    use moduleBasicDataStruct
    integer,intent(in)::a(:) !< array a
    integer,intent(in)::b(:) !< array b
    integer,allocatable::findUnion(:) !< intersection of a and b
    
    allocate(findUnion(size(a)))
    findUnion(:)=a(:)
    do i=1,size(b)
      if(.not.any(a(:)==b(i)))then
        call pushArr(findUnion,b(i))
      end if
    end do
  end function
  
  !> apply to array a intersection operation with b
  pure subroutine applIntersection(a,b)
    use moduleBasicDataStruct
    integer,intent(inout),allocatable::a(:) !< array a
    integer,intent(in)::b(:) !< array b
    integer,allocatable::temp(:)
    
    temp=findIntersection(a,b)
    deallocate(a)
    call move_alloc(temp,a)
  end subroutine
  
  !> apply to array a union operation with b
  pure subroutine applUnion(a,b)
    use moduleBasicDataStruct
    integer,intent(inout),allocatable::a(:) !< array a
    integer,intent(in)::b(:) !< array b
    integer,allocatable::temp(:)
    
    temp=findUnion(a,b)
    deallocate(a)
    call move_alloc(temp,a)
  end subroutine
  
end module