!----------------------------------------------------------------------------- best with 100 columns

!> simple (mostly integer) set operations
module moduleSimpleSetLogic
  private
  
  ! public procedures
  public findIntersection
  public findUnion
  public findComplement
  public applIntersection
  public applUnion
  public applComplement
  
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
    integer,allocatable::findUnion(:) !< union of a and b
    
    allocate(findUnion(size(a)))
    findUnion(:)=a(:)
    do i=1,size(b)
      if(.not.any(a(:)==b(i)))then
        call pushArr(findUnion,b(i))
      end if
    end do
  end function
  
  !> find the complement of b in a
  pure function findComplement(a,b)
    integer,intent(in)::a(:) !< array a
    integer,intent(in)::b(:) !< array b
    integer,allocatable::findComplement(:) !< complement of b in a
    logical mask(size(a))
    
    mask(:)=.false.
    do i=1,size(b)
      mask(:)=mask(:).or.(a(:)==b(i))
    end do
    mask(:)=.not.mask(:)
    allocate(findComplement(count(mask)))
    j=0
    do i=1,size(a)
      if(mask(i))then
        j=j+1
        findComplement(j)=a(i)
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
    
    if(.not.allocated(a))then
      allocate(a(0))
    end if
    temp=findUnion(a,b)
    deallocate(a)
    call move_alloc(temp,a)
  end subroutine
  
  !> apply to array a the complement of b
  pure subroutine applComplement(a,b)
    integer,intent(inout),allocatable::a(:) !< array a
    integer,intent(in)::b(:) !< array b
    integer,allocatable::temp(:)
    
    temp=findComplement(a,b)
    deallocate(a)
    call move_alloc(temp,a)
  end subroutine
  
end module
