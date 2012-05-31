!----------------------------------------------------------------------------- best with 100 columns

!> Basic Data Structures
module moduleBasicDataStruct
  private
  
  !> heterogeneous 1D integer allocatable array
  type,public::typeHtr1DIArr
    integer,allocatable::dat(:) !< 1D integer allocatable array
  contains
    !FIXME:final::purgeHtr1DIArr
  end type
  
  !> reallocate a generic allocatable array
  interface reallocArr
    module procedure::realloc2DIArr
    module procedure::realloc1DDArr
    module procedure::realloc2DDArr
    module procedure::realloc1DHtr1DIArr
  end interface
  public::reallocArr
  
  !> extend a generic allocatable array
  interface extendArr
    module procedure::extendIArr
  end interface
  public::extendArr
  
  !> push a scalar or vector to the back of a generic allocatable array
  interface pushArr
    module procedure::pushIScalArr
    module procedure::pushIVectArr
  end interface
  public::pushArr
  
contains
  
  !> destructor of typeHtr1DIArr
  elemental subroutine purgeHtr1DIArr(this)
    type(typeHtr1DIArr),intent(inout)::this !< this 1D integer allocatable array
    
    if(allocated(this%dat)) deallocate(this%dat)
  end subroutine
  
  !> reallocate the 2D integer allocatable array arr to be m by n
  pure subroutine realloc2DIArr(arr,m,n)
    integer,allocatable,intent(inout)::arr(:,:) !< 2D integer array to be reallocated
    integer,intent(in)::m !< number of rows
    integer,intent(in)::n !< number of columns
    
    if(allocated(arr))then
      if(size(arr,1)/=m.or.size(arr,2)/=n)then
        deallocate(arr)
        allocate(arr(m,n))
      end if
    else
      allocate(arr(m,n))
    end if
  end subroutine
  
  !> reallocate the 1D double allocatable array arr to have m entries
  pure subroutine realloc1DDArr(arr,m)
    double precision,allocatable,intent(inout)::arr(:) !< 1D double array to be reallocated
    integer,intent(in)::m !< number of entries
    
    if(allocated(arr))then
      if(size(arr)/=m)then
        deallocate(arr)
        allocate(arr(m))
      end if
    else
      allocate(arr(m))
    end if
  end subroutine
  
  !> reallocate the 2D double allocatable array arr to be m by n
  pure subroutine realloc2DDArr(arr,m,n)
    double precision,allocatable,intent(inout)::arr(:,:) !< 2D double array to be reallocated
    integer,intent(in)::m !< number of rows
    integer,intent(in)::n !< number of columns
    
    if(allocated(arr))then
      if(size(arr,1)/=m.or.size(arr,2)/=n)then
        deallocate(arr)
        allocate(arr(m,n))
      end if
    else
      allocate(arr(m,n))
    end if
  end subroutine
  
  !> reallocate the 1D allocatable array arr of the heterogeneous 1D integer allocatable array
  !> to have m heterogeneous elements
  pure subroutine realloc1DHtr1DIArr(arr,m)
    type(typeHtr1DIArr),allocatable,intent(inout)::arr(:) !< 1D array to be reallocated
    integer,intent(in)::m !< number of entries
    
    if(allocated(arr))then
      if(size(arr)/=m)then
        deallocate(arr)
        allocate(arr(m))
      end if
    else
      allocate(arr(m))
    end if
  end subroutine
  
  !> extend the integer allocatable array arr by l or by 1
  pure subroutine extendIArr(arr,l)
    integer,allocatable,intent(inout)::arr(:) !< integer array to be extended
    integer,intent(in),optional::l !< length to be extended
    integer,allocatable::temp(:)
    
    if(present(l))then
      k=l
    else
      k=1
    end if
    if(allocated(arr))then
      n=size(arr)
      allocate(temp(n+k))
      temp(1:n)=arr(:)
      call move_alloc(temp,arr)
    else
      allocate(arr(k))
    end if
  end subroutine
  
  !> push scalar v to the back of the integer allocatable array
  pure subroutine pushIScalArr(arr,v)
    integer,allocatable,intent(inout)::arr(:) !< integer array to be extended
    integer,intent(in)::v !< integer scalar to be pushed
    
    call extendArr(arr)
    arr(size(arr))=v
  end subroutine
  
  !> push vector v to the back of the integer allocatable array
  pure subroutine pushIVectArr(arr,v)
    integer,allocatable,intent(inout)::arr(:) !< integer array to be extended
    integer,intent(in)::v(:) !< integer vector to be pushed
    
    call extendArr(arr,size(v))
    arr(size(arr)-size(v)+1:size(arr))=v(:)
  end subroutine
  
end module
