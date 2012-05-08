!----------------------------------------------------------------------------- best with 100 columns

!> Basic Data Structures
module moduleBasicDataStruct
  private
  
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
