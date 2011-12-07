!----------------------------------------------------------------------------- best with 100 columns

!******************************
! miscellaneous data structure
!******************************
module moduleMiscDataStruct
  private
  
  !-------------------
  ! pointers to array
  !-------------------
  type::typePtrScalArray
    double precision,pointer::ptr(:)
  end type
  public::typePtrScalArray
  type::typePtrVectArray
    double precision,pointer::ptr(:,:)
  end type
  public::typePtrVectArray
  type::typePtrTensArray
    double precision,pointer::ptr(:,:,:)
  end type
  public::typePtrTensArray
  
  !------------------------------------
  ! extend a generic allocatable array
  !------------------------------------
  interface extendArray
    module procedure::extendIntegerArray
    module procedure::extendDoubleArray
  end interface
  public::extendArray
  
contains
  
  !-----------------------------------------------
  ! extend the integer allocatable array arr by l
  !-----------------------------------------------
  pure subroutine extendIntegerArray(arr,l)
    integer,allocatable,intent(inout)::arr(:)
    integer,intent(in)::l
    integer,allocatable::temp(:)
    
    if(allocated(arr))then
      n=size(arr)
      allocate(temp(n+l))
      temp(1:n)=arr(:)
      call move_alloc(temp,arr)
    else
      allocate(arr(l))
    end if
  end subroutine
  
  !--------------------------------------------------------
  ! extend the double precision allocatable array arr by l
  !--------------------------------------------------------
  pure subroutine extendDoubleArray(arr,l)
    double precision,allocatable,intent(inout)::arr(:)
    integer,intent(in)::l
    double precision,allocatable::temp(:)
    
    if(allocated(arr))then
      n=size(arr)
      allocate(temp(n+l))
      temp(1:n)=arr(:)
      call move_alloc(temp,arr)
    else
      allocate(arr(l))
    end if
  end subroutine
  
end module
