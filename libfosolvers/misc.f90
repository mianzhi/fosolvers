!----------------------------------------------------------------------------- best with 100 columns

!******************************
! miscellaneous data structure
!******************************
module moduleMiscDataStruct
  private
  
  ! constants
  integer,parameter,private::DATA_NAME_LENGTH=3
  integer,parameter,public::VAL_TYPE=1
  integer,parameter,public::TAB1D_TYPE=2
  
  !-------------------
  ! pointers to array
  !-------------------
  type,public::typePtrScalArray
    double precision,pointer::ptr(:)
  end type
  type,public::typePtrVectArray
    double precision,pointer::ptr(:,:)
  end type
  type,public::typePtrTensArray
    double precision,pointer::ptr(:,:,:)
  end type
  
  !----------
  ! data set
  !----------
  type,public::typeDataItem
    character(DATA_NAME_LENGTH) DataName
    integer DataType
    double precision Val
    double precision,allocatable::Tab1d(:,:)
  contains
    procedure,public::specify=>specifyDataItem
    procedure,public::lookup=>lookupDataItem
    !TODO: final::cleanDataItem
  end type
  
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
  
  !------------------------------------
  ! specify at type for this data item
  !------------------------------------
  elemental subroutine specifyDataItem(this,datatype,n)
    class(typeDataItem),intent(inout)::this
    integer,intent(in)::datatype
    integer,intent(in),optional::n
    integer tablength
    
    this%DataType=datatype
    if(allocated(this%Tab1d))then
      deallocate(this%Tab1d)
    end if
    select case(datatype)
      case(VAL_TYPE)
        ! nothing to do
      case(TAB1D_TYPE)
        if(present(n))then
          tablength=n
        else
          tablength=0
        end if
        allocate(this%Tab1d(tablength,2))
      case default
    end select
  end subroutine
  
  !----------------------------------
  ! lookup value from this data item
  !----------------------------------
  elemental function lookupDataItem(this,x)
    class(typeDataItem),intent(in)::this
    double precision,intent(in),optional::x
    double precision lookupDataItem
    double precision pos
    
    lookupDataItem=0d0
    select case(this%DataType)
      case(VAL_TYPE)
        lookupDataItem=this%Val
      case(TAB1D_TYPE)
        if(present(x))then
          pos=x
        else
          pos=0d0
        end if
        do i=2,size(this%Tab1d,1)
          if(pos<this%Tab1d(i,1))then
            lookupDataItem=this%Tab1d(i-1,2)+(pos-this%Tab1d(i-1,1))&
            &                                *(this%Tab1d(i,2)-this%Tab1d(i-1,2))&
            &                                /(this%Tab1d(i,1)-this%Tab1d(i-1,1))
            exit
          end if
        end do
        if(i>size(this%Tab1d,1))then
          lookupDataItem=this%Tab1d(i-2,2)+(pos-this%Tab1d(i-2,1))&
          &                                *(this%Tab1d(i-1,2)-this%Tab1d(i-2,2))&
          &                                /(this%Tab1d(i-1,1)-this%Tab1d(i-2,1))
        end if
      case default
    end select
  end function
  
  !----------------------------
  ! destructor of typeDataItem
  !----------------------------
  elemental subroutine cleanDataItem(this)
    type(typeDataItem),intent(inout)::this
    
    if(allocated(this%Tab1d))then
      deallocate(this%Tab1d)
    end if
  end subroutine
  
end module
