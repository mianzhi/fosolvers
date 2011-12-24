!----------------------------------------------------------------------------- best with 100 columns

!******************************
! miscellaneous data structure
!******************************
module moduleMiscDataStruct
  private
  
  ! constants
  integer,parameter,public::DATA_NAME_LENGTH=5
  integer,parameter,public::VAL_TYPE=1
  integer,parameter,public::TAB1D_TYPE=2
  
  !------------------------------------
  ! extend a generic allocatable array
  !------------------------------------
  interface extendArray
    module procedure::extendIntegerArray
    module procedure::extendDoubleArray
  end interface
  public::extendArray
  
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
  
  type,public::typeDataSet
    type(typeDataItem),pointer::DataItem(:)=>null()
    type(typeDataItem),pointer::ptrLast
  contains
    procedure,public::extend=>extendDataSet
    procedure,public::push=>pushDataSet
    procedure,public::lookup=>lookupDataSet
    !TODO: final::cleanDataSet
  end type
  
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
  
  !--------------------------------------
  ! extend l data items in this data set
  !--------------------------------------
  elemental subroutine extendDataSet(this,l)
    class(typeDataSet),intent(inout)::this
    integer,intent(in)::l
    type(typeDataItem),pointer::temp(:)
    
    if(associated(this%DataItem))then
      n=size(this%DataItem)
      allocate(temp(n+l))
      temp(1:n)=this%DataItem(:)
      deallocate(this%DataItem)
      this%DataItem=>temp
      nullify(temp)
    else
      allocate(this%DataItem(l))
    end if
    this%ptrLast=>this%DataItem(size(this%DataItem))
  end subroutine
  
  !--------------------------------
  ! push dataitem to this data set
  !--------------------------------
  elemental subroutine pushDataSet(this,dataitem)
    class(typeDataSet),intent(inout)::this
    type(typeDataItem),intent(in),value::dataitem
    
    call this%extend(1)
    this%ptrLast%DataName=dataitem%DataName
    this%ptrLast%DataType=dataitem%DataType
    this%ptrLast%Val=dataitem%Val
    if(allocated(dataitem%Tab1d))then
      allocate(this%ptrLast%Tab1d(size(dataitem%Tab1d,1),2))
      this%ptrLast%Tab1d(:,:)=dataitem%Tab1d(:,:)
    end if
  end subroutine
  
  !------------------------------------------------
  ! lookup value named dataname from this data set
  !------------------------------------------------
  function lookupDataSet(this,dataname,x,stat)
    class(typeDataSet),intent(in)::this
    character(DATA_NAME_LENGTH),intent(in)::dataname
    double precision,optional,intent(in)::x
    integer,optional,intent(out)::stat
    double precision lookupDataSet
    
    lookupDataSet=0d0
    if(present(stat))then
      stat=0
    end if
    if(associated(this%DataItem))then
      do i=1,size(this%DataItem)
        if(this%DataItem(i)%DataName==dataname)then
          if(present(x))then
            lookupDataSet=this%DataItem(i)%lookup(x)
          else
            lookupDataSet=this%DataItem(i)%lookup()
          end if
          exit
        end if
      end do
      if(present(stat).and.i>size(this%DataItem))then
        stat=stat+1 ! data not found
      end if
    else
      if(present(stat))then
        stat=stat+2 ! data not associated
      end if
    end if
  end function
  
  !---------------------------
  ! destructor of typeDataSet
  !---------------------------
  elemental subroutine cleanDataSet(this)
    type(typeDataSet),intent(inout)::this
    
    if(associated(this%DataItem))then
      deallocate(this%DataItem)
    end if
  end subroutine
  
end module
