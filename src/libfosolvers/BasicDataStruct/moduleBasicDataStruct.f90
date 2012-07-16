!----------------------------------------------------------------------------- best with 100 columns

!> Basic Data Structures
module moduleBasicDataStruct
  private
  
  ! constants
  integer,parameter,public::CONST_TYPE=0 !< constant type
  integer,parameter,public::TAB1D_TYPE=1 !< 1d table type
  
  !> heterogeneous 1D integer allocatable array
  type,public::typeHtr1DIArr
    integer,allocatable::dat(:) !< 1D integer allocatable array
  contains
    !FIXME:final::purgeHtr1DIArr
  end type
  
  !> heterogeneous 2D double allocatable array
  type,public::typeHtr2DDArr
    double precision,allocatable::dat(:,:) !< 2D double allocatable array
  contains
    !FIXME:final::purgeHtr2DDArr
  end type
  
  !> general data linked list
  type,public::typeGenDatLlist
    !TODO: change to character(:) (waiting new gcc)
    character(1),allocatable::key(:) !< key of the data item
    type(typeGenDatLlist),pointer::next=>null() !< the next data item in the list
    integer datType !< type of data
    double precision,allocatable::dat(:,:) !< raw data
  contains
    procedure,public::push=>pushGenDatLlist
    procedure,public::test=>testGenDatLlist
    procedure,public::get=>getGenDatLlist
    procedure,public::clear=>clearGenDatLlist
    !FIXME:final::purgeGenDatLlist
  end type
  
  !> reallocate a generic allocatable array
  interface reallocArr
    module procedure::realloc1DIArr
    module procedure::realloc2DIArr
    module procedure::realloc1DDArr
    module procedure::realloc2DDArr
    module procedure::reallocStr
    module procedure::realloc1DHtr1DIArr
    module procedure::realloc1DHtr2DDArr
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
  
  !> destructor of typeHtr2DDArr
  elemental subroutine purgeHtr2DDArr(this)
    type(typeHtr2DDArr),intent(inout)::this !< this 2D double allocatable array
    
    if(allocated(this%dat)) deallocate(this%dat)
  end subroutine
  
  !> push item to general data linked list
  recursive pure subroutine pushGenDatLlist(this,item)
    class(typeGenDatLlist),intent(inout)::this !< this linked list item
    type(typeGenDatLlist),intent(in)::item !< item to be pushed
    !TODO: change to character(:) (waiting new gcc)
    character(:),allocatable::thiskey,itemkey
    if(allocated(this%key))then
      allocate(character(size(this%key))::thiskey)
      forall(i=1:size(this%key))
        thiskey(i:i)=this%key(i)
      end forall
    else
      thiskey=''
    end if
    allocate(character(size(item%key))::itemkey)
    forall(i=1:size(item%key))
      itemkey(i:i)=item%key(i)
    end forall
    
    if(allocated(this%key))then
      if(thiskey/=itemkey)then
        if(.not.associated(this%next))then
          allocate(this%next)
        end if
        call this%next%push(item)
      end if
    else
      allocate(this%key(size(item%key)),source=item%key)
      this%datType=item%datType
      !TODO:remove the allocation bounds (waiting new gcc)
      allocate(this%dat(lbound(item%dat,1):ubound(item%dat,1),lbound(item%dat,2):ubound(item%dat,2)),source=item%dat)
    end if
  end subroutine
  
  !> test if the key exists in the linked list
  recursive pure function testGenDatLlist(this,key)
    class(typeGenDatLlist),intent(in)::this !< this linked list item
    character(*),intent(in)::key !< the key to be tested
    logical testGenDatLlist !< whether the key exists in this linked list
    !TODO: change to character(:) (waiting new gcc)
    character(:),allocatable::thiskey
    allocate(character(size(this%key))::thiskey)
    forall(i=1:size(this%key))
      thiskey(i:i)=this%key(i)
    end forall
    
    testGenDatLlist=.false.
    if((key==thiskey))then
      testGenDatLlist=.true.
    else
      if(associated(this%next))then
        testGenDatLlist=this%next%test(key)
      end if
    end if
  end function
  
  !> get the value according to the key from this linked list
  recursive pure function getGenDatLlist(this,key,arg)
    class(typeGenDatLlist),intent(in)::this !< this linked list item
    character(*),intent(in)::key !< the key to be found
    double precision,intent(in),optional::arg(*) !< optional arguments
    double precision getGenDatLlist !< the value associated with the key
    double precision x
    !TODO: change to character(:) (waiting new gcc)
    character(:),allocatable::thiskey
    allocate(character(size(this%key))::thiskey)
    forall(i=1:size(this%key))
      thiskey(i:i)=this%key(i)
    end forall
    
    getGenDatLlist=0d0
    if(key==thiskey)then
      select case(this%datType)
      case(CONST_TYPE)
        getGenDatLlist=this%dat(1,1)
      case(TAB1D_TYPE)
        if(present(arg))then
          x=arg(1)
        else
          x=0d0
        end if
        do i=2,size(this%dat,1)
          if(x<this%dat(i,0))then
            getGenDatLlist=this%dat(i-1,1)+(x-this%dat(i-1,0))&
            &                              *(this%dat(i,1)-this%dat(i-1,1))&
            &                              /(this%dat(i,0)-this%dat(i-1,0))
            exit
          end if
        end do
        if(i>ubound(this%dat,1))then
          getGenDatLlist=this%dat(i-2,1)+(x-this%dat(i-2,0))&
          &                              *(this%dat(i-1,1)-this%dat(i-2,1))&
          &                              /(this%dat(i-1,0)-this%dat(i-2,0))
        end if
      case default
      end select
    else
      if(associated(this%next))then
        if(present(arg))then
          getGenDatLlist=this%next%get(key,arg)
        else
          getGenDatLlist=this%next%get(key)
        end if
      end if
    end if
  end function
  
  !> clear this linked list
  recursive pure subroutine clearGenDatLlist(this)
    class(typeGenDatLlist),intent(inout)::this !< this linked list item
    
    if(associated(this%next))then
      call this%next%clear()
      nullify(this%next)
    end if
    if(allocated(this%key)) deallocate(this%key)
    if(allocated(this%dat)) deallocate(this%dat)
  end subroutine
  
  !> destructor of GenDatLlist
  elemental subroutine purgeGenDatLlist(this)
    type(typeGenDatLlist),intent(inout)::this !< this linked list item
    
    if(allocated(this%key)) deallocate(this%key)
    if(allocated(this%dat)) deallocate(this%dat)
    if(associated(this%next)) nullify(this%next)
  end subroutine
  
  !> reallocate the 1D integer array arr to have m entries
  pure subroutine realloc1DIArr(arr,m)
    integer,allocatable,intent(inout)::arr(:) !< 1D integer array to be reallocated
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
  
  !> reallocate the string str to have m characters
  pure subroutine reallocStr(str,m)
    character(:),allocatable,intent(inout)::str !< string to be reallocated
    integer,intent(in)::m !< number of characters
    
    if(allocated(str))then
      if(len(str)/=m)then
        deallocate(str)
        allocate(character(m)::str)
      end if
    else
      allocate(character(m)::str)
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
  
  !> reallocate the 1D allocatable array arr of the heterogeneous 2D double allocatable array
  !> to have m heterogeneous elements
  pure subroutine realloc1DHtr2DDArr(arr,m)
    type(typeHtr2DDArr),allocatable,intent(inout)::arr(:) !< 1D array to be reallocated
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
