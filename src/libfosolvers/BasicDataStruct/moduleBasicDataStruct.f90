!----------------------------------------------------------------------------- best with 100 columns

!> Basic Data Structures
module moduleBasicDataStruct
  private
  
  !> list of integer scalars
  type,public::typeListIntegerScal
    integer,allocatable::dat(:) !< data pool
    integer::length !< length of the data pool that is used
  contains
    procedure,public::clear=>clearListIntegerScal
    procedure,public::extend=>extendListIntegerScal
    generic,public::push=>pushListIntegerScalOne,pushListIntegerScalMulti
      procedure::pushListIntegerScalOne
      procedure::pushListIntegerScalMulti
    generic,public::get=>getListIntegerScalOne,getListIntegerScalMulti
      procedure::getListIntegerScalOne
      procedure::getListIntegerScalMulti
    !FIXME:final::cleanListIntegerScal
  end type
  
  !> list of double scalars
  type,public::typeListDoubleScal
    double precision,allocatable::dat(:) !< data pool
    integer::length !< length of the data pool that is used
  contains
    procedure,public::clear=>clearListDoubleScal
    procedure,public::extend=>extendListDoubleScal
    generic,public::push=>pushListDoubleScalOne,pushListDoubleScalMulti
      procedure::pushListDoubleScalOne
      procedure::pushListDoubleScalMulti
    generic,public::get=>getListDoubleScalOne,getListDoubleScalMulti
      procedure::getListDoubleScalOne
      procedure::getListDoubleScalMulti
    !FIXME:final::cleanListDoubleScal
  end type
  
contains
  
  !> clear this ListIntegerScal
  elemental subroutine clearListIntegerScal(this)
    class(typeListIntegerScal),intent(inout)::this
    
    this%length=0
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> extend this ListIntegerScal by doubling size of dat or add n spaces
  elemental subroutine extendListIntegerScal(this,n)
    class(typeListIntegerScal),intent(inout)::this !< this ListIntegerScal
    integer,intent(in),optional::n !< number of entries expected to be extended
    integer,allocatable::newDat(:)
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%length=0
    end if
    if(present(n))then
      allocate(newDat(m+n))
    else
      allocate(newDat(max(1,(m*2))))
    end if
    if(m>0)then
      newDat(1:m)=this%dat(:)
    end if
    call move_alloc(newDat,this%dat)
  end subroutine
  
  !> push one val to the back of this ListIntegerScal
  subroutine pushListIntegerScalOne(this,val)
    class(typeListIntegerScal),intent(inout)::this !< this ListIntegerScal
    integer,intent(in)::val !< scalar val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%length=0
    end if
    do while(this%length+1>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%length+1)=val
    this%length=this%length+1
  end subroutine
  
  !> push multiple val to the back of this ListIntegerScal
  subroutine pushListIntegerScalMulti(this,val)
    class(typeListIntegerScal),intent(inout)::this !< this ListIntegerScal
    integer,intent(in)::val(:) !< vector val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%length=0
    end if
    n=size(val)
    do while(this%length+n>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%length+1:this%length+n)=val(:)
    this%length=this%length+n
  end subroutine
  
  !> get n_th val from this ListIntegerScal
  function getListIntegerScalOne(this,n) result(val)
    class(typeListIntegerScal),intent(in)::this !< this ListIntegerScal
    integer,intent(in)::n !< scalar index to be looked up
    integer val
    
    val=this%dat(n)
  end function
  
  !> get multiple val from this ListIntegerScal
  function getListIntegerScalMulti(this,list) result(val)
    class(typeListIntegerScal),intent(in)::this !< this ListIntegerScal
    integer,intent(in)::list(:) !< vector index to be looked up
    integer val(size(list))
    
    val(:)=this%dat(list(:))
  end function
  
  !> destructor of ListIntegerScal
  elemental subroutine cleanListIntegerScal(this)
    type(typeListIntegerScal),intent(inout)::this !< this ListDoubleScal
    
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> clear this ListDoubleScal
  elemental subroutine clearListDoubleScal(this)
    class(typeListDoubleScal),intent(inout)::this !< this ListDoubleScal
    
    this%length=0
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> extend this ListDoubleScal by doubling size of dat or add n spaces
  elemental subroutine extendListDoubleScal(this,n)
    class(typeListDoubleScal),intent(inout)::this !< this ListDoubleScal
    integer,intent(in),optional::n !< number of entries expected to be extended
    double precision,allocatable::newDat(:)
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%length=0
    end if
    if(present(n))then
      allocate(newDat(m+n))
    else
      allocate(newDat(max(1,(m*2))))
    end if
    if(m>0)then
      newDat(1:m)=this%dat(:)
    end if
    call move_alloc(newDat,this%dat)
  end subroutine
  
  !> push one val to the back of this ListDoubleScal
  subroutine pushListDoubleScalOne(this,val)
    class(typeListDoubleScal),intent(inout)::this !< this ListDoubleScal
    double precision,intent(in)::val !< scalar val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%length=0
    end if
    do while(this%length+1>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%length+1)=val
    this%length=this%length+1
  end subroutine
  
  !> push multiple val to the back of this ListDoubleScal
  subroutine pushListDoubleScalMulti(this,val)
    class(typeListDoubleScal),intent(inout)::this !< this ListDoubleScal
    double precision,intent(in)::val(:) !< vector val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%length=0
    end if
    n=size(val)
    do while(this%length+n>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%length+1:this%length+n)=val(:)
    this%length=this%length+n
  end subroutine
  
  !> get n_th val from this ListDoubleScal
  function getListDoubleScalOne(this,n) result(val)
    class(typeListDoubleScal),intent(in)::this !< this ListDoubleScal
    integer,intent(in)::n !< scalar index to be looked up
    double precision val
    
    val=this%dat(n)
  end function
  
  !> get multiple val from this ListDoubleScal
  function getListDoubleScalMulti(this,list) result(val)
    class(typeListDoubleScal),intent(in)::this !< this ListDoubleScal
    integer,intent(in)::list(:) !< vector index to be looked up
    double precision val(size(list))
    
    val(:)=this%dat(list(:))
  end function
  
  !> destructor of ListDoubleScal
  elemental subroutine cleanListDoubleScal(this)
    type(typeListDoubleScal),intent(inout)::this !< this ListDoubleScal
    
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
end module
