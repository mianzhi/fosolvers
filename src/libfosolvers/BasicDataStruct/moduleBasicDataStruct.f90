!----------------------------------------------------------------------------- best with 100 columns

!> Basic Data Structures
module moduleBasicDataStruct
  private
  
  !> list of integer scalars
  type,public::typeListIScal
    integer,allocatable::dat(:) !< data pool
    integer::len !< length of the data pool that is used
  contains
    procedure,public::clear=>clearListIScal
    procedure,public::extend=>extendListIScal
    generic,public::push=>pushListIScalOne,pushListIScalMulti
      procedure::pushListIScalOne
      procedure::pushListIScalMulti
    generic,public::get=>getListIScalOne,getListIScalMulti
      procedure::getListIScalOne
      procedure::getListIScalMulti
    procedure,public::addto=>addListIScaltoPack
    procedure,public::recover=>recoverListIScalfromPack
    !FIXME:final::purgeListIScal
  end type
  
  !> list of double scalars
  type,public::typeListDScal
    double precision,allocatable::dat(:) !< data pool
    integer::len !< length of the data pool that is used
  contains
    procedure,public::clear=>clearListDScal
    procedure,public::extend=>extendListDScal
    generic,public::push=>pushListDScalOne,pushListDScalMulti
      procedure::pushListDScalOne
      procedure::pushListDScalMulti
    generic,public::get=>getListDScalOne,getListDScalMulti
      procedure::getListDScalOne
      procedure::getListDScalMulti
    procedure,public::addto=>addListDScaltoPack
    procedure,public::recover=>recoverListDScalfromPack
    !FIXME:final::purgeListDScal
  end type
  
  !> list of double vectors
  type,public::typeListDVect
    double precision,allocatable::dat(:,:) !< data pool
    integer::dim !< dimension of the vectors
    integer::len !< length of the data pool that is used
  contains
    procedure,public::clear=>clearListDVect
    procedure,public::specDim=>specListDVectDim
    procedure,public::extend=>extendListDVect
    generic,public::push=>pushListDVectOne,pushListDVectMulti
      procedure::pushListDVectOne
      procedure::pushListDVectMulti
    generic,public::get=>getListDVectOne,getListDVectMulti
      procedure::getListDVectOne
      procedure::getListDVectMulti
    procedure,public::addto=>addListDVecttoPack
    procedure,public::recover=>recoverListDVectfromPack
    !FIXME:final::purgeListDVect
  end type
  
  !> serialized data package
  type,public::typeSerialPack
    type(typeListIScal)::iDat !< integer data pool
    type(typeListDScal)::dDat !< double data pool
    integer iPtr !< pointer for integer data pool
    integer dPtr !< pointer for double data pool
  contains
    procedure,public::clear=>clearSerialPack
    procedure,public::resetPtr=>resetSerialPackPtr
  end type
  
contains
  
  !> clear this ListIScal
  elemental subroutine clearListIScal(this)
    class(typeListIScal),intent(inout)::this
    
    this%len=0
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> extend this ListIScal by doubling size of dat or add n spaces
  elemental subroutine extendListIScal(this,n)
    class(typeListIScal),intent(inout)::this !< this ListIScal
    integer,intent(in),optional::n !< number of entries expected to be extended
    integer,allocatable::newDat(:)
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%len=0
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
  
  !> push one val to the back of this ListIScal
  subroutine pushListIScalOne(this,val)
    class(typeListIScal),intent(inout)::this !< this ListIScal
    integer,intent(in)::val !< scalar val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%len=0
    end if
    do while(this%len+1>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%len+1)=val
    this%len=this%len+1
  end subroutine
  
  !> push multiple val to the back of this ListIScal
  subroutine pushListIScalMulti(this,val)
    class(typeListIScal),intent(inout)::this !< this ListIScal
    integer,intent(in)::val(:) !< vector val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%len=0
    end if
    n=size(val)
    do while(this%len+n>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%len+1:this%len+n)=val(:)
    this%len=this%len+n
  end subroutine
  
  !> get n_th val from this ListIScal
  function getListIScalOne(this,n) result(val)
    class(typeListIScal),intent(in)::this !< this ListIScal
    integer,intent(in)::n !< scalar index to be looked up
    integer val
    
    val=this%dat(n)
  end function
  
  !> get multiple val from this ListIScal
  function getListIScalMulti(this,list) result(val)
    class(typeListIScal),intent(in)::this !< this ListIScal
    integer,intent(in)::list(:) !< vector index to be looked up
    integer val(size(list))
    
    val(:)=this%dat(list(:))
  end function
  
  !> add this ListIScal to serialized data package sp
  subroutine addListIScaltoPack(this,sp)
    class(typeListIScal),intent(in)::this !< this ListIScal
    type(typeSerialPack),intent(inout)::sp !< target package
    
    if(allocated(this%dat))then
      call sp%iDat%push(this%len)
      call sp%iDat%push(this%dat(1:this%len))
    else
      call sp%iDat%push(0)
    end if
  end subroutine
  
  !> recover this ListIScal from serialized data package sp
  subroutine recoverListIScalfromPack(this,sp)
    class(typeListIScal),intent(inout)::this !< this ListIScal
    type(typeSerialPack),intent(inout)::sp !< source package
    
    call this%clear()
    sp%iPtr=sp%iPtr+1
    m=sp%iDat%get(sp%iPtr)
    this%len=m
    if(m>0)then
      call this%extend(m)
      call this%push(sp%iDat%get([(sp%iPtr+i,i=1,m)]))
      sp%iPtr=sp%iPtr+m
    end if
  end subroutine
  
  !> destructor of ListIScal
  elemental subroutine purgeListIScal(this)
    type(typeListIScal),intent(inout)::this !< this ListIScal
    
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> clear this ListDScal
  elemental subroutine clearListDScal(this)
    class(typeListDScal),intent(inout)::this !< this ListDScal
    
    this%len=0
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> extend this ListDScal by doubling size of dat or add n spaces
  elemental subroutine extendListDScal(this,n)
    class(typeListDScal),intent(inout)::this !< this ListDScal
    integer,intent(in),optional::n !< number of entries expected to be extended
    double precision,allocatable::newDat(:)
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%len=0
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
  
  !> push one val to the back of this ListDScal
  subroutine pushListDScalOne(this,val)
    class(typeListDScal),intent(inout)::this !< this ListDScal
    double precision,intent(in)::val !< scalar val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%len=0
    end if
    do while(this%len+1>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%len+1)=val
    this%len=this%len+1
  end subroutine
  
  !> push multiple val to the back of this ListDScal
  subroutine pushListDScalMulti(this,val)
    class(typeListDScal),intent(inout)::this !< this ListDScal
    double precision,intent(in)::val(:) !< vector val to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat)
    else
      m=0
      this%len=0
    end if
    n=size(val)
    do while(this%len+n>m)
      call this%extend()
      m=size(this%dat)
    end do
    this%dat(this%len+1:this%len+n)=val(:)
    this%len=this%len+n
  end subroutine
  
  !> get n_th val from this ListDScal
  function getListDScalOne(this,n) result(val)
    class(typeListDScal),intent(in)::this !< this ListDScal
    integer,intent(in)::n !< scalar index to be looked up
    double precision val
    
    val=this%dat(n)
  end function
  
  !> get multiple val from this ListDScal
  function getListDScalMulti(this,list) result(val)
    class(typeListDScal),intent(in)::this !< this ListDScal
    integer,intent(in)::list(:) !< vector index to be looked up
    double precision val(size(list))
    
    val(:)=this%dat(list(:))
  end function
  
  !> add this ListDScal to serialized data package sp
  subroutine addListDScaltoPack(this,sp)
    class(typeListDScal),intent(in)::this !< this ListDScal
    type(typeSerialPack),intent(inout)::sp !< target package
    
    if(allocated(this%dat))then
      call sp%iDat%push(this%len)
      call sp%dDat%push(this%dat(1:this%len))
    else
      call sp%iDat%push(0)
    end if
  end subroutine
  
  !> recover this ListDScal from serialized data package sp
  subroutine recoverListDScalfromPack(this,sp)
    class(typeListDScal),intent(inout)::this !< this ListDScal
    type(typeSerialPack),intent(inout)::sp !< source package
    
    call this%clear()
    sp%iPtr=sp%iPtr+1
    m=sp%iDat%get(sp%iPtr)
    this%len=m
    if(m>0)then
      call this%extend(m)
      call this%push(sp%dDat%get([(sp%dPtr+i,i=1,m)]))
      sp%dPtr=sp%dPtr+m
    end if
  end subroutine
  
  !> destructor of ListDScal
  elemental subroutine purgeListDScal(this)
    type(typeListDScal),intent(inout)::this !< this ListDScal
    
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> clear this ListDVect
  elemental subroutine clearListDVect(this)
    class(typeListDVect),intent(inout)::this !< this ListDVect
    
    this%dim=0
    this%len=0
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> specify the vectors in this ListDVect to be n dimensional
  elemental subroutine specListDVectDim(this,n)
    class(typeListDVect),intent(inout)::this !< this ListDVect
    integer,intent(in)::n !< vector dimension
    
    this%dim=n
  end subroutine
  
  !> extend this ListDVect by doubling size of dat or add n spaces
  elemental subroutine extendListDVect(this,n)
    class(typeListDVect),intent(inout)::this !< this ListDVect
    integer,intent(in),optional::n !< number of entries expected to be extended
    double precision,allocatable::newDat(:,:)
    
    if(allocated(this%dat))then
      m=size(this%dat,2)
    else
      m=0
      this%len=0
    end if
    if(present(n))then
      allocate(newDat(this%dim,m+n))
    else
      allocate(newDat(this%dim,max(1,(m*2))))
    end if
    if(m>0)then
      newDat(:,1:m)=this%dat(:,:)
    end if
    call move_alloc(newDat,this%dat)
  end subroutine
  
  !> push one vector to the back of this ListDVect
  subroutine pushListDVectOne(this,vector)
    class(typeListDVect),intent(inout)::this !< this ListDVect
    double precision,intent(in)::vector(this%dim) !< vector to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat,2)
    else
      m=0
      this%len=0
    end if
    do while(this%len+1>m)
      call this%extend()
      m=size(this%dat,2)
    end do
    this%dat(:,this%len+1)=vector(:)
    this%len=this%len+1
  end subroutine
  
  !> push multiple vectors to the back of this ListDVect
  subroutine pushListDVectMulti(this,vectors)
    class(typeListDVect),intent(inout)::this !< this ListDVect
    double precision,intent(in)::vectors(:,:) !< vectors to be pushed in
    
    if(allocated(this%dat))then
      m=size(this%dat,2)
    else
      m=0
      this%len=0
    end if
    n=size(vectors,2)
    do while(this%len+n>m)
      call this%extend()
      m=size(this%dat,2)
    end do
    this%dat(:,this%len+1:this%len+n)=vectors(:,:)
    this%len=this%len+n
  end subroutine
  
  !> get n_th vector from this ListDVect
  function getListDVectOne(this,n) result(vector)
    class(typeListDVect),intent(in)::this !< this ListDVect
    integer,intent(in)::n !< vector index to be looked up
    double precision vector(this%dim)
    
    vector(:)=this%dat(:,n)
  end function
  
  !> get multiple vectors from this ListDVect
  function getListDVectMulti(this,list) result(vectors)
    class(typeListDVect),intent(in)::this !< this ListDVect
    integer,intent(in)::list(:) !< vector index to be looked up
    double precision vectors(this%dim,size(list))
    
    vectors(:,:)=this%dat(:,list(:))
  end function
  
  !> add this ListDVect to serialized data package sp
  subroutine addListDVecttoPack(this,sp)
    class(typeListDVect),intent(in)::this !< this ListDVect
    type(typeSerialPack),intent(inout)::sp !< target package
    
    call sp%iDat%push(this%dim)
    if(allocated(this%dat))then
      call sp%iDat%push(this%len)
      do i=1,this%len
        call sp%dDat%push(this%dat(:,i))
      end do
    else
      call sp%iDat%push(0)
    end if
  end subroutine
  
  !> recover this ListDVect from serialized data package sp
  subroutine recoverListDVectfromPack(this,sp)
    class(typeListDVect),intent(inout)::this !< this ListDVect
    type(typeSerialPack),intent(inout)::sp !< source package
    
    call this%clear()
    sp%iPtr=sp%iPtr+1
    k=sp%iDat%get(sp%iPtr)
    this%dim=k
    sp%iPtr=sp%iPtr+1
    m=sp%iDat%get(sp%iPtr)
    this%len=m
    if(m>0)then
      call this%extend(m)
      do i=1,m
        call this%push(sp%dDat%get([(sp%dPtr+j,j=1,k)]))
        sp%dPtr=sp%dPtr+k
      end do
    end if
  end subroutine
  
  !> destructor of ListDVect
  elemental subroutine purgeListDVect(this)
    type(typeListDVect),intent(inout)::this !< this ListDVect
    
    if(allocated(this%dat))then
      deallocate(this%dat)
    end if
  end subroutine
  
  !> clear this SerialPack
  elemental subroutine clearSerialPack(this)
    class(typeSerialPack),intent(inout)::this !< this SerialPack
    
    call this%iDat%clear
    call this%dDat%clear
    call this%resetPtr
  end subroutine
  
  !> reset the pointers of this SerialPack
  elemental subroutine resetSerialPackPtr(this)
    class(typeSerialPack),intent(inout)::this !< this SerialPack
    
    this%iPtr=0
    this%dPtr=0
  end subroutine
  
end module
