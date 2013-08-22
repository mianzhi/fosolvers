!----------------------------------------------------------------------------- best with 100 columns

!> spatial hashing
module moduleSpatialHashing
  use moduleBasicDataStruct
  private
  
  !> constants
  integer,parameter,public::DIMS=3 !< dimensions
  
  !> spatial hash table
  type,public::typeSHT
    double precision SideLength !< side length of the box
    double precision BoundBox(DIMS,2) !< bound box coordinates
    type(typeHtr1DIArr),allocatable::list(:,:,:) !< the position list
  contains
    procedure,public::hash=>hashSHT
    procedure,public::lookup=>lookupSHT
    procedure,public::fill=>fillSHT
    procedure,public::findNeib=>findNeibSHT
    !FIXME:final::purgeSHT
  end type
  
  !> public procedures
  public::hashSHT
  
contains
  
  !> find the hashing result, a 3-integer-combination
  pure function hashSHT(this,pos)
    class(typeSHT),intent(in)::this !< this SHT
    double precision,intent(in)::pos(DIMS) !< the position to be hashed
    integer hashSHT(DIMS) !< the index of the box
    
    hashSHT(:)=ceiling((pos(:)-this%BoundBox(:,1))/this%SideLength)
  end function
  
  !> lookup the list of points with index
  pure function lookupSHT(this,ind)
    class(typeSHT),intent(in)::this !< this SHT
    integer,intent(in)::ind(DIMS) !< the 3-integer-combination
    integer,allocatable::lookupSHT(:) !< the list of points with box having given index
    
    if(allocated(this%list(ind(1),ind(2),ind(3))%dat))then
      call reallocArr(lookupSHT,size(this%list(ind(1),ind(2),ind(3))%dat))
      lookupSHT(:)=this%list(ind(1),ind(2),ind(3))%dat(:)
    end if
  end function
  
  !> fill the SHT with a list of coordinates, use roughly m boxes
  pure subroutine fillSHT(this,pos,m)
    class(typeSHT),intent(inout)::this !< this SHT
    double precision,intent(in)::pos(:,:) !< the list of coordinates
    integer,intent(in)::m !< rough number of boxes
    double precision lwh(DIMS)
    integer ind(DIMS)
    
    this%BoundBox(:,1)=minval(pos,2)
    this%BoundBox(:,2)=maxval(pos,2)
    lwh(:)=this%BoundBox(:,2)-this%BoundBox(:,1)
    this%SideLength=(product(lwh)/dble(m))**(1d0/3d0)
    if(allocated(this%list)) deallocate(this%list)
    allocate(this%list(0:ceiling(lwh(1)/this%SideLength)+1,&
    &                  0:ceiling(lwh(2)/this%SideLength)+1,&
    &                  0:ceiling(lwh(3)/this%SideLength)+1))
    do i=1,size(pos,2)
      ind=this%hash(pos(:,i))
      call pushArr(this%list(ind(1),ind(2),ind(3))%dat,i)
    end do
  end subroutine
  
  !> find a list of neighbor points using SHT
  pure function findNeibSHT(this,pos)
    class(typeSHT),intent(in)::this !< this SHT
    double precision,intent(in)::pos(DIMS) !< the reference position
    integer,allocatable::findNeibSHT(:) !< the list of neighbor points
    integer ind(DIMS)
    
    ind=this%hash(pos(:))
    findNeibSHT=this%lookup(ind)
  end function
  
  !> destructor of SHT
  elemental subroutine purgeSHT(this)
    type(typeSHT),intent(inout)::this !< this SHT
    
    if(allocated(this%list)) deallocate(this%list)
  end subroutine
  
end module
