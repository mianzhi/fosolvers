!----------------------------------------------------------------------------- best with 100 columns

!> polygons and polyhedrons module
module modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> base type of polygons and polyhedrons
  type,public::polyX
    integer::nN !< number of nodes
    integer::nE !< number of elements (facets or cells)
    double precision,allocatable::pN(:,:) !< node position
    integer,allocatable::sE(:) !< shape of elements
    integer,allocatable::nNE(:) !< number of nodes in each element
    integer,allocatable::iNE(:,:) !< connectivity table
    logical::isUp !< if auxiliary data is updated
    double precision,allocatable::p(:,:) !< element position
  contains
    procedure,public::init=>initPolyX
    procedure,public::clear=>clearPolyX
    procedure,public::up=>upPolyX
    final::purgePolyX
  end type
  
contains
  
  !> initialize this polyX
  elemental subroutine initPolyX(this,nN,nE,m)
    class(polyX),intent(inout)::this !< this polyX
    integer,intent(in)::nN !< number of nodes
    integer,intent(in)::nE !< number of elements
    integer,intent(in)::m !< maximum number of nodes per elements
    
    call this%clear()
    this%nN=nN
    this%nE=nE
    allocate(this%pN(DIMS,nN))
    allocate(this%sE(nE))
    allocate(this%nNE(nE))
    allocate(this%iNE(m,nE))
    allocate(this%p(DIMS,nE))
    this%isUp=.false.
  end subroutine
  
  !> clear this polyX
  elemental subroutine clearPolyX(this)
    class(polyX),intent(inout)::this !< this polyX
    
    if(allocated(this%pN)) deallocate(this%pN)
    if(allocated(this%sE)) deallocate(this%sE)
    if(allocated(this%nNE)) deallocate(this%nNE)
    if(allocated(this%iNE)) deallocate(this%iNE)
    if(allocated(this%p)) deallocate(this%p)
  end subroutine
  
  !> update this polyX
  subroutine upPolyX(this)
    use modGeometry
    class(polyX),intent(inout)::this !< this polyX
    
    if(.not.this%isUp)then
      ! calculate element position
      forall(i=1:this%nE)
        this%p(:,i)=sum(this%pN(:,this%iNE(1:this%nNE(i),i)),2)/dble(this%nNE(i))
      end forall
      this%isUp=.true.
    end if
  end subroutine
  
  !> destructor of polyX
  elemental subroutine purgePolyX(this)
    type(polyX),intent(inout)::this !< this polyX
    
    call this%clear()
  end subroutine
  
end module
