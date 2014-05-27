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
  contains
    procedure,public::init=>initPolyX
    procedure,public::clear=>clearPolyX
    procedure,public::p=>pPolyX
    !FIXME:final::purgePolyX
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
  end subroutine
  
  !> clear this polyX
  elemental subroutine clearPolyX(this)
    class(polyX),intent(inout)::this !< this polyX
    
    if(allocated(this%pN)) deallocate(this%pN)
    if(allocated(this%sE)) deallocate(this%sE)
    if(allocated(this%nNE)) deallocate(this%nNE)
    if(allocated(this%iNE)) deallocate(this%iNE)
  end subroutine
  
  !> center position of element k
  pure function pPolyX(this,k)
    class(polyX),intent(in)::this !< this polyX
    integer,intent(in)::k !< element index
    double precision::pPolyX(DIMS)
    
    pPolyX(:)=sum(this%pN(:,this%iNE(:,k)),2)/dble(this%nNE(k))
  end function
  
  !> destructor of polyX
  elemental subroutine purgePolyX(this)
    type(polyX),intent(inout)::this !< this polyX
    
    call this%clear()
  end subroutine
  
end module
