!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron and polygon finite volume grid module
module modPolyFvGrid
  use modPolyGrid
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> polyhedron and polygon finite volume grid type
  type,extends(polyGrid),public::polyFvGrid
    integer::nC !< number of cells
    integer::nP !< number of pairs of elements
    integer,allocatable::iEP(:,:) !< indices of elements of each pair
    integer,allocatable::neib(:,:) !< neighbor list
  contains
    procedure,public::clear=>clearPolyFvGrid
    procedure,public::up=>upPolyFvGrid
    !FIXME:final::purgePolyFvGrid
  end type
  
contains
  
  !> clear this polyFvGrid
  elemental subroutine clearPolyFvGrid(this)
    class(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    call this%polyGrid%clear()
    if(allocated(this%iEP)) deallocate(this%iEP)
    if(allocated(this%neib)) deallocate(this%neib)
  end subroutine
  
  !> update this polyFvGrid
  elemental subroutine upPolyFvGrid(this)
    class(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    if(.not.this%isUp)then
      call this%polyGrid%up()
      call sortPolyFvGrid(this)
      !TODO nP,iEP,neib
    end if
  end subroutine
  
  !> sort elements such that blocks are in the front
  elemental subroutine sortPolyFvGrid(grid)
    class(polyFvGrid),intent(inout)::grid !< the polyFvGrid
    integer,allocatable::sE(:),nNE(:),iNE(:,:),gid(:)
    double precision,allocatable::v(:)
    
    !FIXME: remove the array specification work-around
    allocate(sE(grid%nE),source=grid%sE)
    allocate(nNE(grid%nE),source=grid%nNE)
    allocate(iNE(size(grid%iNE,1),grid%nE),source=grid%iNE)
    allocate(gid(grid%nE),source=grid%gid)
    allocate(v(grid%nE),source=grid%v)
    j=0
    do i=1,grid%nE
      if(sE(i)==TET.or.se(i)==HEX)then
        j=j+1
        grid%sE(j)=sE(i)
        grid%nNE(j)=nNE(i)
        grid%iNE(:,j)=iNE(:,i)
        grid%gid(j)=gid(i)
        grid%v(j)=v(i)
      end if
    end do
    grid%nC=j
    do i=1,grid%nE
      if(sE(i)==TRI.or.se(i)==QUAD)then
        j=j+1
        grid%sE(j)=sE(i)
        grid%nNE(j)=nNE(i)
        grid%iNE(:,j)=iNE(:,i)
        grid%gid(j)=gid(i)
        grid%v(j)=0d0
      end if
    end do
    deallocate(sE)
    deallocate(nNE)
    deallocate(iNE)
    deallocate(gid)
    deallocate(v)
  end subroutine
  
  !> destructor of polyGrid
  elemental subroutine purgePolyFvGrid(this)
    type(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    call this%clear()
  end subroutine
  
end module
