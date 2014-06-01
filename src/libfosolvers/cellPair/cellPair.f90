!----------------------------------------------------------------------------- best with 100 columns

!> cell pair module (for FVM flux)
module modCellPair
  private
  
  !> constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> pairs of cells
  type::cellPair
    integer::nP !< number of pairs
    integer,allocatable::iC(:,:) !< index of cells
  contains
    procedure::initCellPair
    procedure::clear=>clearCellPair
    !FIXME:final::purgeCellPair
  end type
  
  !> pairs of polyGrid cells
  type,extends(cellPair),public::polyPair
    double precision,allocatable::a(:) !< areas of interface
    double precision,allocatable::n(:,:) !< normal vectors
    double precision,allocatable::v(:,:,:) !< vectors cell_1 -> interface -> cell_2 
  contains
    procedure,public::init=>initPolyPair
    procedure,public::clear=>clearPolyPair
    !FIXME:final::purgePolyPair
  end type
  
  !> pairs of otGrid cells
  type,extends(cellPair),public::otPair
    integer,allocatable::dir(:) !< directions of interface
  contains
    procedure,public::init=>initOtPair
    procedure,public::clear=>clearOtPair
    !FIXME:final::purgeOtPair
  end type
  
contains
  
  !> initialize this cellPair
  elemental subroutine initCellPair(this,nP)
    class(cellPair),intent(inout)::this !< this cellPair
    integer,intent(in)::nP !< number of pairs
    
    call this%clear()
    this%nP=nP
    allocate(this%iC(2,nP))
  end subroutine
  
  !> clear this cellPair
  elemental subroutine clearCellPair(this)
    class(cellPair),intent(inout)::this !< this cellPair
    
    if(allocated(this%iC)) deallocate(this%iC)
  end subroutine
  
  !> destructor of cellPair
  elemental subroutine purgeCellPair(this)
    type(cellPair),intent(inout)::this !< this cellPair
    
    call this%clear()
  end subroutine
  
  !> initialize this polyPair
  elemental subroutine initPolyPair(this,grid)
    use modPolyGrid
    class(polyPair),intent(inout)::this !< this polyPair
    type(polyGrid),intent(in)::grid !< the input grid
    
    call this%initCellPair(grid%nE)
    allocate(this%a(grid%nE))
    allocate(this%n(DIMS,grid%nE))
    allocate(this%v(DIMS,2,grid%nE))
    !TODO:fillup data (nP,iC,a,n,v)
  end subroutine
  
  !> clear this polyPair
  elemental subroutine clearPolyPair(this)
    class(polyPair),intent(inout)::this !< this polyPair
    
    call this%cellPair%clear()
    if(allocated(this%a)) deallocate(this%a)
    if(allocated(this%n)) deallocate(this%n)
    if(allocated(this%v)) deallocate(this%v)
  end subroutine
  
  !> destructor of polyPair
  elemental subroutine purgePolyPair(this)
    type(polyPair),intent(inout)::this !< this polyPair
    
    call this%clear()
  end subroutine
  
  !> initialize this otPair
  elemental subroutine initOtPair(this,grid)
    use modOtGrid
    class(otPair),intent(inout)::this !< this otPair
    type(otGrid),intent(in)::grid !< the input grid
    
    call this%initCellPair(grid%nC)
    allocate(this%dir(grid%nC))
    !TODO:fillup data (nP,iC,dir)
  end subroutine
  
  !> clear this otPair
  elemental subroutine clearOtPair(this)
    class(otPair),intent(inout)::this !< this otPair
    
    call this%cellPair%clear()
    if(allocated(this%dir)) deallocate(this%dir)
  end subroutine
  
  !> destructor of otPair
  elemental subroutine purgeOtPair(this)
    type(otPair),intent(inout)::this !< this otPair
    
    call this%clear()
  end subroutine
  
end module
