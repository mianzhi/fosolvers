!----------------------------------------------------------------------------- best with 100 columns

!> 1-dimensional grid
module moduleGrid1D
  private
  
  !> 1-dimensional grid data and procedures
  type,public::typeGrid1D
    integer::nNode !< number of nodes
    double precision,allocatable::NodePos(:) !< node position
    integer::nCell !< number of cells
    double precision,allocatable::CellPos(:) !< cell position
  contains
    procedure,public::init=>initGrid1D
    procedure,public::clear=>clearGrid1D
    !FIXME:final::purgeGrid1D
    procedure,public::genUniform=>genUniformGrid1D
  end type
  
contains
  
  !> initialize this Grid1D
  elemental subroutine initGrid1D(this)
    class(typeGrid1D),intent(inout)::this !< this grid
    
    this%nNode=0
    this%nCell=0
    call this%clear()
  end subroutine
  
  !> clear this Grid1D
  elemental subroutine clearGrid1D(this)
    class(typeGrid1D),intent(inout)::this !< this grid
    
    if(allocated(this%NodePos)) deallocate(this%NodePos)
    if(allocated(this%CellPos)) deallocate(this%CellPos)
  end subroutine
  
  !> destructor of typeGrid1D
  elemental subroutine purgeGrid1D(this)
    type(typeGrid1D),intent(inout)::this !< this grid
    
    call this%clear()
  end subroutine
  
  !> generate uniform grid
  elemental subroutine genUniformGrid1D(this,boundL,boundR,nCell)
    class(typeGrid1D),intent(inout)::this !< this grid
    double precision,intent(in)::boundL !< left bound
    double precision,intent(in)::boundR !< right bound
    integer,intent(in)::nCell !< number of cells
    double precision h
    
    call this%clear()
    this%nNode=nCell+1
    this%nCell=nCell
    h=(boundR-boundL)/dble(nCell)
    allocate(this%NodePos(this%nNode))
    allocate(this%CellPos(this%nCell))
    this%NodePos(1)=boundL
    do i=1,this%nCell
      this%NodePos(i+1)=this%NodePos(i)+h
      this%CellPos(i)=this%NodePos(i)+h/2d0
    end do
  end subroutine
  
end module
