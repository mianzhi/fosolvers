!----------------------------------------------------------------------------- best with 100 columns

!> 1-dimensional grid
module moduleGrid1D
  private
  
  ! constants
  integer,parameter,public::BIND_NODE=1 !< bind with node
  integer,parameter,public::BIND_CELL=2 !< bind with cell
  
  !> 1-dimensional grid data and procedures
  type,public::typeGrid1D
    integer::nNode !< number of nodes
    double precision,allocatable::NodePos(:) !< node position
    integer::nCell !< number of cells
    double precision,allocatable::CellPos(:) !< cell position
    double precision,allocatable::CellWidth(:) !< cell width
    double precision,allocatable::NodeWidth(:) !< node width
  contains
    procedure,public::init=>initGrid1D
    procedure,public::clear=>clearGrid1D
    !FIXME:final::purgeGrid1D
    procedure,public::genUniform
    procedure,public::update=>updateGrid1D
  end type
  
  ! public procedures
  public::ClN
  public::CrN
  public::ClC
  public::CrC
  public::NlC
  public::NrC
  public::NlN
  public::NrN
  
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
    if(allocated(this%CellWidth)) deallocate(this%CellWidth)
    if(allocated(this%NodeWidth)) deallocate(this%NodeWidth)
  end subroutine
  
  !> destructor of typeGrid1D
  elemental subroutine purgeGrid1D(this)
    type(typeGrid1D),intent(inout)::this !< this grid
    
    call this%clear()
  end subroutine
  
  !> generate uniform grid
  elemental subroutine genUniform(this,boundL,boundR,nCell)
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
    this%NodePos(1)=boundL
    do i=1,this%nCell
      this%NodePos(i+1)=this%NodePos(i)+h
    end do
    call this%update()
  end subroutine
  
  !> update the auxiliary data of this grid
  elemental subroutine updateGrid1D(this)
    use moduleBasicDataStruct
    class(typeGrid1D),intent(inout)::this !< this grid
    
    call reallocArr(this%CellPos,this%nCell)
    call reallocArr(this%CellWidth,this%nCell)
    call reallocArr(this%NodeWidth,this%nNode)
    this%NodeWidth(:)=0d0
    do i=1,this%nCell
      this%CellPos(i)=(this%NodePos(ClN(i))+this%NodePos(CrN(i)))/2d0
      this%CellWidth(i)=abs(this%NodePos(CrN(i))-this%NodePos(ClN(i)))
      this%NodeWidth(ClN(i))=this%NodeWidth(ClN(i))+this%CellWidth(i)/2d0
      this%NodeWidth(CrN(i))=this%NodeWidth(CrN(i))+this%CellWidth(i)/2d0
    end do
  end subroutine
  
  !> cell left node
  elemental function ClN(n)
    integer,intent(in)::n !< cell index
    integer ClN !< result
    
    ClN=n
  end function
  
  !> cell right node
  elemental function CrN(n)
    integer,intent(in)::n !< cell index
    integer CrN !< result
    
    CrN=n+1
  end function
  
  !> cell left cell
  elemental function ClC(n)
    integer,intent(in)::n !< cell index
    integer ClC !< result
    
    ClC=n-1
  end function
  
  !> cell right cell
  elemental function CrC(n)
    integer,intent(in)::n !< cell index
    integer CrC !< result
    
    CrC=n+1
  end function
  
  !> node left cell
  elemental function NlC(n)
    integer,intent(in)::n !< node index
    integer NlC !< result
    
    NlC=n-1
  end function
  
  !> node right cell
  elemental function NrC(n)
    integer,intent(in)::n !< node index
    integer NrC !< result
    
    NrC=n
  end function
  
  !> node left node
  elemental function NlN(n)
    integer,intent(in)::n !< node index
    integer NlN !< result
    
    NlN=n-1
  end function
  
  !> node right node
  elemental function NrN(n)
    integer,intent(in)::n !< node index
    integer NrN !< result
    
    NrN=n+1
  end function
  
end module
