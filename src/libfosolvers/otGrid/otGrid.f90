!----------------------------------------------------------------------------- best with 100 columns

!> octree grid module
module modOtGrid
  use iso_c_binding,only:c_ptr,c_long
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  integer,parameter::N_FACE=6 !< number of faces per cell
  
  ! public constants
  integer,parameter,public::OT_OUT=0
  integer,parameter,public::OT_IN=1
  integer,parameter,public::OT_CUT=2
  
  !> octree grid
  type,public::otGrid
    ! basic data
    double precision::p0(DIMS) !< minimum coordinate
    double precision::h0 !< edge length of the root
    type(c_ptr),private::otree !< pointer to the c octree root node
    integer::nC !< number of cells
    integer(c_long),allocatable::oid(:) !< octal ids
    integer,allocatable::lvl(:) !< tree levels
    integer,allocatable::neib(:,:) !< neighbor lists
    integer,allocatable::flag(:) !< flag (0-out, 1-in, 2-cut)
  contains
    procedure,public::init=>initOtGrid
    procedure,public::clear=>clearOtGrid
    procedure,public::p=>pOtGrid
    procedure,public::h=>hOtGrid
    procedure,public::a=>aOtGrid
    procedure,public::v=>vOtGrid
    !FIXME:final::purgeOtGrid
  end type
  
  ! interface to the otree code in c
  interface
    !> initialize tree
    function ot_init() bind(c,name='init')
      use iso_c_binding
      type(c_ptr)::ot_init !< root pointer
    end function
    
    !> clear tree recursively
    subroutine ot_clear(root) bind(c,name='clear')
      use iso_c_binding
      type(c_ptr),value::root !< root pointer
    end subroutine
    
    !> grow tree by given level
    subroutine ot_grow(root,lvl) bind(c,name='grow')
      use iso_c_binding
      type(c_ptr),value::root !< root pointer
      integer(c_int),value::lvl !< levels to grow
    end subroutine
    
    !> split node
    subroutine ot_split(root,i) bind(c,name='splitId')
      use iso_c_binding
      type(c_ptr),value::root !< root pointer
      integer(c_long),value::i !< octal id
    end subroutine
    
    !> trim node
    subroutine ot_trim(root,i) bind(c,name='trimId')
      use iso_c_binding
      type(c_ptr),value::root !< root pointer
      integer(c_long),value::i !< octal id
    end subroutine
    
    !> number of leaves
    function ot_nLf(root) bind(c,name='nLf')
      use iso_c_binding
      type(c_ptr),value::root !< root pointer
      integer(c_int)::ot_nLf !< number of leaves
    end function
    
    !> generate connectivity tables
    subroutine ot_genConn(root,id,lvl,neib) bind(c,name='genConn')
      use iso_c_binding
      type(c_ptr),value::root !< root pointer
      integer(c_long)::id(*) !< octal id
      integer(c_int)::lvl(*) !< level
      integer(c_int)::neib(*) !< neighbor list
    end subroutine
    
    !> number of binary digits of a given integer
    function ot_nBinDigits(i) bind(c,name='nBinDigits')
      use iso_c_binding
      integer(c_long),value::i !< given integer
      integer(c_int)::ot_nBinDigits !< number of binary digits
    end function
    
    !> physical position of an octal id
    subroutine ot_getPos(i,p) bind(c,name='getPos')
      use iso_c_binding
      integer(c_long),value::i !< octal id
      real(c_double)::p(3) !< position
    end subroutine
  end interface
  
contains
  
  !> initialize this otGrid
  subroutine initOtGrid(this,p,h,lvl)
    class(otGrid),intent(inout)::this !< this otGrid
    double precision,intent(in)::p(DIMS) !< minimum coordinate
    double precision,intent(in)::h !< edge length of the root
    integer,optional::lvl !< initial tree depth
    
    if(present(lvl))then
      l=lvl
    else
      l=3
    end if
    call this%clear()
    this%p0(:)=p(:)
    this%h0=h
    this%otree=ot_init()
    call ot_grow(this%otree,l)
    this%nC=ot_nLf(this%otree)
    allocate(this%oid(this%nC))
    allocate(this%lvl(this%nC))
    allocate(this%neib(N_FACE,this%nC))
    allocate(this%flag(this%nC))
    call ot_genConn(this%otree,this%oid,this%lvl,this%neib)
    this%flag(:)=OT_IN
  end subroutine
  
  !> clear this otGrid
  elemental subroutine clearOtGrid(this)
    class(otGrid),intent(inout)::this !< this otGrid
    
    if(allocated(this%oid)) deallocate(this%oid)
    if(allocated(this%lvl)) deallocate(this%lvl)
    if(allocated(this%neib)) deallocate(this%neib)
    if(allocated(this%flag)) deallocate(this%flag)
  end subroutine
  
  !> position of cell k
  function pOtGrid(this,k)
    class(otGrid),intent(in)::this !< this otGrid
    integer,intent(in)::k !< cell index
    double precision::pOtGrid(DIMS) !< result
    double precision::p(DIMS)
    
    call ot_getPos(this%oid(k),p)
    pOtGrid(:)=this%p0(:)+this%h0*p(:)
  end function
  
  !> edge length of cell k
  elemental function hOtGrid(this,k)
    class(otGrid),intent(in)::this !< this otGrid
    integer,intent(in)::k !< cell index
    double precision::hOtGrid !< result
    
    hOtGrid=this%h0/dble(2**this%lvl(k))
  end function
  
  !> face area of cell k
  elemental function aOtGrid(this,k)
    class(otGrid),intent(in)::this !< this otGrid
    integer,intent(in)::k !< cell index
    double precision::aOtGrid !< result
    
    aOtGrid=(this%h0/dble(2**this%lvl(k)))**2
  end function
  
  !> volume of cell k
  elemental function vOtGrid(this,k)
    class(otGrid),intent(in)::this !< this otGrid
    integer,intent(in)::k !< cell index
    double precision::vOtGrid !< result
    
    vOtGrid=(this%h0/dble(2**this%lvl(k)))**3
  end function
  
  !> destructor of otGrid
  elemental subroutine purgeOtGrid(this)
    type(otGrid),intent(inout)::this !< this otGrid
    
    call this%clear()
  end subroutine
  
end module
