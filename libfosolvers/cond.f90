!----------------------------------------------------------------------------- best with 100 columns

!***********************************
! simulation condition related data
!***********************************
module moduleCond
  private
  
  integer,parameter,public::COND_NAME_LEN=2
  
  !----------------------------------
  ! simulation condition information
  !----------------------------------
  type,public::typeCond
    ! type of condition can be stored in 'what'
    !   e.g. what='Dr' ~ Dirichlet
    character(COND_NAME_LEN),public::what
    double precision,public,allocatable::Val(:)
    integer,public,allocatable::Tab(:)
  end type
  
  type,public::typeCondList
    type(typeCond),public,allocatable::Cond(:)
  contains
    procedure,public::getSpace=>getCondSpace
  end type
  type(typeCondList),public,allocatable,save::CondNode(:),CondFacet(:),CondEle(:)
  
contains
  
  function getCondSpace(this)
    class(typeCondList),intent(inout)::this
    integer getCondSpace
    type(typeCond),allocatable::tempCond(:)
    
    getCondSpace=0
    
    if(allocated(this%Cond))then
      getCondSpace=size(this%Cond)+1
      allocate(tempCond(getCondSpace))
      tempCond(1:getCondSpace-1)=this%Cond(:)
      call move_alloc(tempCond,this%Cond)
    else
      allocate(this%Cond(1))
      getCondSpace=1
    end if
  end function
end module
