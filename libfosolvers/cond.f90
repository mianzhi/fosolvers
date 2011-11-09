!----------------------------------------------------------------------------- best with 100 columns

!***********************************
! simulation condition related data
!***********************************
module moduleCond
  
  !----------------------------------
  ! simulation condition information
  !----------------------------------
  type::typeCond
    ! type of condition can be stored in 'what'
    !   e.g. what='Dr' ~ Dirichlet
    character(2) what
    double precision,allocatable::Val(:)
    integer,allocatable::Tab(:)
  end type
  
  type::typeCondList
    type(typeCond),allocatable::Cond(:)
  end type
  type(typeCondList),allocatable,save::CondNode(:),CondFacet(:),CondEle(:)
  
end module
