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
  end type
  type(typeCondList),public,allocatable,save::CondNode(:),CondFacet(:),CondEle(:)
  
end module
