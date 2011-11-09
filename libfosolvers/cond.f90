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
  
  !---------------------
  ! 1-dimensional table
  !---------------------
  type::typeTab1d
    integer length
    double precision,pointer::x(:)
    double precision,pointer::y(:)
  contains
    procedure::lookup=>lookupTab1d
  end type
  type(typeTab1d),allocatable,save::dataTab1d(:)
  
contains
  
  !-----------------------------------------------------------
  ! look up 1-dimensional table with independent variable pos
  !-----------------------------------------------------------
  function lookupTab1d(this,pos)
    class(typeTab1d),intent(in)::this
    double precision,intent(in)::pos
    double precision lookupTab1d
    lookupTab1d=0d0
    do i=2,this%length-1
      if(pos<this%x(i))then
        lookupTab1d=(pos-this%x(i-1))/(this%x(i)-this%x(i-1))*(this%y(i)-this%y(i-1))+this%y(i-1)
        exit
      end if
    end do
    if(i==this%length)then
      lookupTab1d=(pos-this%x(i-1))/(this%x(i)-this%x(i-1))*(this%y(i)-this%y(i-1))+this%y(i-1)
    end if
  end function
end module
