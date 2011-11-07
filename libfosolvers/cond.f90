!----------------------------------------------------------------------------- best with 100 columns

!***********************************
! simulation condition related data
!***********************************
module moduleCond
  
  !----------------------------------
  ! simulation condition information
  !----------------------------------
  type::typeCond
    ! geometric entity where the condition is applied
    integer GeoEnti
    ! type of condition can be stored in what
    !   e.g. what='Dr' ~ Dirichlet
    character(2) what
    ! four double float numbers is available for parameters of the condition
    !   e.g. val=1d0 ~ Dirichlet value is 1d0
    double precision val
    double precision val2
    double precision val3
    double precision val4
    ! index of 1-dimensional tables associated with the condition
    !   e.g. tab=10 ~ source value which is a function of t is stored in table 10
    integer tab
    integer tab2
  end type
  type(typeCond),allocatable,save::Conditions(:)
  
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
