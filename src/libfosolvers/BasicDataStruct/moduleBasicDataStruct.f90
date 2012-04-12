!----------------------------------------------------------------------------- best with 100 columns

!*******************
! Basic Data Struct
!*******************
module moduleBasicDataStruct
  private
  
  !-------------------------
  ! list of integer scalars
  !-------------------------
  type,public::typeListIntegerScal
    integer,allocatable,private::dat(:)
    integer,public::length
  contains
  end type
!  interface typeListIntegerScal
!    module procedure::constructListIntegerScal
!  end interface
  
contains
  
  !------------------------------------
  ! constructor of typeListIntegerScal
  !------------------------------------
  elemental function constructListIntegerScal(initSpace)
    integer,intent(in),optional::initSpace
    type(typeListIntegerScal)::constructListIntegerScal
    
    constructListIntegerScal%length=0
    if(present(initSpace))then
      allocate(constructListIntegerScal%dat(initSpace))
    end if
  end function
  
end module
