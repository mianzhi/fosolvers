!----------------------------------------------------------------------------- best with 100 columns

!***********************
! simulation conditions
!***********************
module moduleCond
  use moduleMiscDataStruct
  private
  
  ! conditions associated with nodes, facets and blocks
  type(typeDataSet),public,allocatable,save::condNode(:)
  type(typeDataSet),public,allocatable,save::condFacet(:)
  type(typeDataSet),public,allocatable,save::condBlock(:)
  
contains
end module
