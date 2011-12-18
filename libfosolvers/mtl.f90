!----------------------------------------------------------------------------- best with 100 columns

!***************
! material data
!***************
module moduleMtl
  use moduleMiscDataStruct
  private
  
  ! constants
  
  ! material data
  type(typeDataSet),public,allocatable,save::Mtl(:)

contains
end module
