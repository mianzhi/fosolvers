!----------------------------------------------------------------------------- best with 100 columns

!***************
! material data
!***************
module moduleMtl
  use moduleMiscDataStruct
  private
  
  ! material data
  type(typeDataSet),public,allocatable,save::Mtl(:)

contains
end module
