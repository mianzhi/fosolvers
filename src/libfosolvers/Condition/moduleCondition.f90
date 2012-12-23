!----------------------------------------------------------------------------- best with 100 columns

!> simulation condition
module moduleCondition
  use moduleBasicDataStruct
  private
  
  !> condition
  type,public::typeCondition
    integer Ent !< geometric entity associated with the condition
    integer bind !< bind with node, line, facet or block
    type(typeGenDatLlist)::dat !< condition data
  contains
    !FIXME:final::purgeCondition
  end type
  
contains
  
  !> destructor of typeCondition
  elemental subroutine purgeCondition(this)
    type(typeCondition),intent(inout)::this !< this condition
    
    call this%dat%clear()
  end subroutine
  
end module
