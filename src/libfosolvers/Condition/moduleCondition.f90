!----------------------------------------------------------------------------- best with 100 columns

!> simulation condition
module moduleCondition
  use moduleBasicDataStruct
  private
  
  !> condition
  type,public::typeCondition
    integer Ent !< geometric entity associated with the condition
    !TODO: change to character(:) (waiting new gcc)
    character(1),allocatable::key(:) !< key of the data item
    type(typeGenDatLlist)::dat !< condition data
  contains
    !FIXME:final::purgeCondition
  end type
  
  ! procedures
  public::findCondition
  
contains
  
  !> destructor of typeCondition
  elemental subroutine purgeCondition(this)
    type(typeCondition),intent(inout)::this !< this condition
    
    call this%dat%clear()
  end subroutine
  
  !> find the index of a condition having name str and associated with entity ent from cond
  pure function findCondition(cond,ent,str)
    type(typeCondition),intent(in)::cond(:) !< the condition list to find from
    integer,intent(in)::ent !< entity
    character(*),intent(in)::str !< the name of the condition
    integer findCondition !< the index of the condition in cond
    !TODO:remove tempkey
    character(1),allocatable::tempstr(:)
    
    findCondition=0
    !TODO:remove from here
    allocate(tempstr(len(str)))
    forall(i=1:len(str))
      tempstr(i)=str(i:i)
    end forall
    !TODO:to here
    do i=1,size(cond)
      if(cond(i)%Ent==ent.and.all(cond(i)%key(:)==tempstr(:)))then
        findCondition=i
        exit
      end if
    end do
  end function
  
end module
