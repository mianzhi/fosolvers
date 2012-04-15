!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testList()  
end program

subroutine testList()
  use utilityTest
  use moduleBasicDataStruct
  type(typeListIntegerScal)::ListIntegerScal(2)
  type(typeListDoubleScal)::ListDoubleScal(2)
  
  call ListIntegerScal(1)%extend(2)
  call ListIntegerScal(1)%push(1)
  call ListIntegerScal(2)%push([1,2,3,4,5,6])
  if(ListIntegerScal(1)%get(1)==1.and.&
  &  all(ListIntegerScal(2)%get([1,3,6])==[1,3,6]).and.&
  &  size(ListIntegerScal(1)%dat)==2.and.&
  &  size(ListIntegerScal(2)%dat)==8.and.&
  &  ListIntegerScal(2)%length==6)then
    call showPass('ListIntegerScal push and get')
  else
    call showFail('ListIntegerScal push and get')
  end if
  call ListIntegerScal(1)%clear()
  if(ListIntegerScal(1)%length==0.and..not.allocated(ListIntegerScal(1)%dat))then
    call showPass('ListIntegerScal clear')
  end if
  
  call ListDoubleScal(1)%extend(2)
  call ListDoubleScal(1)%push(1d0)
  call ListDoubleScal(2)%push([1d0,2d0,3d0,4d0,5d0,6d0])
  if(abs(ListDoubleScal(1)%get(1)-1d0)<TOLERANCE.and.&
  &  norm2(ListDoubleScal(2)%get([1,3,6])-[1d0,3d0,6d0])<TOLERANCE.and.&
  &  size(ListDoubleScal(1)%dat)==2.and.&
  &  size(ListDoubleScal(2)%dat)==8.and.&
  &  ListDoubleScal(2)%length==6)then
    call showPass('ListDoubleScal push and get')
  else
    call showFail('ListDoubleScal push and get')
  end if
  call ListDoubleScal(1)%clear()
  if(ListDoubleScal(1)%length==0.and..not.allocated(ListDoubleScal(1)%dat))then
    call showPass('ListDoubleScal clear')
  end if
end subroutine
