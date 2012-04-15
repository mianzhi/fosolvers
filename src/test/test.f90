!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testList()
  call testSerialPack()
end program

subroutine testList()
  use utilityTest
  use moduleBasicDataStruct
  type(typeListIScal)::ListIScal(2)
  type(typeListDScal)::ListDScal(2)
  
  call ListIScal(1)%extend(2)
  call ListIScal(1)%push(1)
  call ListIScal(2)%push([1,2,3,4,5,6])
  if(ListIScal(1)%get(1)==1.and.&
  &  all(ListIScal(2)%get([1,3,6])==[1,3,6]).and.&
  &  size(ListIScal(1)%dat)==2.and.&
  &  size(ListIScal(2)%dat)==8.and.&
  &  ListIScal(2)%length==6)then
    call showPass('ListIScal push and get')
  else
    call showFail('ListIScal push and get')
  end if
  call ListIScal(1)%clear()
  if(ListIScal(1)%length==0.and..not.allocated(ListIScal(1)%dat))then
    call showPass('ListIScal clear')
  end if
  
  call ListDScal(1)%extend(2)
  call ListDScal(1)%push(1d0)
  call ListDScal(2)%push([1d0,2d0,3d0,4d0,5d0,6d0])
  if(abs(ListDScal(1)%get(1)-1d0)<TOLERANCE.and.&
  &  norm2(ListDScal(2)%get([1,3,6])-[1d0,3d0,6d0])<TOLERANCE.and.&
  &  size(ListDScal(1)%dat)==2.and.&
  &  size(ListDScal(2)%dat)==8.and.&
  &  ListDScal(2)%length==6)then
    call showPass('ListDScal push and get')
  else
    call showFail('ListDScal push and get')
  end if
  call ListDScal(1)%clear()
  if(ListDScal(1)%length==0.and..not.allocated(ListDScal(1)%dat))then
    call showPass('ListDScal clear')
  end if
end subroutine

subroutine testSerialPack()
  use utilityTest
  use moduleBasicDataStruct
  type(typeSerialPack)::sp
  type(typeListIScal)::list(3)
  
  call sp%clear()
  call list(1)%push([3,4,5])
  call list(1)%addto(sp)
  call list(2)%addto(sp)
  call list(1)%addto(sp)
  call sp%resetPtr()
  call list(3)%recover(sp)
  call list(1)%recover(sp)
  call list(2)%recover(sp)
  if(all(list(3)%get([(i,i=1,list(3)%length)])==[3,4,5]).and.&
  &  list(1)%length==0.and.&
  &  all(list(2)%get([(i,i=1,list(2)%length)])==[3,4,5]))then
    call showPass('SerialPack store ListIScal')
  else
    call showFail('SerialPack store ListIScal')
  end if
  
  call sp%clear()
  if(sp%iPtr==0.and.sp%dPtr==0)then
    call showPass('SerialPack clear')
  else
    call showFail('SerialPack clear')
  end if
end subroutine
