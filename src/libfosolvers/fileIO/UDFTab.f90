!----------------------------------------------------------------------------- best with 100 columns

!> read UDFs from fid into udf
subroutine readUDFTab(fid,udf)
  use modUDF
  use iso_c_binding
  integer,intent(in)::fid !< file id
  type(UDFTab),intent(inout)::udf !< result
  integer,parameter::STR_LEN=400
  character(STR_LEN)::str
  integer(C_LONG)::evaluator_create_
  
  read(fid,*)n
  call udf%init(n)
  do i=1,n
    str=''
    read(fid,'(a)')str
    udf%ptr(i)=evaluator_create_(trim(adjustl(str)))
  end do
end subroutine
