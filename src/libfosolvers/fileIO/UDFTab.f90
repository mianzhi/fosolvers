!----------------------------------------------------------------------------- best with 100 columns

!> read UDFs from fid into udf
subroutine readUDFTab(fid,udf)
  use modUDF
  use iso_c_binding
  integer,intent(in)::fid !< file id
  type(UDFTab),intent(inout)::udf !< result
  integer,parameter::STR_LEN=400
  character(len=STR_LEN,kind=C_CHAR)::str
  
  read(fid,*)n
  call udf%init(n)
  do i=1,n
    str=''
    read(fid,'(a)')str
    udf%ptr(i)=evaluator_create(trim(adjustl(str))//C_NULL_CHAR)
  end do
end subroutine
