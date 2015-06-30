!----------------------------------------------------------------------------- best with 100 columns

function UDF1()
  use modFileIO
  use modUDF
  integer UDF1
  type(UDFTab)::udf
  
  UDF1=0
  open(10,file='data/udf',action='read')
  call readUDFTab(10,udf)
  if(abs(udf%eval(1,[1d0,2d0,3d0],4d0)-1234d0)>1d-8)then
    UDF1=UDF1+1
  end if
  if(abs(udf%eval(2,[1d0,2d0,3d0],4d0)-4d0)>1d-8)then
    UDF1=UDF1+10
  end if
  if(abs(udf%eval(3,[1d0,2d0,3d0],4d0)-1d0)>1d-8)then
    UDF1=UDF1+100
  end if
  call udf%clear()
  close(10)
end function
