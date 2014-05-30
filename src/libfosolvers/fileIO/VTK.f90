!----------------------------------------------------------------------------- best with 100 columns

!> read VTK from fid into poly
subroutine readVTKPolyX(fid,poly)
  use modPolyX
  integer,intent(in)::fid !< file id
  class(polyX),intent(inout)::poly !< result
  integer::ierr
  character(400)::tmpStr
  double precision,allocatable::p(:)
  integer,allocatable::c(:),t(:),np(:)
  
  rewind(fid)
  ierr=0
  do while(ierr==0)
    tmpStr=''
    read(fid,'(a)',iostat=ierr),tmpStr
    if(tmpStr(1:6)=='POINTS')then
      read(tmpStr(7:),*),n
      allocate(p(n*3))
      read(fid,*),p(:)
    end if
    if(tmpStr(1:5)=='CELLS')then
      read(tmpStr(6:),*),m,l
      allocate(c(l))
      read(fid,*),c(:)
    end if
    if(tmpStr(1:10)=='CELL_TYPES')then
      read(tmpStr(11:),*),m
      allocate(t(m))
      read(fid,*),t(:)
    end if
  end do
  allocate(np(m))
  do i=1,m
    call transShape(t(i),np(i))
  end do
  call poly%init(n,m,maxval(np))
  poly%pN(:,:)=reshape(p,shape(poly%pN))
  poly%sE(:)=t(:)
  poly%nNE(:)=np(:)
  j=1
  do i=1,m
    poly%iNE(1:c(j),i)=c(j+1:j+c(j))+1
    j=j+c(j)+1
  end do
  deallocate(p)
  deallocate(c)
  deallocate(t)
  deallocate(np)

contains
  
  !> translate the shape type of VTK
  subroutine transShape(t,n)
    use modPolyMesh,only:TRI,TRI_N,QUAD,QUAD_N
    use modPolyGrid,only:TET,TET_N,HEX,HEX_N
    integer,intent(inout)::t !< shape type
    integer,intent(out)::n !< number of nodes
    
    select case(t)
    case(5)
      t=TRI
      n=TRI_N
    case(9)
      t=QUAD
      n=QUAD_N
    case(10)
      t=TET
      n=TET_N
    case(12)
      t=HEX
      n=HEX_N
    case default
      t=0
      n=0
    end select
  end subroutine
  
end subroutine
