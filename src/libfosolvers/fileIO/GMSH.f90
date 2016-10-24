!----------------------------------------------------------------------------- best with 100 columns

!> read GMSH from fid into polyGrid
subroutine readGMSHPolyGrid(fid,grid)
  use modPolyGrid
  integer,intent(in)::fid !< file id
  class(polyGrid),intent(inout)::grid !< result
  integer::ierr
  character(400)::tmpStr
  double precision,allocatable::p(:,:)
  integer,allocatable::c(:,:),t(:),np(:),g(:)
  integer::tmpIArr(20)
  integer,parameter::MAX_N_PER_E=27
  
  rewind(fid)
  ierr=0
  do while(ierr==0)
    tmpStr=''
    read(fid,'(a)',iostat=ierr)tmpStr
    if(tmpStr(1:6)=='$Nodes')then
      read(fid,*)n
      allocate(p(3,n))
      do i=1,n
        read(fid,*)j,p(:,i)
      end do
    end if
    if(tmpStr(1:9)=='$Elements')then
      read(fid,*)m
      allocate(c(MAX_N_PER_E,m))
      allocate(t(m))
      allocate(np(m))
      allocate(g(m))
      c(:,:)=0
      do i=1,m
        tmpStr=''
        read(fid,'(a)')tmpStr
        read(tmpStr,*)j,t(i),l
        call transShape(t(i),np(i))
        read(tmpStr,*)j,k,tmpIArr(1:l+1),c(1:np(i),i)
        g(i)=tmpIArr(3)
      end do
      m=count(np(:)>0) ! exclude un-supported elements
      call grid%init(n,m,maxval(np))
      if(allocated(p))then
        grid%pN(:,:)=p(:,1:grid%nN)
      else
        write(*,'(a)')"[E] readGMSHPolyGrid(): reading elements before nodes"
      end if
      j=0
      do i=1,size(np)
        if(np(i)>0)then
          j=j+1
          grid%sE(j)=t(i)
          grid%nNE(j)=np(i)
          grid%iNE(1:grid%nNE(j),j)=c(1:grid%nNE(j),i)
          grid%gid(j)=g(i)
        end if
      end do
    end if
  end do
  deallocate(p)
  deallocate(c)
  deallocate(t)
  deallocate(np)
  deallocate(g)

contains
  
  !> translate the shape type of GMSH
  pure subroutine transShape(t,n)
    use modPolyGrid,only:TRI,TRI_N,QUAD,QUAD_N,TET,TET_N,HEX,HEX_N
    integer,intent(inout)::t !< shape type
    integer,intent(out)::n !< number of nodes
    
    select case(t)
    case(2)
      t=TRI
      n=TRI_N
    case(3)
      t=QUAD
      n=QUAD_N
    case(4)
      t=TET
      n=TET_N
    case(5)
      t=HEX
      n=HEX_N
    case default
      t=0
      n=0
    end select
  end subroutine
  
end subroutine
