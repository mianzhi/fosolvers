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
  pure subroutine transShape(t,n)
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

!> write VTK into fid from poly
subroutine writeVTKPolyX(fid,poly)
  use modPolyX
  integer,intent(in)::fid !< file id
  class(polyX),intent(in)::poly !< source
  
  write(fid,'(a)'),'# vtk DataFile Version 3.0'
  write(fid,'(a)'),'Generated by FOSolverS'
  write(fid,'(a)'),'ASCII'
  write(fid,'(a)'),''
  write(fid,'(a)'),'DATASET UNSTRUCTURED_GRID'
  write(fid,'(a,i10,a)'),'POINTS ',poly%nN,' double'
  do i=1,poly%nN
    write(fid,*),poly%pN(:,i)
  end do
  write(fid,'(a)'),''
  write(fid,'(a,i10,i11)'),'CELLS',poly%nE,poly%nE+sum(poly%nNE)
  do i=1,poly%nE
    write(fid,*),poly%nNE(i),poly%iNE(:,i)-1
  end do
  write(fid,'(a)'),''
  write(fid,'(a,i10)'),'CELL_TYPES',poly%nE
  do i=1,poly%nE
    write(fid,*),transShape(poly%sE(i))
  end do
  
contains
  
  !> translate to the shape type of VTK
  pure function transShape(s)
    use modPolyMesh,only:TRI,QUAD
    use modPolyGrid,only:TET,HEX
    integer,intent(in)::s !< shape type
    integer::transShape !< shape type in VTK
    
    select case(s)
    case(TRI)
      transShape=5
    case(QUAD)
      transShape=9
    case(TET)
      transShape=10
    case(HEX)
      transShape=12
    case default
      transShape=0
    end select
  end function
  
end subroutine

!> write VTK data header into fid
subroutine writeVTKHead(fid,poly,k)
  use modPolyX
  integer,intent(in)::fid !< file id
  class(polyX),intent(in)::poly !< polyX
  integer,intent(in)::k !< header switch
  integer,parameter::N_DATA=1 !< node data
  integer,parameter::E_DATA=2 !< element data
  
  write(fid,'(a)'),''
  select case(k)
  case(N_DATA)
    write(fid,'(a,i10)'),'POINT_DATA',poly%nN
  case(E_DATA)
    write(fid,'(a,i10)'),'CELL_DATA',poly%nE ! cell means element in VTK
  case default
  end select
end subroutine

!> write VTK into fid from nodal scalar field
subroutine writeVTKScal(fid,key,a)
  integer,intent(in)::fid !< file id
  character(*),intent(in)::key !< data name
  double precision,intent(in)::a(:) !< data
  
  write(fid,'(a)'),''
  write(fid,'(a)'),'SCALARS '//trim(adjustl(key))//' double'
  write(fid,'(a)'),'LOOKUP_TABLE default'
  do i=1,size(a)
    write(fid,*),a(i)
  end do
end subroutine

!> write VTK into fid from nodal vector field
subroutine writeVTKVect(fid,key,a)
  integer,intent(in)::fid !< file id
  character(*),intent(in)::key !< data name
  double precision,intent(in)::a(:,:) !< data
  
  write(fid,'(a)'),''
  write(fid,'(a)'),'VECTORS '//trim(adjustl(key))//' double'
  do i=1,size(a,2)
    write(fid,*),a(:,i)
  end do
end subroutine
