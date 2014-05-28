!----------------------------------------------------------------------------- best with 100 columns

!> read GTS from fid into mesh
subroutine readGTS(fid,mesh)
  use modpolyMesh
  integer,intent(in)::fid !< file id
  type(polyMesh),intent(inout)::mesh !< result
  integer::nPoint,nEdge,nTri
  integer,allocatable::lEdge(:,:),lTri(:,:)
  
  read(fid,*),nPoint,nEdge,nTri
  allocate(lEdge(2,nEdge))
  allocate(lTri(3,nTri))
  call mesh%init(nPoint,nTri,TRI_N)
  mesh%sE(:)=TRI
  mesh%nNE(:)=TRI_N
  mesh%iNE(:,:)=0
  do i=1,nPoint
    read(fid,*),mesh%pN(:,i)
  end do
  do i=1,nEdge
    read(fid,*),lEdge(:,i)
  end do
  do i=1,nTri
    read(fid,*),lTri(:,i)
  end do
  do i=1,nTri
    if(all(lEdge(:,lTri(2,i))/=lEdge(1,lTri(1,i))))then
      mesh%iNE(1:2,i)=lEdge(:,lTri(1,i))
    else
      mesh%iNE(1:2,i)=[lEdge(2,lTri(1,i)),lEdge(1,lTri(1,i))]
    end if
    n1=lEdge(1,lTri(2,i))
    n2=lEdge(2,lTri(2,i))
    mesh%iNE(3,i)=merge(n1,n2,all(lEdge(:,lTri(1,i))/=n1))
  end do
  deallocate(lEdge)
  deallocate(lTri)
end subroutine
