!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  use moduleFileIO
  use moduleMiscDataStruct
  
!  class(typeDataSet),pointer::vp(3)
!  type(typeDataSetVal)::v
!  v%Val=2
!  allocate(vp,source=v)
!  write(*,*),vp(1)%lookup()
  
  double precision,target,allocatable::v(:),vp(:,:)
  
  call readmsh('grid.msh',50,verbose=.true.)
  
  allocate(v(nFacet))
  allocate(vp(nNode,DIMS))
  call addWrite(v)
  call addWrite(vp)
  do i=1,nFacet
    v(i)=norm2(Facet(i)%PC)
  end do
  do i=1,nNode
    vp(i,:)=Node(i)%Pos(:)
  end do
  
  call writerst('rst.msh',55)
end program
