!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  use moduleFileIO
  use moduleMiscDataStruct
  
  type(typeDataItem)::v
  call v%specify(TAB1D_TYPE,5)
  v%Tab1d(:,1)=[1d0,2d0,3d0,4d0,5d0]
  v%Tab1d(:,2)=[1d1,2d1,3d1,4d1,5d1]
  write(*,*),v%lookup(-2d0)
  
!  double precision,target,allocatable::v(:),vp(:,:)
!  
!  call readmsh('grid.msh',50,verbose=.true.)
!  
!  allocate(v(nFacet))
!  allocate(vp(nNode,DIMS))
!  call addWrite(v)
!  call addWrite(vp)
!  do i=1,nFacet
!    v(i)=norm2(Facet(i)%PC)
!  end do
!  do i=1,nNode
!    vp(i,:)=Node(i)%Pos(:)
!  end do
!  
!  call writerst('rst.msh',55)
end program
