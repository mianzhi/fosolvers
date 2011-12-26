!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  use moduleFileIO
  use moduleMtl
  use moduleCond
  
  call readMsh('grid.msh',50,verbose=.true.)
  call readCond('conditions.cod',52)
  write(*,*),condNode(512)%lookup('vname'),condNode(512)%lookup('tname',5d-1)
  write(*,*),condFacet(322)%lookup('fffc1'),condFacet(322)%lookup('fffc2',0.6d0)
  
!  call readMtl('materials.mtl',51)
!  write(*,*),size(Mtl),size(Mtl(1)%DataItem)
!  write(*,*),Mtl(1)%lookup('Denst')
!  write(*,*),Mtl(1)%lookup('SpcHt',475d0)
!  write(*,*),Mtl(1)%lookup('ThrmCd',5d2)
!  write(*,*),Mtl(1)%lookup('YounM')
!  write(*,*),Mtl(1)%lookup('PoisR')
!  write(*,*),Mtl(1)%lookup('Stren')
  
!  type(typeDataSet)::v
!  type(typeDataItem)::ve
!  
!  write(*,*),v%lookup('name3',stat=ierr)
!  write(*,*),ierr
!  
!  call v%extend(1)
!  call v%ptrLast%specify(TAB1D_TYPE,5)
!  v%ptrLast%DataName='name1'
!  v%ptrLast%Tab1d(:,1)=[1d0,2d0,3d0,4d0,5d0]
!  v%ptrLast%Tab1d(:,2)=[1d1,2d1,3d1,4d1,5d1]
!  
!  call v%extend(1)
!  call v%ptrLast%specify(VAL_TYPE)
!  v%ptrLast%DataName='name2'
!  v%ptrLast%Val=250d0
!  
!  call ve%specify(TAB1D_TYPE,5)
!  ve%DataName='name4'
!  ve%Tab1d(:,1)=[1d0,2d0,3d0,4d0,5d0]
!  ve%Tab1d(:,2)=[1d2,2d2,3d2,4d2,5d2]
!  call v%push(ve)
!  v%ptrLast%DataName='name4'
!  
!  write(*,*),v%lookup('name1',-2d0)
!  write(*,*),v%lookup('name2')
!  write(*,*),v%lookup('name4',3d0)
  
!  double precision,target,allocatable::v(:),vp(:,:)
!  
!  call readMsh('grid.msh',50,verbose=.true.)
!  
!  allocate(v(nFacet+100))
!  allocate(vp(nNode,DIMS))
!  call addWrite(v(1:nFacet),binding=RST_BIND_FACET)
!  call addWrite(vp,binding=RST_BIND_NODE)
!  do i=1,nFacet
!    v(i)=norm2(Facet(i)%PC)
!  end do
!  do i=1,nNode
!    vp(i,:)=Node(i)%Pos(:)
!  end do
!  
!  call writeRst('rst.msh',55)
end program
