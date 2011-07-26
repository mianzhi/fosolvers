!author: Wang, Mianzhi
!subroutines:
! writedata()

subroutine writedata()
!write data file rst.msh for post-process
  use globalvar
  integer ifile,ncell,nvertice
  character(200) tempstring
  ncell=size(element,1) !number of cells
  nvertice=size(vertice,1) !number of vertices
  ifile=10
  open(ifile,file='rst.msh',status='replace')
  !write header
  write(ifile,'(a)'),'$MeshFormat'
  write(ifile,'(a)'),'2.2 0 8'
  write(ifile,'(a)'),'$EndMeshFormat'
  !write nodes
  write(ifile,'(a)'),'$Nodes'
  write(tempstring,*),nvertice
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nvertice
    write(tempstring,'(i8,1x,g14.7,1x,g14.7,1x,g14.7)'),i,px(i),py(i),pz(i)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndNodes'
  !write elements
  write(ifile,'(a)'),'$Elements'
  write(tempstring,*),ncell
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,ncell
    write(tempstring,'(i8,1x,i1,1x,i1,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)'),i,5,0,&
                       &element(i,5),element(i,1),element(i,2),element(i,6),&
                       &element(i,8),element(i,4),element(i,3),element(i,7)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElements'
  !write velocities
  write(ifile,'(a)'),'$NodeData'
  write(ifile,'(i1)'),1 !1 string tag
  write(ifile,'(a)'),'"velocity"' !name of the data-set
  write(ifile,'(i1)'),1 !1 real tag
  write(ifile,'(f6.4)'),0.0 !time
  write(ifile,'(i1)'),3 !3 integer tags
  write(ifile,'(i1)'),0 !time step index
  write(ifile,'(i1)'),3 !3 components vector
  write(tempstring,*),nvertice
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nvertice
    write(tempstring,'(i8,1x,g14.7,1x,g14.7,1x,g14.7)'),i,vu(i),vv(i),vw(i)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndNodeData'
  !write pressure
  write(ifile,'(a)'),'$ElementData'
  write(ifile,'(i1)'),1 !1 string tag
  write(ifile,'(a)'),'"pressure"' !name of the data-set
  write(ifile,'(i1)'),1 !1 real tag
  write(ifile,'(f6.4)'),0.0 !time
  write(ifile,'(i1)'),3 !3 integer tags
  write(ifile,'(i1)'),0 !time step index
  write(ifile,'(i1)'),1 !1 component scaler
  write(tempstring,*),ncell
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,ncell
    write(tempstring,'(i8,1x,g14.7)'),i,Pr(i)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElementData'
  !write density
  write(ifile,'(a)'),'$ElementData'
  write(ifile,'(i1)'),1 !1 string tag
  write(ifile,'(a)'),'"density"' !name of the data-set
  write(ifile,'(i1)'),1 !1 real tag
  write(ifile,'(f6.4)'),0.0 !time
  write(ifile,'(i1)'),3 !3 integer tags
  write(ifile,'(i1)'),0 !time step index
  write(ifile,'(i1)'),1 !1 component scaler
  write(tempstring,*),ncell
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,ncell
    write(tempstring,'(i8,1x,g14.7)'),i,rr(i)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElementData'
  write(ifile,'(a)'),'$EndElementData'
  !write temperature
  write(ifile,'(a)'),'$ElementData'
  write(ifile,'(i1)'),1 !1 string tag
  write(ifile,'(a)'),'"temperature"' !name of the data-set
  write(ifile,'(i1)'),1 !1 real tag
  write(ifile,'(f6.4)'),0.0 !time
  write(ifile,'(i1)'),3 !3 integer tags
  write(ifile,'(i1)'),0 !time step index
  write(ifile,'(i1)'),1 !1 component scaler
  write(tempstring,*),ncell
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,ncell
    write(tempstring,'(i8,1x,g14.7)'),i,KT(i)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElementData'
  close(ifile)
end subroutine
