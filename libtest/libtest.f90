program libtest
  use moduleGrid
  double precision,allocatable::value(:)
  double precision g(3)
  call readmsh()
  call sortEle()
  call updateFacetPara()
  call updateElePara()
  allocate(value(nNode))
  do i=1,nNode
    value(i)=10d0*Node(i)%Pos(1)+5d0*Node(i)%Pos(2)-5d0*Node(i)%Pos(3)
  end do
  do i=1,nNode
    call findNodeGradScal(i,value,g)
    write(*,*),i,g
  end do
end program
