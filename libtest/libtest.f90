program libtest
  use moduleGrid
  !double precision,allocatable::value(:,:)
  !double precision g(3,3)
  call readmsh()
  call sortEle()
  call updateFacetPara()
  call updateElePara()
  !allocate(value(nEle,3))
  !do i=1,nEle
  !  value(i,:)=[10d0*Ele(i)%PC(1)+5d0*Ele(i)%PC(2),5d0*Ele(i)%PC(2),-5d0*Ele(i)%PC(3)]
  !end do
  do i=1,nFacet
    write(*,*),i,Facet(i)%NodeInd(:)
  !  call findEleGradVect(i,value,g)
  !  write(*,*),i
  !  write(*,*),g(1,:)
  !  write(*,*),g(2,:)
  !  write(*,*),g(3,:)
  end do
end program
