program libtest
  use moduleGrid
  double precision,allocatable::value(:)
  double precision g(3)
  call readmsh()
  call updateFacetPara()
  !call sortEle()
  call updateElePara()
  allocate(value(nEle))
  do i=1,nEle
    value(i)=10d0*Ele(i)%PC(1)+5d0*Ele(i)%PC(2)-5d0*Ele(i)%PC(3)
  end do
  do i=1,nEle
    call findEleGradScal(i,value,g)
    write(*,*),i,g(:)
  end do
end program
