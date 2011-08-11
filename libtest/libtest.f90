program libtest
  use moduleGrid
  call readmsh()
  call sortEle()
  do i=1,nEle
    write(*,*),i,Ele(i)%findPC()
  end do
end program
