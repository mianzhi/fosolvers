program libtest
  use moduleGrid
  call readmsh()
  call updateNeib()
  do i=1,nEle
    do j=1,6
      k=Ele(i)%Neib(j)
      write(*,*),i,j,k
    end do
  end do
end program
