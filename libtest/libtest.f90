program libtest
  use moduleGrid
  call readmsh()
  do i=1,nPoint
    write(*,*),i,Point(i)%findPos()
  end do
end program

