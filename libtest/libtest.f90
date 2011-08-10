program libtest
  use moduleGrid
  call readmsh()
  do i=1,nEle
    do j=1,Ele(i)%SurfNum
      write(*,*),i,j,Ele(i)%findSurfNorm(j)
    end do
  end do
end program
