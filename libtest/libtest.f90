program libtest
  use moduleGrid
  call readmsh()
  call updateElePara()
  do i=1,nEle
    write(*,*),i,Ele(i)%SurfNorm(1,:)
  end do
end program
