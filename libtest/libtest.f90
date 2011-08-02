program libtest
  use moduleGrid
  double precision range(2,3)
  call readmsh()
  call findBound(range)
  write(*,*),range(1,:)
  write(*,*),range(2,:)
end program

