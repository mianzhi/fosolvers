program libtest
  use moduleGrid
  call readmsh()
  write(*,*),find3PNorm(Node(1)%Pos(:),Node(2)%Pos(:),Node(3)%Pos(:))
end program

