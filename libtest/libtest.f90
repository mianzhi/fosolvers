program libtest
  use moduleGrid
  double precision rst1(3),rst2(3)
  allocate(Node(8))
  Node(1)%Pos=[0,0,0]
  Node(2)%Pos=[1,0,0]
  Node(3)%Pos=[1,1,0]
  Node(4)%Pos=[0,1,0]
  Node(5)%Pos=[0,0,1]
  Node(6)%Pos=[1,0,1]
  Node(7)%Pos=[1,1,1]
  Node(8)%Pos=[0,1,1]
  allocate(Tet(1))
  Tet(1)%NodeInd(:)=[1,2,5,8]
  Tet(1)%GeoEnti=10
  allocate(Hex(1))
  Hex(1)%NodeInd(:)=[1,2,3,4,5,6,7,8]
  Hex(1)%GeoEnti=11
  allocate(Ele(2))
  Ele(1)%ShapeType=5
  Ele(1)%ShapeInd=1
  Ele(2)%ShapeType=4
  Ele(2)%ShapeInd=1
  call Ele(1)%findPC(rst1)
  call Ele(2)%findPC(rst2)
  write(*,*),rst1,rst2
  write(*,*),Ele(1)%findVol(),Ele(2)%findVol()
  write(*,*),Ele(1)%findGeoEnti(),Ele(2)%findGeoEnti()
end program

