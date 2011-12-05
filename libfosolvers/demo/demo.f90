!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  
  allocate(Node(3))
  Node(1)%Pos=[0d0,0d0,0d0]
  Node(2)%Pos=[1d0,1d0,0d0]
  Node(3)%Pos=[0d0,1d0,0d0]
  allocate(Facet(1))
  call Facet(1)%specify(TRI_TYPE)
  Facet(1)%NodeInd=[1,2,3]
  call Facet(1)%updateNorm()
  write(*,*),Facet(1)%Norm
end program
