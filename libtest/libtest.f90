program libtest
  use moduleGrid
  call readmsh()
  do i=1,nFacet
    write(*,*),i,Facet(i)%findGeoEnti(),Facet(i)%ShapeType,Facet(i)%findArea()
  end do
  write(*,*),size(Facet)
end program

