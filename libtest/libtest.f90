program libtest
  use moduleGrid
  call readmsh()
  call updateFacetPara()
  call updateElePara()
  do i=1,nFacet
    write(*,*),i,Facet(i)%PC(:),Facet(i)%Area,Facet(i)%Norm(:)
  end do
end program
