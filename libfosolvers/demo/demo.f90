!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  use moduleFileIO
  
  double precision,target,allocatable::v(:)
  
  call readmsh('grid.msh',50,verbose=.true.)
  call initWriteEnv(0,0,0,1,0,0,0,0,0)
  
  allocate(v(nFacet))
  rstFacetScal(1)%ptr=>v
  do i=1,nFacet
    v(i)=norm2(Facet(i)%PC)
  end do
  
  call writerst('rst.msh',55)
end program
