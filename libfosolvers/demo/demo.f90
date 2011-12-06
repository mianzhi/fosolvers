!----------------------------------------------------------------------------- best with 100 columns

program demo
  use moduleGrid
  use moduleFileIO
  
  call readmsh('grid.msh',50,verbose=.true.)
end program
