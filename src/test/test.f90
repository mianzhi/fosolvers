!----------------------------------------------------------------------------- best with 100 columns

program test
  use modOtGrid
  
  type(otGrid)::grid
  call grid%init([0d0,0d0,0d0],1d0)
  call grid%clear()
end program
