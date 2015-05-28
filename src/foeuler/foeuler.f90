!----------------------------------------------------------------------------- best with 100 columns

!> environment for Euler problem
module modEuler
  use modPolyFvGrid
  use iso_c_binding
  
  public
  
  type(polyFvGrid)::grid
end module

!> Euler solver
program foeuler
  use modEuler
end program
