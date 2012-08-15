!----------------------------------------------------------------------------- best with 100 columns

!> grid inspection
module moduleGridInspection
  use moduleGrid
  private
  
  ! public procedures
  public::findBoundBox
  
contains
  
  !> find the bounding box of the grid
  pure function findBoundBox(grid)
    type(typeGrid),intent(in)::grid !< the grid
    double precision findBoundBox(DIMS,2) !< the two points defining the bounding box
    
    findBoundBox(:,1)=minval(grid%NodePos,2)
    findBoundBox(:,2)=maxval(grid%NodePos,2)
  end function
  
end module
