!----------------------------------------------------------------------------- best with 100 columns

!> grid operations
module moduleGridOperation
  use moduleGrid
  private
  
  ! public procedures
  public::splitGridPrt
  
contains
  
  !> split the original grid into resulting grids with respect to partition
  subroutine splitGridPrt(origin,rst)
    type(typeGrid),intent(in)::origin !< the original grid
    type(typeGrid),intent(inout),allocatable::rst(:) !< the resulting grids
    
    if(allocated(rst)) deallocate(rst)
    allocate(rst(origin%nPrt))
  end subroutine
  
end module
