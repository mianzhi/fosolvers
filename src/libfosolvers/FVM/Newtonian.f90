!----------------------------------------------------------------------------- best with 100 columns

!> Newtonian fluid constitutive relations for FVM
module modNewtonian
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding viscose force
  interface findViscForce
    module procedure::findViscForcePoly
  end interface
  public::findViscForce
  
contains
  
  !> find Newtonian viscous force on polyFvGrid
  !> \f[ \int_A \hat{n} \cdot \mathbf{\tau} dA \f]
  subroutine findViscForcePoly(grid,u,visc,dRhou)
    use modPolyFvGrid
    use modDiffusion
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::u(:,:) !< state variables
    double precision,intent(in)::visc(:) !< viscosity
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision,allocatable::tmpVisc(:,:)
    
    call grid%up()
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    allocate(tmpVisc(DIMS,size(visc)))
    forall(i=1:size(visc))
      tmpVisc(:,i)=visc(i)
    end forall
    call findDiff(grid,u,tmpVisc,dRhou)
    ! TODO complete the rest part of div(tau)
    deallocate(tmpVisc)
  end subroutine
  
end module
