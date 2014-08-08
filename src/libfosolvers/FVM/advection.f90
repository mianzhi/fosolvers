!----------------------------------------------------------------------------- best with 100 columns

!> advection for FVM
module modAdvection
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding advection
  interface findAdv
    module procedure::findAdvPolyVect
    module procedure::findAdvPolyScal
  end interface
  public::findAdv
  
contains
  
  !> find advection rate of vector s by velocity u on polyFvGrid
  !> \f[ \int_A \mathbf{s} (\mathbf{u} \cdot \hat{n}) dA \f]
  !> NOTE: s, u are reconstructed data
  subroutine findAdvPolyVect(grid,s,u,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< reconstructed input data
    double precision,intent(in)::u(:,:) !< reconstructed velocity field
    double precision,allocatable,intent(inout)::adv(:,:) !< advection output
    double precision::F(size(s,1))
    
    call grid%up()
    if(.not.(allocated(adv)))then
      allocate(adv(size(s,1),grid%nC))
    end if
    adv(:,:)=0d0
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        F(:)=-grid%aP(i)*dot_product(u(:,i),grid%normP(:,i))*s(:,i)
        adv(:,m)=adv(:,m)+F(:)
        adv(:,n)=adv(:,n)-F(:)
      end if
    end do
  end subroutine
  
  !> find advection rate of scalar s by velocity u on polyFvGrid
  !> \f[ \int_A s (\mathbf{u} \cdot \hat{n}) dA \f]
  !> NOTE: s, u are reconstructed data
  subroutine findAdvPolyScal(grid,s,u,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< reconstructed input data
    double precision,intent(in)::u(:,:) !< reconstructed velocity field
    double precision,allocatable,intent(inout)::adv(:) !< advection output
    double precision,allocatable::sv(:,:),advv(:,:)
    
    allocate(sv(1,size(s)))
    if(allocated(adv))then
      allocate(advv(1,size(adv)))
    end if
    sv(1,:)=s(:)
    call findAdvPolyVect(grid,sv,u,advv)
    if(.not.allocated(adv))then
      allocate(adv(size(advv,2)),source=advv(1,:))!FIXME:remove work-around
    else
      adv(:)=advv(1,:)
    end if
    deallocate(sv)
    deallocate(advv)
  end subroutine
  
end module
