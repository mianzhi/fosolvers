!----------------------------------------------------------------------------- best with 100 columns

!> pressure force for FVM
module modPressure
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding pressure force
  interface findPresForce
    module procedure::findPresForcePoly
    module procedure::findPresForcePolyNoCorrect
  end interface
  public::findPresForce
  
contains
  
  !> find pressure force on polyFvGrid
  !> \f[ \int_A p \hat{n} dA \f]
  subroutine findPresForcePoly(grid,p,gradP,dRhou)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::p(:) !< pressure
    double precision,intent(in)::gradP(:,:) !< gradient of p
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision::flow(DIMS),pf,gradPf(DIMS),vIP(DIMS),vMP(DIMS),vPN(DIMS),eps
    
    call grid%up()
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(p).and.n<=size(p))then
        if(n<=grid%nC)then ! internal pairs
          vMP(:)=grid%pP(:,i)-grid%p(:,m)
          vPN(:)=grid%p(:,n)-grid%pP(:,i)
          eps=norm2(vPN)/(norm2(vMP)+norm2(vPN))
          pf=eps*p(m)+(1d0-eps)*p(n)
          gradPf(:)=eps*gradP(:,m)+(1d0-eps)*gradP(:,n)
          vIP(:)=grid%pP(:,i)-(eps*grid%p(:,m)+(1d0-eps)*grid%p(:,n))
          pf=pf+dot_product(gradPf(:),vIP(:))
          flow(:)=-pf*grid%aP(i)*grid%normP(:,i)
          dRhou(:,m)=dRhou(:,m)+flow(:)
          dRhou(:,n)=dRhou(:,n)-flow(:)
        else ! boundary pairs
          pf=0.5d0*(p(m)+p(n))
          vIP(:)=grid%pP(:,i)-grid%p(:,m)
          vIP(:)=vIP(:)-grid%normP(:,i)*dot_product(grid%normP(:,i),vIP(:))
          pf=pf+dot_product(gradP(:,m),vIP(:))
          flow(:)=-pf*grid%aP(i)*grid%normP(:,i)
          dRhou(:,m)=dRhou(:,m)+flow(:)
        end if
      end if
    end do
  end subroutine
  
  !> find pressure force on polyFvGrid without skewness correction by pressure gradient
  !> \f[ \int_A p \hat{n} dA \f]
  subroutine findPresForcePolyNoCorrect(grid,p,dRhou)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::p(:) !< pressure
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision,allocatable::gradP(:,:) !< gradient of p
    
    allocate(gradP(DIMS,grid%nC))
    gradP(:,:)=0d0
    call findPresForcePoly(grid,p,gradP,dRhou)
    deallocate(gradP)
  end subroutine
  
end module
