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
    double precision::flow(DIMS),pf,vMP(DIMS),vNP(DIMS),vM1P(DIMS),vN1P(DIMS),&
    &                 vMM1(DIMS),vNN1(DIMS),eps
    
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
          vNP(:)=grid%pP(:,i)-grid%p(:,n)
          vM1P(:)=grid%normP(:,i)*dot_product(vMP(:),grid%normP(:,i))
          vN1P(:)=grid%normP(:,i)*dot_product(vNP(:),grid%normP(:,i))
          vMM1(:)=vMP(:)-vM1P(:)
          vNN1(:)=vNP(:)-vN1P(:)
          eps=norm2(vN1P)/(norm2(vM1P)+norm2(vN1P))
          pf=eps*(p(m)+dot_product(gradP(:,m),vMM1(:)))&
          &  +(1d0-eps)*(p(n)+dot_product(gradP(:,n),vNN1(:)))
          flow(:)=-pf*grid%aP(i)*grid%normP(:,i)
          dRhou(:,m)=dRhou(:,m)+flow(:)
          dRhou(:,n)=dRhou(:,n)-flow(:)
        else ! boundary pairs
          vMP(:)=grid%pP(:,i)-grid%p(:,m)
          vM1P(:)=grid%normP(:,i)*dot_product(vMP(:),grid%normP(:,i))
          pf=0.5d0*(p(m)+p(n))+dot_product(gradP(:,m),vMP(:)-vM1P(:))
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
