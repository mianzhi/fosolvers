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
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::u(:,:) !< state variables
    double precision,intent(in)::visc(:) !< viscosity
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision::fPF,tauF(DIMS,DIMS),flow(DIMS),utP(DIMS),utF(DIMS)
    double precision,allocatable::tau(:,:,:)
    
    call grid%up()
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    allocate(tau(DIMS,DIMS,grid%nC))
    call findGrad(grid,u,tau)
    forall(i=1:grid%nC)
      tau(:,:,i)=visc(i)*((tau(:,:,i)+transpose(tau(:,:,i)))&
      &                   -2d0/3d0*(tau(1,1,i)+tau(2,2,i)+tau(3,3,i))&
      &                    *reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[DIMS,DIMS]))
    end forall
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(u,2).and.n<=size(u,2))then
        if(n<=grid%nC)then ! internal pairs
          fPF=norm2(grid%pP(:,i)-grid%p(:,n))&
          &   /(norm2(grid%pP(:,i)-grid%p(:,m))+norm2(grid%pP(:,i)-grid%p(:,n)))
          tauF(:,:)=fPF*tau(:,:,m)+(1d0-fPF)*tau(:,:,n)
          flow(:)=grid%aP(i)*matmul(grid%normP(:,i),tauF(:,:))
          dRhou(:,m)=dRhou(:,m)+flow(:)
          dRhou(:,n)=dRhou(:,n)-flow(:)
        else ! boundary pairs
          utP(:)=u(:,m)-grid%normP(:,i)*dot_product(u(:,m),grid%normP(:,i))
          utF(:)=u(:,n)-grid%normP(:,i)*dot_product(u(:,n),grid%normP(:,i))
          flow(:)=grid%aP(i)*visc(m)*(utF(:)-utP(:))&
          &       /(2d0*dot_product(grid%p(:,n)-grid%p(:,m),grid%normP(:,i)))
          dRhou(:,m)=dRhou(:,m)+flow(:)
        end if
      end if
    end do
    deallocate(tau)
  end subroutine
  
end module
