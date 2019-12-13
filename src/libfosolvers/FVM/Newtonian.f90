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
  subroutine findViscForcePoly(grid,u,gradU,visc,dRhou)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::u(:,:) !< state variables
    double precision,intent(in)::gradU(:,:,:) !< gradient of u
    double precision,intent(in)::visc(:) !< viscosity
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision::sf(DIMS),tf(DIMS),rf(DIMS),dPF,Afs,flow(DIMS),fPF,viscF,&
    &                 gradUF(DIMS,DIMS)
    
    call grid%up()
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    !$omp parallel do default(shared)&
    !$omp& private(m,n,fPF,viscF,gradUF,sf,dPF,k,l,tf,rf,Afs,flow)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(n<=grid%nC)then
        fPF=norm2(grid%pP(:,i)-grid%p(:,n))&
        &   /(norm2(grid%pP(:,i)-grid%p(:,m))+norm2(grid%pP(:,i)-grid%p(:,n)))
        viscF=fPF*visc(m)+(1d0-fPF)*visc(n)
        gradUF(:,:)=fPF*gradU(:,:,m)+(1d0-fPF)*gradU(:,:,n)
      else
        viscF=visc(m)
        gradUF(:,:)=gradU(:,:,m)
      end if
      sf(:)=grid%p(:,n)-grid%p(:,m)
      dPF=norm2(sf)
      sf(:)=sf(:)/dPF
      k=maxloc(abs(grid%normP(:,i)),dim=1)
      l=merge(1,k+1,k==3)
      tf(:)=0d0
      tf(l)=grid%normP(k,i)
      tf(k)=-grid%normP(l,i)
      tf(:)=tf(:)/norm2(tf)
      rf(1)=grid%normP(2,i)*tf(3)-grid%normP(3,i)*tf(2)
      rf(2)=-grid%normP(1,i)*tf(3)+grid%normP(3,i)*tf(1)
      rf(3)=grid%normP(1,i)*tf(2)-grid%normP(2,i)*tf(1)
      Afs=grid%aP(i)/dot_product(sf,grid%normP(:,i))
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        flow(:)=viscF*Afs*((u(:,n)-u(:,m))/dPF&
        &                  -matmul(transpose(gradUF(:,:)),&
        &                          dot_product(sf,tf)*tf+dot_product(sf,rf)*rf))
      else ! boundary pairs
        flow(:)=viscF*Afs*(u(:,n)-u(:,m))/(2d0*dPF)
      end if
      flow(:)=flow(:)+viscF*grid%aP(i)*[dot_product(gradUF(1,:),grid%normP(:,i)),&
      &                                 dot_product(gradUF(2,:),grid%normP(:,i)),&
      &                                 dot_product(gradUF(3,:),grid%normP(:,i))]&
      &      -2d0/3d0*viscF*grid%aP(i)*(gradUF(1,1)+gradUF(2,2)+gradUF(3,3))*grid%normP(:,i)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        do j=1,DIMS
          !$omp atomic
          dRhou(j,m)=dRhou(j,m)+flow(j)
          !$omp atomic
          dRhou(j,n)=dRhou(j,n)-flow(j)
        end do
      else ! boundary pairs
        do j=1,DIMS
          !$omp atomic
          dRhou(j,m)=dRhou(j,m)+flow(j)
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine
  
end module
