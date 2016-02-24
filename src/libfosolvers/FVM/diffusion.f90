!----------------------------------------------------------------------------- best with 100 columns

!> diffusion for FVM
module modDiffusion
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding diffusion
  interface findDiff
    module procedure::findDiffPolyVect
    module procedure::findDiffPolyScal
  end interface
  public::findDiff
  
  !> generic finding viscose force
  interface findViscForce
    module procedure::findViscForcePoly
  end interface
  public::findViscForce
  
contains
  
  !> find diffusion due to gradient of vector s on polyFvGrid
  !> \f[ \int_A D_i \nabla \mathbf{s} \cdot \hat{n} dA \f]
  subroutine findDiffPolyVect(grid,s,d,diff)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< state variables
    double precision,intent(in)::d(:,:) !< diffusivities
    double precision,allocatable,intent(inout)::diff(:,:) !< diffusion output
    double precision::sf(DIMS),tf(DIMS),rf(DIMS),dPF,Afs,flow(size(s,1)),fPF,dF(size(s,1)),&
    &                 gradSF(DIMS,size(s,1))
    double precision,allocatable::gradS(:,:,:)
    
    call grid%up()
    if(.not.(allocated(diff)))then
      allocate(diff(size(s,1),grid%nC))
    end if
    diff(:,:)=0d0
    allocate(gradS(DIMS,size(s,1),grid%nC))
    call findGrad(grid,s(:,1:grid%nC),gradS)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(s,2).and.n<=size(s,2))then
        if(n<=grid%nC)then
          fPF=norm2(grid%pP(:,i)-grid%p(:,n))&
          &   /(norm2(grid%pP(:,i)-grid%p(:,m))+norm2(grid%pP(:,i)-grid%p(:,n)))
          dF(:)=fPF*d(:,m)+(1d0-fPF)*d(:,n)
          gradSF(:,:)=fPF*gradS(:,:,m)+(1d0-fPF)*gradS(:,:,n)
        else
          dF(:)=d(:,m)
          gradSF(:,:)=gradS(:,:,m)
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
          flow(:)=dF(:)*Afs*((s(:,n)-s(:,m))/dPF&
          &                  -matmul(transpose(gradSF(:,:)),&
          &                          dot_product(sf,tf)*tf+dot_product(sf,rf)*rf))
          diff(:,m)=diff(:,m)+flow(:)
          diff(:,n)=diff(:,n)-flow(:)
        else ! boundary pairs
          flow(:)=dF(:)*Afs*(s(:,n)-s(:,m))/(2d0*dPF)
          diff(:,m)=diff(:,m)+flow(:)
        end if
      end if
    end do
    deallocate(gradS)
  end subroutine
  
  !> find diffusion due to gradient of scalar s on polyFvGrid
  !> \f[ \int_A D \nabla \mathbf{s} \cdot \hat{n} dA \f]
  subroutine findDiffPolyScal(grid,s,d,diff)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< state variable
    double precision,intent(in)::d(:) !< diffusivity
    double precision,allocatable,intent(inout)::diff(:) !< diffusion output
    double precision,allocatable::sv(:,:),dv(:,:),diffv(:,:)
    
    allocate(sv(1,size(s)))
    allocate(dv(1,size(s)))
    if(allocated(diff))then
      allocate(diffv(1,size(diff)))
    end if
    sv(1,:)=s(:)
    dv(1,:)=d(:)
    call findDiffPolyVect(grid,sv,dv,diffv)
    if(.not.allocated(diff))then
      allocate(diff(size(diffv,2)),source=diffv(1,:))!FIXME:remove work-around
    else
      diff(:)=diffv(1,:)
    end if
    deallocate(sv)
    deallocate(dv)
    deallocate(diffv)
  end subroutine
  
  !> find Newtonian viscous force on polyFvGrid
  !> \f[ \int_A \hat{n} \cdot \mathbf{\tao} dA \f]
  subroutine findViscForcePoly(grid,u,visc,dRhou)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::u(:,:) !< state variables
    double precision,intent(in)::visc(:) !< viscosity
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision::fPF,taoF(DIMS,DIMS),flow(DIMS),utP(DIMS),utF(DIMS)
    double precision,allocatable::tao(:,:,:)
    
    call grid%up()
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    allocate(tao(DIMS,DIMS,grid%nC))
    call findGrad(grid,u(:,1:grid%nC),tao)
    forall(i=1:grid%nC)
      tao(:,:,i)=visc(i)*((tao(:,:,i)+transpose(tao(:,:,i)))&
      &                   -2d0/3d0*reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[DIMS,DIMS]))
    end forall
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(u,2).and.n<=size(u,2))then
        if(n<=grid%nC)then ! internal pairs
          fPF=norm2(grid%pP(:,i)-grid%p(:,n))&
          &   /(norm2(grid%pP(:,i)-grid%p(:,m))+norm2(grid%pP(:,i)-grid%p(:,n)))
          taoF(:,:)=fPF*tao(:,:,m)+(1d0-fPF)*tao(:,:,n)
          flow(:)=grid%aP(i)*matmul(grid%normP(:,i),taoF(:,:))
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
    deallocate(tao)
  end subroutine
  
end module
