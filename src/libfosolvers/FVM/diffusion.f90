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
  
contains
  
  !> find diffusion due to gradient of vector s on polyFvGrid
  !> \f[ \int_A D_i \nabla \mathbf{s} \cdot \hat{n} dA \f]
  subroutine findDiffPolyVect(grid,s,gradS,d,diff)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< state variables
    double precision,intent(in)::gradS(:,:,:) !< gradient of s
    double precision,intent(in)::d(:,:) !< diffusivities
    double precision,allocatable,intent(inout)::diff(:,:) !< diffusion output
    double precision::sf(DIMS),dPF,Afs,flow(size(s,1)),fPF,dF(size(s,1)),gradSF(DIMS,size(s,1))
    
    call grid%up()
    if(.not.(allocated(diff)))then
      allocate(diff(size(s,1),grid%nC))
    end if
    diff(:,:)=0d0
    !$omp parallel do default(shared)&
    !$omp& private(m,n,fPF,dF,gradSF,sf,dPF,Afs,flow,j)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
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
      Afs=grid%aP(i)/dot_product(sf,grid%normP(:,i))
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        flow(:)=dF(:)*Afs*((s(:,n)-s(:,m))/dPF&
        &                  -matmul(transpose(gradSF(:,:)),&
        &                          sf-grid%normP(:,i)*dot_product(sf,grid%normP(:,i))))
        do j=1,size(flow)
          !$omp atomic
          diff(j,m)=diff(j,m)+flow(j)
          !$omp atomic
          diff(j,n)=diff(j,n)-flow(j)
        end do
      else ! boundary pairs
        flow(:)=dF(:)*Afs*(s(:,n)-s(:,m))/(2d0*dPF)
        do j=1,size(flow)
          !$omp atomic
          diff(j,m)=diff(j,m)+flow(j)
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine
  
  !> find diffusion due to gradient of scalar s on polyFvGrid
  !> \f[ \int_A D \nabla \mathbf{s} \cdot \hat{n} dA \f]
  subroutine findDiffPolyScal(grid,s,gradS,d,diff)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< state variable
    double precision,allocatable::gradS(:,:) !< gradient of s
    double precision,intent(in)::d(:) !< diffusivity
    double precision,allocatable,intent(inout)::diff(:) !< diffusion output
    double precision::sf(DIMS),dPF,Afs,flow,fPF,dF,gradSF(DIMS)
    
    call grid%up()
    if(.not.(allocated(diff)))then
      allocate(diff(grid%nC))
    end if
    diff(:)=0d0
    !$omp parallel do default(shared)&
    !$omp& private(m,n,fPF,dF,gradSF,sf,dPF,Afs,flow)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(n<=grid%nC)then
        fPF=norm2(grid%pP(:,i)-grid%p(:,n))&
        &   /(norm2(grid%pP(:,i)-grid%p(:,m))+norm2(grid%pP(:,i)-grid%p(:,n)))
        dF=fPF*d(m)+(1d0-fPF)*d(n)
        gradSF(:)=fPF*gradS(:,m)+(1d0-fPF)*gradS(:,n)
      else
        dF=d(m)
        gradSF(:)=gradS(:,m)
      end if
      sf(:)=grid%p(:,n)-grid%p(:,m)
      dPF=norm2(sf)
      sf(:)=sf(:)/dPF
      Afs=grid%aP(i)/dot_product(sf,grid%normP(:,i))
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        flow=dF*Afs*((s(n)-s(m))/dPF&
        &            -dot_product(gradSF(:),sf-grid%normP(:,i)*dot_product(sf,grid%normP(:,i))))
        !$omp atomic
        diff(m)=diff(m)+flow
        !$omp atomic
        diff(n)=diff(n)-flow
      else ! boundary pairs
        flow=dF*Afs*(s(n)-s(m))/(2d0*dPF)
        !$omp atomic
        diff(m)=diff(m)+flow
      end if
    end do
    !$omp end parallel do
  end subroutine
  
end module
