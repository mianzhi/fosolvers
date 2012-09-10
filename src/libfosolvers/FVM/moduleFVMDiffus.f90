!----------------------------------------------------------------------------- best with 100 columns

!> FVM diffusion schemes
module moduleFVMDiffus
  private
  
  !> find diffusion through interface
  interface findDiffus
    module procedure findDiffusORTHScal
    module procedure findDiffusORTHVect
    module procedure findDiffusSDScal
    module procedure findDiffusSDVect
  end interface
  public findDiffus
  
contains
  
  !> find diffusion driven by vector v through interface using orthogonal scheme
  !> \f[ \int_{intf} \Gamma \nabla \mathbf{v} \cdot \hat{n} dA \f]
  function findDiffusORTHVect(gamm,gammbind,v,grid)
    use moduleGrid
    use moduleInterpolation
    double precision,intent(in)::gamm(:) !< the diffusivity
    integer,intent(in)::gammbind !< bind diffusivity with node/block/interface
    double precision,intent(in)::v(:,:) !< the vector which drives the diffusion
    type(typeGrid),intent(inout)::grid !< the grid on which v is defined
    double precision findDiffusORTHVect(size(v,1),grid%nBlock) !< increment due to diffusion
    double precision flowRate(size(v,1))
    double precision,allocatable::gammIntf(:)
    
    call grid%updateIntfPos()
    call grid%updateIntfArea()
    call grid%updateBlockPos()
    findDiffusORTHVect(:,:)=0d0
    allocate(gammIntf(grid%nIntf))
    select case(gammbind)
    case(BIND_NODE)
      gammIntf(:)=itplNode2Intf(gamm,grid)
    case(BIND_BLOCK)
      gammIntf(:)=itplBlock2Intf(gamm,grid)
    case(BIND_INTF)
      gammIntf(:)=gamm(:)
    case default
    end select
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      flowRate(:)=gammIntf(i)*(v(:,n)-v(:,m))*grid%IntfArea(i)&
      &           /norm2(grid%BlockPos(:,n)-grid%BlockPos(:,m))
      findDiffusORTHVect(:,m)=findDiffusORTHVect(:,m)+flowRate(:)
      findDiffusORTHVect(:,n)=findDiffusORTHVect(:,n)-flowRate(:)
    end do
    deallocate(gammIntf)
  end function
  
  !> find diffusion driven by scalar v through interface using orthogonal scheme
  !> \f[ \int_{intf} \Gamma \nabla v \cdot \hat{n} dA \f]
  function findDiffusORTHScal(gamm,gammbind,v,grid)
    use moduleGrid
    double precision,intent(in)::gamm(:) !< the diffusivity
    integer,intent(in)::gammbind !< bind diffusivity with node/block/interface
    double precision,intent(in)::v(:) !< the scalar which drives the diffusion
    type(typeGrid),intent(inout)::grid !< the grid on which v is defined
    double precision findDiffusORTHScal(grid%nBlock) !< increment due to diffusion
    double precision vv(1,size(v)),vrst(1,grid%nBlock)
    
    vv(1,:)=v(:)
    vrst=findDiffusORTHVect(gamm,gammbind,vv,grid)
    findDiffusORTHScal(:)=vrst(1,:)
  end function
  
  !> find diffusion driven by vector v through interface using surface decomposition scheme
  !> \f[ \int_{intf} \Gamma \nabla \mathbf{v} \cdot \hat{n} dA \f]
  function findDiffusSDVect(gamm,gammbind,v,grid,grad)
    use moduleGrid
    use moduleInterpolation
    double precision,intent(in)::gamm(:) !< the diffusivity
    integer,intent(in)::gammbind !< bind diffusivity with node/block/interface
    double precision,intent(in)::v(:,:) !< the vector which drives the diffusion
    type(typeGrid),intent(inout)::grid !< the grid on which v is defined
    double precision,intent(in)::grad(DIMS,size(v,1),size(v,2)) !< the gradient of v
    double precision findDiffusSDVect(size(v,1),grid%nBlock) !< increment due to diffusion
    double precision flowRate(size(v,1)),dPF,sf(DIMS),tf(DIMS),rf(DIMS),Afs
    double precision,allocatable::gradIntf(:,:,:)
    double precision,allocatable::gammIntf(:)
    
    call grid%updateIntfPos()
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    call grid%updateBlockPos()
    findDiffusSDVect(:,:)=0d0
    allocate(gradIntf(size(grad,1),size(grad,2),grid%nIntf))
    gradIntf=reshape(itplBlock2Intf(reshape(grad,[size(grad,1)*size(grad,2),size(grad,3)]),grid),&
    &                [size(gradIntf,1),size(gradIntf,2),size(gradIntf,3)])
    allocate(gammIntf(grid%nIntf))
    select case(gammbind)
    case(BIND_NODE)
      gammIntf(:)=itplNode2Intf(gamm,grid)
    case(BIND_BLOCK)
      gammIntf(:)=itplBlock2Intf(gamm,grid)
    case(BIND_INTF)
      gammIntf(:)=gamm(:)
    case default
    end select
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      sf(:)=grid%BlockPos(:,n)-grid%BlockPos(:,m)
      dPF=norm2(sf)
      sf(:)=sf(:)/dPF
      tf(:)=grid%NodePos(:,grid%Intf(i)%iNode(1))-grid%NodePos(:,grid%Intf(i)%iNode(2))
      tf(:)=tf(:)/norm2(tf)
      rf(1)=grid%IntfNorm(2,i)*tf(3)-grid%IntfNorm(3,i)*tf(2)
      rf(2)=-grid%IntfNorm(1,i)*tf(3)+grid%IntfNorm(3,i)*tf(1)
      rf(3)=grid%IntfNorm(1,i)*tf(2)-grid%IntfNorm(2,i)*tf(1)
      Afs=grid%IntfArea(i)/dot_product(sf,grid%IntfNorm(:,i))
      flowRate(:)=gammIntf(i)*(Afs*((v(:,n)-v(:,m))/dPF&
      &                        -matmul(transpose(gradIntf(:,:,i)),&
      &                                dot_product(sf,tf)*tf+dot_product(sf,rf)*rf)))
      findDiffusSDVect(:,m)=findDiffusSDVect(:,m)+flowRate(:)
      findDiffusSDVect(:,n)=findDiffusSDVect(:,n)-flowRate(:)
    end do
    deallocate(gradIntf)
    deallocate(gammIntf)
  end function
  
  !> find diffusion driven by scalar v through interface using surface decomposition scheme
  !> \f[ \int_{intf} \Gamma \nabla v \cdot \hat{n} dA \f]
  function findDiffusSDScal(gamm,gammbind,v,grid,grad)
    use moduleGrid
    double precision,intent(in)::gamm(:) !< the diffusivity
    integer,intent(in)::gammbind !< bind diffusivity with node/block/interface
    double precision,intent(in)::v(:) !< the scalar which drives the diffusion
    type(typeGrid),intent(inout)::grid !< the grid on which v is defined
    double precision,intent(in)::grad(DIMS,size(v)) !< the gradient of v
    double precision findDiffusSDScal(grid%nBlock) !< increment due to diffusion
    double precision vv(1,size(v)),vgrad(DIMS,1,size(v)),vrst(1,grid%nBlock)
    
    vv(1,:)=v(:)
    vgrad(:,1,:)=grad(:,:)
    vrst=findDiffusSDVect(gamm,gammbind,vv,grid,grad)
    findDiffusSDScal(:)=vrst(1,:)
  end function
  
end module
