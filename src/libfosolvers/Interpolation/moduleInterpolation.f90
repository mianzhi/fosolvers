!----------------------------------------------------------------------------- best with 100 columns

!> interpolation schemes
module moduleInterpolation
  private
  
  !> interpolate from block to interface
  interface itplBlock2Intf
    module procedure::itplBCDVect
    module procedure::itplBCDScal
    module procedure::itplBCDCVect
    module procedure::itplBCDCScal
  end interface
  public itplBlock2Intf
  
  !> interpolate from node to interface
  interface itplNode2Intf
    module procedure::itplNode2IntfVect
    module procedure::itplNode2IntfScal
  end interface
  public itplNode2Intf
  
contains
  
  !> interpolate vector v within grid from block to interface using block center direction scheme
  function itplBCDVect(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:) !< block data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBCDVect(:,:) !< interpolated data on interface
    double precision PF(DIMS),PS(DIMS),alpha
    
    call grid%updateIntfPos()
    call grid%updateBlockPos()
    call reallocArr(itplBCDVect,size(v,1),grid%nIntf)
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      PF(:)=grid%BlockPos(:,n)-grid%BlockPos(:,m)
      PS(:)=grid%IntfPos(:,i)-grid%BlockPos(:,m)
      alpha=dot_product(PS,PF)/dot_product(PF,PF)
      itplBCDVect(:,i)=(1d0-alpha)*v(:,m)+alpha*v(:,n)
    end do
  end function
  
  !> interpolate scalar v within grid from block to interface using block center direction scheme
  function itplBCDScal(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:) !< block data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBCDScal(:) !< interpolated data on interface
    double precision vv(1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vrst=itplBCDVect(vv,grid)
    call reallocArr(itplBCDScal,size(vrst,2))
    itplBCDScal(:)=vrst(1,:)
  end function
  
  !> interpolate vector v within grid from block to interface using BCD scheme corrected using grad
  function itplBCDCVect(v,grad,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:) !< block data to be interpolated
    double precision,intent(in)::grad(:,:,:) !< gradient of v
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBCDCVect(:,:) !< interpolated data on interface
    double precision dPS,dSF,alpha,Saux(DIMS)
    
    call grid%updateIntfPos()
    call grid%updateBlockPos()
    call reallocArr(itplBCDCVect,size(v,1),grid%nIntf)
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      dPS=norm2(grid%IntfPos(:,i)-grid%BlockPos(:,m))
      dSF=norm2(grid%BlockPos(:,n)-grid%IntfPos(:,i))
      alpha=dPS/(dPS+dSF)
      Saux(:)=grid%BlockPos(:,m)+alpha*(grid%BlockPos(:,n)-grid%BlockPos(:,m))
      itplBCDCVect(:,i)=(1d0-alpha)*v(:,m)+alpha*v(:,n)&
      &                 +matmul(transpose(grad(:,:,m)),grid%IntfPos(:,i)-Saux(:))
    end do
  end function
  
  !> interpolate scalar v within grid from block to interface using BCD scheme corrected using grad
  function itplBCDCScal(v,grad,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:) !< block data to be interpolated
    double precision,intent(in)::grad(:,:) !< gradient of v
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBCDCScal(:) !< interpolated data on interface
    double precision vv(1,size(v)),vgrad(DIMS,1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vgrad(:,1,:)=grad(:,:)
    vrst=itplBCDCVect(vv,vgrad,grid)
    call reallocArr(itplBCDCScal,size(vrst,2))
    itplBCDCScal(:)=vrst(1,:)
  end function
  
  !> interpolate vector v within grid from node to interface
  function itplNode2IntfVect(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2IntfVect(:,:) !< interpolated data on interface
    
    call grid%updateIntf()
    call reallocArr(itplNode2IntfVect,size(v,1),grid%nIntf)
    forall(i=1:grid%nIntf)
      itplNode2IntfVect(:,i)=sum(v(:,grid%Intf(i)%iNode(:)),2)/dble(grid%Intf(i)%nNode)
    end forall
  end function
  
  !> interpolate scaler v within grid from node to interface
  function itplNode2IntfScal(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2IntfScal(:) !< interpolated data on interface
    double precision vv(1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vrst=itplNode2IntfVect(vv,grid)
    call reallocArr(itplNode2IntfScal,size(vrst,2))
    itplNode2IntfScal(:)=vrst(1,:)
  end function
  
end module
