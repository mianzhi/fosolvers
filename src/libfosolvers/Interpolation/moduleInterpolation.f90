!----------------------------------------------------------------------------- best with 100 columns

!> interpolation schemes
module moduleInterpolation
  private
  
  !> interpolate from block to node
  interface itplBlock2Node
    module procedure::itplBlock2NodeVect
    module procedure::itplBlock2NodeScal
  end interface
  public itplBlock2Node
  
  !> interpolate from block to interface
  interface itplBlock2Intf
    module procedure::itplBCDVect
    module procedure::itplBCDScal
    module procedure::itplBCDCVect
    module procedure::itplBCDCScal
  end interface
  public itplBlock2Intf
  
  !> interpolate from node to block
  interface itplNode2Block
    module procedure::itplNode2BlockMat
    module procedure::itplNode2BlockVect
    module procedure::itplNode2BlockScal
  end interface
  public itplNode2Block
  
  !> interpolate from node to edge
  interface itplNode2Edge
    module procedure::itplNode2EdgeVect
    module procedure::itplNode2EdgeScal
  end interface
  public itplNode2Edge
  
  !> interpolate from node to interface
  interface itplNode2Intf
    module procedure::itplNode2IntfVect
    module procedure::itplNode2IntfScal
  end interface
  public itplNode2Intf
  
contains
  
  !> interpolate vector v within grid from block to node
  function itplBlock2NodeVect(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:) !< block data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBlock2NodeVect(:,:) !< interpolated data on node
    double precision weight,weightTot,temp(size(v,1))
    
    call grid%updateBlockPos()
    call grid%updateNodeNeib()
    call reallocArr(itplBlock2NodeVect,size(v,1),grid%nNode)
    do i=1,grid%nNode
      weightTot=0d0
      temp(:)=0d0
      do j=1,size(grid%NodeNeibBlock(i)%dat)
        k=grid%NodeNeibBlock(i)%dat(j)
        weight=1d0/dot_product(grid%BlockPos(:,k)-grid%NodePos(:,i),&
        &                      grid%BlockPos(:,k)-grid%NodePos(:,i))
        weightTot=weightTot+weight
        temp(:)=temp(:)+weight*v(:,k)
      end do
      itplBlock2NodeVect(:,i)=temp(:)/weightTot
    end do
  end function
  
  !> interpolate scalar v within grid from block to node
  function itplBlock2NodeScal(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:) !< block data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBlock2NodeScal(:) !< interpolated data on node
    double precision vv(1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vrst=itplBlock2NodeVect(vv,grid)
    call reallocArr(itplBlock2NodeScal,size(vrst,2))
    itplBlock2NodeScal(:)=vrst(1,:)
  end function
  
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
  
  !> interpolate matrix v within grid from node to block
  function itplNode2BlockMat(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:,:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2BlockMat(:,:,:) !< interpolated data on block
    
    call reallocArr(itplNode2BlockMat,size(v,1),size(v,2),grid%nBlock)
    forall(i=1:grid%nBlock)
      itplNode2BlockMat(:,:,i)=sum(v(:,:,grid%Block(i)%iNode(:)),3)/dble(grid%Block(i)%nNode)
    end forall
  end function
  
  !> interpolate vector v within grid from node to block
  function itplNode2BlockVect(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2BlockVect(:,:) !< interpolated data on block
    
    call reallocArr(itplNode2BlockVect,size(v,1),grid%nBlock)
    forall(i=1:grid%nBlock)
      itplNode2BlockVect(:,i)=sum(v(:,grid%Block(i)%iNode(:)),2)/dble(grid%Block(i)%nNode)
    end forall
  end function
  
  !> interpolate scaler v within grid from node to block
  function itplNode2BlockScal(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2BlockScal(:) !< interpolated data on block
    double precision vv(1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vrst=itplNode2BlockVect(vv,grid)
    call reallocArr(itplNode2BlockScal,size(vrst,2))
    itplNode2BlockScal(:)=vrst(1,:)
  end function
  
  !> interpolate vector v within grid from node to edge
  function itplNode2EdgeVect(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:,:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2EdgeVect(:,:) !< interpolated data on edge
    
    call grid%updateEdge()
    call reallocArr(itplNode2EdgeVect,size(v,1),grid%nEdge)
    forall(i=1:grid%nEdge)
      itplNode2EdgeVect(:,i)=sum(v(:,grid%Edge(i)%iNode(:)),2)/dble(grid%Edge(i)%nNode)
    end forall
  end function
  
  !> interpolate scaler v within grid from node to edge
  function itplNode2EdgeScal(v,grid)
    use moduleGrid
    use moduleBasicDataStruct
    double precision,intent(in)::v(:) !< node data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplNode2EdgeScal(:) !< interpolated data on edge
    double precision vv(1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vrst=itplNode2EdgeVect(vv,grid)
    call reallocArr(itplNode2EdgeScal,size(vrst,2))
    itplNode2EdgeScal(:)=vrst(1,:)
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
