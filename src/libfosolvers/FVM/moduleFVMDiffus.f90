!----------------------------------------------------------------------------- best with 100 columns

!> FVM diffusion schemes
module moduleFVMDiffus
  private
  
  !> find diffusion through interface
  interface findDiffus
    module procedure findDiffusORTHScal
    module procedure findDiffusORTHVect
  end interface
  public findDiffus
  
contains
  
  !> find diffusion driven by vector v through interface using orthogonal scheme
  !> \f[ \int_{intf} \nabla \mathbf{v} \cdot \hat{n} dA \f]
  function findDiffusORTHVect(v,grid)
    use moduleGrid
    double precision,intent(in)::v(:,:) !< the vector which drives the diffusion
    type(typeGrid),intent(inout)::grid !< the grid on which v is defined
    double precision findDiffusORTHVect(size(v,1),grid%nBlock)
    double precision flowRate(size(v,1))
    
    call grid%updateIntfPos()
    call grid%updateIntfArea()
    call grid%updateBlockPos()
    findDiffusORTHVect(:,:)=0d0
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      flowRate(:)=(v(:,n)-v(:,m))*grid%IntfArea(i)/norm2(grid%BlockPos(:,n)-grid%BlockPos(:,m))
      findDiffusORTHVect(:,m)=findDiffusORTHVect(:,m)+flowRate(:)
      findDiffusORTHVect(:,n)=findDiffusORTHVect(:,n)-flowRate(:)
    end do
  end function
  
  !> find diffusion driven by scalar v through interface using orthogonal scheme
  !> \f[ \int_{intf} \nabla v \cdot \hat{n} dA \f]
  function findDiffusORTHScal(v,grid)
    use moduleGrid
    double precision,intent(in)::v(:) !< the scalar which drives the diffusion
    type(typeGrid),intent(inout)::grid !< the grid on which v is defined
    double precision findDiffusORTHScal(grid%nBlock)
    double precision vv(1,size(v)),vrst(1,grid%nBlock)
    
    vv(1,:)=v(:)
    vrst=findDiffusORTHVect(vv,grid)
    findDiffusORTHScal(:)=vrst(1,:)
  end function
  
end module
