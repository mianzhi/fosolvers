!----------------------------------------------------------------------------- best with 100 columns

!> interpolation schemes
module moduleInterpolation
  private
  
  !> interpolate from block to interface
  interface itplBlock2Intf
    module procedure::itplBCDVect
    module procedure::itplBCDScal
  end interface
  public itplBlock2Intf
  
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
    double precision,intent(in)::v(:) !< block data to be interpolated
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,allocatable::itplBCDScal(:) !< interpolated data on interface
    double precision vv(1,size(v))
    double precision,allocatable::vrst(:,:)
    
    vv(1,:)=v(:)
    vrst=itplBCDVect(vv,grid)
    allocate(itplBCDScal(size(vrst,2)))
    itplBCDScal(:)=vrst(1,:)
  end function
  
end module
