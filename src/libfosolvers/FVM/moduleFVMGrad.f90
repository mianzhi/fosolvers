!----------------------------------------------------------------------------- best with 100 columns

!> gradient schemes for FVM
module moduleFVMGrad
  private
  
  ! constants
  integer,parameter::GRAD_MIN_NEIB=4 !< minimum number of neighbor data sets to find gradient
  
  !> generic find gradient
  interface findGrad
    module procedure::findGradScal
    module procedure::findGradVect
  end interface
  public findGrad
  
contains
  
  !> find gradient of vector field v on grid
  function findGradVect(v,grid,bind)
    use moduleGrid
    use moduleSimpleSetLogic
    double precision,intent(in)::v(:,:) !< vector data
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    integer,intent(in)::bind !< bind with node/block
    double precision::findGradVect(DIMS,size(v,1),size(v,2)) !< the gradient of v
    integer,allocatable::lNeib(:)
    integer nEq,rank,iwork(50),er
    integer,parameter::lwork=900 ! according to returned work(1) with lwork=-1
    double precision work(lwork)
    double precision,allocatable::dx(:,:),dv(:,:),s(:)
    
    m=size(v,1)
    
    select case(bind)
    case(BIND_NODE)
      call grid%updateNodeNeib()
      do i=1,grid%nNode
        do j=1,size(grid%NodeNeibBlock(i)%dat)
          call applUnion(lNeib,grid%Block(grid%NodeNeibBlock(i)%dat(j))%iNode(:))
        end do
        call applComplement(lNeib,[i])
        nEq=size(lNeib)
        allocate(dx(nEq,DIMS))
        allocate(dv(max(nEq,DIMS),m))
        allocate(s(min(nEq,DIMS)))
        forall(j=1:nEq)
          dx(j,:)=grid%NodePos(:,lNeib(j))-grid%NodePos(:,i)
          dv(j,:)=v(:,lNeib(j))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        call DGELSD(nEq,DIMS,m,dx,nEq,dv,max(nEq,DIMS),s,-1d0,rank,work,lwork,iwork,er)
        findGradVect(:,:,i)=dv(1:DIMS,:)
        deallocate(lNeib)
        deallocate(dx)
        deallocate(dv)
        deallocate(s)
      end do
    case(BIND_BLOCK)
      call grid%updateBlockNeib()
      call grid%updateBlockPos()
      do i=1,grid%nBlock
        call applUnion(lNeib,grid%BlockNeibBlock(i)%dat(:))
        call applComplement(lNeib,[i,0])
        nEq=size(lNeib)
        do while(nEq<GRAD_MIN_NEIB) ! include more data
          do j=1,nEq
            call applUnion(lNeib,grid%BlockNeibBlock(lNeib(j))%dat(:))
          end do
          call applComplement(lNeib,[i,0])
          nEq=size(lNeib)
        end do
        allocate(dx(nEq,DIMS))
        allocate(dv(max(nEq,DIMS),m))
        allocate(s(min(nEq,DIMS)))
        forall(j=1:nEq)
          dx(j,:)=grid%BlockPos(:,lNeib(j))-grid%BlockPos(:,i)
          dv(j,:)=v(:,lNeib(j))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        call DGELSD(nEq,DIMS,m,dx,nEq,dv,max(nEq,DIMS),s,-1d0,rank,work,lwork,iwork,er)
        findGradVect(:,:,i)=dv(1:DIMS,:)
        deallocate(lNeib)
        deallocate(dx)
        deallocate(dv)
        deallocate(s)
      end do
    case default
    end select
  end function
  
  !> find gradient of scalar field v on grid
  function findGradScal(v,grid,bind)
    use moduleGrid
    double precision,intent(in)::v(:) !< scalar data
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    integer,intent(in)::bind !< bind with node/block
    double precision::findGradScal(DIMS,size(v)) !< the gradient of v
    double precision vv(1,size(v)),vrst(DIMS,1,size(v))
    
    vv(1,:)=v(:)
    vrst=findGradVect(vv,grid,bind)
    findGradScal(:,:)=vrst(:,1,:)
  end function
  
end module
