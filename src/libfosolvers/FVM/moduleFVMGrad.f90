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
  function findGradVect(v,bind,grid,FacetVal)
    use moduleGrid
    use moduleBasicDataStruct
    use moduleSimpleSetLogic
    double precision,intent(in)::v(:,:) !< vector data
    integer,intent(in)::bind !< bind with node/block
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,intent(in),optional::FacetVal(size(v,1),grid%nFacet) !< vector data on facet
    double precision::findGradVect(DIMS,size(v,1),size(v,2)) !< the gradient of v
    integer,allocatable::lNeib(:),lFacet(:)
    integer nEq,rank,iwork(50),er
    integer,parameter::lwork=900 ! according to returned work(1) with lwork=-1
    double precision work(lwork)
    double precision,allocatable::dx(:,:),dv(:,:),s(:)
    
    m=size(v,1)
    
    select case(bind)
    case(BIND_NODE)
      call grid%updateNodeNeib()
      call grid%updateFacetPos()
      do i=1,grid%nNode
        do j=1,size(grid%NodeNeibBlock(i)%dat)
          call applUnion(lNeib,grid%Block(grid%NodeNeibBlock(i)%dat(j))%iNode(:))
        end do
        call applComplement(lNeib,[i])
        allocate(lFacet(0))
        if(present(FacetVal))then
          if(allocated(grid%NodeNeibFacet(i)%dat))then
            call pushArr(lFacet,grid%NodeNeibFacet(i)%dat)
          end if
        end if
        nEq=size(lNeib)+size(lFacet)
        allocate(dx(nEq,DIMS))
        allocate(dv(max(nEq,DIMS),m))
        allocate(s(min(nEq,DIMS)))
        forall(j=1:size(lNeib))
          dx(j,:)=grid%NodePos(:,lNeib(j))-grid%NodePos(:,i)
          dv(j,:)=v(:,lNeib(j))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        forall(j=size(lNeib)+1:nEq)
          dx(j,:)=grid%FacetPos(:,lFacet(j-size(lNeib)))-grid%NodePos(:,i)
          dv(j,:)=FacetVal(:,lFacet(j-size(lNeib)))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        call DGELSD(nEq,DIMS,m,dx,nEq,dv,max(nEq,DIMS),s,-1d0,rank,work,lwork,iwork,er)
        findGradVect(:,:,i)=dv(1:DIMS,:)
        deallocate(lNeib)
        deallocate(lFacet)
        deallocate(dx)
        deallocate(dv)
        deallocate(s)
      end do
    case(BIND_BLOCK)
      call grid%updateBlockNeib()
      call grid%updateBlockPos()
      call grid%updateFacetPos()
      do i=1,grid%nBlock
        call applUnion(lNeib,grid%BlockNeibBlock(i)%dat(:))
        call applComplement(lNeib,[i,0])
        allocate(lFacet(0))
        if(present(FacetVal))then
          call applUnion(lFacet,grid%BlockNeibFacet(i)%dat(:))
          call applComplement(lFacet,[0])
        end if
        nEq=size(lNeib)+size(lFacet)
        do while(nEq<GRAD_MIN_NEIB) ! include more data
          do j=1,nEq
            call applUnion(lNeib,grid%BlockNeibBlock(lNeib(j))%dat(:))
          end do
          call applComplement(lNeib,[i,0])
          nEq=size(lNeib)+size(lFacet)
        end do
        allocate(dx(nEq,DIMS))
        allocate(dv(max(nEq,DIMS),m))
        allocate(s(min(nEq,DIMS)))
        forall(j=1:size(lNeib))
          dx(j,:)=grid%BlockPos(:,lNeib(j))-grid%BlockPos(:,i)
          dv(j,:)=v(:,lNeib(j))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        forall(j=size(lNeib)+1:nEq)
          dx(j,:)=grid%FacetPos(:,lFacet(j-size(lNeib)))-grid%BlockPos(:,i)
          dv(j,:)=FacetVal(:,lFacet(j-size(lNeib)))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        call DGELSD(nEq,DIMS,m,dx,nEq,dv,max(nEq,DIMS),s,-1d0,rank,work,lwork,iwork,er)
        findGradVect(:,:,i)=dv(1:DIMS,:)
        deallocate(lNeib)
        deallocate(lFacet)
        deallocate(dx)
        deallocate(dv)
        deallocate(s)
      end do
    case default
    end select
  end function
  
  !> find gradient of scalar field v on grid
  function findGradScal(v,bind,grid,FacetVal)
    use moduleGrid
    double precision,intent(in)::v(:) !< scalar data
    integer,intent(in)::bind !< bind with node/block
    type(typeGrid),intent(inout)::grid !< grid on which v is defined
    double precision,intent(in),optional::FacetVal(grid%nFacet) !< scalar data on facet
    double precision::findGradScal(DIMS,size(v)) !< the gradient of v
    double precision vv(1,size(v)),vrst(DIMS,1,size(v))
    double precision,allocatable::vFacetVal(:,:)
    
    vv(1,:)=v(:)
    if(present(FacetVal))then
      allocate(vFacetVal(1,size(FacetVal)))
      vFacetVal(1,:)=FacetVal(:)
      vrst=findGradVect(vv,bind,grid,vFacetVal)
      deallocate(vFacetVal)
    else
      vrst=findGradVect(vv,bind,grid)
    end if
    findGradScal(:,:)=vrst(:,1,:)
  end function
  
end module
