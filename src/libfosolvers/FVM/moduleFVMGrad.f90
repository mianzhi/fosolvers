!----------------------------------------------------------------------------- best with 100 columns

!> gradient schemes for FVM
module moduleFVMGrad
  private
  
  ! constants
  integer,parameter::GRAD_MIN_NEIB=6 !< minimum number of neighbor data sets to find gradient
  
  !> generic find gradient
  interface findGrad
    ! schemes for 3-D unstructured grid
    module procedure::findGradScal
    module procedure::findGradVect
    ! schemes for 1-D grid
    module procedure::findGrad1DScal
    module procedure::findGrad1DVect
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
          n=size(lNeib)
          do j=1,n
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
    double precision vv(1,size(v)),vFacetVal(1,grid%nFacet),vrst(DIMS,1,size(v))
    
    vv(1,:)=v(:)
    if(present(FacetVal))then
      vFacetVal(1,:)=FacetVal(:)
      vrst=findGradVect(vv,bind,grid,vFacetVal)
    else
      vrst=findGradVect(vv,bind,grid)
    end if
    findGradScal(:,:)=vrst(:,1,:)
  end function
  
  !> find derivative of vector v along x
  function findGrad1DVect(v,bind,grid)
    use moduleGrid1D
    double precision,intent(in)::v(:,:) !< vector data
    integer,intent(in)::bind !< bind with node/cell
    type(typeGrid1D),intent(inout)::grid !< grid on which v is defined
    double precision findGrad1DVect(size(v,1),size(v,2)) !< the result
    
    select case(bind)
    case(BIND_NODE)
      forall(i=2:grid%nNode-1)
        findGrad1DVect(:,i)=((grid%NodePos(i)-grid%NodePos(i-1))&
        &                    *(v(:,i+1)-v(:,i))/(grid%NodePos(i+1)-grid%NodePos(i))&
        &                   +(grid%NodePos(i+1)-grid%NodePos(i))&
        &                    *(v(:,i)-v(:,i-1))/(grid%NodePos(i)-grid%NodePos(i-1)))&
        &                   /(grid%NodePos(i+1)-grid%NodePos(i-1))
      end forall
      findGrad1DVect(:,1)=(v(:,2)-v(:,1))/(grid%NodePos(2)-grid%NodePos(1))
      findGrad1DVect(:,grid%nNode)=(v(:,grid%nNode)-v(:,grid%nNode-1))&
      &                            /(grid%NodePos(grid%nNode)-grid%NodePos(grid%nNode-1))
    case(BIND_CELL)
      forall(i=2:grid%nCell-1)
        findGrad1DVect(:,i)=((grid%CellPos(i)-grid%CellPos(i-1))&
        &                    *(v(:,i+1)-v(:,i))/(grid%CellPos(i+1)-grid%CellPos(i))&
        &                   +(grid%CellPos(i+1)-grid%CellPos(i))&
        &                    *(v(:,i)-v(:,i-1))/(grid%CellPos(i)-grid%CellPos(i-1)))&
        &                   /(grid%CellPos(i+1)-grid%CellPos(i-1))
      end forall
      findGrad1DVect(:,1)=(v(:,2)-v(:,1))/(grid%CellPos(2)-grid%CellPos(1))
      findGrad1DVect(:,grid%nCell)=(v(:,grid%nCell)-v(:,grid%nCell-1))&
      &                            /(grid%CellPos(grid%nCell)-grid%CellPos(grid%nCell-1))
    case default
    end select
  end function
  
  !> find derivative of scalar v along x
  function findGrad1DScal(v,bind,grid)
    use moduleGrid1D
    double precision,intent(in)::v(:) !< scalar data
    integer,intent(in)::bind !< bind with node/cell
    type(typeGrid1D),intent(inout)::grid !< grid on which v is defined
    double precision findGrad1DScal(size(v)) !< the result
    double precision vv(1,size(v)),vrst(1,size(v))
    
    vv(1,:)=v(:)
    vrst=findGrad1DVect(vv,bind,grid)
    findGrad1DScal(:)=vrst(1,:)
  end function
  
end module
