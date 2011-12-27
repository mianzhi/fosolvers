!----------------------------------------------------------------------------- best with 100 columns

!***********************
! finite voluem schemes
!***********************
module moduleFVM
  private
  
  ! generic find gradient
  interface findGrad
    module procedure::findGradScal
    module procedure::findGradVect
  end interface
  public findGrad
  
contains
  
  !------------------------------------------------------------
  ! find the gradient of a vector field at the k_th node/block
  !------------------------------------------------------------
  function findGradVect(k,v,binding)
    use moduleGrid
    use moduleMiscDataStruct
    integer,intent(in)::k
    double precision,intent(in)::v(:,:)
    integer,intent(in),optional::binding
    double precision findGradVect(size(v,2),DIMS)
    integer finalBinding
    integer,allocatable::listNeib(:)
    double precision listAlpha(Block(k)%SurfNum),vs,direction(DIMS)
    ! lapack related variables
    double precision,allocatable::dx(:,:),dv(:,:),s(:)
    integer,parameter::lwork=900 ! according to returned work(1) with lwork=-1
    double precision work(lwork)
    integer m,nEq,rank,iwork(50),er
    
    m=size(v,2)
    nEq=0
    
    if(present(binding))then
      finalBinding=binding
    else
      if(size(v,1)==nNode)then
        finalBinding=BIND_NODE
      else if(size(v,1)==nFacet)then
        finalBinding=BIND_FACET
      else if(size(v,1)==nBlock)then
        finalBinding=BIND_BLOCK
      else
        finalBinding=0
      end if
    end if
    
    select case(finalBinding)
      case(BIND_NODE) ! v is the node-bind data
      case(BIND_BLOCK) ! v is the block-bind data
        ! construct the neighbour list
        do i=1,Block(k)%SurfNum
          if(Block(k)%Neib(i)>0)then
            call extendArray(listNeib,1)
            nEq=nEq+1
            listNeib(nEq)=Block(k)%Neib(i)
          end if
        end do
        if(nEq<DIMS)then ! include more data
          do i=1,nEq
            do j=1,Block(listNeib(i))%SurfNum
              if(Block(listNeib(i))%Neib(j)>0.and.Block(listNeib(i))%Neib(j)/=k&
              &  .and.all(listNeib(:)/=Block(listNeib(i))%Neib(j)))then
                call extendArray(listNeib,1)
                nEq=nEq+1
                listNeib(nEq)=Block(listNeib(i))%Neib(j)
              end if
            end do
          end do
        end if
        ! weighted least-squares gradient evaluation
        allocate(dx(nEq,DIMS))
        allocate(dv(max(nEq,DIMS),m))
        allocate(s(min(nEq,DIMS)))
        forall(i=1:nEq)
          dx(i,:)=Block(listNeib(i))%PC(:)-Block(k)%PC(:)
          dv(i,:)=v(listNeib(i),:)-v(k,:)
          dv(i,:)=dv(i,:)/dot_product(dx(i,:),dx(i,:))
          dx(i,:)=dx(i,:)/dot_product(dx(i,:),dx(i,:))
        end forall
        call DGELSD(nEq,DIMS,m,dx,nEq,dv,max(nEq,DIMS),s,-1d0,rank,work,lwork,iwork,er)
        findGradVect(:,:)=transpose(dv(1:DIMS,:))
        ! regulate the gradient
        do i=1,m
          listAlpha(:)=1d0
          do j=1,Block(k)%SurfNum
            if(Block(k)%Neib(j)>0)then
              direction(:)=Block(Block(k)%Neib(j))%PC(:)-Block(k)%PC(:)
              direction(:)=direction(:)/norm2(direction)
              if(abs(dot_product(findGradVect(m,:)/norm2(findGradVect(m,:)),&
              &                  direction(:)))>0.5d0)then
                vs=v(k,m)+dot_product(Block(k)%SurfPC(j,:)-Block(k)%PC(:),findGradVect(m,:))
                if(vs>v(k,m))then
                  listAlpha(j)=min(1d0,(max(v(k,m),v(Block(k)%Neib(j),m))-v(k,m))/(vs-v(k,m)))
                else
                  listAlpha(j)=min(1d0,(min(v(k,m),v(Block(k)%Neib(j),m))-v(k,m))/(vs-v(k,m)))
                end if
              end if
            end if
          end do
          findGradVect(i,:)=minval(listAlpha)*findGradVect(i,:)
        end do
      case default
    end select
    ! clean up
    if(allocated(listNeib))then
      deallocate(listNeib)
    end if
    if(allocated(dx))then
      deallocate(dx)
    end if
    if(allocated(dv))then
      deallocate(dv)
    end if
    if(allocated(s))then
      deallocate(s)
    end if
  end function
  
  !------------------------------------------------------------
  ! find the gradient of a scaler field at the k_th node/block
  !------------------------------------------------------------
  function findGradScal(k,v,binding)
    use moduleGrid
    integer,intent(in)::k
    double precision,intent(in)::v(:)
    integer,intent(in),optional::binding
    double precision findGradScal(DIMS)
    double precision vv(size(v),1),vrst(1,DIMS)
    
    vv(:,1)=v(:)
    if(present(binding))then
      vrst=findGradVect(k,vv,binding)
    else
      vrst=findGradVect(k,vv)
    end if
    findGradScal(:)=vrst(1,:)
  end function
  
end module
