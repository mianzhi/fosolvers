!----------------------------------------------------------------------------- best with 100 columns

!***********************
! finite volume schemes
!***********************
module moduleFVM
  private
  
  ! constants
  integer,parameter::GRAD_MIN_NEIB=5
  
  ! generic find gradient
  interface findGrad
    module procedure::findGradScal
    module procedure::findGradVect
  end interface
  public findGrad
  
  ! generic block surface interpolation (inner surface)
  interface itplBCD
    module procedure::itplBCDScal
    module procedure::itplBCDVect
  end interface
  public itplBCD
  interface itplCBCD
    module procedure::itplCBCDScal
    module procedure::itplCBCDVect
  end interface
  public itplCBCD
  
  ! generic block surface interpolation (boundary surface)
  interface itplBS
    module procedure::itplBSScal
    module procedure::itplBSVect
  end interface
  public itplBS
  interface itplBSAP
    module procedure::itplBSAPScal
    module procedure::itplBSAPVect
  end interface
  public itplBSAP
  
  ! generic diffusion flux evaluation (inner surface)
  interface diffuseORTH
    module procedure::diffuseORTHScal
    module procedure::diffuseORTHVect
  end interface
  public diffuseORTH
  interface diffuseSD
    module procedure::diffuseSDScal
    module procedure::diffuseSDVect
  end interface
  public diffuseSD
  
  ! generic diffusion flux evaluation (boundary surface)
  interface diffuseBSORTH
    module procedure::diffuseBSORTHScal
    module procedure::diffuseBSORTHVect
  end interface
  public diffuseBSORTH
  
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
        ! construct the neighbour list
        do i=1,size(Node(k)%BlockInd)
          do j=1,Block(Node(k)%BlockInd(i))%NodeNum
            if(Block(Node(k)%BlockInd(i))%NodeInd(j)/=k)then
              call extendArray(listNeib,1)
              nEq=nEq+1
              listNeib(nEq)=Block(Node(k)%BlockInd(i))%NodeInd(j)
            end if
          end do
        end do
        ! weighted least-squares gradient evaluation
        allocate(dx(nEq,DIMS))
        allocate(dv(max(nEq,DIMS),m))
        allocate(s(min(nEq,DIMS)))
        forall(i=1:nEq)
          dx(i,:)=Node(listNeib(i))%Pos(:)-Node(k)%Pos(:)
          dv(i,:)=v(listNeib(i),:)-v(k,:)
          dv(i,:)=dv(i,:)/dot_product(dx(i,:),dx(i,:))
          dx(i,:)=dx(i,:)/dot_product(dx(i,:),dx(i,:))
        end forall
        call DGELSD(nEq,DIMS,m,dx,nEq,dv,max(nEq,DIMS),s,-1d0,rank,work,lwork,iwork,er)
        findGradVect(:,:)=transpose(dv(1:DIMS,:))
        
      case(BIND_BLOCK) ! v is the block-bind data
        ! construct the neighbour list
        do i=1,Block(k)%SurfNum
          if(Block(k)%Neib(i)>0)then
            call extendArray(listNeib,1)
            nEq=nEq+1
            listNeib(nEq)=Block(k)%Neib(i)
          end if
        end do
        do while(nEq<GRAD_MIN_NEIB) ! include more data
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
        end do
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
  
  !----------------------------------------------------------------------------------------------
  ! interpolate block vector v to n_th surface of m_th block using block-centre-direction scheme
  !----------------------------------------------------------------------------------------------
  function itplBCDVect(m,n,v)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision itplBCDVect(size(v,2))
    double precision PF(DIMS),PS(DIMS),alpha
    
    itplBCDVect(:)=0d0
    if(Block(m)%Neib(n)>0)then
      PF(:)=Block(Block(m)%Neib(n))%PC(:)-Block(m)%PC(:)
      PS(:)=Block(m)%SurfPC(n,:)-Block(m)%PC(:)
      alpha=dot_product(PS,PF)/dot_product(PF,PF)
      itplBCDVect(:)=(1d0-alpha)*v(m,:)+alpha*v(Block(m)%Neib(n),:)
    else
      call showError('invalid inner surface number.')
    end if
  end function
  
  !----------------------------------------------------------------------------------------------
  ! interpolate block scalar v to n_th surface of m_th block using block-centre-direction scheme
  !----------------------------------------------------------------------------------------------
  function itplBCDScal(m,n,v)
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision itplBCDScal
    double precision vv(size(v),1),vrst(1)
    
    vv(:,1)=v(:)
    vrst=itplBCDVect(m,n,vv)
    itplBCDScal=vrst(1)
  end function
  
  !----------------------------------------------------------
  ! interpolate block vector v to n_th surface of m_th block
  ! using corrected block-centre-direction scheme
  !----------------------------------------------------------
  function itplCBCDVect(m,n,v,grad)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision,intent(in),optional::grad(:,:,:)
    double precision itplCBCDVect(size(v,2))
    double precision gradP(size(v,2),DIMS),Saux(DIMS),alpha,dPS,dSF
    
    itplCBCDVect(:)=0d0
    if(Block(m)%Neib(n)>0)then
      if(present(grad))then
        gradP(:,:)=grad(m,:,:)
      else
        gradP=findGradVect(m,v,binding=BIND_BLOCK)
      end if
      dPS=norm2(Block(m)%SurfPC(n,:)-Block(m)%PC(:))
      dSF=norm2(Block(Block(m)%Neib(n))%PC(:)-Block(m)%SurfPC(n,:))
      alpha=dPS/(dPS+dSF)
      Saux(:)=Block(m)%PC(:)+alpha*(Block(Block(m)%Neib(n))%PC(:)-Block(m)%PC(:))
      itplCBCDVect(:)=(1d0-alpha)*v(m,:)+alpha*v(Block(m)%Neib(n),:)&
      &               +matmul(gradP(:,:),Block(m)%SurfPC(n,:)-Saux(:))
    else
      call showError('invalid inner surface number.')
    end if
  end function
  
  !----------------------------------------------------------
  ! interpolate block scalar v to n_th surface of m_th block
  ! using corrected block-centre-direction scheme
  !----------------------------------------------------------
  function itplCBCDScal(m,n,v,grad)
    use moduleGrid
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision,intent(in),optional::grad(:,:)
    double precision itplCBCDScal
    double precision vv(size(v),1),vrst(1),vgrad(size(v),1,DIMS)
    
    vv(:,1)=v(:)
    if(present(grad))then
      vgrad(:,1,:)=grad(:,:)
      vrst=itplCBCDVect(m,n,vv,grad=vgrad)
    else
      vrst=itplCBCDVect(m,n,vv)
    end if
    itplCBCDScal=vrst(1)
  end function
  
  !-----------------------------------------------------------------------------
  ! interpolate block vector v to n_th surface (boundary surface) of m_th block
  !-----------------------------------------------------------------------------
  function itplBSVect(m,n,v,ghostVal,grad)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision,intent(in),optional::ghostVal(size(v,2))
    double precision,intent(in),optional::grad(:,:,:)
    double precision itplBSVect(size(v,2))
    double precision gradP(size(v,2),DIMS)
    
    itplBSVect(:)=0d0
    if(Block(m)%Neib(n)<0)then
      if(present(ghostVal))then
        itplBSVect(:)=(v(m,:)+ghostVal(:))/2d0
      else
        if(present(grad))then
          gradP(:,:)=grad(m,:,:)
        else
          gradP=findGradVect(m,v,binding=BIND_BLOCK)
        end if
        itplBSVect(:)=v(m,:)+matmul(gradP(:,:),Block(m)%PC(:)-Block(m)%SurfPC(n,:))
      end if
    else
      call showError('invalid boundary surface number.')
    end if
  end function
  
  !-----------------------------------------------------------------------------
  ! interpolate block scalar v to n_th surface (boundary surface) of m_th block
  !-----------------------------------------------------------------------------
  function itplBSScal(m,n,v,ghostVal,grad)
    use moduleGrid
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision,intent(in),optional::ghostVal
    double precision,intent(in),optional::grad(:,:)
    double precision itplBSScal
    double precision vv(size(v),1),vrst(1),vghostVal(1),vgrad(size(v),1,DIMS)
    
    vv(:,1)=v(:)
    if(present(ghostVal))then
      vghostVal(1)=ghostVal
      vrst=itplBSVect(m,n,vv,ghostVal=vghostVal)
    else
      vv(:,1)=v(:)
      if(present(grad))then
        vgrad(:,1,:)=grad(:,:)
        vrst=itplBSVect(m,n,vv,grad=vgrad)
      else
        vrst=itplBSVect(m,n,vv)
      end if
    end if
    itplBSScal=vrst(1)
  end function
  
  !-----------------------------------------------------------------------------
  ! interpolate block vector v to n_th surface (boundary surface) of m_th block
  ! using auxiliary-point scheme
  !-----------------------------------------------------------------------------
  function itplBSAPVect(m,n,v,ghostVal,grad)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision,intent(in)::ghostVal(size(v,2))
    double precision,intent(in),optional::grad(:,:,:)
    double precision itplBSAPVect(size(v,2))
    double precision gradP(size(v,2),DIMS),distA(DIMS)
    
    itplBSAPVect(:)=0d0
    if(Block(m)%Neib(n)<0)then
      if(present(grad))then
        gradP(:,:)=grad(m,:,:)
      else
        gradP=findGradVect(m,v,binding=BIND_BLOCK)
      end if
      distA(:)=Block(m)%SurfPC(n,:)&
      &       -Block(m)%SurfNorm(n,:)&
      &        *dot_product(Block(m)%SurfNorm(n,:),Block(m)%SurfPC(n,:)-Block(m)%PC(:))&
      &       -Block(m)%PC(:)
      itplBSAPVect(:)=(v(m,:)+matmul(gradP(:,:),distA(:))+ghostVal(:))/2d0
    else
      call showError('invalid boundary surface number.')
    end if
  end function
  
  !-----------------------------------------------------------------------------
  ! interpolate block scalar v to n_th surface (boundary surface) of m_th block
  ! using auxiliary-point scheme
  !-----------------------------------------------------------------------------
  function itplBSAPScal(m,n,v,ghostVal,grad)
    use moduleGrid
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision,intent(in)::ghostVal
    double precision,intent(in),optional::grad(:,:)
    double precision itplBSAPScal
    double precision vv(size(v),1),vrst(1),vghostVal(1),vgrad(size(v),1,DIMS)
    
    vv(:,1)=v(:)
    vghostVal(1)=ghostVal
    if(present(grad))then
      vgrad(:,1,:)=grad(:,:)
      vrst=itplBSAPVect(m,n,vv,vghostVal,grad=vgrad)
    else
      vrst=itplBSAPVect(m,n,vv,vghostVal)
    end if
    itplBSAPScal=vrst(1)
  end function
  
  !--------------------------------------------------------------------------
  ! evaluate the diffusion flux into the m_th block through its n_th surface
  ! driven by vector v using orthogonal scheme
  !
  !  / /
  !  | | __      ^
  !  | | \/[v] * n dA
  !  | |
  !  / /
  !  Surf n of Block m
  !--------------------------------------------------------------------------
  function diffuseORTHVect(m,n,v)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision diffuseORTHVect(size(v,2))
    
    diffuseORTHVect(:)=0d0
    if(Block(m)%Neib(n)>0)then
      diffuseORTHVect(:)=(v(Block(m)%Neib(n),:)-v(m,:))*Block(m)%SurfArea(n)&
      &                  /norm2(Block(Block(m)%Neib(n))%PC(:)-Block(m)%PC(:))
    else
      call showError('invalid inner surface number.')
    end if
  end function
  
  !--------------------------------------------------------------------------
  ! evaluate the diffusion flux into the m_th block through its n_th surface
  ! driven by scalar v using orthogonal scheme
  !
  !  / /
  !  | | __    ^
  !  | | \/v * n dA
  !  | |
  !  / /
  !  Surf n of Block m
  !--------------------------------------------------------------------------
  function diffuseORTHScal(m,n,v)
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision diffuseORTHScal
    double precision vv(size(v),1),vrst(1)
    
    vv(:,1)=v(:)
    vrst=diffuseORTHVect(m,n,vv)
    diffuseORTHScal=vrst(1)
  end function
  
  !---------------------------------------------------------------------------------------------
  ! evaluate the diffusion flux into the m_th block through its n_th surface (boundary surface)
  ! driven by vector v using orthogonal scheme
  !
  !  / /
  !  | | __      ^
  !  | | \/[v] * n dA
  !  | |
  !  / /
  !  Surf n of Block m
  !---------------------------------------------------------------------------------------------
  function diffuseBSORTHVect(m,n,v,ghostVal)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision,intent(in)::ghostVal(size(v,2))
    double precision diffuseBSORTHVect(size(v,2))
    
    diffuseBSORTHVect(:)=0d0
    if(Block(m)%Neib(n)<0)then
      diffuseBSORTHVect(:)=(ghostVal(:)-v(m,:))*Block(m)%SurfArea(n)/2d0&
      &                    /norm2(Block(m)%SurfPC(n,:)-Block(m)%PC(:))
    else
      call showError('invalid boundary surface number.')
    end if
  end function
  
  !---------------------------------------------------------------------------------------------
  ! evaluate the diffusion flux into the m_th block through its n_th surface (boundary surface)
  ! driven by scalar v using orthogonal scheme
  !
  !  / /
  !  | | __    ^
  !  | | \/v * n dA
  !  | |
  !  / /
  !  Surf n of Block m
  !---------------------------------------------------------------------------------------------
  function diffuseBSORTHScal(m,n,v,ghostVal)
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision,intent(in)::ghostVal
    double precision diffuseBSORTHScal
    double precision vv(size(v),1),vghostVal(1),vrst(1)
    
    vv(:,1)=v(:)
    vghostVal(1)=ghostVal
    vrst=diffuseBSORTHVect(m,n,vv,vghostVal)
    diffuseBSORTHScal=vrst(1)
  end function
  
  !--------------------------------------------------------------------------
  ! evaluate the diffusion flux into the m_th block through its n_th surface
  ! driven by vector v using surface decomposition scheme
  !
  !  / /
  !  | | __      ^
  !  | | \/[v] * n dA
  !  | |
  !  / /
  !  Surf n of Block m
  !--------------------------------------------------------------------------
  function diffuseSDVect(m,n,v,grad)
    use moduleGrid
    use moduleUtility
    integer,intent(in)::m,n
    double precision,intent(in)::v(:,:)
    double precision,intent(in),optional::grad(size(v,1),size(v,2),DIMS)
    double precision diffuseSDVect(size(v,2))
    double precision gradField(size(v,1),size(v,2),DIMS),gradAvg(size(v,2),DIMS),&
    &                dPF,sf(DIMS),tf(DIMS),rf(DIMS),Afs
    
    diffuseSDVect(:)=0d0
    if(Block(m)%Neib(n)>0)then
      if(present(grad))then
        do i=1,size(v,2)
          gradAvg(i,:)=itplBCDVect(m,n,grad(:,i,:))
        end do
      else
        gradField(m,:,:)=findGradVect(m,v,binding=BIND_BLOCK)
        gradField(Block(m)%Neib(n),:,:)=findGradVect(Block(m)%Neib(n),v,binding=BIND_BLOCK)
        do i=1,size(v,2)
          gradAvg(i,:)=itplBCDVect(m,n,gradField(:,i,:))
        end do
      end if
      sf(:)=Block(Block(m)%Neib(n))%PC(:)-Block(m)%PC(:)
      dPF=norm2(sf)
      sf(:)=sf(:)/dPF
      tf(:)=Node(Block(m)%SurfNodeInd(n,1))%Pos(:)-Node(Block(m)%SurfNodeInd(n,2))%Pos(:)
      tf(:)=tf(:)/norm2(tf)
      rf(1)=Block(m)%SurfNorm(n,2)*tf(3)-Block(m)%SurfNorm(n,3)*tf(2)
      rf(2)=-Block(m)%SurfNorm(n,1)*tf(3)+Block(m)%SurfNorm(n,3)*tf(1)
      rf(3)=Block(m)%SurfNorm(n,1)*tf(2)-Block(m)%SurfNorm(n,2)*tf(1)
      Afs=Block(m)%SurfArea(n)/dot_product(sf,Block(m)%SurfNorm(n,:))
      diffuseSDVect(:)=Afs*((v(Block(m)%Neib(n),:)-v(m,:))/dPF&
      &                     -matmul(gradAvg,dot_product(sf,tf)*tf+dot_product(sf,rf)*rf))
    else
      call showError('invalid inner surface number.')
    end if
  end function
  
  !--------------------------------------------------------------------------
  ! evaluate the diffusion flux into the m_th block through its n_th surface
  ! driven by scalar v using surface decomposition scheme
  !
  !  / /
  !  | | __    ^
  !  | | \/v * n dA
  !  | |
  !  / /
  !  Surf n of Block m
  !--------------------------------------------------------------------------
  function diffuseSDScal(m,n,v,grad)
    use moduleGrid
    integer,intent(in)::m,n
    double precision,intent(in)::v(:)
    double precision,intent(in),optional::grad(:,:)
    double precision diffuseSDScal
    double precision vv(size(v),1),vrst(1),vgrad(size(v),1,DIMS)
    
    vv(:,1)=v(:)
    if(present(grad))then
      vgrad(:,1,:)=grad(:,:)
      vrst=diffuseSDVect(m,n,vv,grad=vgrad)
    else
      vrst=diffuseSDVect(m,n,vv)
    end if
    diffuseSDScal=vrst(1)
  end function
    
end module
