!----------------------------------------------------------------------------- best with 100 columns

!> FVM convection schemes
module moduleFVMConvect
  private
  
  !> find convection through interface
  interface findConvect
    module procedure findConvectUpWindScal
    module procedure findConvectUpWindVect
    module procedure findConvectTVDScal
    module procedure findConvectTVDVect
  end interface
  public findConvect
  
  !> find displacement convection (for ALE rezoning)
  interface findDispConvect
    module procedure findDispConvectUpWindScal
    module procedure findDispConvectUpWindVect
    module procedure findDispConvectTVDScal
    module procedure findDispConvectTVDVect
  end interface
  public findDispConvect
  
  ! flux limiters
  interface
    pure function modelLimiter(r)
      double precision,intent(in)::r
      double precision modelLimiter
    end function
  end interface
  procedure(modelLimiter),pointer::vanLeer=>vanLeerLimiter
  procedure(modelLimiter),pointer::minmod=>minmodLimiter
  public vanLeer
  public minmod
  
contains
  
  !> find the convection of vector phi driven by velocity u through interface using upwind scheme
  !> \f[ \int_{intf} \mathbf{\Phi} (\mathbf{u} \cdot \hat{n}) dA \f]
  function findConvectUpWindVect(phi,u,ubind,grid)
    use moduleGrid
    double precision,intent(in)::phi(:,:) !< variable to be convected
    double precision,intent(in)::u(:,:) !< velocity which drives the convection
    integer,intent(in)::ubind !< bind velocity with node/block/interface
    type(typeGrid),intent(inout)::grid !< the grid
    double precision findConvectUpWindVect(size(phi,1),grid%nBlock) !< increment due to convection
    double precision F,flowRate(size(phi,1))
    double precision,allocatable::uIntf(:,:)
    
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    findConvectUpWindVect(:,:)=0d0
    allocate(uIntf(DIMS,grid%nIntf))
    select case(ubind)
    case(BIND_NODE)
      !TODO:interpolation from node to interface
    case(BIND_BLOCK)
      !TODO:interpolation for co-located u
    case(BIND_INTF)
      uIntf(:,:)=u(:,:)
    case default
    end select
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      F=dot_product(grid%IntfNorm(:,i),uIntf(:,i))*grid%IntfArea(i)
      if(F>0d0)then
        flowRate(:)=-F*phi(:,m)
      else
        flowRate(:)=-F*phi(:,n)
      end if
      findConvectUpWindVect(:,m)=findConvectUpWindVect(:,m)+flowRate(:)
      findConvectUpWindVect(:,n)=findConvectUpWindVect(:,n)-flowRate(:)
    end do
  end function
  
  !> find the convection of scalar phi driven by velocity u through interface using upwind scheme
  !> \f[ \int_{intf} \Phi (\mathbf{u} \cdot \hat{n}) dA \f]
  function findConvectUpWindScal(phi,u,ubind,grid)
    use moduleGrid
    double precision,intent(in)::phi(:) !< variable to be convected
    double precision,intent(in)::u(:,:) !< velocity which drives the convection
    integer,intent(in)::ubind !< bind velocity with node/block/interface
    type(typeGrid),intent(inout)::grid !< the grid
    double precision findConvectUpWindScal(grid%nBlock) !< increment due to convection
    double precision vphi(1,size(phi)),vrst(1,grid%nBlock)
    
    vphi(1,:)=phi(:)
    vrst=findConvectUpWindVect(vphi,u,ubind,grid)
    findConvectUpWindScal(:)=vrst(1,:)
  end function
  
  !> van Leer flux limiter
  pure function vanLeerLimiter(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision vanLeerLimiter !< limit function
    
    vanLeerLimiter=(r+abs(r))/(1d0+abs(r))
  end function
  
  !> minmod flux limiter
  pure function minmodLimiter(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision minmodLimiter !< limit function
    
    minmodLimiter=max(0d0,min(r,1d0))
  end function
  
  !> find the convection of vector phi driven by velocity u through interface using TVD scheme
  !> \f[ \int_{intf} \mathbf{\Phi} (\mathbf{u} \cdot \hat{n}) dA \f]
  function findConvectTVDVect(phi,u,ubind,grid,grad,limiter)
    use moduleGrid
    double precision,intent(in)::phi(:,:) !< variable to be convected
    double precision,intent(in)::u(:,:) !< velocity which drives the convection
    integer,intent(in)::ubind !< bind velocity with node/block/interface
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::grad(DIMS,size(phi,1),size(phi,2)) !< gradient of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findConvectTVDVect(size(phi,1),grid%nBlock) !< increment due to convection
    double precision F,phiU(size(phi,1)),phiD(size(phi,1)),r,flowRate
    double precision,allocatable::uIntf(:,:)
    procedure(modelLimiter),pointer::lim
    
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    call grid%updateBlockPos()
    findConvectTVDVect(:,:)=0d0
    allocate(uIntf(DIMS,grid%nIntf))
    select case(ubind)
    case(BIND_NODE)
      !TODO:interpolation from node to interface
    case(BIND_BLOCK)
      !TODO:interpolation for co-located u
    case(BIND_INTF)
      uIntf(:,:)=u(:,:)
    case default
    end select
    if(present(limiter))then
      lim=>limiter
    else
      lim=>vanLeerLimiter
    end if
    do i=1,grid%nIntf
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      F=dot_product(grid%IntfNorm(:,i),uIntf(:,i))*grid%IntfArea(i)
      if(F>0d0)then
        phiU(:)=phi(:,m)
        phiD(:)=phi(:,n)
      else
        phiU(:)=phi(:,n)
        phiD(:)=phi(:,m)
      end if
      do j=1,size(phi,1)
        if(abs(phiD(j)-phiU(j))>tiny(1d0))then
          r=2d0*dot_product(grad(:,j,m),grid%BlockPos(:,n)-grid%BlockPos(:,m))&
          &    /(phiD(j)-phiU(j))-1d0
          flowRate=-F*(phiU(j)+lim(r)*(phiD(j)-phiU(j))/2d0)
        else
          flowRate=-F*phiU(j)
        end if
        findConvectTVDVect(j,m)=findConvectTVDVect(j,m)+flowRate
        findConvectTVDVect(j,n)=findConvectTVDVect(j,n)-flowRate
      end do
    end do
  end function
  
  !> find the convection of scalar phi driven by velocity u through interface using TVD scheme
  !> \f[ \int_{intf} \Phi (\mathbf{u} \cdot \hat{n}) dA \f]
  function findConvectTVDScal(phi,u,ubind,grid,grad,limiter)
    use moduleGrid
    double precision,intent(in)::phi(:) !< variable to be convected
    double precision,intent(in)::u(:,:) !< velocity which drives the convection
    integer,intent(in)::ubind !< bind velocity with node/block/interface
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::grad(DIMS,size(phi)) !< gradient of phi
    procedure(modelLimiter),pointer,intent(in),optional::limiter !< the flux limiter
    double precision findConvectTVDScal(grid%nBlock) !< increment due to convection
    double precision vphi(1,size(phi)),vgrad(DIMS,1,size(phi)),vrst(1,grid%nBlock)
    
    vphi(1,:)=phi(:)
    vgrad(:,1,:)=grad(:,:)
    if(present(limiter))then
      vrst=findConvectTVDVect(vphi,u,ubind,grid,vgrad,limiter=limiter)
    else
      vrst=findConvectTVDVect(vphi,u,ubind,grid,vgrad)
    end if
    findConvectTVDScal(:)=vrst(1,:)
  end function
  
  !> find the convection of vector phi due to the displacement of block surface using upwind scheme
  !> \f[ \int_{V^{block}_{disp}} \mathbf{\Phi} dV \f]
  function findDispConvectUpWindVect(phi,disp,grid)
    use moduleGrid
    use moduleInterpolation
    use moduleGridOperation
    double precision,intent(in)::phi(:,:) !< variable to be convected
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision findDispConvectUpWindVect(size(phi,1),grid%nBlock) !< increment of phi
    double precision dVolFracMax,flowRate(size(phi,1)),phiTemp(size(phi,1),grid%nBlock)
    double precision,allocatable::dVol(:),dispIntf(:,:)
    double precision,parameter::LIM_DVOL_FRAC=0.2d0
    
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    call grid%updateBlockVol()
    findDispConvectUpWindVect(:,:)=0d0
    allocate(dVol(grid%nIntf))
    dispIntf=itplNode2Intf(disp,grid)
    forall(i=1:grid%nIntf)
      dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
    end forall
    dVolFracMax=0d0
    do i=1,grid%nIntf
      dVolFracMax=max(dVolFracMax,abs(dVol(i)/minval(grid%BlockVol(grid%IntfNeibBlock(:,i)))))
    end do
    k=ceiling(dVolFracMax/LIM_DVOL_FRAC) ! number of sub-cycle steps
    dispIntf(:,:)=dispIntf(:,:)/dble(k)
    phiTemp(:,:)=phi(:,:)
    do l=1,k
      call grid%updateIntfArea()
      call grid%updateIntfNorm()
      call grid%updateBlockVol()
      forall(i=1:grid%nIntf)
        dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
      end forall
      do i=1,grid%nIntf
        m=grid%IntfNeibBlock(1,i)
        n=grid%IntfNeibBlock(2,i)
        if(dVol(i)>0d0)then
          flowRate(:)=dVol(i)*phiTemp(:,n)
        else
          flowRate(:)=dVol(i)*phiTemp(:,m)
        end if
        findDispConvectUpWindVect(:,m)=findDispConvectUpWindVect(:,m)+flowRate(:)
        findDispConvectUpWindVect(:,n)=findDispConvectUpWindVect(:,n)-flowRate(:)
        phiTemp(:,m)=phiTemp(:,m)+flowRate(:)/grid%BlockVol(m)
        phiTemp(:,n)=phiTemp(:,n)-flowRate(:)/grid%BlockVol(n)
      end do
      call mvGrid(grid,disp/dble(k))
    end do
    call mvGrid(grid,-disp)
    deallocate(dVol)
    deallocate(dispIntf)
  end function
  
  !> find the convection of scalar phi due to the displacement of block surface using upwind scheme
  !> \f[ \int_{V^{block}_{disp}} \Phi dV \f]
  function findDispConvectUpWindScal(phi,disp,grid)
    use moduleGrid
    double precision,intent(in)::phi(:) !< variable to be convected
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision findDispConvectUpWindScal(grid%nBlock) !< increment of phi
    double precision vphi(1,size(phi)),vrst(1,grid%nBlock)
    
    vphi(1,:)=phi(:)
    vrst=findDispConvectUpWindVect(vphi,disp,grid)
    findDispConvectUpWindScal(:)=vrst(1,:)
  end function
  
  !> find the convection of vector phi due to the displacement of block surface using TVD scheme
  !> \f[ \int_{V^{block}_{disp}} \mathbf{\Phi} dV \f]
  function findDispConvectTVDVect(phi,disp,grid,grad,limiter)
    use moduleGrid
    use moduleGridOperation
    use moduleInterpolation
    use moduleFVMGrad
    double precision,intent(in)::phi(:,:) !< variable to be convected
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision,intent(in)::grad(DIMS,size(phi,1),size(phi,2)) !< gradient of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findDispConvectTVDVect(size(phi,1),grid%nBlock) !< increment of phi
    double precision dVolFracMax,phiTemp(size(phi,1),grid%nBlock),&
    &                gradTemp(DIMS,size(phi,1),grid%nBlock),phiU(size(phi,1)),phiD(size(phi,1)),r,&
    &                flowRate
    double precision,allocatable::dVol(:),dispIntf(:,:)
    double precision,parameter::LIM_DVOL_FRAC=0.2d0
    procedure(modelLimiter),pointer::lim
    
    call grid%updateIntfArea()
    call grid%updateIntfNorm()
    call grid%updateBlockPos()
    call grid%updateBlockVol()
    findDispConvectTVDVect(:,:)=0d0
    if(present(limiter))then
      lim=>limiter
    else
      lim=>vanLeerLimiter
    end if
    allocate(dVol(grid%nIntf))
    dispIntf=itplNode2Intf(disp,grid)
    forall(i=1:grid%nIntf)
      dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
    end forall
    dVolFracMax=0d0
    do i=1,grid%nIntf
      dVolFracMax=max(dVolFracMax,abs(dVol(i)/minval(grid%BlockVol(grid%IntfNeibBlock(:,i)))))
    end do
    k=ceiling(dVolFracMax/LIM_DVOL_FRAC) ! number of sub-cycle steps
    dispIntf(:,:)=dispIntf(:,:)/dble(k)
    phiTemp(:,:)=phi(:,:)
    gradTemp(:,:,:)=grad(:,:,:)
    do l=1,k
      call grid%updateIntfArea()
      call grid%updateIntfNorm()
      call grid%updateBlockPos()
      call grid%updateBlockVol()
      if(l/=1)then
        gradTemp=findGrad(phiTemp,grid,BIND_BLOCK)
      end if
      forall(i=1:grid%nIntf)
        dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
      end forall
      do i=1,grid%nIntf
        m=grid%IntfNeibBlock(1,i)
        n=grid%IntfNeibBlock(2,i)
        if(dVol(i)>0d0)then
          phiU(:)=phiTemp(:,n)
          phiD(:)=phiTemp(:,m)
        else
          phiU(:)=phiTemp(:,m)
          phiD(:)=phiTemp(:,n)
        end if
        do j=1,size(phi,1)
          if(abs(phiD(j)-phiU(j))>tiny(1d0))then
            r=2d0*dot_product(gradTemp(:,j,n),grid%BlockPos(:,m)-grid%BlockPos(:,n))&
            &    /(phiD(j)-phiU(j))-1d0
            flowRate=dVol(i)*(phiU(j)+lim(r)*(phiD(j)-phiU(j))/2d0)
          else
            flowRate=dVol(i)*phiU(j)
          end if
          findDispConvectTVDVect(j,m)=findDispConvectTVDVect(j,m)+flowRate
          findDispConvectTVDVect(j,n)=findDispConvectTVDVect(j,n)-flowRate
          phiTemp(j,m)=phiTemp(j,m)+flowRate/grid%BlockVol(m)
          phiTemp(j,n)=phiTemp(j,n)-flowRate/grid%BlockVol(n)
        end do
      end do
      call mvGrid(grid,disp/dble(k))
    end do
    call mvGrid(grid,-disp)
    deallocate(dVol)
    deallocate(dispIntf)
  end function
  
  !> find the convection of scalar phi due to the displacement of block surface using TVD scheme
  !> \f[ \int_{V^{block}_{disp}} \Phi dV \f]
  function findDispConvectTVDScal(phi,disp,grid,grad,limiter)
    use moduleGrid
    double precision,intent(in)::phi(:) !< variable to be convected
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision,intent(in)::grad(DIMS,size(phi)) !< gradient of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findDispConvectTVDScal(grid%nBlock) !< increment of phi
    double precision vphi(1,size(phi)),vgrad(DIMS,1,size(phi)),vrst(1,grid%nBlock)
    
    vphi(1,:)=phi(:)
    vgrad(:,1,:)=grad(:,:)
    if(present(limiter))then
      vrst=findDispConvectTVDVect(vphi,disp,grid,vgrad,limiter=limiter)
    else
      vrst=findDispConvectTVDVect(vphi,disp,grid,vgrad)
    end if
    findDispConvectTVDScal(:)=vrst(1,:)
  end function
  
end module
