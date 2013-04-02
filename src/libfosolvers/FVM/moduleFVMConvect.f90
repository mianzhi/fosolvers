!----------------------------------------------------------------------------- best with 100 columns

!> FVM convection schemes
module moduleFVMConvect
  private
  
  !> find convection through interface
  interface findConvect
    ! schemes for 3-D unstructured grid
    module procedure findConvectUpWindScal
    module procedure findConvectUpWindVect
    module procedure findConvectTVDScal
    module procedure findConvectTVDVect
    ! schemes for 1-D grid
    module procedure findConvect1DTVDScal
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
    double precision F,r,flowRate
    double precision,allocatable::uIntf(:,:)
    integer up,dn
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
      F=dot_product(grid%IntfNorm(:,i),uIntf(:,i))*grid%IntfArea(i)
      m=grid%IntfNeibBlock(1,i)
      n=grid%IntfNeibBlock(2,i)
      if(F>0d0)then
        up=m
        dn=n
      else
        up=n
        dn=m
      end if
      do j=1,size(phi,1)
        if(abs(phi(j,up)-phi(j,dn))>tiny(1d0))then
          r=2d0*dot_product(grad(:,j,up),grid%BlockPos(:,dn)-grid%BlockPos(:,up))&
          &    /(phi(j,dn)-phi(j,up))-1d0
          flowRate=-F*(phi(j,up)+lim(r)*(phi(j,dn)-phi(j,up))/2d0)
        else
          flowRate=-F*phi(j,up)
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
  
  !> find the convection of vector phi driven by velocity u through node using TVD scheme
  !> \f[ u\mathbf{\Phi} \f]
  function findConvect1DTVDVect(phi,u,ubind,grid,phix,limiter)
    use moduleGrid1D
    use moduleInterpolation
    double precision,intent(in)::phi(:,:) !< variable to be convected
    double precision,intent(in)::u(:) !< velocity which drives the convection
    integer,intent(in)::ubind !< bind velocity with node/cell
    type(typeGrid1D),intent(inout)::grid !< the grid
    double precision,intent(in)::phix(size(phi,1),size(phi,2)) !< x derivative of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findConvect1DTVDVect(size(phi,1),grid%nCell) !< increment due to convection
    double precision r,flowRate
    double precision,allocatable::uNode(:)
    integer up,dn
    procedure(modelLimiter),pointer::lim
    
    findConvect1DTVDVect(:,:)=0d0
    allocate(uNode(grid%nNode))
    select case(ubind)
    case(BIND_NODE)
      uNode(:)=u(:)
    case(BIND_CELL)
      uNode=itplCell2Node(u,grid)
    case default
    end select
    if(present(limiter))then
      lim=>limiter
    else
      lim=>vanLeerLimiter
    end if
    do i=2,grid%nNode-1
      m=NlC(i)
      n=NrC(i)
      if(uNode(i)>0d0)then
        up=m
        dn=n
      else
        up=n
        dn=m
      end if
      do j=1,size(phi,1)
        if(abs(phi(j,up)-phi(j,dn))>tiny(1d0))then
          r=2d0*phix(j,up)*(grid%CellPos(dn)-grid%CellPos(up))&
          &    /(phi(j,dn)-phi(j,up))-1d0
          flowRate=-uNode(i)*(phi(j,up)+lim(r)*(phi(j,dn)-phi(j,up))/2d0)
        else
          flowRate=-uNode(i)*phi(j,up)
        end if
        findConvect1DTVDVect(j,m)=findConvect1DTVDVect(j,m)+flowRate
        findConvect1DTVDVect(j,n)=findConvect1DTVDVect(j,n)-flowRate
      end do
    end do
  end function
  
  !> find the convection of scalar phi driven by velocity u through node using TVD scheme
  !> \f[ u\Phi \f]
  function findConvect1DTVDScal(phi,u,ubind,grid,phix,limiter)
    use moduleGrid1D
    double precision,intent(in)::phi(:) !< variable to be convected
    double precision,intent(in)::u(:) !< velocity which drives the convection
    integer,intent(in)::ubind !< bind velocity with node/cell
    type(typeGrid1D),intent(inout)::grid !< the grid
    double precision,intent(in)::phix(size(phi)) !< x derivative of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findConvect1DTVDScal(grid%nCell) !< increment due to convection
    double precision vphi(1,size(phi)),vphix(1,size(phi)),vrst(1,grid%nCell)
    
    vphi(1,:)=phi(:)
    vphix(1,:)=phix(:)
    if(present(limiter))then
      vrst=findConvect1DTVDVect(vphi,u,ubind,grid,vphix,limiter=limiter)
    else
      vrst=findConvect1DTVDVect(vphi,u,ubind,grid,vphix)
    end if
    findConvect1DTVDScal(:)=vrst(1,:)
  end function
  
  !> find the convection of vector phi due to the displacement of CV surface using upwind scheme
  !> \f[ \int_{V^{CV}_{disp}} \mathbf{\Phi} dV \f]
  function findDispConvectUpWindVect(phi,bind,disp,grid)
    use moduleGrid
    use moduleInterpolation
    use moduleGridOperation
    double precision,intent(in)::phi(:,:) !< variable to be convected
    integer,intent(in)::bind !< phi bind with node/block
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision findDispConvectUpWindVect(size(phi,1),size(phi,2)) !< increment of phi
    double precision dVolFracMax,flowRate(size(phi,1)),phiTemp1(size(phi,1),size(phi,2)),&
    &                phiTemp2(size(phi,1),size(phi,2))
    double precision,allocatable::dVol(:),dispIntf(:,:)
    double precision,parameter::LIM_DVOL_FRAC=0.2d0
    
    findDispConvectUpWindVect(:,:)=0d0
    dVolFracMax=0d0
    select case(bind)
    case(BIND_NODE)
      call grid%updateDualBlock()
      allocate(dVol(grid%nEdge))
      dispIntf=itplNode2Edge(disp,grid)
      forall(i=1:grid%nEdge)
        dVol(i)=dot_product(grid%EAreaVect(:,i),dispIntf(:,i))
      end forall
      do i=1,grid%nEdge
        dVolFracMax=max(dVolFracMax,abs(dVol(i)/minval(grid%NodeVol(grid%Edge(i)%iNode(:)))))
      end do
      k=ceiling(dVolFracMax/LIM_DVOL_FRAC) ! number of sub-cycle steps
      dispIntf(:,:)=dispIntf(:,:)/dble(k)
      forall(i=1:grid%nNode)
        phiTemp1(:,i)=phi(:,i)*grid%NodeVol(i)
      end forall
      do l=1,k
        call grid%updateDualBlock()
        forall(i=1:grid%nEdge)
          dVol(i)=dot_product(grid%EAreaVect(:,i),dispIntf(:,i))
        end forall
        forall(i=1:grid%nNode)
          phiTemp2(:,i)=phiTemp1(:,i)/grid%NodeVol(i)
        end forall
        do i=1,grid%nEdge
          m=grid%Edge(i)%iNode(1)
          n=grid%Edge(i)%iNode(2)
          if(dVol(i)>0d0)then
            flowRate(:)=dVol(i)*phiTemp2(:,n)
          else
            flowRate(:)=dVol(i)*phiTemp2(:,m)
          end if
          findDispConvectUpWindVect(:,m)=findDispConvectUpWindVect(:,m)+flowRate(:)
          findDispConvectUpWindVect(:,n)=findDispConvectUpWindVect(:,n)-flowRate(:)
          phiTemp1(:,m)=phiTemp1(:,m)+flowRate(:)
          phiTemp1(:,n)=phiTemp1(:,n)-flowRate(:)
        end do
        call mvGrid(grid,disp/dble(k))
      end do
      call mvGrid(grid,-disp)
      deallocate(dVol)
      deallocate(dispIntf)
    case(BIND_BLOCK)
      call grid%updateIntfArea()
      call grid%updateIntfNorm()
      call grid%updateBlockVol()
      allocate(dVol(grid%nIntf))
      dispIntf=itplNode2Intf(disp,grid)
      forall(i=1:grid%nIntf)
        dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
      end forall
      do i=1,grid%nIntf
        dVolFracMax=max(dVolFracMax,abs(dVol(i)/minval(grid%BlockVol(grid%IntfNeibBlock(:,i)))))
      end do
      k=ceiling(dVolFracMax/LIM_DVOL_FRAC) ! number of sub-cycle steps
      dispIntf(:,:)=dispIntf(:,:)/dble(k)
      forall(i=1:grid%nBlock)
        phiTemp1(:,i)=phi(:,i)*grid%BlockVol(i)
      end forall
      do l=1,k
        call grid%updateIntfArea()
        call grid%updateIntfNorm()
        call grid%updateBlockVol()
        forall(i=1:grid%nIntf)
          dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
        end forall
        forall(i=1:grid%nBlock)
          phiTemp2(:,i)=phiTemp1(:,i)/grid%BlockVol(i)
        end forall
        do i=1,grid%nIntf
          m=grid%IntfNeibBlock(1,i)
          n=grid%IntfNeibBlock(2,i)
          if(dVol(i)>0d0)then
            flowRate(:)=dVol(i)*phiTemp2(:,n)
          else
            flowRate(:)=dVol(i)*phiTemp2(:,m)
          end if
          findDispConvectUpWindVect(:,m)=findDispConvectUpWindVect(:,m)+flowRate(:)
          findDispConvectUpWindVect(:,n)=findDispConvectUpWindVect(:,n)-flowRate(:)
          phiTemp1(:,m)=phiTemp1(:,m)+flowRate(:)
          phiTemp1(:,n)=phiTemp1(:,n)-flowRate(:)
        end do
        call mvGrid(grid,disp/dble(k))
      end do
      call mvGrid(grid,-disp)
      deallocate(dVol)
      deallocate(dispIntf)
    case default
    end select
  end function
  
  !> find the convection of scalar phi due to the displacement of CV surface using upwind scheme
  !> \f[ \int_{V^{CV}_{disp}} \Phi dV \f]
  function findDispConvectUpWindScal(phi,bind,disp,grid)
    use moduleGrid
    double precision,intent(in)::phi(:) !< variable to be convected
    integer,intent(in)::bind !< phi bind with node/block
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision findDispConvectUpWindScal(size(phi)) !< increment of phi
    double precision vphi(1,size(phi)),vrst(1,size(phi))
    
    vphi(1,:)=phi(:)
    vrst=findDispConvectUpWindVect(vphi,bind,disp,grid)
    findDispConvectUpWindScal(:)=vrst(1,:)
  end function
  
  !> find the convection of vector phi due to the displacement of CV surface using TVD scheme
  !> \f[ \int_{V^{CV}_{disp}} \mathbf{\Phi} dV \f]
  function findDispConvectTVDVect(phi,bind,disp,grid,grad,limiter)
    use moduleGrid
    use moduleGridOperation
    use moduleInterpolation
    use moduleFVMGrad
    double precision,intent(in)::phi(:,:) !< variable to be convected
    integer,intent(in)::bind !< phi bind with node/block
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision,intent(in)::grad(DIMS,size(phi,1),size(phi,2)) !< gradient of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findDispConvectTVDVect(size(phi,1),size(phi,2)) !< increment of phi
    double precision dVolFracMax,r,flowRate,gradTemp(DIMS,size(phi,1),size(phi,2)),&
    &                phiTemp1(size(phi,1),size(phi,2)),phiTemp2(size(phi,1),size(phi,2))
    double precision,allocatable::dVol(:),dispIntf(:,:)
    double precision,parameter::LIM_DVOL_FRAC=0.2d0
    procedure(modelLimiter),pointer::lim
    integer up,dn
    
    findDispConvectTVDVect(:,:)=0d0
    dVolFracMax=0d0
    if(present(limiter))then
      lim=>limiter
    else
      lim=>vanLeerLimiter
    end if
    select case(bind)
    case(BIND_NODE)
      call grid%updateDualBlock()
      allocate(dVol(grid%nEdge))
      dispIntf=itplNode2Edge(disp,grid)
      forall(i=1:grid%nEdge)
        dVol(i)=dot_product(grid%EAreaVect(:,i),dispIntf(:,i))
      end forall
      do i=1,grid%nEdge
        dVolFracMax=max(dVolFracMax,abs(dVol(i)/minval(grid%NodeVol(grid%Edge(i)%iNode(:)))))
      end do
      k=ceiling(dVolFracMax/LIM_DVOL_FRAC) ! number of sub-cycle steps
      dispIntf(:,:)=dispIntf(:,:)/dble(k)
      forall(i=1:grid%nNode)
        phiTemp1(:,i)=phi(:,i)*grid%NodeVol(i)
      end forall
      gradTemp(:,:,:)=grad(:,:,:)
      do l=1,k
        call grid%updateDualBlock()
        forall(i=1:grid%nEdge)
          dVol(i)=dot_product(grid%EAreaVect(:,i),dispIntf(:,i))
        end forall
        forall(i=1:grid%nNode)
          phiTemp2(:,i)=phiTemp1(:,i)/grid%NodeVol(i)
        end forall
        if(l/=1)then
          gradTemp=findGrad(phiTemp2,BIND_NODE,grid)
        end if
        do i=1,grid%nEdge
          m=grid%Edge(i)%iNode(1)
          n=grid%Edge(i)%iNode(2)
          if(dVol(i)>0d0)then
            up=n
            dn=m
          else
            up=m
            dn=n
          end if
          do j=1,size(phi,1)
            if(abs(phiTemp2(j,up)-phiTemp2(j,dn))>tiny(1d0))then
              r=2d0*dot_product(gradTemp(:,j,up),grid%NodePos(:,dn)-grid%NodePos(:,up))&
              &    /(phiTemp2(j,dn)-phiTemp2(j,up))-1d0
              flowRate=dVol(i)*(phiTemp2(j,up)+lim(r)*(phiTemp2(j,dn)-phiTemp2(j,up))/2d0)
            else
              flowRate=dVol(i)*phiTemp2(j,up)
            end if
            findDispConvectTVDVect(j,m)=findDispConvectTVDVect(j,m)+flowRate
            findDispConvectTVDVect(j,n)=findDispConvectTVDVect(j,n)-flowRate
            phiTemp1(j,m)=phiTemp1(j,m)+flowRate
            phiTemp1(j,n)=phiTemp1(j,n)-flowRate
          end do
        end do
        call mvGrid(grid,disp/dble(k))
      end do
      call mvGrid(grid,-disp)
      deallocate(dVol)
      deallocate(dispIntf)
    case(BIND_BLOCK)
      call grid%updateIntfArea()
      call grid%updateIntfNorm()
      call grid%updateBlockPos()
      call grid%updateBlockVol()
      allocate(dVol(grid%nIntf))
      dispIntf=itplNode2Intf(disp,grid)
      forall(i=1:grid%nIntf)
        dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
      end forall
      do i=1,grid%nIntf
        dVolFracMax=max(dVolFracMax,abs(dVol(i)/minval(grid%BlockVol(grid%IntfNeibBlock(:,i)))))
      end do
      k=ceiling(dVolFracMax/LIM_DVOL_FRAC) ! number of sub-cycle steps
      dispIntf(:,:)=dispIntf(:,:)/dble(k)
      forall(i=1:grid%nBlock)
        phiTemp1(:,i)=phi(:,i)*grid%BlockVol(i)
      end forall
      gradTemp(:,:,:)=grad(:,:,:)
      do l=1,k
        call grid%updateIntfArea()
        call grid%updateIntfNorm()
        call grid%updateBlockPos()
        call grid%updateBlockVol()
        forall(i=1:grid%nIntf)
          dVol(i)=grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),dispIntf(:,i))
        end forall
        forall(i=1:grid%nBlock)
          phiTemp2(:,i)=phiTemp1(:,i)/grid%BlockVol(i)
        end forall
        if(l/=1)then
          gradTemp=findGrad(phiTemp2,BIND_BLOCK,grid)
        end if
        do i=1,grid%nIntf
          m=grid%IntfNeibBlock(1,i)
          n=grid%IntfNeibBlock(2,i)
          if(dVol(i)>0d0)then
            up=n
            dn=m
          else
            up=m
            dn=n
          end if
          do j=1,size(phi,1)
            if(abs(phiTemp2(j,up)-phiTemp2(j,dn))>tiny(1d0))then
              r=2d0*dot_product(gradTemp(:,j,up),grid%BlockPos(:,dn)-grid%BlockPos(:,up))&
              &    /(phiTemp2(j,dn)-phiTemp2(j,up))-1d0
              flowRate=dVol(i)*(phiTemp2(j,up)+lim(r)*(phiTemp2(j,dn)-phiTemp2(j,up))/2d0)
            else
              flowRate=dVol(i)*phiTemp2(j,up)
            end if
            findDispConvectTVDVect(j,m)=findDispConvectTVDVect(j,m)+flowRate
            findDispConvectTVDVect(j,n)=findDispConvectTVDVect(j,n)-flowRate
            phiTemp1(j,m)=phiTemp1(j,m)+flowRate
            phiTemp1(j,n)=phiTemp1(j,n)-flowRate
          end do
        end do
        call mvGrid(grid,disp/dble(k))
      end do
      call mvGrid(grid,-disp)
      deallocate(dVol)
      deallocate(dispIntf)
    case default
    end select
  end function
  
  !> find the convection of scalar phi due to the displacement of CV surface using TVD scheme
  !> \f[ \int_{V^{CV}_{disp}} \Phi dV \f]
  function findDispConvectTVDScal(phi,bind,disp,grid,grad,limiter)
    use moduleGrid
    double precision,intent(in)::phi(:) !< variable to be convected
    integer,intent(in)::bind !< phi bind with node/block
    type(typeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::disp(DIMS,grid%nNode) !< node displacement
    double precision,intent(in)::grad(DIMS,size(phi)) !< gradient of phi
    procedure(modelLimiter),pointer,optional::limiter !< the flux limiter
    double precision findDispConvectTVDScal(size(phi)) !< increment of phi
    double precision vphi(1,size(phi)),vgrad(DIMS,1,size(phi)),vrst(1,size(phi))
    
    vphi(1,:)=phi(:)
    vgrad(:,1,:)=grad(:,:)
    if(present(limiter))then
      vrst=findDispConvectTVDVect(vphi,bind,disp,grid,vgrad,limiter=limiter)
    else
      vrst=findDispConvectTVDVect(vphi,bind,disp,grid,vgrad)
    end if
    findDispConvectTVDScal(:)=vrst(1,:)
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
  
end module
