!----------------------------------------------------------------------------- best with 100 columns

!> advection for FVM
module modAdvection
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding mass flow rate through pairs
  interface findMassFlow
    module procedure::findMassFlowMagicPoly
  end interface
  public::findMassFlow
  
  !> generic finding variable flow rate through pairs using mass flow rate
  interface findVarFlow
    module procedure::findVarFlowPolyVect
    module procedure::findVarFlowPolyScal
  end interface
  public::findVarFlow
  
  !> generic finding advection rate into each cell
  interface findAdv
    module procedure::findAdvFlowPolyVect
    module procedure::findAdvFlowPolyScal
    module procedure::findAdvPolyVect
    module procedure::findAdvPolyScal
    module procedure::findAdvPolyEuler
    module procedure::findAdvOtVect
    module procedure::findAdvOtScal
  end interface
  public::findAdv
  
  !> generic finding advection Jacobian
  interface findAdvJac
    module procedure::findAdvJacPolyEuler
  end interface
  public::findAdvJac
  
contains
  
  !> find mass flow rate on polyFvGrid with Rhie-Chow u interpolation and TVD rho reconstruction
  !> \f[ \int_A \rho\mathbf{u} \cdot \hat{n} dA \f]
  subroutine findMassFlowMagicPoly(grid,rho,u,p,presF,dt,flow)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::rho(:) !< density
    double precision,intent(in)::u(:,:) !< velocity
    double precision,intent(in)::p(:) !< pressure
    double precision,intent(in)::presF(:,:) !< pressure force
    double precision,intent(in)::dt !< time step size
    double precision,allocatable,intent(inout)::flow(:) !< mass flow rate output
    integer::up,dn
    double precision::vMN(DIMS),vMP(DIMS),vPN(DIMS),flux(DIMS),eps,rhoUp,rhoDn,dRho,r,rhof
    double precision,allocatable::gradRho(:,:)
    
    call grid%up()
    if(.not.allocated(flow))then
      allocate(flow(grid%nP))
    end if
    flow(:)=0d0
    allocate(gradRho(DIMS,grid%nC))
    call findGrad(grid,rho,gradRho)
    !$omp parallel do default(shared)&
    !$omp& private(m,n,up,dn,rhoUp,rhoDn,dRho,r,rhof,vMN,vMP,vPN,eps,flux)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(dot_product(u(:,m)+u(:,n),grid%normP(:,i))>=0d0)then ! upwind by velocity
        up=m
        dn=n
      else
        up=n
        dn=m
      end if
      rhoUp=rho(up)
      rhoDn=rho(dn)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        dRho=dot_product(grid%p(:,dn)-grid%p(:,up),gradRho(:,up))
        r=merge(2d0*dRho/(rhoDn-rhoUp)-1d0,0d0,abs(rhoDn-rhoUp)>tiny(1d0))
        rhof=rhoUp+0.5d0*vanAlbada(r)*(rhoDn-rhoUp)
        vMN(:)=grid%p(:,n)-grid%p(:,m)
        vMP(:)=grid%pP(:,i)-grid%p(:,m)
        vPN(:)=grid%p(:,n)-grid%pP(:,i)
        eps=norm2(vPN)/(norm2(vMP)+norm2(vPN))
        flux(:)=rhof*(eps*u(:,m)+(1d0-eps)*u(:,n))&
        &       -dt/(1d0/rho(m)+1d0/rho(n))*&
        &        (presF(:,m)/grid%v(m)/rho(m)+presF(:,n)/grid%v(n)/rho(n))&
        &       +dt*vMN(:)*(p(m)-p(n))/dot_product(vMN,vMN)
      else ! boundary pairs
        flux(:)=rho(up)*0.5d0*(u(:,m)+u(:,n))
        if(abs(dot_product(grid%normP(:,i),flux(:)))>tiny(1d0))then
          vMP(:)=grid%pP(:,i)-grid%p(:,m)
          flux(:)=flux(:)-0.5d0*dt*presF(:,m)/grid%v(m)&
          &       +0.5d0*dt*grid%normP(:,i)*(p(m)-p(n))/(2d0*dot_product(vMP,grid%normP(:,i)))
        end if ! Rhie-Chow does not initiate penetration through boundary pair
      end if
      flow(i)=grid%aP(i)*dot_product(grid%normP(:,i),flux(:))
    end do
    !$omp end parallel do
    deallocate(gradRho)
  end subroutine
  
  !> find flow rate of vector v on polyFvGrid using mass flow rate
  !> \f[ \int_A \rho\mathbf{v}(\mathbf{u} \cdot \hat{n}) dA \f]
  subroutine findVarFlowPolyVect(grid,v,mFlow,flow)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:,:) !< transported variables
    double precision,intent(in)::mFlow(:) !< mass flow through pairs
    double precision,allocatable,intent(inout)::flow(:,:) !< v flow rate output
    double precision,allocatable::gradV(:,:,:)
    double precision::vf(size(v,1)),dV(size(v,1)),r(size(v,1))
    integer::up,dn
    
    call grid%up()
    if(.not.allocated(flow))then
      allocate(flow(size(v,1),grid%nP))
    end if
    flow(:,:)=0d0
    allocate(gradV(DIMS,size(v,1),grid%nC))
    call findGrad(grid,v,gradV)
    !$omp parallel do default(shared)&
    !$omp& private(m,n,up,dn,dV,r,vf)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(v,2).and.n<=size(v,2))then
        if(abs(mFlow(i))<=tiny(1d0))then ! no mass flux
          cycle
        else ! upwinding by mass flux
          if(mFlow(i)>=0d0)then
            up=m
            dn=n
          else
            up=n
            dn=m
          end if
          if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
            dV(:)=matmul(grid%p(:,dn)-grid%p(:,up),gradV(:,:,up))
            r(:)=merge(2d0*dV/(v(:,dn)-v(:,up))-1d0,0d0,abs(v(:,dn)-v(:,up))>tiny(1d0))
            vf(:)=v(:,up)+0.5d0*vanAlbada(r)*(v(:,dn)-v(:,up))
          else ! boundary pairs
            vf(:)=v(:,up)
          end if
          flow(:,i)=vf(:)*mFlow(i)
        end if
      end if
    end do
    !$omp end parallel do
    deallocate(gradV)
  end subroutine
  
  !> find flow rate of scalar v on polyFvGrid using mass flow rate
  !> \f[ \int_A \rho v(\mathbf{u} \cdot \hat{n}) dA \f]
  subroutine findVarFlowPolyScal(grid,v,mFlow,flow)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:) !< transported variable
    double precision,intent(in)::mFlow(:) !< mass flow through pairs
    double precision,allocatable,intent(inout)::flow(:) !< v flow rate output
    double precision,allocatable::vv(:,:),flowv(:,:)
    
    allocate(vv(1,size(v)))
    if(allocated(flow))then
      allocate(flowv(1,size(flow)))
    end if
    vv(1,:)=v(:)
    call findVarFlowPolyVect(grid,vv,mFlow,flowv)
    if(.not.allocated(flow))then
      allocate(flow(size(flowv,2)),source=flowv(1,:))!FIXME:remove work-around
    else
      flow(:)=flowv(1,:)
    end if
    deallocate(vv)
    deallocate(flowv)
  end subroutine
  
  !> find advection according to the vector flow rate through pairs on polyFvGrid
  subroutine findAdvFlowPolyVect(grid,flow,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::flow(:,:) !< flow rates
    double precision,allocatable,intent(inout)::adv(:,:) !< advection output
    
    call grid%up()
    if(.not.allocated(adv))then
      allocate(adv(size(flow,1),grid%nC))
    end if
    adv(:,:)=0d0
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        adv(:,m)=adv(:,m)-flow(:,i)
        adv(:,n)=adv(:,n)+flow(:,i)
      else ! boundary pairs
        adv(:,m)=adv(:,m)-flow(:,i)
      end if
    end do
  end subroutine
  
  !> find advection according to the scalar flow rate through pairs on polyFvGrid
  subroutine findAdvFlowPolyScal(grid,flow,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::flow(:) !< flow rates
    double precision,allocatable,intent(inout)::adv(:) !< advection output
    
    call grid%up()
    if(.not.allocated(adv))then
      allocate(adv(grid%nC))
    end if
    adv(:)=0d0
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        adv(m)=adv(m)-flow(i)
        adv(n)=adv(n)+flow(i)
      else ! boundary pairs
        adv(m)=adv(m)-flow(i)
      end if
    end do
  end subroutine
  
  !> find advection due to flux f depending on vector s on polyFvGrid
  !> \f[ \int_A \mathbf{f}(\mathbf{s}) \cdot \hat{n} dA \f]
  subroutine findAdvPolyVect(grid,s,f,adv)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< state variables
    double precision,intent(in)::f(:,:,:) !< fluxes
    double precision,allocatable,intent(inout)::adv(:,:) !< advection output
    integer::up,dn
    double precision::flow,fUp,fDn,df,r
    double precision,allocatable::gradF(:,:,:)
    
    call grid%up()
    if(.not.allocated(adv))then
      allocate(adv(size(s,1),grid%nC))
    end if
    adv(:,:)=0d0
    allocate(gradF(DIMS,DIMS*size(s,1),grid%nC))
    call findGrad(grid,reshape(f(:,:,1:grid%nC),&
    &                          [size(f(:,:,1:grid%nC),1)*size(f(:,:,1:grid%nC),2),&
    &                           size(f(:,:,1:grid%nC),3)]),gradF)
    !$omp parallel do default(shared)&
    !$omp& private(m,n,up,dn,fUp,fDn,df,r,flow)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(s,2).and.m<=size(f,3).and.n<=size(s,2).and.n<=size(f,3))then
        do j=1,size(s,1)
          if(abs(dot_product(f(:,j,m)+f(:,j,n),grid%normP(:,i)))<=tiny(1d0))then ! canceling flux
            cycle
          else if(abs(s(j,m))>tiny(1d0).and.abs(s(j,n))>tiny(1d0))then ! upwinding by velocity
            if(dot_product(f(:,j,m)/s(j,m)+f(:,j,n)/s(j,n),grid%normP(:,i))>=0d0)then
              up=m
              dn=n
            else
              up=n
              dn=m
            end if
          else ! upwinding by flux
            if(dot_product(f(:,j,m)+f(:,j,n),grid%normP(:,i))>=0d0)then
              up=m
              dn=n
            else
              up=n
              dn=m
            end if
          end if
          fUp=dot_product(f(:,j,up),grid%normP(:,i))
          fDn=dot_product(f(:,j,dn),grid%normP(:,i))
          !$omp critical
          if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
            df=dot_product(matmul(grid%p(:,dn)-grid%p(:,up),gradF(:,j*DIMS-(DIMS-1):j*DIMS,up)),&
            &              grid%normP(:,i))
            r=merge(2d0*df/(fDn-fUp)-1d0,0d0,abs(fDn-fUp)>tiny(1d0))
            flow=-grid%aP(i)*(fUp+0.5d0*vanAlbada(r)*(fDn-fUp))
            adv(j,m)=adv(j,m)+flow
            adv(j,n)=adv(j,n)-flow
          else ! boundary pairs
            flow=-grid%aP(i)*fUp
            adv(j,m)=adv(j,m)+flow
          end if
          !$omp end critical
        end do
      end if
    end do
    !$omp end parallel do
    deallocate(gradF)
  end subroutine
  
  !> find advection due to flux f depending on scalar s on polyFvGrid
  !> \f[ \int_A \mathbf{f}(s) \cdot \hat{n} dA \f]
  subroutine findAdvPolyScal(grid,s,f,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< state variable
    double precision,intent(in)::f(:,:) !< flux
    double precision,allocatable,intent(inout)::adv(:) !< advection output
    double precision,allocatable::sv(:,:),fv(:,:,:),advv(:,:)
    
    allocate(sv(1,size(s)))
    allocate(fv(DIMS,1,size(s)))
    if(allocated(adv))then
      allocate(advv(1,size(adv)))
    end if
    sv(1,:)=s(:)
    fv(:,1,:)=f(:,:)
    call findAdvPolyVect(grid,sv,fv,advv)
    if(.not.allocated(adv))then
      allocate(adv(size(advv,2)),source=advv(1,:))!FIXME:remove work-around
    else
      adv(:)=advv(1,:)
    end if
    deallocate(sv)
    deallocate(fv)
    deallocate(advv)
  end subroutine
  
  !> find advection of Euler system on polyFvGrid
  subroutine findAdvPolyEuler(grid,rho,rhou,rhoE,p,gamm,dRho,dRhou,dRhoE)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::rho(:) !< cell-averaged rho
    double precision,intent(in)::rhou(:,:) !< cell-averaged rhou
    double precision,intent(in)::rhoE(:) !< cell-averaged rhoE
    double precision,intent(in)::p(:) !< cell-averaged pressure
    double precision,intent(in)::gamm !< ratio of heat capacities
    double precision,allocatable,intent(inout)::dRho(:) !< net flux of rho
    double precision,allocatable,intent(inout)::dRhou(:,:) !< net flux of rhou
    double precision,allocatable,intent(inout)::dRhoE(:) !< net flux of rhoE
    double precision,allocatable,save::u(:,:),H(:)
    double precision,allocatable,save::f1c(:,:),f2c(:,:,:),f3c(:,:)
    double precision::rhoAvg,uAvg(DIMS),HAvg,cAvg,uNormAvg,rhoJump,pJump,uJump(DIMS),uNormJump
    double precision::dFEntropy(DIMS+2),dFAcoustic1(DIMS+2),dFAcoustic2(DIMS+2),flow(DIMS+2)
    
    call grid%up()
    if(.not.allocated(dRho))then
      allocate(dRho(grid%nC))
    end if
    dRho(:)=0d0
    if(.not.allocated(dRhou))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    if(.not.allocated(dRhoE))then
      allocate(dRhoE(grid%nC))
    end if
    dRhoE(:)=0d0
    ! find auxiliary state and flux vectors in cell
    k=minval([size(rho),size(rhou,2),size(rhoE),size(p)])
    if(.not.allocated(u))then
      allocate(u(DIMS,k))
      allocate(H(k))
      allocate(f1c(DIMS,k))
      allocate(f2c(DIMS,DIMS,k))
      allocate(f3c(DIMS,k))
    else if(size(u,2)<k)then
      deallocate(u,H)
      deallocate(f1c,f2c,f3c)
      allocate(u(DIMS,k))
      allocate(H(k))
      allocate(f1c(DIMS,k))
      allocate(f2c(DIMS,DIMS,k))
      allocate(f3c(DIMS,k))
    end if
    !$omp parallel do default(shared)
    do i=1,k
      u(:,i)=rhou(:,i)/rho(i)
      H(i)=(rhoE(i)+p(i))/rho(i)
      f1c(:,i)=rhou(:,i)
      f2c(:,1,i)=rhou(1,i)*u(:,i)+[p(i),0d0,0d0]
      f2c(:,2,i)=rhou(2,i)*u(:,i)+[0d0,p(i),0d0]
      f2c(:,3,i)=rhou(3,i)*u(:,i)+[0d0,0d0,p(i)]
      f3c(:,i)=(rhoE(i)+p(i))*u(:,i)
    end do
    !$omp end parallel do
    ! Roe flux difference splitting
    !$omp parallel do default(shared)&
    !$omp& private(m,n,rhoAvg,uAvg,Havg,cAvg,uNormAvg,rhoJump,pJump,uJump,uNormJump,&
    !$omp&         dFEntropy,dFAcoustic1,dFAcoustic2,flow)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=k.and.n<=k)then
        rhoAvg=sqrt(rho(m)*rho(n))
        uAvg(:)=(sqrt(rho(m))*u(:,m)+sqrt(rho(n))*u(:,n))/(sqrt(rho(m))+sqrt(rho(n)))
        uNormAvg=dot_product(uAvg,grid%normP(:,i))
        HAvg=(sqrt(rho(m))*H(m)+sqrt(rho(n))*H(n))/(sqrt(rho(m))+sqrt(rho(n)))
        cAvg=sqrt((gamm-1d0)*(HAvg-0.5d0*(dot_product(uAvg,uAvg))))
        rhoJump=rho(n)-rho(m)
        pJump=p(n)-p(m)
        uJump(:)=u(:,n)-u(:,m)
        uNormJump=dot_product(uJump,grid%normP(:,i))
        dFEntropy(:)=abs(uNormAvg)*&
        &            ((rhoJump-pJump/cAvg**2)*&
        &             [1d0,uAvg(1),uAvg(2),uAvg(3),0.5d0*dot_product(uAvg,uAvg)]&
        &             +rhoAvg*&
        &             [0d0,&
        &              uJump(1)-grid%normP(1,i)*uNormJump,&
        &              uJump(2)-grid%normP(2,i)*uNormJump,&
        &              uJump(3)-grid%normP(3,i)*uNormJump,&
        &              dot_product(uAvg,uJump)-uNormAvg*uNormJump])
        dFAcoustic1(:)=abs(uNormAvg-cAvg)*(0.5d0*(pJump-rhoAvg*cAvg*uNormJump)/cAvg**2)*&
        &              [1d0,&
        &               uAvg(1)-grid%normP(1,i)*cAvg,&
        &               uAvg(2)-grid%normP(2,i)*cAvg,&
        &               uAvg(3)-grid%normP(3,i)*cAvg,&
        &               HAvg-uNormAvg*cAvg]
        dFAcoustic2(:)=abs(uNormAvg+cAvg)*(0.5d0*(pJump+rhoAvg*cAvg*uNormJump)/cAvg**2)*&
        &              [1d0,&
        &               uAvg(1)+grid%normP(1,i)*cAvg,&
        &               uAvg(2)+grid%normP(2,i)*cAvg,&
        &               uAvg(3)+grid%normP(3,i)*cAvg,&
        &               HAvg+uNormAvg*cAvg]
        flow(:)=grid%aP(i)*&
        &       ([dot_product(0.5d0*(f1c(:,m)+f1c(:,n)),grid%normP(:,i)),&
        &         dot_product(0.5d0*(f2c(:,1,m)+f2c(:,1,n)),grid%normP(:,i)),&
        &         dot_product(0.5d0*(f2c(:,2,m)+f2c(:,2,n)),grid%normP(:,i)),&
        &         dot_product(0.5d0*(f2c(:,3,m)+f2c(:,3,n)),grid%normP(:,i)),&
        &         dot_product(0.5d0*(f3c(:,m)+f3c(:,n)),grid%normP(:,i))]&
        &        -0.5d0*(dFEntropy(:)+dFAcoustic1(:)+dFAcoustic2(:)))
        !$omp critical
        if(n<=grid%nC)then
          dRho(m)=dRho(m)-flow(1)
          dRho(n)=dRho(n)+flow(1)
          dRhou(:,m)=dRhou(:,m)-flow(2:4)
          dRhou(:,n)=dRhou(:,n)+flow(2:4)
          dRhoE(m)=dRhoE(m)-flow(5)
          dRhoE(n)=dRhoE(n)+flow(5)
        else
          dRho(m)=dRho(m)-flow(1)
          dRhou(:,m)=dRhou(:,m)-flow(2:4)
          dRhoE(m)=dRhoE(m)-flow(5)
        end if
        !$omp end critical
      end if
    end do
    !$end omp parallel do
  end subroutine
  
  !> find advection Jacobian of Euler system on polyFvGrid
  subroutine findAdvJacPolyEuler(grid,rho,rhou,rhoE,p,gamm,JacP,JacC)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::rho(:) !< cell-averaged rho
    double precision,intent(in)::rhou(:,:) !< cell-averaged rhou
    double precision,intent(in)::rhoE(:) !< cell-averaged rhoE
    double precision,intent(in)::p(:) !< cell-averaged pressure
    double precision,intent(in)::gamm !< ratio of heat capacities
    double precision,allocatable,intent(inout)::JacP(:,:,:) !< Jacobian for each pair
    double precision,allocatable,intent(inout)::JacC(:,:,:) !< Jacobian for each cell
    double precision,allocatable,save::u(:,:),H(:),c(:)
    double precision::uAvg(DIMS),HAvg,cAvg
    
    call grid%up()
    if(.not.allocated(JacP))then
      allocate(JacP(10,10,grid%nP))
    end if
    if(.not.allocated(JacC))then
      allocate(JacC(5,size(grid%neib,1)*5+5,grid%nC))
    end if
    JacP(:,:,:)=0d0
    JacC(:,:,:)=0d0
    ! find auxiliary state
    k=minval([size(rho),size(rhou,2),size(rhoE),size(p)])
    if(.not.allocated(u))then
      allocate(u(DIMS,k))
      allocate(H(k))
      allocate(c(k))
    else if(size(u,2)<k)then
      deallocate(u,H,c)
      allocate(u(DIMS,k))
      allocate(H(k))
      allocate(c(k))
    end if
    !$omp parallel do default(shared)
    do i=1,k
      u(:,i)=rhou(:,i)/rho(i)
      H(i)=(rhoE(i)+p(i))/rho(i)
      c(i)=sqrt((gamm-1d0)*(H(i)-0.5d0*(dot_product(u(:,i),u(:,i)))))
    end do
    !$omp end parallel do
    ! find approximate Jacobian with quasi-constant absolute flux Jacobian
    !$omp parallel do default(shared)&
    !$omp& private(m,n,uAvg,Havg,cAvg)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=k.and.n<=k)then
        uAvg(:)=(sqrt(rho(m))*u(:,m)+sqrt(rho(n))*u(:,n))/(sqrt(rho(m))+sqrt(rho(n)))
        HAvg=(sqrt(rho(m))*H(m)+sqrt(rho(n))*H(n))/(sqrt(rho(m))+sqrt(rho(n)))
        cAvg=sqrt((gamm-1d0)*(HAvg-0.5d0*(dot_product(uAvg,uAvg))))
        if(n<=grid%nC)then
          JacP(1:5,1:5,i)=-0.5d0*grid%aP(i)*(eulerNFJ(grid%normP(:,i),u(:,m),H(m),gamm)&
          &                                  +eulerANFJ(grid%normP(:,i),uAvg(:),HAvg,cAvg,gamm))
          JacP(1:5,6:10,i)=-0.5d0*grid%aP(i)*(eulerNFJ(grid%normP(:,i),u(:,n),H(n),gamm)&
          &                                   -eulerANFJ(grid%normP(:,i),uAvg(:),HAvg,cAvg,gamm))
          JacP(6:10,1:10,i)=-JacP(1:5,1:10,i)
        else
          JacP(1:5,1:5,i)=-0.5d0*grid%aP(i)*(eulerNFJ(grid%normP(:,i),u(:,m),H(m),gamm)&
          &                                 +eulerANFJ(grid%normP(:,i),uAvg(:),HAvg,cAvg,gamm))
        end if
      end if
    end do
    !$end omp parallel do
    ! copy Jacobian to cell
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC)then
        JacC(:,1:5,m)=JacC(:,1:5,m)+JacP(1:5,1:5,i)
      end if
      if(n<=grid%nC)then
        JacC(:,1:5,n)=JacC(:,1:5,n)+JacP(6:10,6:10,i)
      end if
      if(m<=grid%nC.and.n<=grid%nC)then
        do j=1,size(grid%neib,1)
          if(grid%neib(j,m)==n)then
            JacC(:,j*5+1:j*5+5,m)=JacC(:,j*5+1:j*5+5,m)+JacP(1:5,6:10,i)
            exit
          end if
        end do
        do j=1,size(grid%neib,1)
          if(grid%neib(j,n)==m)then
            JacC(:,j*5+1:j*5+5,n)=JacC(:,j*5+1:j*5+5,n)+JacP(6:10,1:5,i)
            exit
          end if
        end do
      end if
    end do
  end subroutine
  
  !> find advection due to flux f depending on vector s on otGrid
  !> \f[ \int_A \mathbf{f}(\mathbf{s}) \cdot \hat{n} dA \f]
  subroutine findAdvOtVect(grid,s,f,adv)
    use modOtGrid
    use modGradient
    class(otGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< state variables
    double precision,intent(in)::f(:,:,:) !< fluxes
    double precision,allocatable,intent(inout)::adv(:,:) !< advection output
    
    if(.not.allocated(adv))then
      allocate(adv(size(s,2),grid%nC))
    end if
    adv(:,:)=0d0
    !TODO
  end subroutine
  
  !> find advection due to flux f depending on scalar s on otGrid
  !> \f[ \int_A \mathbf{f}(s) \cdot \hat{n} dA \f]
  subroutine findAdvOtScal(grid,s,f,adv)
    use modOtGrid
    class(otGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< state variable
    double precision,intent(in)::f(:,:) !< flux
    double precision,allocatable,intent(inout)::adv(:) !< advection output
    double precision,allocatable::sv(:,:),fv(:,:,:),advv(:,:)
    
    allocate(sv(1,size(s)))
    allocate(fv(DIMS,1,size(s)))
    if(allocated(adv))then
      allocate(advv(1,size(adv)))
    end if
    sv(1,:)=s(:)
    fv(:,1,:)=f(:,:)
    call findAdvOtVect(grid,sv,fv,advv)
    if(.not.allocated(adv))then
      allocate(adv(size(advv,2)),source=advv(1,:))!FIXME:remove work-around
    else
      adv(:)=advv(1,:)
    end if
    deallocate(sv)
    deallocate(fv)
    deallocate(advv)
  end subroutine
  
  !> van Albada flux limiter
  elemental function vanAlbada(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision::vanAlbada !< limit function
    double precision,parameter::R_LMT=1d20
    
    vanAlbada=merge((r+r**2)/(1d0+r**2),1d0,abs(r)<R_LMT)
  end function
  
  !> Euler normal flux Jacobian
  function eulerNFJ(norm,u,H,gamm)
    double precision,intent(in)::norm(DIMS) !< normal vector
    double precision,intent(in)::u(DIMS) !< velocity
    double precision,intent(in)::H !< total mass-specific enthalpy
    double precision,intent(in)::gamm !< ratio of heat capacity
    double precision::eulerNFJ(5,5) !< result
    double precision::v,vn,tmp(DIMS,DIMS)
    
    v=norm2(u)
    vn=dot_product(u,norm)
    tmp(:,:)=0d0
    call dger(DIMS,DIMS,1d0,u,1,norm,1,tmp,DIMS)
    call dger(DIMS,DIMS,1d0-gamm,norm,1,u,1,tmp,DIMS)
    forall(i=1:DIMS)
      tmp(i,i)=tmp(i,i)+vn
    end forall
    eulerNFJ(:,:)=0d0
    eulerNFJ(2:4,1)=0.5d0*(gamm-1d0)*v**2*norm(:)-u(:)*vn
    eulerNFJ(5,1)=(0.5d0*(gamm-1d0)*v**2-H)*vn
    eulerNFJ(1,2:4)=norm(:)
    eulerNFJ(2:4,2:4)=tmp(:,:)
    eulerNFJ(5,2:4)=H*norm(:)-(gamm-1d0)*u(:)*vn
    eulerNFJ(2:4,5)=(gamm-1d0)*norm(:)
    eulerNFJ(5,5)=gamm*vn
  end function
  
  !> Euler absolute normal flux Jacobian
  function eulerANFJ(norm,u,H,c,gamm)
    double precision,intent(in)::norm(DIMS) !< normal vector
    double precision,intent(in)::u(DIMS) !< velocity
    double precision,intent(in)::H !< total mass-specific enthalpy
    double precision,intent(in)::c !< speed of sound
    double precision,intent(in)::gamm !< ratio of heat capacity
    double precision::eulerANFJ(5,5) !< result
    double precision::v,vn,M,Mn,tmp1(DIMS,DIMS),tmp2(5,5)
    
    v=norm2(u)
    vn=dot_product(u,norm)
    M=v/c
    Mn=vn/c
    eulerANFJ(:,:)=0d0
    tmp1(:,:)=0d0
    call dger(DIMS,DIMS,(gamm-1d0)/c**2,u,1,u,1,tmp1,DIMS)
    call dger(DIMS,DIMS,-1d0,norm,1,norm,1,tmp1,DIMS)
    forall(i=1:DIMS)
      tmp1(i,i)=tmp1(i,i)+1d0
    end forall
    tmp2(1,1)=1d0-0.5d0*(gamm-1d0)*M**2
    tmp2(2:4,1)=-(1d0+0.5d0*(gamm-1d0)*M**2)*u(:)+vn*norm(:)
    tmp2(5,1)=vn**2-0.5d0*v**2*(1d0+0.5d0*(gamm-1d0)*M**2)
    tmp2(1,2:4)=(gamm-1d0)/c**2*u(:)
    tmp2(2:4,2:4)=tmp1(:,:)
    tmp2(5,2:4)=(1d0+0.5d0*(gamm-1d0)*M**2)*u(:)-vn*norm(:)
    tmp2(1,5)=-(gamm-1d0)/c**2
    tmp2(2:4,5)=-(gamm-1d0)/c**2*u(:)
    tmp2(5,5)=-0.5d0*(gamm-1d0)*M**2
    eulerANFJ(:,:)=abs(vn)*tmp2(:,:)
    tmp1(:,:)=0d0
    call dger(DIMS,DIMS,1d0,u-c*norm,1,-0.5d0*(gamm-1d0)/c**2*u-0.5d0*norm/c,1,tmp1,DIMS)
    tmp2(1,1)=0.25d0*(gamm-1d0)*M**2+0.5d0*Mn
    tmp2(2:4,1)=(u(:)-c*norm(:))*tmp2(1,1)
    tmp2(5,1)=(H-vn*c)*tmp2(1,1)
    tmp2(1,2:4)=-0.5d0*(gamm-1d0)/c**2*u(:)-0.5d0*norm(:)/c
    tmp2(2:4,2:4)=tmp1(:,:)
    tmp2(5,2:4)=(H-vn*c)*tmp2(1,2:4)
    tmp2(1,5)=0.5d0*(gamm-1d0)/c**2
    tmp2(2:4,5)=(u(:)-c*norm(:))*tmp2(1,5)
    tmp2(5,5)=(H-vn*c)*tmp2(1,5)
    eulerANFJ(:,:)=eulerANFJ(:,:)+abs(vn-c)*tmp2(:,:)
    tmp1(:,:)=0d0
    call dger(DIMS,DIMS,1d0,u+c*norm,1,-0.5d0*(gamm-1d0)/c**2*u+0.5d0*norm/c,1,tmp1,DIMS)
    tmp2(1,1)=0.25d0*(gamm-1d0)*M**2-0.5d0*Mn
    tmp2(2:4,1)=(u(:)+c*norm(:))*tmp2(1,1)
    tmp2(5,1)=(H+vn*c)*tmp2(1,1)
    tmp2(1,2:4)=-0.5d0*(gamm-1d0)/c**2*u(:)+0.5d0*norm(:)/c
    tmp2(2:4,2:4)=tmp1(:,:)
    tmp2(5,2:4)=(H+vn*c)*tmp2(1,2:4)
    tmp2(1,5)=0.5d0*(gamm-1d0)/c**2
    tmp2(2:4,5)=(u(:)+c*norm(:))*tmp2(1,5)
    tmp2(5,5)=(H+vn*c)*tmp2(1,5)
    eulerANFJ(:,:)=eulerANFJ(:,:)+abs(vn+c)*tmp2(:,:)
  end function
end module
