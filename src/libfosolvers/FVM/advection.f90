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
    module procedure::findAdvOtVect
    module procedure::findAdvOtScal
  end interface
  public::findAdv
  
contains
  
  !> find mass flow rate on polyFvGrid with Rhie-Chow u interpolation and TVD rho reconstruction
  !> \f[ \int_A \rho\mathbf{u} \cdot \hat{n} dA \f]
  subroutine findMassFlowMagicPoly(grid,rho,u,p,presF,dt,flow,gradRho,sensP)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::rho(:) !< density
    double precision,intent(in)::u(:,:) !< velocity
    double precision,intent(in)::p(:) !< pressure
    double precision,intent(in)::presF(:,:) !< pressure force
    double precision,intent(in)::dt !< time step size
    double precision,allocatable,intent(inout)::flow(:) !< mass flow rate output
    double precision,intent(in),optional::gradRho(:,:) !< gradient of rho
    double precision,allocatable,intent(inout),optional::sensP(:,:) !< flow's sensitivity on p
    integer::up,dn
    double precision::vMN(DIMS),vMP(DIMS),vPN(DIMS),flux(DIMS),eps,rhoUp,rhoDn,dRho,r,rhof
    logical::sameDirection
    
    call grid%up()
    if(.not.allocated(flow))then
      allocate(flow(grid%nP))
    end if
    flow(:)=0d0
    if(present(sensP))then
      if(.not.allocated(sensP))then
        allocate(sensP(2,grid%nP))
      end if
      sensP(:,:)=0d0
    end if
    !$omp parallel do default(shared)&
    !$omp& private(m,n,up,dn,rhoUp,rhoDn,dRho,r,rhof,vMN,vMP,vPN,eps,flux)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      sameDirection=dot_product(u(:,m)+u(:,n),grid%normP(:,i))>=0d0 ! upwind by velocity
      up=merge(m,n,sameDirection)
      dn=merge(m,n,.not.sameDirection)
      rhoUp=rho(up)
      rhoDn=rho(dn)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        if(present(gradRho))then
          dRho=dot_product(grid%p(:,dn)-grid%p(:,up),gradRho(:,up))
          r=merge(2d0*dRho/(rhoDn-rhoUp)-1d0,0d0,abs(rhoDn-rhoUp)>tiny(1d0))
          rhof=rhoUp+0.5d0*vanAlbada(r)*(rhoDn-rhoUp)
        else
          rhof=rhoUp
        end if
        vMN(:)=grid%p(:,n)-grid%p(:,m)
        vMP(:)=grid%pP(:,i)-grid%p(:,m)
        vPN(:)=grid%p(:,n)-grid%pP(:,i)
        eps=norm2(vPN)/(norm2(vMP)+norm2(vPN))
        flux(:)=rhof*(eps*u(:,m)+(1d0-eps)*u(:,n))&
        &       -dt/(1d0/rho(m)+1d0/rho(n))*&
        &        (presF(:,m)/grid%v(m)/rho(m)+presF(:,n)/grid%v(n)/rho(n))&
        &       +dt*grid%normP(:,i)*(p(m)-p(n))/dot_product(vMN,grid%normP(:,i))
        if(present(sensP))then
          sensP(1,i)=grid%aP(i)*dt/dot_product(vMN,grid%normP(:,i))
          sensP(2,i)=-grid%aP(i)*dt/dot_product(vMN,grid%normP(:,i))
        end if
      else ! boundary pairs
        flux(:)=rhoUp*0.5d0*(u(:,m)+u(:,n))
        if(abs(dot_product(grid%normP(:,i),flux(:)))>tiny(1d0))then
          vMP(:)=grid%pP(:,i)-grid%p(:,m)
          flux(:)=flux(:)-0.5d0*dt*presF(:,m)/grid%v(m)&
          &       +0.5d0*dt*grid%normP(:,i)*(p(m)-p(n))/(2d0*dot_product(vMP,grid%normP(:,i)))
          if(present(sensP))then
            sensP(1,i)=grid%aP(i)*0.5d0*dt/(2d0*dot_product(vMP,grid%normP(:,i)))
            sensP(2,i)=-grid%aP(i)*0.5d0*dt/(2d0*dot_product(vMP,grid%normP(:,i)))
          end if
        end if ! Rhie-Chow does not initiate penetration through boundary pair
      end if
      flow(i)=grid%aP(i)*dot_product(grid%normP(:,i),flux(:))
    end do
    !$omp end parallel do
  end subroutine
  
  !> find flow rate of vector v on polyFvGrid using mass flow rate
  !> \f[ \int_A \rho\mathbf{v}(\mathbf{u} \cdot \hat{n}) dA \f]
  subroutine findVarFlowPolyVect(grid,v,mFlow,flow,gradV)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:,:) !< transported variables
    double precision,intent(in)::mFlow(:) !< mass flow through pairs
    double precision,allocatable,intent(inout)::flow(:,:) !< v flow rate output
    double precision,intent(in),optional::gradV(:,:,:) !< gradient of v
    double precision::vf(size(v,1)),dV(size(v,1)),r(size(v,1))
    integer::up,dn
    logical::sameDirection
    
    call grid%up()
    if(.not.allocated(flow))then
      allocate(flow(size(v,1),grid%nP))
    end if
    flow(:,:)=0d0
    !$omp parallel do default(shared)&
    !$omp& private(m,n,up,dn,dV,r,vf)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(abs(mFlow(i))<=tiny(1d0))then ! no mass flux
        cycle
      else ! upwinding by mass flux
        sameDirection=mFlow(i)>=0d0
        up=merge(m,n,sameDirection)
        dn=merge(m,n,.not.sameDirection)
        vf(:)=v(:,up)
        if(m<=grid%nC.and.n<=grid%nC.and.present(gradV))then ! internal pairs and gradV available
          dV(:)=matmul(grid%p(:,dn)-grid%p(:,up),gradV(:,:,up))
          r(:)=merge(2d0*dV/(v(:,dn)-v(:,up))-1d0,0d0,abs(v(:,dn)-v(:,up))>tiny(1d0))
          vf(:)=vf(:)+0.5d0*vanAlbada(r)*(v(:,dn)-v(:,up))
        end if
        flow(:,i)=vf(:)*mFlow(i)
      end if
    end do
    !$omp end parallel do
  end subroutine
  
  !> find flow rate of scalar v on polyFvGrid using mass flow rate
  !> \f[ \int_A \rho v(\mathbf{u} \cdot \hat{n}) dA \f]
  subroutine findVarFlowPolyScal(grid,v,mFlow,flow,gradV)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:) !< transported variable
    double precision,intent(in)::mFlow(:) !< mass flow through pairs
    double precision,allocatable,intent(inout)::flow(:) !< v flow rate output
    double precision,intent(in),optional::gradV(:,:) !< gradient of v
    double precision::vf,dV,r
    integer::up,dn
    logical::sameDirection
    
    call grid%up()
    if(.not.allocated(flow))then
      allocate(flow(grid%nP))
    end if
    flow(:)=0d0
    !$omp parallel do default(shared)&
    !$omp& private(m,n,up,dn,dV,r,vf)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(abs(mFlow(i))<=tiny(1d0))then ! no mass flux
        cycle
      else ! upwinding by mass flux
        sameDirection=mFlow(i)>=0d0
        up=merge(m,n,sameDirection)
        dn=merge(m,n,.not.sameDirection)
        vf=v(up)
        if(m<=grid%nC.and.n<=grid%nC.and.present(gradV))then ! internal pairs and gradV available
          dV=dot_product(grid%p(:,dn)-grid%p(:,up),gradV(:,up))
          r=merge(2d0*dV/(v(dn)-v(up))-1d0,0d0,abs(v(dn)-v(up))>tiny(1d0))
          vf=vf+0.5d0*vanAlbada(r)*(v(dn)-v(up))
        end if
        flow(i)=vf*mFlow(i)
      end if
    end do
    !$omp end parallel do
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
    
    if(r<0d0)then
      vanAlbada=0d0
    else if(r<R_LMT)then
      vanAlbada=(r+r**2)/(1d0+r**2)
    else
      vanAlbada=1d0
    end if
  end function
  
end module
