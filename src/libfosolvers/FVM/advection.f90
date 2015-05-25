!----------------------------------------------------------------------------- best with 100 columns

!> advection for FVM
module modAdvection
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding advection
  interface findAdv
    module procedure::findAdvPolyVect
    module procedure::findAdvPolyScal
    module procedure::findAdvPolyEuler
    module procedure::findAdvOtVect
    module procedure::findAdvOtScal
  end interface
  public::findAdv
  
contains
  
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
    if(.not.(allocated(adv)))then
      allocate(adv(size(s,2),grid%nC))
    end if
    adv(:,:)=0d0
    allocate(gradF(DIMS,DIMS*size(s,1),grid%nC))
    call findGrad(grid,reshape(f(:,:,1:grid%nC),&
    &                          [size(f(:,:,1:grid%nC),1)*size(f(:,:,1:grid%nC),2),&
    &                           size(f(:,:,1:grid%nC),3)]),gradF)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(s,2).and.m<=size(f,3).and.n<=size(s,2).and.n<=size(f,3))then
        do j=1,size(s,1)
          if(abs(s(j,m)-s(j,n))<=tiny(1d0))then
            up=m
            dn=n
          else if(dot_product((f(:,j,m)-f(:,j,n))/(s(j,m)-s(j,n)),grid%normP(:,i))>=0d0)then
            up=m
            dn=n
          else
            up=n
            dn=m
          end if
          fUp=dot_product(f(:,j,up),grid%normP(:,i))
          fDn=dot_product(f(:,j,dn),grid%normP(:,i))
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
        end do
      end if
    end do
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
  subroutine findAdvPolyEuler(grid,rho,rhou,rhoE,p,dRho,dRhou,dRhoE)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::rho(:) !< cell-averaged rho
    double precision,intent(in)::rhou(:,:) !< cell-averaged rhou
    double precision,intent(in)::rhoE(:) !< cell-averaged rhoE
    double precision,intent(in)::p(:) !< cell-averaged pressure
    double precision,allocatable,intent(inout)::dRho(:) !< net flux of rho
    double precision,allocatable,intent(inout)::dRhou(:,:) !< net flux of rhou
    double precision,allocatable,intent(inout)::dRhoE(:) !< net flux of rhoE
    double precision,allocatable::u(:,:),H(:)
    double precision,allocatable::f1c(:,:),f2c(:,:,:),f3c(:,:)
    double precision::rhoAvg,uAvg(DIMS),HAvg,cAvg,uNormAvg,rhoJump,pJump,uJump(DIMS),uNormJump
    double precision::dFEntropy(DIMS+2),dFAcoustic1(DIMS+2),dFAcoustic2(DIMS+2),flow(DIMS+2)
    double precision,parameter::GAMM=1.4d0 ! TODO other gamma and real gas
    
    call grid%up()
    if(.not.(allocated(dRho)))then
      allocate(dRho(grid%nC))
    end if
    dRho(:)=0d0
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    if(.not.(allocated(dRhoE)))then
      allocate(dRhoE(grid%nC))
    end if
    dRhoE(:)=0d0
    ! find auxiliary state and flux vectors in cell
    allocate(u(DIMS,grid%nC))
    allocate(H(grid%nC))
    allocate(f1c(DIMS,grid%nC))
    allocate(f2c(DIMS,DIMS,grid%nC))
    allocate(f3c(DIMS,grid%nC))
    forall(i=1:grid%nC)
      u(:,i)=rhou(:,i)/rho(i)
      H(i)=(rhoE(i)+p(i))/rho(i)
      f1c(:,i)=rhou(:,i)
      f2c(:,1,i)=rhou(1,i)*u(:,i)+[p(i),0d0,0d0]
      f2c(:,2,i)=rhou(2,i)*u(:,i)+[0d0,p(i),0d0]
      f2c(:,3,i)=rhou(3,i)*u(:,i)+[0d0,0d0,p(i)]
      f3c(:,i)=(rhoE(i)+p(i))*u(:,i)
    end forall
    ! Roe flux difference splitting
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        rhoAvg=sqrt(rho(m)*rho(n))
        uAvg(:)=(sqrt(rho(m))*u(:,m)+sqrt(rho(n))*u(:,n))/(sqrt(rho(m))+sqrt(rho(n)))
        uNormAvg=dot_product(uAvg,grid%normP(:,i))
        HAvg=(sqrt(rho(m))*H(m)+sqrt(rho(n))*H(n))/(sqrt(rho(m))+sqrt(rho(n)))
        cAvg=sqrt((GAMM-1d0)*(HAvg-0.5d0*(dot_product(uAvg,uAvg)))) ! TODO arbitrary gamma; real gas
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
        dRho(m)=dRho(m)-flow(1)
        dRho(n)=dRho(n)+flow(1)
        dRhou(:,m)=dRhou(:,m)-flow(2:4)
        dRhou(:,n)=dRhou(:,n)+flow(2:4)
        dRhoE(m)=dRhoE(m)-flow(5)
        dRhoE(n)=dRhoE(n)+flow(5)
      end if
    end do
    deallocate(u,H)
    deallocate(f1c,f2c,f3c)
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
    
    if(.not.(allocated(adv)))then
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
  
end module
