!----------------------------------------------------------------------------- best with 100 columns

!> Euler system for FVM
module modEuler
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding Euler advection rate into each cell
  interface findEuler
    module procedure::findEulerPoly
  end interface
  public::findEuler
  
  !> generic finding Euler advection Jacobian
  interface findEulerJac
    module procedure::findEulerJacPoly
  end interface
  public::findEulerJac
  
contains
  
  !> find advection of Euler system on polyFvGrid
  subroutine findEulerPoly(grid,rho,rhou,rhoE,p,gamm,dRho,dRhou,dRhoE)
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
    !$omp parallel default(shared)
    !$omp do
    do i=1,k
      u(:,i)=rhou(:,i)/rho(i)
      H(i)=(rhoE(i)+p(i))/rho(i)
      f1c(:,i)=rhou(:,i)
      f2c(:,1,i)=rhou(1,i)*u(:,i)+[p(i),0d0,0d0]
      f2c(:,2,i)=rhou(2,i)*u(:,i)+[0d0,p(i),0d0]
      f2c(:,3,i)=rhou(3,i)*u(:,i)+[0d0,0d0,p(i)]
      f3c(:,i)=(rhoE(i)+p(i))*u(:,i)
    end do
    !$omp end do
    ! Roe flux difference splitting
    !$omp do private(m,n,rhoAvg,uAvg,Havg,cAvg,uNormAvg,rhoJump,pJump,uJump,uNormJump,&
    !$omp&           dFEntropy,dFAcoustic1,dFAcoustic2,flow)
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
        if(n<=grid%nC)then
          !$omp atomic
          dRho(m)=dRho(m)-flow(1)
          !$omp atomic
          dRho(n)=dRho(n)+flow(1)
          do j=1,DIMS
            !$omp atomic
            dRhou(j,m)=dRhou(j,m)-flow(j+1)
            !$omp atomic
            dRhou(j,n)=dRhou(j,n)+flow(j+1)
          end do
          !$omp atomic
          dRhoE(m)=dRhoE(m)-flow(5)
          !$omp atomic
          dRhoE(n)=dRhoE(n)+flow(5)
        else
          !$omp atomic
          dRho(m)=dRho(m)-flow(1)
          do j=1,DIMS
            !$omp atomic
            dRhou(j,m)=dRhou(j,m)-flow(j+1)
          end do
          !$omp atomic
          dRhoE(m)=dRhoE(m)-flow(5)
        end if
      end if
    end do
    !$omp end do
    !$omp end parallel
  end subroutine
  
  !> find advection Jacobian of Euler system on polyFvGrid
  subroutine findEulerJacPoly(grid,rho,rhou,rhoE,p,gamm,JacP,JacC)
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
    !$omp parallel default(shared)
    !$omp do
    do i=1,k
      u(:,i)=rhou(:,i)/rho(i)
      H(i)=(rhoE(i)+p(i))/rho(i)
      c(i)=sqrt((gamm-1d0)*(H(i)-0.5d0*(dot_product(u(:,i),u(:,i)))))
    end do
    !$omp end do
    ! find approximate Jacobian with quasi-constant absolute flux Jacobian
    !$omp do private(m,n,uAvg,Havg,cAvg)
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
    !$omp end do
    !$omp end parallel
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
