!----------------------------------------------------------------------------- best with 100 columns

!> set the preconditioner
subroutine fkpset(s,ss,r,rs,tmp1,tmp2,ier)
  use modPbc
  double precision::s(*) !< state
  double precision::ss(*) !< scale of state
  double precision::r(*) !< residual
  double precision::rs(*) !< scale of residual
  double precision::tmp1(*),tmp2(*) !< work space
  integer::ier !< error indicator
  
  call extractVar(s(1:nEq),p,u,temp)
  call setBC()
  call recoverState(p,u,temp,Y,rho,rhou,rhoE)
  do i=1,grid%nC
    ! FIXME calculation based on cantera
    tmp1(i)=dot_product(u(:,i),u(:,i))
    tmp2(i)=287.058d0/(1.4d0-1d0)-tmp1(i)/temp(i)
    localPrec(:,1,i)=[287.058d0*temp(i)+287.058d0*tmp1(i)/tmp2(i),&
    &                 -u(1,i)/rho(i)-u(1,i)*tmp1(i)/rho(i)/tmp2(i)/temp(i),&
    &                 -u(2,i)/rho(i)-u(2,i)*tmp1(i)/rho(i)/tmp2(i)/temp(i),&
    &                 -u(3,i)/rho(i)-u(3,i)*tmp1(i)/rho(i)/tmp2(i)/temp(i),&
    &                 tmp1(i)/rho(i)/tmp2(i)]
    localPrec(:,2,i)=[-u(1,i)*287.058d0/tmp2(i),&
    &                 1d0/rho(i)+u(1,i)**2/rho(i)/tmp2(i)/temp(i),&
    &                 u(1,i)*u(2,i)/rho(i)/tmp2(i)/temp(i),&
    &                 u(1,i)*u(3,i)/rho(i)/tmp2(i)/temp(i),&
    &                 -u(1,i)/rho(i)/tmp2(i)]
    localPrec(:,3,i)=[-u(2,i)*287.058d0/tmp2(i),&
    &                 u(1,i)*u(2,i)/rho(i)/tmp2(i)/temp(i),&
    &                 1d0/rho(i)+u(2,i)**2/rho(i)/tmp2(i)/temp(i),&
    &                 u(2,i)*u(3,i)/rho(i)/tmp2(i)/temp(i),&
    &                 -u(2,i)/rho(i)/tmp2(i)]
    localPrec(:,4,i)=[-u(3,i)*287.058d0/tmp2(i),&
    &                 u(1,i)*u(3,i)/rho(i)/tmp2(i)/temp(i),&
    &                 u(2,i)*u(3,i)/rho(i)/tmp2(i)/temp(i),&
    &                 1d0/rho(i)+u(3,i)**2/rho(i)/tmp2(i)/temp(i),&
    &                 -u(3,i)/rho(i)/tmp2(i)]
    localPrec(:,5,i)=[p(i)/rho0(i)/tmp2(i)/temp(i),&
    &                 -u(1,i)/rho0(i)/tmp2(i)/temp(i),&
    &                 -u(2,i)/rho0(i)/tmp2(i)/temp(i),&
    &                 -u(3,i)/rho0(i)/tmp2(i)/temp(i),&
    &                 1d0/rho0(i)/tmp2(i)]
    j=(i-1)*5
    forall(m=1:5,n=1:5)
      localPrec(m,n,i)=localPrec(m,n,i)*ss(j+m)/rs(j+n)
    end forall
  end do
  ier=0
end subroutine

!> solve the preconditioning problem
subroutine fkpsol(s,ss,r,rs,v,tmp,ier)
  use modPbc
  double precision::s(*) !< state
  double precision::ss(*) !< scale of state
  double precision::r(*) !< residual
  double precision::rs(*) !< scale of residual
  double precision::v(*) !< RHS and solution of the preconditioning problem
  double precision::tmp(*) !< work space
  integer::ier !< error indicator
  
  tmp(1:nEq)=v(1:nEq)
  do i=1,grid%nC
    j=(i-1)*5
    v(j+1:j+5)=matmul(localPrec(:,:,i),tmp(j+1:j+5))
  end do
  ier=0
end subroutine
