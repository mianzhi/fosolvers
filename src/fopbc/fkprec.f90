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
  
  tmp(1:nEq)=v(1:nEq)*rscale(1:nEq)
  do i=1,grid%nC
    j=(i-1)*5
    v(j+2:j+4)=tmp(j+2:j+4)/rho(i)
    ! FIXME calculation based on cantera
    tmp(j+5)=tmp(j+5)-rho0(i)*dot_product(u0(:,i),v(j+2:j+4))
    v(j+5)=tmp(j+5)/(rho0(i)*287.058d0/(1.4d0-1d0))
    tmp(j+1)=tmp(j+1)+p(i)/287.058d0/temp(i)**2*v(j+5)
    v(j+1)=tmp(j+1)*287.058d0*temp(i)
  end do
  !v(1:nEq)=v(1:nEq)*ss(1:nEq)
  ier=0
end subroutine
