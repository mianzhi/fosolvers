!----------------------------------------------------------------------------- best with 100 columns

!> the nonlinear residual function to be solved
subroutine fkfun(s,r,ier)
  use modPbc
  use modPolyFvGrid
  use modGradient
  double precision::s(*)
  double precision::r(*)
  integer::ier
  
  call extractVar(s(1:nEq),p,u,temp)
  call recoverState(p,u,temp,Y,rho,rhou,rhoE)
  call findGrad(grid,p(1:grid%nC),gradP)
  ier=0
end subroutine
