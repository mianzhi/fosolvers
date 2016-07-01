!----------------------------------------------------------------------------- best with 100 columns

!> the nonlinear residual function to be solved
subroutine fkfun(s,r,ier)
  use modPbc
  use modPolyFvGrid
  use modGradient
  use modAdvection
  use modDiffusion
  use modNewtonian
  use modPressure
  use modRhieChow
  double precision::s(*)
  double precision::r(*)
  integer::ier
  
  call extractVar(s(1:nEq),p,u,temp)
  call setBC()
  call recoverState(p,u,temp,Y,rho,rhou,rhoE)
  call findGrad(grid,p,gradP)
  forall(i=1:grid%nE)
    carrier(:,i)=[rho(i),rhou(:,i),rhoE(i)+p(i)]
    forall(j=1:5)
      flux(:,j,i)=carrier(j,i)*u(:,i)
    end forall
  end forall
  call findAdv(grid,carrier,flux,flow)
  call addRhieChow(grid,carrier,p,gradP,rho,dt,flow)
  call findViscForce(grid,u,visc,viscF)
  call findPresForce(grid,p,presF)
  call findDiff(grid,temp,cond,condQ)
  do i=1,grid%nC
    j=(i-1)*5
    r(j+1)=rho(i)-rho0(i)-dt/grid%v(i)*flow(1,i)
    r(j+2:j+4)=rhou(:,i)-rhou0(:,i)-dt/grid%v(i)*(flow(2:4,i)+presF(:,i)+viscF(:,i))
    r(j+5)=rhoE(i)-rhoE0(i)-dt/grid%v(i)*(flow(5,i)+condQ(i))
  end do
  r(1:nEq)=r(1:nEq)/rscale(1:nEq)
  ier=0
end subroutine
