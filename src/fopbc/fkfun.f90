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
  double precision,save,allocatable::gradP(:,:),viscF(:,:),presF(:,:),condQ(:),&
  &                                  tmp1(:,:),tmp2(:,:,:),tmp3(:,:)
  integer::ier
  
  if(.not.allocated(gradP))then
    allocate(gradP(DIMS,grid%nC))
  else if(size(gradP,2)/=grid%nC)then
    deallocate(gradP)
    allocate(gradP(DIMS,grid%nC))
  end if
  if(.not.allocated(viscF))then
    allocate(viscF(DIMS,grid%nC))
  else if(size(viscF,2)/=grid%nC)then
    deallocate(viscF)
    allocate(viscF(DIMS,grid%nC))
  end if
  if(.not.allocated(presF))then
    allocate(presF(DIMS,grid%nC))
  else if(size(presF,2)/=grid%nC)then
    deallocate(presF)
    allocate(presF(DIMS,grid%nC))
  end if
  if(.not.allocated(condQ))then
    allocate(condQ(grid%nC))
  else if(size(condQ)/=grid%nC)then
    deallocate(condQ)
    allocate(condQ(grid%nC))
  end if
  ! {rho,rhou,rhoH+rhoKE}
  if(.not.allocated(tmp1))then
    allocate(tmp1(5,grid%nE))
  else if(size(tmp1,2)/=grid%nE)then
    deallocate(tmp1)
    allocate(tmp1(5,grid%nE))
  end if
  ! flux of {rho,rhou,rhoH+rhoKE}
  if(.not.allocated(tmp2))then
    allocate(tmp2(DIMS,5,grid%nE))
  else if(size(tmp2,3)/=grid%nE)then
    deallocate(tmp2)
    allocate(tmp2(DIMS,5,grid%nE))
  end if
  ! flow rate of {rho,rhou,rhoH+rhoKE}
  if(.not.allocated(tmp3))then
    allocate(tmp3(5,grid%nC))
  else if(size(tmp3,2)/=grid%nC)then
    deallocate(tmp3)
    allocate(tmp3(5,grid%nC))
  end if
  
  call extractVar(s(1:nEq),p,u,temp)
  ! FIXME apply BC on p,u,temp here
  call recoverState(p,u,temp,Y,rho,rhou,rhoE)
  call findGrad(grid,p(1:grid%nC),gradP)
  forall(i=1:grid%nE)
    tmp1(:,i)=[rho(i),rhou(:,i),rhoE(i)+p(i)+0.5d0*rho(i)*dot_product(u(:,i),u(:,i))]
    forall(j=1:5)
      tmp2(:,j,i)=tmp1(j,i)*u(:,i)
    end forall
  end forall
  call findAdv(grid,tmp1,tmp2,tmp3)
  call addRhieChow(grid,tmp1,p,gradP,rho,dt,tmp3)
  call findViscForce(grid,u,visc,viscF)
  call findPresForce(grid,p,presF)
  call findDiff(grid,temp,cond,condQ)
  do i=1,grid%nC
    j=(i-1)*5
    r(j+1)=(rho(i)-rho0(i))-dt/grid%v(i)*tmp3(1,i)
    r(j+2:j+4)=rhou(:,i)-rhou0(:,i)-dt/grid%v(i)*(tmp3(2:4,i)+presF(:,i)+viscF(:,i))
    r(j+5)=rhoE(i)-rhoE0(i)-dt/grid%v(i)*(tmp3(5,i)+condQ(i))
  end do
  r(1:nEq)=r(1:nEq)*rscale(1:nEq)
  ier=0
end subroutine
