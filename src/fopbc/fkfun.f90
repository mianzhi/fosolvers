!----------------------------------------------------------------------------- best with 100 columns

!> the nonlinear residual function to be solved
subroutine fkfun(s,r,ier)
  use modPbc
  use modPolyFvGrid
  use modGradient
  use modAdvection
  use modRhieChow
  double precision::s(*)
  double precision::r(*)
  double precision,save,allocatable::gradP(:,:),tmp1(:,:),tmp2(:,:,:),tmp3(:,:)
  integer::ier
  
  if(.not.allocated(gradP))then
    allocate(gradP(DIMS,grid%nC))
  else if(size(gradP,2)/=grid%nC)then
    deallocate(gradP)
    allocate(gradP(DIMS,grid%nC))
  end if
  ! {rho,rhou,rhoE,rhoKE}
  if(.not.allocated(tmp1))then
    allocate(tmp1(6,grid%nE))
  else if(size(tmp1,2)/=grid%nE)then
    deallocate(tmp1)
    allocate(tmp1(6,grid%nE))
  end if
  ! flux of {rho,rhou,rhoE,rhoKE}
  if(.not.allocated(tmp2))then
    allocate(tmp2(DIMS,6,grid%nE))
  else if(size(tmp2,3)/=grid%nE)then
    deallocate(tmp2)
    allocate(tmp2(DIMS,6,grid%nE))
  end if
  ! flow rate of {rho,rhou,rhoE,rhoKE}
  if(.not.allocated(tmp3))then
    allocate(tmp3(6,grid%nC))
  else if(size(tmp3,2)/=grid%nC)then
    deallocate(tmp3)
    allocate(tmp3(6,grid%nC))
  end if
  
  call extractVar(s(1:nEq),p,u,temp)
  ! FIXME apply BC on p,u,temp here
  call recoverState(p,u,temp,Y,rho,rhou,rhoE)
  call findGrad(grid,p(1:grid%nC),gradP)
  forall(i=1:grid%nE)
    tmp1(:,i)=[rho(i),rhou(:,i),rhoE(i),0.5d0*rho(i)*dot_product(u(:,i),u(:,i))]
    forall(j=1:6)
      tmp2(:,j,i)=tmp1(j,i)*u(:,i)
    end forall
  end forall
  call findAdv(grid,tmp1,tmp2,tmp3)
  call addRhieChow(grid,tmp1,p,gradP,rho,dt,tmp3)
  ier=0
end subroutine
