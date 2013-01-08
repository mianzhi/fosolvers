!----------------------------------------------------------------------------- best with 100 columns

!> the momentum residual function (not coupled with pressure)
function resMom(testU1d)
  use moduleGrid
  use moduleInterpolation
  use moduleCondition
  use miscNS
  double precision testU1d(:) !< the test velocity (unwrapped as 1d array)
  double precision resMom(size(testU1d)) !< the momentum residual function
  double precision testU(DIMS,size(testU1d)/DIMS),testMom(DIMS,size(testU1d)/DIMS),&
  &                resMomWrapped(DIMS,size(testU1d)/DIMS)
  
  call grid%updateDualBlock()
  testU=reshape(testU1d,[DIMS,grid%nNode])
  tao=findTao(testU)
  forall(i=1:grid%nNode)
    testMom(:,i)=testU(:,i)*rhoNode(i)*grid%NodeVol(i)
  end forall
  ! pressure force and viscous force
  resMomWrapped(:,:)=-testMom(:,:)+Mom(:,:)
  do i=1,grid%nBlock
    do j=1,grid%Block(i)%nNode
      n=grid%Block(i)%iNode(j)
      do k=1,size(grid%NodeNeibBlock(n)%dat)
        if(grid%NodeNeibBlock(n)%dat(k)==i)then
          resMomWrapped(:,n)=resMomWrapped(:,n)&
          &                  -dt*grid%NBAreaVect(n)%dat(:,k)*p(i)&
          &                  +dt*matmul(grid%NBAreaVect(n)%dat(:,k),tao(:,:,i))
          exit
        end if
      end do
    end do
  end do
  ! boundary conditions
  do i=1,grid%nFacet
    l=findCondition(condition,grid%Facet(i)%Ent,'Wall')
    if(l>0)then
      do j=1,grid%Facet(i)%nNode
        k=grid%Facet(i)%iNode(j)
        resMomWrapped(:,k)=-testMom(:,k)+[0d0,0d0,0d0]*rhoNode(k)*grid%NodeVol(k)
      end do
    end if
  end do
  resMom=reshape(resMomWrapped,[DIMS*grid%nNode])
end function

!> the energy residual function (not includes pressure work)
function resEnergy(testT)
  use moduleGrid
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleCondition
  use miscNS
  double precision testT(:) !< the test temperature
  double precision resEnergy(size(testT)) !< the energy residual function
  double precision testEnergy(size(testT)),viscWork
  
  call grid%updateBlockVol()
  call grid%updateIntfArea()
  call grid%updateIntfNorm()
  call grid%updateFacetNeib()
  call grid%updateFacetPos()
  call grid%updateFacetArea()
  call grid%updateFacetNorm()
  forall(i=1:grid%nBlock)
    testEnergy(i)=(200d0/(gamm-1)*testT(i)+dot_product(uBlock(:,i),uBlock(:,i))/2d0)& !TODO:IE=IE(p,T)
    &             *rho(i)*grid%BlockVol(i)
  end forall
  ! work done by viscous force
  resEnergy(:)=-testEnergy(:)+Energy(:)
  do i=1,grid%nIntf
    m=grid%IntfNeibBlock(1,i)
    n=grid%IntfNeibBlock(2,i)
    viscWork=dt*grid%IntfArea(i)*dot_product(matmul(taoIntf(:,:,i),grid%IntfNorm(:,i)),uIntf(:,i))
    resEnergy(m)=resEnergy(m)+viscWork
    resEnergy(n)=resEnergy(n)-viscWork
  end do
  ! heat conduction
  gradT=findGrad(testT,grid,BIND_BLOCK)
  resEnergy=resEnergy+dt*findDiffus(thermK,BIND_BLOCK,testT,grid,gradT)
  ! boundary conditions
  do i=1,grid%nFacet
    l=findCondition(condition,grid%Facet(i)%Ent,'Wall')
    if(l>0)then
      m=maxval(grid%FacetNeibBlock(:,i))
      if(condition(l)%dat%test('Wall_Heat_Flux'))then
        resEnergy(m)=resEnergy(m)+dt*grid%FacetArea(i)*condition(l)%dat%get('Wall_Heat_Flux',&
        &                         [grid%FacetPos(1,i),grid%FacetPos(2,i),grid%FacetPos(3,i),t])
      end if
      if(condition(l)%dat%test('Wall_Temperature'))then
        !TODO: wall temperature condition
      end if
    end if
  end do
end function

!> the pressure residual function
function resPressure(testP)
  use moduleGrid
  use moduleGridOperation
  use miscNS
  double precision testP(:) !< the test pressure
  double precision resPressure(size(testP)) !< the pressure residual function
  double precision,allocatable::testU(:,:),tempVol(:)
  
  allocate(testU(DIMS,grid%nNode))
  allocate(tempVol(grid%nBlock))
  tempVol(:)=grid%BlockVol(:)
  ! resolve velocity using momentum equation with pressure being tested
  testU(:,:)=u(:,:)
  do i=1,grid%nNode
    do j=1,size(grid%NodeNeibBlock(i)%dat)
      testU(:,i)=testU(:,i)-dt*grid%NBAreaVect(i)%dat(:,j)*p(grid%NodeNeibBlock(i)%dat(j))&
      &                     /grid%NodeVol(i)/rhoNode(i)
    end do
  end do
  do i=1,grid%nFacet
    if(findCondition(condition,grid%Facet(i)%Ent,'Wall')>0)then
      do j=1,grid%Facet(i)%nNode
        k=grid%Facet(i)%iNode(j)
        testU(:,k)=0d0
      end do
    end if
  end do
  ! find pressure residual with moved grid
  call mvGrid(grid,dt*testU)
  call grid%updateBlockVol()
  resPressure(:)=-testP(:)+(gamm-1d0)*(IEnergy(:)+oldP(:)*(tempVol(:)-grid%BlockVol(:)))/tempVol(:) !TODO:p=p(IE,v)
  call mvGrid(grid,-dt*testU)
  grid%BlockVol(:)=tempVol(:)
  deallocate(testU)
  deallocate(tempVol)
end function
