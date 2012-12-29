!----------------------------------------------------------------------------- best with 100 columns

!> the momentum residual function (not coupled with pressure)
function resMom(testU1d)
  use moduleGrid
  use moduleFVMGrad
  use moduleInterpolation
  use moduleCondition
  use miscNS
  double precision testU1d(:) !< the test velocity (unwrapped as 1d array)
  double precision resMom(size(testU1d)) !< the momentum residual function
  double precision testU(DIMS,size(testU1d)/DIMS),gradU(DIMS,DIMS,size(testU1d)/DIMS),&
  &                tao(DIMS,DIMS,size(testU1d)/DIMS),testMom(DIMS,size(testU1d)/DIMS),&
  &                resMomWrapped(DIMS,size(testU1d)/DIMS)
  !FIXME:replace the local viscosities
  double precision visc,viscRate
  visc=1d0
  viscRate=-2d0/3d0
  
  call grid%updateDualBlock()
  testU=reshape(testU1d,[DIMS,grid%nNode])
  gradU=findGrad(testU,grid,BIND_NODE)
  forall(l=1:grid%nNode)
    tao(:,:,l)=visc*(gradU(:,:,l)+transpose(gradU(:,:,l)))
    forall(i=1:DIMS)
      tao(i,i,l)=tao(i,i,l)+viscRate*visc*sum([(gradU(j,j,l),j=1,DIMS)])
    end forall
  end forall
  forall(i=1:grid%nNode)
    testMom(:,i)=testU(:,i)*rhoNode(i)*grid%NodeVol(i)
  end forall
  resMomWrapped(:,:)=-testMom(:,:)+Mom(:,:)
  do i=1,grid%nNode
    do j=1,size(grid%NodeNeibBlock(i)%dat)
      resMomWrapped(:,i)=resMomWrapped(:,i)&
      &                  -dt*grid%NBAreaVect(i)%dat(:,j)*p(grid%NodeNeibBlock(i)%dat(j))&
      &                  +dt*matmul(grid%NBAreaVect(i)%dat(:,j),tao(:,:,i))
    end do
  end do
  do i=1,grid%nFacet
    if(findCondition(condition,grid%Facet(i)%Ent,'Wall')>0)then
      do j=1,grid%Facet(i)%nNode
        k=grid%Facet(i)%iNode(j)
        resMomWrapped(:,k)=-testMom(:,k)+[0d0,0d0,0d0]*rhoNode(k)*grid%NodeVol(k)
      end do
    end if
  end do
  resMom=reshape(resMomWrapped,[DIMS*grid%nNode])
end function
