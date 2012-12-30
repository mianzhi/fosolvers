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
  double precision testU(DIMS,size(testU1d)/DIMS),testMom(DIMS,size(testU1d)/DIMS),&
  &                resMomWrapped(DIMS,size(testU1d)/DIMS)
  
  call grid%updateDualBlock()
  testU=reshape(testU1d,[DIMS,grid%nNode])
  tao=findTao(testU)
  forall(i=1:grid%nNode)
    testMom(:,i)=testU(:,i)*rhoNode(i)*grid%NodeVol(i)
  end forall
  resMomWrapped(:,:)=-testMom(:,:)+Mom(:,:)
  do i=1,grid%nBlock
    do j=1,grid%Block(i)%nNode
      n=grid%Block(i)%iNode(j)
      do k=1,size(grid%NodeNeibBlock(n)%dat)
        if(grid%NodeNeibBlock(n)%dat(k)==i)then
          resMomWrapped(:,n)=resMomWrapped(:,n)&
          &                  -dt*grid%NBAreaVect(n)%dat(:,k)*p(i)&
          &                  +dt*matmul(grid%NBAreaVect(n)%dat(:,k),tao(:,:,i))
        end if
      end do
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

!> the energy residual function (not includes pressure work)
function resEnergy(testT)
  use moduleGrid
  use moduleCondition
  use miscNS
  double precision testT(:) !< the test temperature
  double precision resEnergy(size(testT)) !< the energy residual function
  
  resEnergy(:)=0d0
end function
