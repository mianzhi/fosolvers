!----------------------------------------------------------------------------- best with 100 columns

!> global variables for foeuler
module gVarEuler
  use moduleGrid
  public
  
  type(typeGrid)::grid !< the grid
  double precision,allocatable::rho(:) !< density
  double precision,allocatable::u(:,:) !< velocity
  double precision,allocatable::p(:) !< pressure
  double precision,allocatable::E(:) !< internal_energy+kinetic_energy
  double precision,allocatable::rhou(:,:) !< momentum per unit volume
  double precision,allocatable::rhoE(:) !< rho*E
  double precision gamm !< gamma=c_p/c_v
  double precision t !< current time
  double precision dt !< time step size
  double precision tFinal !< final time
  integer iStep !< current time step number
  integer nStep !< total number of time steps
end module

!> foeuler main program
program foeuler
  use moduleBasicDataStruct
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  use moduleInterpolation
  use moduleFVMGrad
  use moduleFVMConvect
  use moduleMPIComm
  use gVarEuler
  double precision,allocatable::rhoNode(:) !< density at node
  double precision,allocatable::uBlock(:,:) !< velocity at block
  double precision,allocatable::uIntf(:,:) !< velocity at interface
  double precision,allocatable::pIntf(:) !< pressure at interface
  double precision,allocatable::tempMass(:) !< temporary block mass
  double precision,allocatable::tempMom(:,:) !< temporary extensive node momentum
  double precision,allocatable::tempEnergy(:) !< temporary extensive block energy
  double precision pWork !< pressure work done on block surface
  
  call initMPI()
  if(pidMPI==0)then
    ! read grid
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    call grid%updateIntf()
    close(12)
    ! allocate storage
    allocate(rho(grid%nBlock))
    allocate(u(DIMS,grid%nNode))
    allocate(p(grid%nBlock))
    allocate(E(grid%nBlock))
    allocate(rhou(DIMS,grid%nNode))
    allocate(rhoE(grid%nBlock))
    allocate(rhoNode(grid%nNode))
    allocate(uBlock(DIMS,grid%nBlock))
    allocate(uIntf(DIMS,grid%nIntf))
    allocate(pIntf(grid%nIntf))
    allocate(tempMass(grid%nBlock))
    allocate(tempMom(DIMS,grid%nNode))
    allocate(tempEnergy(grid%nBlock))
    ! simulation control
    dt=1d-6
    dFinal=1d-4
    nStep=ceiling(dFinal/dt)
    ! initial value of variables
    call grid%updateBlockPos()
    gamm=1.4d0
    u(:,:)=0d0
    rhou(:,:)=0d0
    forall(i=1:grid%nBlock)
      rho(i)=merge(1d0,0.125d0,grid%BlockPos(1,i)<0.5d0)
      p(i)=merge(1d5,1d4,grid%BlockPos(1,i)<0.5d0)
      E(i)=p(i)/rho(i)/(gamm-1d0)
      rhoE(i)=rho(i)*E(i)
    end forall
    t=0d0
    ! write initial states
    open(13,file='rstEuler.msh',status='replace')
    call writeGMSH(13,grid)
    call writeGMSH(13,rho,grid,BIND_BLOCK,'rho',0,t)
    call writeGMSH(13,u,grid,BIND_NODE,'u',0,t)
    call writeGMSH(13,p,grid,BIND_BLOCK,'p',0,t)
    ! advance in time
    do iStep=1,nStep
      call grid%updateDualBlock()
      call grid%updateBlockVol()
      call grid%updateFacetNorm()
      call grid%updateIntfArea()
      call grid%updateIntfNorm()
      ! Lagrangian step
      tempMass(:)=rho(:)*grid%BlockVol(:)
      forall(i=1:grid%nNode)
        tempMom(:,i)=rhou(:,i)*grid%NodeVol(i)
      end forall
      tempEnergy(:)=rhoE(:)*grid%BlockVol(:)
      do i=1,grid%nNode
        do j=1,size(grid%NodeNeibBlock(i)%dat)
          tempMom(:,i)=tempMom(:,i)&
          &            -dt*grid%NBAreaVect(i)%dat(:,j)*p(grid%NodeNeibBlock(i)%dat(j))
        end do
      end do
      uIntf=itplNode2Intf(u,grid)
      pIntf=itplBlock2Intf(p,grid)
      do i=1,grid%nIntf
        m=grid%IntfNeibBlock(1,i)
        n=grid%IntfNeibBlock(2,i)
        pWork=dt*pIntf(i)*grid%IntfArea(i)*dot_product(grid%IntfNorm(:,i),uIntf(:,i))
        tempEnergy(m)=tempEnergy(m)-pWork
        tempEnergy(n)=tempEnergy(n)+pWork
      end do
      ! apply BC
      do i=1,grid%nFacet
        do j=1,grid%Facet(i)%nNode
          k=grid%Facet(i)%iNode(j)
          tempMom(:,k)=tempMom(:,k)&
          &            -grid%FacetNorm(:,i)*dot_product(grid%FacetNorm(:,i),tempMom(:,k))
        end do
      end do
      ! recover intensive state
      rhoNode=1d0!TODO:itplBlock2Node(rho,grid)
      forall(i=1:grid%nNode)
        rhou(:,i)=tempMom(:,i)/grid%NodeVol(i)
        u(:,i)=rhou(:,i)/rhoNode(i)
      end forall
      uBlock=0d0!TODO:itplNode2Block(u,grid)
      forall(i=1:grid%nBlock)
        rhoE(i)=tempEnergy(i)/grid%BlockVol(i)
      end forall
      rho(:)=tempMass(:)/grid%BlockVol(:)
      call mvGrid(grid,dt*u)
      ! Euler rezoning
      call mvGrid(grid,-dt*u)
      ! recover state
      rhoNode=1d0!TODO:itplBlock2Node(rho,grid)
      forall(i=1:grid%nNode)
        rhou(:,i)=tempMom(:,i)/grid%NodeVol(i)
        u(:,i)=rhou(:,i)/rhoNode(i)
      end forall
      uBlock=0d0!TODO:itplNode2Block(u,grid)
      forall(i=1:grid%nBlock)
        rhoE(i)=tempEnergy(i)/grid%BlockVol(i)
        E(i)=rhoE(i)/rho(i)
        p(i)=(E(i)-dot_product(uBlock(:,i),uBlock(:,i))/2d0)*rho(i)*(gamm-1d0)
      end forall
      rho(:)=tempMass(:)/grid%BlockVol(:)
      t=t+dt
      ! write results
      call writeGMSH(13,rho,grid,BIND_BLOCK,'rho',iStep,t)
      call writeGMSH(13,u,grid,BIND_NODE,'u',iStep,t)
      call writeGMSH(13,p,grid,BIND_BLOCK,'p',iStep,t)
    end do
    deallocate(rhoNode)
    deallocate(uBlock)
    deallocate(uIntf)
    deallocate(pIntf)
    deallocate(tempMass)
    deallocate(tempMom)
    deallocate(tempEnergy)
    close(13)
  else
  end if
  call finalMPI()
end program
