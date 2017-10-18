!----------------------------------------------------------------------------- best with 100 columns

!> gradient for FVM
module modGradient
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding gradient
  interface findGrad
    module procedure::findGradPolyVect
    module procedure::findGradPolyScal
    module procedure::findGradOtVect
    module procedure::findGradOtScal
  end interface
  public::findGrad
  
contains
  
  !> find gradient of vector v on polyFvGrid
  subroutine findGradPolyVect(grid,v,gradv)
    use modPolyFvGrid
    use modPolyGrid
    use modGeometry
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:,:) !< input data
    double precision,allocatable,intent(inout)::gradv(:,:,:) !< gradient output
    logical::useCellOnly
    integer,parameter::MAX_N_NEAR=25
    integer,parameter::MAX_L_WORK=2000
    integer,parameter::L_IWORK=200
    integer::nNear,iNear(MAX_N_NEAR),lwork,iwork(L_IWORK),rank,ier
    double precision::dx(MAX_N_NEAR,DIMS),dv(MAX_N_NEAR,size(v,1)),stat(DIMS),work(MAX_L_WORK),rcond
    double precision::tmp,eps,tmpP(DIMS)
    
    call grid%up()
    m=size(v,1) ! number of components
    ! handling non-cell elements
    useCellOnly=size(v,2)<grid%nE
    if(.not.allocated(gradv))then
      allocate(gradv(DIMS,m,grid%nC))
    end if
    ! find gradient
    !$omp parallel do default(shared)&
    !$omp& private(nNear,iNear,dx,dv,vF,lwork,work,iwork,stat,rcond,rank,ier,j)
    do i=1,grid%nC
      ! process the list of nearby points
      nNear=0
      do j=1,size(grid%near,1)
        if(grid%near(j,i)>0.and.grid%near(j,i)<=merge(grid%nC,grid%nE,useCellOnly))then
          nNear=nNear+1
          iNear(nNear)=grid%near(j,i)
        else if(grid%near(j,i)==0)then
          exit
        end if
      end do
      ! LSQ method
      if(.not.useCellOnly)then ! use gradient mapping (constructed by DEGSDD)
        forall(j=1:nNear)
          dv(j,:)=v(:,iNear(j))-v(:,i)
        end forall
        gradv(1:DIMS,1:m,i)=matmul(grid%gradMap(1:DIMS,1:nNear,i),dv(1:nNear,1:m))
      else ! use DGELSD
        forall(j=1:nNear)
          dx(j,:)=grid%p(:,iNear(j))-grid%p(:,i)
          dv(j,:)=v(:,iNear(j))-v(:,i)
          dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
          dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
        end forall
        rcond=-1d0
        lwork=-1
        call DGELSD(nNear,DIMS,m,dx,MAX_N_NEAR,dv,MAX_N_NEAR,stat,rcond,rank,work,lwork,iwork,ier)
        lwork=min(MAX_L_WORK,int(work(1)))
        call DGELSD(nNear,DIMS,m,dx,MAX_N_NEAR,dv,MAX_N_NEAR,stat,rcond,rank,work,lwork,iwork,ier)
        gradv(:,:,i)=dv(1:DIMS,:)
      end if
    end do
    !$end omp parallel do
    ! limit the gradient
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      eps=norm2(grid%pP(:,i)-grid%p(:,n))&
      &   /(norm2(grid%pP(:,i)-grid%p(:,m))+norm2(grid%pP(:,i)-grid%p(:,n)))
      tmpP(:)=eps*grid%p(:,m)+(1d0-eps)*grid%p(:,n)
      if(m<=grid%nC.and.n<=grid%nC)then
        do j=1,size(gradv,2)
          tmp=v(j,m)+dot_product(gradv(:,j,m),tmpP(:)-grid%p(:,m))
          if(tmp>max(v(j,m),v(j,n)))then
            gradv(:,j,m)=gradv(:,j,m)&
            &            -(tmpP(:)-grid%p(:,m))*(tmp-max(v(j,m),v(j,n)))&
            &             /dot_product(tmpP(:)-grid%p(:,m),tmpP(:)-grid%p(:,m))
          else if(tmp<min(v(j,m),v(j,n)))then
            gradv(:,j,m)=gradv(:,j,m)&
            &            -(tmpP(:)-grid%p(:,m))*(tmp-min(v(j,m),v(j,n)))&
            &             /dot_product(tmpP(:)-grid%p(:,m),tmpP(:)-grid%p(:,m))
          end if
          tmp=v(j,n)+dot_product(gradv(:,j,n),tmpP(:)-grid%p(:,n))
          if(tmp>max(v(j,m),v(j,n)))then
            gradv(:,j,n)=gradv(:,j,n)&
            &            -(tmpP(:)-grid%p(:,n))*(tmp-max(v(j,m),v(j,n)))&
            &             /dot_product(tmpP(:)-grid%p(:,n),tmpP(:)-grid%p(:,n))
          else if(tmp<min(v(j,m),v(j,n)))then
            gradv(:,j,n)=gradv(:,j,n)&
            &            -(tmpP(:)-grid%p(:,n))*(tmp-min(v(j,m),v(j,n)))&
            &             /dot_product(tmpP(:)-grid%p(:,n),tmpP(:)-grid%p(:,n))
          end if
        end do
      end if
    end do
  end subroutine
  
  !> find gradient of scalar v on polyFvGrid
  subroutine findGradPolyScal(grid,v,gradv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:) !< input data
    double precision,allocatable,intent(inout)::gradv(:,:) !< gradient output
    double precision,allocatable::vv(:,:),gradvv(:,:,:)
    
    allocate(vv(1,size(v)))
    if(allocated(gradv))then
      allocate(gradvv(DIMS,1,size(gradv,2)))
    end if
    vv(1,:)=v(:)
    call findGradPolyVect(grid,vv,gradvv)
    if(.not.allocated(gradv))then
      allocate(gradv(DIMS,size(gradvv,3)),source=gradvv(:,1,:))!FIXME:remove work-around
    else
      gradv(:,:)=gradvv(:,1,:)
    end if
    deallocate(vv)
    deallocate(gradvv)
  end subroutine
  
  !> find gradient of vector v on otGrid
  subroutine findGradOtVect(grid,v,gradv)
    use modOtGrid
    class(otGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:,:) !< input data
    double precision,allocatable,intent(inout)::gradv(:,:,:) !< gradient output
    
    m=size(v,1) ! number of components
    if(.not.allocated(gradv))then
      allocate(gradv(DIMS,m,grid%nC))
    end if
    gradv(:,:,:)=0d0
    do i=1,grid%nC
      if(grid%flag(i)/=OT_OUT)then
        do l=1,DIMS
          if(grid%neib(l*2-1,i)>0.and.grid%neib(l*2,i)>0)then
            !TODO: different level?
            gradv(l,:,i)=(v(:,grid%neib(l*2,i))-v(:,grid%neib(l*2-1,i)))/2d0/grid%h(i)
          else if(grid%neib(l*2-1,i)>0)then
            !TODO: different level?
            gradv(l,:,i)=(v(:,i)-v(:,grid%neib(l*2-1,i)))/grid%h(i)
          else if(grid%neib(l*2,i)>0)then
            !TODO: different level?
            gradv(l,:,i)=(v(:,grid%neib(l*2,i))-v(:,i))/grid%h(i)
          end if
        end do
      end if
    end do
  end subroutine
  
  !> find gradient of scalar v on otGrid
  subroutine findGradOtScal(grid,v,gradv)
    use modOtGrid
    class(otGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:) !< input data
    double precision,allocatable,intent(inout)::gradv(:,:) !< gradient output
    double precision,allocatable::vv(:,:),gradvv(:,:,:)
    
    allocate(vv(1,size(v)))
    if(allocated(gradv))then
      allocate(gradvv(DIMS,1,size(gradv,2)))
    end if
    vv(1,:)=v(:)
    call findGradOtVect(grid,vv,gradvv)
    if(.not.allocated(gradv))then
      allocate(gradv(DIMS,size(gradvv,3)),source=gradvv(:,1,:))!FIXME:remove work-around
    else
      gradv(:,:)=gradvv(:,1,:)
    end if
    deallocate(vv)
    deallocate(gradvv)
  end subroutine
  
end module
