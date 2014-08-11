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
  end interface
  public::findGrad
  
contains
  
  !> find gradient of vector v on polyFvGrid
  subroutine findGradPolyVect(grid,v,gradv)
    use modPolyFvGrid
    use modPolyGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(:,:) !< input data
    double precision,allocatable,intent(inout)::gradv(:,:,:) !< gradient output
    logical::useCellOnly,findCellOnly
    integer,parameter::MIN_N_NEIB=8
    integer,parameter::MAX_N_NEIB=25
    integer,parameter::MAX_L_WORK=2000
    integer,parameter::L_IWORK=200
    integer::nNeib,iNeib(MAX_N_NEIB),lwork,iwork(L_IWORK),ier
    double precision::dx(MAX_N_NEIB,DIMS),dv(MAX_N_NEIB,size(v,1)),stat(DIMS),work(MAX_L_WORK),rcond
    double precision::tmp
    
    call grid%up()
    m=size(v,1) ! number of components
    ! handling non-cell elements
    useCellOnly=.not.(size(v,2)>=grid%nE)
    findCellOnly=.true.
    if(allocated(gradv))then
      if(size(gradv,3)>=grid%nE)then
        findCellOnly=.false.
      end if
    else
      findCellOnly=useCellOnly
      allocate(gradv(DIMS,m,merge(grid%nC,grid%nE,findCellOnly)))
    end if
    ! find gradient
    do i=1,grid%nE
      if(i>grid%nC.and.findCellOnly)then
        exit
      end if
      ! list of neighbors
      if(i>grid%nC)then ! at non-cell
        if(findCellOnly)then
          exit
        end if
        write(*,*),'findGradPolyVect: support cell element only'
        !TODO:nNeib,lNeib at non-cell
        !TODO check useCellOnly
      else ! at cell
        nNeib=0
        do j=1,nF(grid%sE(i))
          if(nNeib<MAX_N_NEIB.and.grid%neib(j,i)>0&
          &  .and.(grid%neib(j,i)<=grid%nC.or..not.useCellOnly))then
            nNeib=nNeib+1
            iNeib(nNeib)=grid%neib(j,i)
          end if
        end do
        ! find neighbors if needed
        do while(nNeib<MIN_N_NEIB)
          n=nNeib
          do k=1,n
            l=iNeib(k)
            do j=1,nF(grid%sE(l))
              if(nNeib<MAX_N_NEIB.and.grid%neib(j,l)>0&
              &  .and.(grid%neib(j,l)<=grid%nC.or..not.useCellOnly)&
              &  .and.all(iNeib(1:nNeib)/=grid%neib(j,l))&
              &  .and.grid%neib(j,l)/=i)then
                nNeib=nNeib+1
                iNeib(nNeib)=grid%neib(j,l)
              end if
            end do
          end do
        end do
      end if
      ! LSQ method
      forall(j=1:nNeib)
        dx(j,:)=grid%p(iNeib(j))-grid%p(i)
        dv(j,:)=v(:,iNeib(j))-v(:,i)
        dv(j,:)=dv(j,:)/dot_product(dx(j,:),dx(j,:))
        dx(j,:)=dx(j,:)/dot_product(dx(j,:),dx(j,:))
      end forall
      rcond=-1d0
      lwork=-1
      call DGELSD(nNeib,DIMS,m,dx,MAX_N_NEIB,dv,MAX_N_NEIB,stat,rcond,rank,work,lwork,iwork,ier)
      lwork=min(MAX_L_WORK,int(work(1)))
      call DGELSD(nNeib,DIMS,m,dx,MAX_N_NEIB,dv,MAX_N_NEIB,stat,rcond,rank,work,lwork,iwork,ier)
      gradv(:,:,i)=dv(1:DIMS,:)
    end do
    ! limit the gradient
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        do j=1,size(gradv,2)
          tmp=v(j,m)+dot_product(gradv(:,j,m),grid%pP(:,i)-grid%p(m))
          if(tmp>max(v(j,m),v(j,n)))then
            gradv(:,j,m)=gradv(:,j,m)&
            &            -(grid%pP(:,i)-grid%p(m))*(tmp-max(v(j,m),v(j,n)))&
            &             /dot_product(grid%pP(:,i)-grid%p(m),grid%pP(:,i)-grid%p(m))
          else if(tmp<min(v(j,m),v(j,n)))then
            gradv(:,j,m)=gradv(:,j,m)&
            &            -(grid%pP(:,i)-grid%p(m))*(tmp-min(v(j,m),v(j,n)))&
            &             /dot_product(grid%pP(:,i)-grid%p(m),grid%pP(:,i)-grid%p(m))
          end if
          tmp=v(j,n)+dot_product(gradv(:,j,n),grid%pP(:,i)-grid%p(n))
          if(tmp>max(v(j,m),v(j,n)))then
            gradv(:,j,n)=gradv(:,j,n)&
            &            -(grid%pP(:,i)-grid%p(n))*(tmp-max(v(j,m),v(j,n)))&
            &             /dot_product(grid%pP(:,i)-grid%p(n),grid%pP(:,i)-grid%p(n))
          else if(tmp<min(v(j,m),v(j,n)))then
            gradv(:,j,n)=gradv(:,j,n)&
            &            -(grid%pP(:,i)-grid%p(n))*(tmp-min(v(j,m),v(j,n)))&
            &             /dot_product(grid%pP(:,i)-grid%p(n),grid%pP(:,i)-grid%p(n))
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
  
end module
