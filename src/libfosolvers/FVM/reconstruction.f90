!----------------------------------------------------------------------------- best with 100 columns

!> interface reconstruction for FVM
module modReconstruction
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic average reconstruction
  interface reconAvg
    module procedure::reconAvgPolyVect
    module procedure::reconAvgPolyScal
  end interface
  public::reconAvg
  
  !> generic upwind reconstruction
  interface reconUW
    module procedure::reconUWPolyVect
    module procedure::reconUWPolyScal
  end interface
  public::reconUW
  
  !> generic TVD-limited reconstruction
  interface reconLtd
    module procedure::reconLtdPolyVect
    module procedure::reconLtdPolyScal
  end interface
  public::reconLtd
  
contains
  
  !> average reconstruction of vector s on PolyFvGrid
  subroutine reconAvgPolyVect(grid,s,grads,sr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< input data
    double precision,intent(in)::grads(:,:,:) !< gradient of s
    double precision,allocatable,intent(inout)::sr(:,:) !< output of s reconstructed on interfaces
    double precision::sm(size(s,1)),sn(size(s,1)),vm(DIMS),vn(DIMS)
    
    call grid%up()
    if(.not.(allocated(sr)))then
      allocate(sr(size(s,1),grid%nP))
    end if
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      vm(:)=grid%pP(:,i)-grid%p(m)
      vn(:)=grid%pP(:,i)-grid%p(n)
      if(m<=grid%nC)then
        sm(:)=s(:,m)+matmul(vm(:),grads(:,:,m))
      end if
      if(n<=grid%nC)then
        sn(:)=s(:,n)+matmul(vn(:),grads(:,:,n))
      end if
      if(m<=grid%nC.and.n<=grid%nC)then
        sr(:,i)=(norm2(vn)*sm(:)+norm2(vm)*sn(:))/(norm2(vm)+norm2(vn))
      else if(m<=grid%nC)then
        sr(:,i)=sm(:)
      else if(n<=grid%nC)then
        sr(:,i)=sn(:)
      else
        write(*,'(a)'),'[E] reconAvgPolyVect: must reconstruct from cell element'
      end if
    end do
  end subroutine
  
  !> average reconstruction of scalar s on PolyFvGrid
  subroutine reconAvgPolyScal(grid,s,grads,sr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< input data
    double precision,intent(in)::grads(:,:) !< gradient of s
    double precision,allocatable,intent(inout)::sr(:) !< output of s reconstructed on interfaces
    double precision,allocatable::sv(:,:),gradsv(:,:,:),srv(:,:)
    
    allocate(sv(1,size(s)))
    allocate(gradsv(DIMS,1,size(grads,2)))
    if(allocated(sr))then
      allocate(srv(1,size(sr)))
    end if
    sv(1,:)=s(:)
    gradsv(:,1,:)=grads(:,:)
    call reconAvgPolyVect(grid,sv,gradsv,srv)
    if(.not.allocated(sr))then
      allocate(sr(size(srv,2)),source=srv(1,:))!FIXME:remove work-around
    else
      sr(:)=srv(1,:)
    end if
    deallocate(sv)
    deallocate(gradsv)
    deallocate(srv)
  end subroutine
  
  !> upwind reconstruction of vector s on PolyFvGrid
  subroutine reconUWPolyVect(grid,s,ur,sr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< input data
    double precision,intent(in)::ur(:,:) !< reconstructed velocity
    double precision,allocatable,intent(inout)::sr(:,:) !< output of s reconstructed on interfaces
    
    call grid%up()
    if(.not.(allocated(sr)))then
      allocate(sr(size(s,1),grid%nP))
    end if
    do i=1,grid%nP
      m=grid%iEP(merge(1,2,dot_product(ur(:,i),grid%normP(:,i))>=0d0),i)
      if(m>grid%nC)then
        m=grid%iEP(1,i)
      end if
      sr(:,i)=s(:,m)
    end do
  end subroutine
  
  !> upwind reconstruction of scalar s on PolyFvGrid
  subroutine reconUWPolyScal(grid,s,ur,sr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< input data
    double precision,intent(in)::ur(:,:) !< reconstructed velocity
    double precision,allocatable,intent(inout)::sr(:) !< output of s reconstructed on interfaces
    double precision,allocatable::sv(:,:),srv(:,:)
    
    allocate(sv(1,size(s)))
    if(allocated(sr))then
      allocate(srv(1,size(sr)))
    end if
    sv(1,:)=s(:)
    call reconUWPolyVect(grid,sv,ur,srv)
    if(.not.allocated(sr))then
      allocate(sr(size(srv,2)),source=srv(1,:))!FIXME:remove work-around
    else
      sr(:)=srv(1,:)
    end if
    deallocate(sv)
    deallocate(srv)
  end subroutine
  
  !> TVD-limited reconstruction of vector s on PolyFvGrid
  subroutine reconLtdPolyVect(grid,s,grads,ur,sr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< input data
    double precision,intent(in)::grads(:,:,:) !< gradient of s
    double precision,intent(in)::ur(:,:) !< reconstructed velocity
    double precision,allocatable,intent(inout)::sr(:,:) !< output of s reconstructed on interfaces
    double precision::r(size(s,1))
    integer::up,dn
    
    call grid%up()
    if(.not.(allocated(sr)))then
      allocate(sr(size(s,1),grid%nP))
    end if
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        if(dot_product(ur(:,i),grid%normP(:,i))>=0d0)then
          up=m
          dn=n
        else
          up=n
          dn=m
        end if
        r(:)=0d0
        forall(j=1:size(r),abs(s(j,dn)-s(j,up))>tiny(1d0))
          r(j)=2d0*dot_product(grid%p(dn)-grid%p(up),grads(:,j,up))/(s(j,dn)-s(j,up))-1d0
        end forall
        sr(:,i)=s(:,up)+vanLeer(r)*(s(:,dn)-s(:,up))/2d0
      else
        sr(:,i)=0d0
      end if
    end do
  end subroutine
  
  !> TVD-limited reconstruction of scalar s on PolyFvGrid
  subroutine reconLtdPolyScal(grid,s,grads,ur,sr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< input data
    double precision,intent(in)::grads(:,:) !< gradient of s
    double precision,intent(in)::ur(:,:) !< reconstructed velocity
    double precision,allocatable,intent(inout)::sr(:) !< output of s reconstructed on interfaces
    double precision,allocatable::sv(:,:),gradsv(:,:,:),srv(:,:)
    
    allocate(sv(1,size(s)))
    allocate(gradsv(DIMS,1,size(grads,2)))
    if(allocated(sr))then
      allocate(srv(1,size(sr)))
    end if
    sv(1,:)=s(:)
    gradsv(:,1,:)=grads(:,:)
    call reconLtdPolyVect(grid,sv,gradsv,ur,srv)
    if(.not.allocated(sr))then
      allocate(sr(size(srv,2)),source=srv(1,:))!FIXME:remove work-around
    else
      sr(:)=srv(1,:)
    end if
    deallocate(sv)
    deallocate(gradsv)
    deallocate(srv)
  end subroutine
  
  !> van Leer flux limiter
  elemental function vanLeer(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision::vanLeer !< limit function
    
    vanLeer=(r+abs(r))/(1d0+abs(r))
  end function
  
  !> minmod flux limiter
  elemental function minmod(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision minmod !< limit function
    
    minmod=max(0d0,min(r,1d0))
  end function
  
end module
