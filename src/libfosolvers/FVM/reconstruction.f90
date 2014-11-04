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
  
  !> generic upwind reconstruction for advection
  interface reconUW
    module procedure::reconUWPolyVect
    module procedure::reconUWPolyScal
  end interface
  public::reconUW
  
  !> generic TVD-limited reconstruction for advection
  interface reconLtd
    module procedure::reconLtdPolyVect
    module procedure::reconLtdPolyScal
  end interface
  public::reconLtd
  
  !> generic surface decomposition reconstruction for diffusion
  interface reconSD
    module procedure::reconSDPolyVect
    module procedure::reconSDPolyScal
  end interface
  public::reconSD
  
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
        sr(:,i)=s(:,up)+vanAlbada(r)*(s(:,dn)-s(:,up))/2d0
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
  
  !> surface decomposition reconstruction of gradient of vector s on PolyFvGrid
  !> NOTE: gradients from SD reconstruction are in the normal direction
  subroutine reconSDPolyVect(grid,s,grads,gradsr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< input data
    double precision,intent(in)::grads(:,:,:) !< gradient of s
    double precision,allocatable,intent(inout)::gradsr(:,:,:) !< grads reconstructed on interfaces
    double precision::dMN,sf(DIMS),tf(DIMS),rf(DIMS),gradavg(DIMS,size(s,1)),tmp(size(s,1)),&
    &                 vm(DIMS),vn(DIMS)
    
    call grid%up()
    if(.not.(allocated(gradsr)))then
      allocate(gradsr(DIMS,size(s,1),grid%nP))
    end if
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        vm(:)=grid%pP(:,i)-grid%p(m)
        vn(:)=grid%pP(:,i)-grid%p(n)
        gradavg(:,:)=(norm2(vn)*grads(:,:,m)+norm2(vm)*grads(:,:,n))/(norm2(vm)+norm2(vn))
        sf(:)=grid%p(n)-grid%p(m)
        dMN=norm2(sf)
        sf(:)=sf(:)/dMN
        k=minloc(abs(sf(:)),1)
        select case(k)
        case(1)
          tf(:)=[0d0,grid%pP(3,i),-grid%pP(2,i)]
        case(2)
          tf(:)=[grid%pP(3,i),0d0,-grid%pP(1,i)]
        case(3)
          tf(:)=[grid%pP(2,i),-grid%pP(1,i),0d0]
        end select
        tf(:)=tf(:)/norm2(tf)
        rf(1)=grid%normP(2,i)*tf(3)-grid%normP(3,i)*tf(2)
        rf(2)=-grid%normP(1,i)*tf(3)+grid%normP(3,i)*tf(1)
        rf(3)=grid%normP(1,i)*tf(2)-grid%normP(2,i)*tf(1)
        tmp(:)=((s(:,n)-s(:,m))/dMN-matmul(dot_product(sf,tf)*tf+dot_product(sf,rf)*rf,gradavg))
        forall(j=1:size(s,1))
          gradsr(:,j,i)=tmp(j)*grid%normP(:,i)/dot_product(sf(:),grid%normP(:,i))
        end forall
      else
        gradsr(:,:,i)=0d0
      end if
    end do
  end subroutine
  
  !> surface decomposition reconstruction of gradient of scalar s on PolyFvGrid
  !> NOTE: gradients from SD reconstruction are in the normal direction
  subroutine reconSDPolyScal(grid,s,grads,gradsr)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< input data
    double precision,intent(in)::grads(:,:) !< gradient of s
    double precision,allocatable,intent(inout)::gradsr(:,:) !< grads reconstructed on interfaces
    double precision,allocatable::sv(:,:),gradsv(:,:,:),gradsrv(:,:,:)
    
    allocate(sv(1,size(s)))
    allocate(gradsv(DIMS,1,size(grads,2)))
    if(allocated(gradsr))then
      allocate(gradsrv(DIMS,1,size(gradsr,2)))
    end if
    sv(1,:)=s(:)
    gradsv(:,1,:)=grads(:,:)
    call reconSDPolyVect(grid,sv,gradsv,gradsrv)
    if(.not.allocated(gradsr))then
      allocate(gradsr(DIMS,size(gradsrv,3)),source=gradsrv(:,1,:))!FIXME:remove work-around
    else
      gradsr(:,:)=gradsrv(:,1,:)
    end if
    deallocate(sv)
    deallocate(gradsv)
    deallocate(gradsrv)
  end subroutine
  
  !> van Leer flux limiter
  elemental function vanLeer(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision::vanLeer !< limit function
    
    vanLeer=(r+abs(r))/(1d0+abs(r))
  end function
  
  !> van Albada flux limiter
  elemental function vanAlbada(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision::vanAlbada !< limit function
    
    vanAlbada=(r+r**2)/(1d0+r**2)
  end function
  
end module
