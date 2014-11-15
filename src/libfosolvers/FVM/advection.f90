!----------------------------------------------------------------------------- best with 100 columns

!> advection for FVM
module modAdvection
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding advection
  interface findAdv
    module procedure::findAdvPolyVect
    module procedure::findAdvPolyScal
  end interface
  public::findAdv
  
contains
  
  !> find advection of vector s by velocity u on polyFvGrid
  !> \f[ \int_A \mathbf{s} (\mathbf{u} \cdot \hat{n}) dA \f]
  subroutine findAdvPolyVect(grid,s,u,adv)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< advected values
    double precision,intent(in)::u(:,:) !< velocity
    double precision,allocatable,intent(inout)::adv(:,:) !< advection output
    integer::up,dn
    double precision::flow(size(s,1)),fUp,fDn,df,r
    double precision,allocatable::flux(:,:),gradFlux(:,:,:)
    
    call grid%up()
    if(.not.(allocated(adv)))then
      allocate(adv(size(s,2),grid%nC))
    end if
    adv(:,:)=0d0
    allocate(flux(DIMS*size(s,1),grid%nC))
    allocate(gradFlux(DIMS,DIMS*size(s,1),grid%nC))
    forall(i=1:grid%nC)
      forall(j=1:size(s,1))
        flux(j*DIMS-(DIMS-1):j*DIMS,i)=s(j,i)*u(:,i)
      end forall
    end forall
    call findGrad(grid,flux,gradFlux)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        if(dot_product(u(:,m)+u(:,n),grid%normP(:,i))>=0d0)then
          up=m
          dn=n
        else
          up=n
          dn=m
        end if
        do j=1,size(s,1)
          fUp=dot_product(flux(j*DIMS-(DIMS-1):j*DIMS,up),grid%normP(:,i))
          fDn=dot_product(flux(j*DIMS-(DIMS-1):j*DIMS,dn),grid%normP(:,i))
          df=dot_product(matmul(grid%p(dn)-grid%p(up),gradFlux(:,j*DIMS-(DIMS-1):j*DIMS,up)),&
          &              grid%normP(:,i))
          r=merge(2d0*df/(fDn-fUp)-1d0,fUp,abs(fDn-fUp)>tiny(1d0))
          flow(j)=-grid%aP(i)*(fUp+0.5d0*vanAlbada(r)*(fDn-fUp))
        end do
        adv(:,m)=adv(:,m)+flow(:)
        adv(:,n)=adv(:,n)-flow(:)
      end if
    end do
    deallocate(flux)
    deallocate(gradFlux)
  end subroutine
  
  !> find advection of scalar s by velocity u on polyFvGrid
  !> \f[ \int_A s (\mathbf{u} \cdot \hat{n}) dA \f]
  subroutine findAdvPolyScal(grid,s,u,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< advected value
    double precision,intent(in)::u(:,:) !< velocity
    double precision,allocatable,intent(inout)::adv(:) !< advection output
    double precision,allocatable::sv(:,:),advv(:,:)
    
    allocate(sv(1,size(s)))
    if(allocated(adv))then
      allocate(advv(1,size(adv)))
    end if
    sv(1,:)=s(:)
    call findAdvPolyVect(grid,sv,u,advv)
    if(.not.allocated(adv))then
      allocate(adv(size(advv,2)),source=advv(1,:))!FIXME:remove work-around
    else
      adv(:)=advv(1,:)
    end if
    deallocate(sv)
    deallocate(advv)
  end subroutine
  
  !> van Albada flux limiter
  elemental function vanAlbada(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision::vanAlbada !< limit function
    
    vanAlbada=(r+r**2)/(1d0+r**2)
  end function
  
end module
