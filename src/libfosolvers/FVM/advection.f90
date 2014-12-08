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
  
  !> find advection due to flux f depending on vector s on polyFvGrid
  !> \f[ \int_A \mathbf{f}(\mathbf{s}) \cdot \hat{n} dA \f]
  subroutine findAdvPolyVect(grid,s,f,adv)
    use modPolyFvGrid
    use modGradient
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< state variables
    double precision,intent(in)::f(:,:,:) !< fluxes
    double precision,allocatable,intent(inout)::adv(:,:) !< advection output
    integer::up,dn
    double precision::flow(size(s,1)),fUp,fDn,df,r
    double precision,allocatable::gradF(:,:,:)
    
    call grid%up()
    if(.not.(allocated(adv)))then
      allocate(adv(size(s,2),grid%nC))
    end if
    adv(:,:)=0d0
    allocate(gradF(DIMS,DIMS*size(s,1),grid%nC))
    call findGrad(grid,reshape(f,[size(f,1)*size(f,2),size(f,3)]),gradF)
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        do j=1,size(s,1)
          if(abs(s(j,m)-s(j,n))<=tiny(1d0))then
            up=m
            dn=n
          else if(dot_product((f(:,j,m)-f(:,j,n))/(s(j,m)-s(j,n)),grid%normP(:,i))>=0d0)then
            up=m
            dn=n
          else
            up=n
            dn=m
          end if
          fUp=dot_product(f(:,j,up),grid%normP(:,i))
          fDn=dot_product(f(:,j,dn),grid%normP(:,i))
          df=dot_product(matmul(grid%p(:,dn)-grid%p(:,up),gradF(:,j*DIMS-(DIMS-1):j*DIMS,up)),&
          &              grid%normP(:,i))
          r=merge(2d0*df/(fDn-fUp)-1d0,0d0,abs(fDn-fUp)>tiny(1d0))
          flow(j)=-grid%aP(i)*(fUp+0.5d0*vanAlbada(r)*(fDn-fUp))
        end do
        adv(:,m)=adv(:,m)+flow(:)
        adv(:,n)=adv(:,n)-flow(:)
      end if
    end do
    deallocate(gradF)
  end subroutine
  
  !> find advection due to flux f depending on scalar s on polyFvGrid
  !> \f[ \int_A \mathbf{f}(s) \cdot \hat{n} dA \f]
  subroutine findAdvPolyScal(grid,s,f,adv)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< state variable
    double precision,intent(in)::f(:,:) !< flux
    double precision,allocatable,intent(inout)::adv(:) !< advection output
    double precision,allocatable::sv(:,:),fv(:,:,:),advv(:,:)
    
    allocate(sv(1,size(s)))
    allocate(fv(DIMS,1,size(s)))
    if(allocated(adv))then
      allocate(advv(1,size(adv)))
    end if
    sv(1,:)=s(:)
    fv(:,1,:)=f(:,:)
    call findAdvPolyVect(grid,sv,fv,advv)
    if(.not.allocated(adv))then
      allocate(adv(size(advv,2)),source=advv(1,:))!FIXME:remove work-around
    else
      adv(:)=advv(1,:)
    end if
    deallocate(sv)
    deallocate(fv)
    deallocate(advv)
  end subroutine
  
  !> van Albada flux limiter
  elemental function vanAlbada(r)
    double precision,intent(in)::r !< successive gradient ratio
    double precision::vanAlbada !< limit function
    double precision,parameter::R_LMT=1d20
    
    vanAlbada=merge((r+r**2)/(1d0+r**2),1d0,abs(r)<R_LMT)
  end function
  
end module
