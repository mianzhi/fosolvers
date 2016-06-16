!----------------------------------------------------------------------------- best with 100 columns

!> pressure force for FVM
module modPressure
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding pressure force
  interface findPresForce
    module procedure::findPresForcePoly
  end interface
  public::findPresForce
  
contains
  
  !> find pressure force on polyFvGrid
  !> \f[ \int_A p \hat{n} dA \f]
  subroutine findPresForcePoly(grid,p,dRhou)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::p(:) !< pressure
    double precision,allocatable,intent(inout)::dRhou(:,:) !< force (net flow of rhou)
    double precision::flow(DIMS),pf,vMP(DIMS),vPN(DIMS),eps
    
    call grid%up()
    if(.not.(allocated(dRhou)))then
      allocate(dRhou(DIMS,grid%nC))
    end if
    dRhou(:,:)=0d0
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=size(p).and.n<=size(p))then
        if(n<=grid%nC)then ! internal pairs
          vMP(:)=grid%pP(:,i)-grid%p(:,m)
          vPN(:)=grid%p(:,n)-grid%pP(:,i)
          eps=norm2(vPN)/(norm2(vMP)+norm2(vPN))
          pf=eps*p(m)+(1d0-eps)*p(n)
          flow(:)=-pf*grid%aP(i)*grid%normP(:,i)
          dRhou(:,m)=dRhou(:,m)+flow(:)
          dRhou(:,n)=dRhou(:,n)-flow(:)
        else ! boundary pairs
          pf=0.5d0*(p(m)+p(n))
          flow(:)=-pf*grid%aP(i)*grid%normP(:,i)
          dRhou(:,m)=dRhou(:,m)+flow(:)
        end if
      end if
    end do
  end subroutine
  
end module
