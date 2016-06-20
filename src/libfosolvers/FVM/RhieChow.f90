!----------------------------------------------------------------------------- best with 100 columns

!> Rhie-Chow (pressure-induced) flow for FVM
module modRhieChow
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic adding Rhie-Chow flow
  interface addRhieChow
    module procedure::addRhieChowPolyVect
    module procedure::addRhieChowPolyScal
  end interface
  public::addRhieChow
  
contains
  
  !> add Rhie-Chow flow of vector s to r on polyFvGrid
  subroutine addRhieChowPolyVect(grid,s,p,gradP,rho,dt,r)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:,:) !< state variables
    double precision,intent(in)::p(:) !< pressure
    double precision,intent(in)::gradP(:,:) !< pressure gradient
    double precision,intent(in)::rho(:) !< density
    double precision,intent(in)::dt !< time step size
    double precision,intent(inout)::r(:,:) !< target of addition
    double precision::flow(size(s,1)),vMN(DIMS),vMP(DIMS),vPN(DIMS),eps
    
    call grid%up()
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        vMN(:)=grid%p(:,n)-grid%p(:,m)
        vMP(:)=grid%pP(:,i)-grid%p(:,m)
        vPN(:)=grid%p(:,n)-grid%pP(:,i)
        eps=norm2(vPN)/(norm2(vMP)+norm2(vPN))
        flow(:)=-(eps*s(:,m)+(1d0-eps)*s(:,n))/(eps*rho(m)+(1d0-eps)*rho(n))*dt*grid%aP(i)&
        &        *dot_product(grid%normP(:,i),eps*gradP(:,m)+(1d0-eps)*gradP(:,n)&
        &                                     -vMN(:)*(p(n)-p(m))/dot_product(vMN,vMN))
        r(:,m)=r(:,m)+flow(:)
        r(:,n)=r(:,n)-flow(:)
      end if
    end do
  end subroutine
  
  !> add Rhie-Chow flow of scalar s to r on polyFvGrid
  subroutine addRhieChowPolyScal(grid,s,p,gradP,rho,dt,r)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::s(:) !< state variable
    double precision,intent(in)::p(:) !< pressure
    double precision,intent(in)::gradP(:,:) !< pressure gradient
    double precision,intent(in)::rho(:) !< density
    double precision,intent(in)::dt !< time step size
    double precision,intent(inout)::r(:) !< target of addition
    double precision::flow,vMN(DIMS),vMP(DIMS),vPN(DIMS),eps
    
    call grid%up()
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then ! internal pairs
        vMN(:)=grid%p(:,n)-grid%p(:,m)
        vMP(:)=grid%pP(:,i)-grid%p(:,m)
        vPN(:)=grid%p(:,n)-grid%pP(:,i)
        eps=norm2(vPN)/(norm2(vMP)+norm2(vPN))
        flow=-(eps*s(m)+(1d0-eps)*s(n))/(eps*rho(m)+(1d0-eps)*rho(n))*dt*grid%aP(i)&
        &     *dot_product(grid%normP(:,i),eps*gradP(:,m)+(1d0-eps)*gradP(:,n)&
        &                                  -vMN(:)*(p(n)-p(m))/norm2(vMN))
        r(m)=r(m)+flow
        r(n)=r(n)-flow
      end if
    end do
  end subroutine
  
end module
