!----------------------------------------------------------------------------- best with 100 columns

!> diffusion for FVM
module modDiffusion
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> generic finding diffusion
  interface findDfs
    module procedure::findDfsPolyVect
    module procedure::findDfsPolyScal
  end interface
  public::findDfs
  
contains

  !> find diffusion rate of vector s by diffusivity gamma on polyFvGrid
  !> \f[ \int_A (\mathbf{I} \mathbf{\Gamma}) (\nabla \mathbf{\mathbf{s}} \cdot \hat{n}) dA \f]
  !> NOTE: gamma, grads are reconstructed data
  subroutine findDfsPolyVect(grid,gam,grads,dfs)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::gam(:,:) !< reconstructed input data
    double precision,intent(in)::grads(:,:,:) !< reconstructed velocity field
    double precision,allocatable,intent(inout)::dfs(:,:) !< diffusion output
    double precision::F(size(grads,2))
    
    call grid%up()
    if(.not.(allocated(dfs)))then
      allocate(dfs(size(grads,2),grid%nC))
    end if
    dfs(:,:)=0d0
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(m<=grid%nC.and.n<=grid%nC)then
        F(:)=grid%aP(i)*gam(:,i)*matmul(grid%normP(:,i),grads(:,:,i))
        dfs(:,m)=dfs(:,m)+F(:)
        dfs(:,n)=dfs(:,n)-F(:)
      end if
    end do
  end subroutine
  
  !> find diffusion rate of scalar s by diffusivity gamma on polyFvGrid
  !> \f[ \int_A \Gamma (\nabla \mathbf{s} \cdot \hat{n}) dA \f]
  !> NOTE: gamma, grads are reconstructed data
  subroutine findDfsPolyScal(grid,gam,grads,dfs)
    use modPolyFvGrid
    class(polyFvGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::gam(:) !< reconstructed input data
    double precision,intent(in)::grads(:,:) !< reconstructed velocity field
    double precision,allocatable,intent(inout)::dfs(:) !< diffusion output
    double precision,allocatable::gamv(:,:),gradsv(:,:,:),dfsv(:,:)
    
    allocate(gamv(1,size(gam)))
    allocate(gradsv(DIMS,1,size(grads,2)))
    if(allocated(dfs))then
      allocate(dfsv(1,size(dfs)))
    end if
    gamv(1,:)=gam(:)
    gradsv(:,1,:)=grads(:,:)
    call findDfsPolyVect(grid,gamv,gradsv,dfsv)
    if(.not.allocated(dfs))then
      allocate(dfs(size(dfsv,2)),source=dfsv(1,:))!FIXME:remove work-around
    else
      dfs(:)=dfsv(1,:)
    end if
    deallocate(gamv)
    deallocate(gradsv)
    deallocate(dfsv)
  end subroutine
  
end module
