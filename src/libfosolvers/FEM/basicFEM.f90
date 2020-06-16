!----------------------------------------------------------------------------- best with 100 columns

!> finite element discretizations
module modBasicFEM
  public
  
contains
  
  !> find Laplacian operator
  !> NOTE: A must be pre-initialized to contain the temporary CSR matrix with duplications
  subroutine findLaplacian(grid,A,isDirichlet)
    use modPolyFeGrid
    use modPolyGrid
    use modSMQ
    use modSparse
    class(polyFeGrid),intent(inout)::grid !< the grid
    class(linEq),intent(inout)::A !< the discretized linear system
    logical,intent(in),optional::isDirichlet(grid%nN) !< whether each node is Dirichlet node 
    integer::nnz,maxNNZ
    integer,allocatable::iA(:),jA(:)
    double precision,allocatable::vA(:)
    
    call grid%up()
    maxNNZ=size(grid%iNE,1)**2*grid%nC
    allocate(iA(maxNNZ)) ! COO temporary matrix
    allocate(jA(maxNNZ))
    allocate(vA(maxNNZ))
    nnz=0
    do i=1,grid%nC
      select case(grid%sE(i))
      case(TET)
        do j=1,TET_N
          if(present(isDirichlet))then
            if(isDirichlet(grid%iNE(j,i)))then
              cycle
            end if
          end if
          do k=1,TET_N
            nnz=nnz+1
            iA(nnz)=grid%iNE(j,i)
            jA(nnz)=grid%iNE(k,i)
            vA(nnz)=0d0
            do l=1,size(TET_QW)
              vA(nnz)=vA(nnz)+dot_product(matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,j,l)),&
              &                           matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,k,l)))*&
              &               grid%detJ(l,i)*TET_QW(l)
            end do
          end do
        end do
      case(TET10)
        do j=1,TET10_N
          if(present(isDirichlet))then
            if(isDirichlet(grid%iNE(j,i)))then
              cycle
            end if
          end if
          do k=1,TET10_N
            nnz=nnz+1
            iA(nnz)=grid%iNE(j,i)
            jA(nnz)=grid%iNE(k,i)
            vA(nnz)=0d0
            do l=1,size(TET10_QW)
              vA(nnz)=vA(nnz)+dot_product(matmul(grid%invJ(:,:,l,i),TET10_GRAD_QP(:,j,l)),&
              &                           matmul(grid%invJ(:,:,l,i),TET10_GRAD_QP(:,k,l)))*&
              &               grid%detJ(l,i)*TET10_QW(l)
            end do
          end do
        end do
      case default
      end select
    end do
    if(present(isDirichlet))then
      do i=1,grid%nN
        if(isDirichlet(i))then
          nnz=nnz+1
          iA(nnz)=i
          jA(nnz)=i
          vA(nnz)=1d0
        end if
      end do
    end if
    call A%setCOO(iA,jA,vA,nnz,job=CSR_CLEAN_SORT)
    deallocate(iA,jA,vA)
  end subroutine
  
  !> find volumetric source discretization factor
  subroutine findVolSrc(grid,src)
    use modPolyFeGrid
    use modPolyGrid
    use modSMQ
    class(polyFeGrid),intent(inout)::grid !< the grid
    double precision,allocatable,intent(inout)::src(:) !< the discretized source factor
    
    call grid%up()
    if(.not.allocated(src))then
      allocate(src(grid%nN))
    end if
    src(:)=0d0
    do i=1,grid%nC
      select case(grid%sE(i))
      case(TET)
        do j=1,TET_N
          do l=1,size(TET_QW)
            src(grid%iNE(j,i))=src(grid%iNE(j,i))+TET_SHAPE_QP(j,l)*grid%detJ(l,i)*TET_QW(l)
          end do
        end do
      case(TET10)
        do j=1,TET10_N
          do l=1,size(TET10_QW)
            src(grid%iNE(j,i))=src(grid%iNE(j,i))+TET10_SHAPE_QP(j,l)*grid%detJ(l,i)*TET10_QW(l)
          end do
        end do
      case default
      end select
    end do
  end subroutine
  
  !> find surface source (i.e. Neumann BC) discretization factor
  subroutine findSurfSrc(grid,src)
    use modPolyFeGrid
    use modPolyGrid
    use modSMQ
    class(polyFeGrid),intent(inout)::grid !< the grid
    double precision,allocatable,intent(inout)::src(:) !< the discretized source factor
    
    call grid%up()
    if(.not.allocated(src))then
      allocate(src(grid%nN))
    end if
    src(:)=0d0
    do i=grid%nC+1,grid%nE
      select case(grid%sE(i))
      case(TRI)
        do j=1,TRI_N
          do l=1,size(TRI_QW)
            src(grid%iNE(j,i))=src(grid%iNE(j,i))+TRI_SHAPE_QP(j,l)*grid%detJ(l,i)*TRI_QW(l)
          end do
        end do
      case(TRI6)
        do j=1,TRI6_N
          do l=1,size(TRI6_QW)
            src(grid%iNE(j,i))=src(grid%iNE(j,i))+TRI6_SHAPE_QP(j,l)*grid%detJ(l,i)*TRI6_QW(l)
          end do
        end do
      case default
      end select
    end do
  end subroutine
  
end module
