!----------------------------------------------------------------------------- best with 100 columns

!> edge-based finite element discretizations
module modEdgeFEM
  public
  
  integer,parameter,private::DIMS=3
  
contains
  
  !> find cell-average curl of vp
  subroutine findCellCurl(grid,vp,rst)
    use modPolyEdgeFeGrid
    use modPolyGrid
    use modSMQ
    class(polyEdgeFeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::vp(grid%nEdge) !< vector potential
    double precision,allocatable,intent(inout)::rst(:,:) !< curl of vp at cells
    double precision::gradNa(DIMS),gradNb(DIMS),curlNj(DIMS)
    
    call grid%up()
    if(.not.allocated(rst))then
      allocate(rst(DIMS,grid%nC))
    end if
    rst(:,:)=0d0
    do i=1,grid%nC
      select case(grid%sE(i))
      case(TET)
        do j=1,TET_N_EDGE
          l=1 ! only 1 quadrature point for TET
          gradNa=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(1,j),l))
          gradNb=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(2,j),l))
          curlNj(1)=2*(gradNa(2)*gradNb(3)-gradNa(3)*gradNb(2))
          curlNj(2)=2*(gradNa(3)*gradNb(1)-gradNa(1)*gradNb(3))
          curlNj(3)=2*(gradNa(1)*gradNb(2)-gradNa(2)*gradNb(1))
          curlNj=curlNj*merge(1,-1,grid%sameDir(j,i))
          rst(:,i)=rst(:,i)+curlNj*vp(grid%iEdgeE(j,i))
        end do
      case default
      end select
    end do
  end subroutine
  
  !> find curl of curl operator
  !> NOTE: A must be pre-initialized to contain the temporary CSR matrix with duplications
  subroutine findCurlCurl(grid,A,isInactive)
    use modPolyEdgeFeGrid
    use modPolyGrid
    use modSMQ
    use modSparse
    class(polyEdgeFeGrid),intent(inout)::grid !< the grid
    class(linEq),intent(inout)::A !< the discretized linear system
    logical,intent(in),optional::isInactive(grid%nEdge) !< whether inactive (gauged and Dirichlet)
    integer::nnz,maxNNZ
    integer,allocatable::iA(:),jA(:)
    double precision,allocatable::vA(:)
    double precision::curlNj(DIMS),curlNk(DIMS),gradNa(DIMS),gradNb(DIMS)
    
    call grid%up()
    maxNNZ=size(grid%iEdgeE,1)**2*grid%nC
    allocate(iA(maxNNZ)) ! COO temporary matrix
    allocate(jA(maxNNZ))
    allocate(vA(maxNNZ))
    nnz=0
    do i=1,grid%nC
      select case(grid%sE(i))
      case(TET)
        do j=1,TET_N_EDGE
          if(present(isInactive))then
            if(isInactive(grid%iEdgeE(j,i)))then
              cycle
            end if
          end if
          do k=1,TET_N_EDGE
            nnz=nnz+1
            iA(nnz)=grid%iEdgeE(j,i)
            jA(nnz)=grid%iEdgeE(k,i)
            vA(nnz)=0d0
            do l=1,size(TET_QW)
              gradNa=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(1,j),l))
              gradNb=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(2,j),l))
              curlNj(1)=2*(gradNa(2)*gradNb(3)-gradNa(3)*gradNb(2))
              curlNj(2)=2*(gradNa(3)*gradNb(1)-gradNa(1)*gradNb(3))
              curlNj(3)=2*(gradNa(1)*gradNb(2)-gradNa(2)*gradNb(1))
              curlNj=curlNj*merge(1,-1,grid%sameDir(j,i))
              gradNa=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(1,k),l))
              gradNb=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(2,k),l))
              curlNk(1)=2*(gradNa(2)*gradNb(3)-gradNa(3)*gradNb(2))
              curlNk(2)=2*(gradNa(3)*gradNb(1)-gradNa(1)*gradNb(3))
              curlNk(3)=2*(gradNa(1)*gradNb(2)-gradNa(2)*gradNb(1))
              curlNk=curlNk*merge(1,-1,grid%sameDir(k,i))
              vA(nnz)=vA(nnz)+dot_product(curlNj,curlNk)*grid%detJ(l,i)*TET_QW(l)
            end do
          end do
        end do
      case default
      end select
    end do
    if(present(isInactive))then
      do i=1,grid%nEdge
        if(isInactive(i))then
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
  
  !> find volumetric source discretization of cell vector v
  subroutine findVolVectSrc(grid,v,rst)
    use modPolyEdgeFeGrid
    use modPolyGrid
    use modSMQ
    class(polyEdgeFeGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::v(DIMS,grid%nC) !< vector in cell (source density)
    double precision,allocatable,intent(inout)::rst(:) !< source discretization on edge
    double precision::Na,Nb,gradNa(DIMS),gradNb(DIMS),Ne(DIMS)
    
    call grid%up()
    if(.not.allocated(rst))then
      allocate(rst(grid%nEdge))
    end if
    rst(:)=0d0
    do i=1,grid%nC
      select case(grid%sE(i))
      case(TET)
        do j=1,TET_N_EDGE
          do l=1,size(TET_QW)
            gradNa=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(1,j),l))
            gradNb=matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,TET_NI_EDGE(2,j),l))
            Na=TET_SHAPE_QP(TET_NI_EDGE(1,j),l)
            Nb=TET_SHAPE_QP(TET_NI_EDGE(2,j),l)
            Ne=(Na*gradNb-Nb*gradNa)*merge(1,-1,grid%sameDir(j,i))
            rst(grid%iEdgeE(j,i))=rst(grid%iEdgeE(j,i))&
            &                     +dot_product(Ne,v(:,i))*grid%detJ(l,i)*TET_QW(l)
          end do
        end do
      case default
      end select
    end do
  end subroutine
  
end module
