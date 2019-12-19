!----------------------------------------------------------------------------- best with 100 columns

!> polygons and polyhedrons module
module modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> base type of polygons and polyhedrons
  type,public::polyX
    integer::nN !< number of nodes
    integer::nE !< number of elements (facets or cells)
    double precision,allocatable::pN(:,:) !< node position
    integer,allocatable::sE(:) !< shape of elements
    integer,allocatable::nNE(:) !< number of nodes in each element
    integer,allocatable::iNE(:,:) !< connectivity table
    logical::isUp !< if auxiliary data is updated
    double precision,allocatable::p(:,:) !< element position
  contains
    procedure,public::init=>initPolyX
    procedure,public::clear=>clearPolyX
    procedure,public::up=>upPolyX
    final::purgePolyX
  end type
  
contains
  
  !> initialize this polyX
  elemental subroutine initPolyX(this,nN,nE,m)
    class(polyX),intent(inout)::this !< this polyX
    integer,intent(in)::nN !< number of nodes
    integer,intent(in)::nE !< number of elements
    integer,intent(in)::m !< maximum number of nodes per elements
    
    call this%clear()
    this%nN=nN
    this%nE=nE
    allocate(this%pN(DIMS,nN))
    allocate(this%sE(nE))
    allocate(this%nNE(nE))
    allocate(this%iNE(m,nE))
    allocate(this%p(DIMS,nE))
    this%isUp=.false.
  end subroutine
  
  !> clear this polyX
  elemental subroutine clearPolyX(this)
    class(polyX),intent(inout)::this !< this polyX
    
    if(allocated(this%pN)) deallocate(this%pN)
    if(allocated(this%sE)) deallocate(this%sE)
    if(allocated(this%nNE)) deallocate(this%nNE)
    if(allocated(this%iNE)) deallocate(this%iNE)
    if(allocated(this%p)) deallocate(this%p)
  end subroutine
  
  !> update this polyX
  subroutine upPolyX(this)
    use modGeometry
    class(polyX),intent(inout)::this !< this polyX
    
    if(.not.this%isUp)then
      call sortPolyX(this)
      ! calculate element position
      forall(i=1:this%nE)
        this%p(:,i)=sum(this%pN(:,this%iNE(1:this%nNE(i),i)),2)/dble(this%nNE(i))
      end forall
      this%isUp=.true.
    end if
  end subroutine
  
  !> sort polyX nodes and elements
  subroutine sortPolyX(grid)
    use metis_interface
    use modSort
    use iso_c_binding
    class(polyX),intent(inout)::grid !< the polyX
    integer(kind=idx_t)::options(0:METIS_NOPTIONS-1)
    integer(kind=idx_t),allocatable::eptr(:),eind(:),xadj(:),adjncy(:),perm(:),iperm(:)
    type(C_PTR)::xadjptr,adjncyptr
    integer(kind=idx_t),pointer::fxadj(:),fadjncy(:)
    integer(C_INT)::ierr
    double precision,allocatable::pN(:,:)
    integer,allocatable::a(:),permE(:),sE(:),nNE(:),iNE(:,:)
    
    ! reorder nodes with METIS (improves performance of FEM-type applications)
    allocate(eptr(grid%nE+1))
    allocate(eind(sum(grid%nNE(:))))
    eptr(1)=1
    do i=1,grid%nE
      eptr(i+1)=eptr(i)+grid%nNE(i)
      eind(eptr(i):eptr(i+1)-1)=grid%iNE(1:grid%nNE(i),i)
    end do
    ierr=METIS_MeshToNodal(int(grid%nE,kind=idx_t),int(grid%nN,kind=idx_t),eptr,eind,&
    &                      int(1,kind=idx_t),xadjptr,adjncyptr)
    call c_f_pointer(xadjptr,fxadj,shape=[grid%nN+1])
    call c_f_pointer(adjncyptr,fadjncy,shape=[fxadj(grid%nN+1)])
    allocate(xadj(grid%nN+1))
    allocate(adjncy(fxadj(grid%nN+1)))
    allocate(perm(grid%nN))
    allocate(iperm(grid%nN))
    xadj(:)=fxadj(1:grid%nN+1)
    adjncy(:)=fadjncy(1:fxadj(grid%nN+1))
    ierr=METIS_SetDefaultOptions(options)
    options(METIS_OPTION_NUMBERING)=1
    ierr=METIS_NodeND(int(grid%nN,kind=idx_t),xadj,adjncy,options=options,perm=perm,iperm=iperm)
    allocate(pN(DIMS,grid%nN))
    pN=grid%pN(:,perm(:))
    grid%pN(:,:)=pN(:,:)
    forall(i=1:grid%nE)
      grid%iNE(1:grid%nNE(i),i)=iperm(grid%iNE(1:grid%nNE(i),i))
    end forall
    deallocate(eptr,eind,xadj,adjncy,perm,iperm,pN)
    ierr=METIS_Free(xadjptr)
    ierr=METIS_Free(adjncyptr)
    ! reorder elements for data locality
    allocate(a(grid%nE))
    allocate(permE(grid%nE))
    a(:)=minval(grid%iNE(:,:),dim=1)
    call sort(a,perm=permE)
    allocate(sE(grid%nE))
    allocate(nNE(grid%nE))
    allocate(iNE(size(grid%iNE,1),grid%nE))
    forall(i=1:grid%nE)
      sE(i)=grid%sE(permE(i))
      nNE(i)=grid%nNE(permE(i))
      iNE(:,i)=grid%iNE(:,permE(i))
    end forall
    grid%sE(:)=sE(:)
    grid%nNE(:)=nNE(:)
    grid%iNE(:,:)=iNE(:,:)
    deallocate(a,permE,sE,nNE,iNE)
  end subroutine
  
  !> destructor of polyX
  elemental subroutine purgePolyX(this)
    type(polyX),intent(inout)::this !< this polyX
    
    call this%clear()
  end subroutine
  
end module
