!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! multi-processing related data (level 1)
!*****************************************
module moduleMP1
  private
  
  ! constants
  integer,public,parameter::ROOT_PID=0
  
  ! variables
  integer,public,save::errMPI,pidMPI,sizeMPI
  integer,allocatable,public,save::statMPI(:)
  
  !-------------------
  ! send generic data
  !-------------------
  interface sendData
    module procedure::sendIntegerScal
    module procedure::sendIntegerVect
    module procedure::sendIntegerMat
    module procedure::sendDoubleScal
    module procedure::sendDoubleVect
    module procedure::sendDoubleMat
    module procedure::sendDataItemScal
    module procedure::sendDataItemVect
    module procedure::sendDataSetScal
    module procedure::sendDataSetVect
    module procedure::sendNodeScal
    module procedure::sendNodeVect
    module procedure::sendPointScal
    module procedure::sendPointVect
    module procedure::sendLineScal
    module procedure::sendLineVect
    module procedure::sendFacetScal
    module procedure::sendFacetVect
    module procedure::sendBlockScal
    module procedure::sendBlockVect
  end interface
  public::sendData
  
  !----------------------
  ! receive generic data
  !----------------------
  interface recvData
    module procedure::recvIntegerScal
    module procedure::recvIntegerVect
    module procedure::recvIntegerVectRealloc
    module procedure::recvIntegerMat
    module procedure::recvIntegerMatRealloc
    module procedure::recvDoubleScal
    module procedure::recvDoubleVect
    module procedure::recvDoubleVectRealloc
    module procedure::recvDoubleMat
    module procedure::recvDoubleMatRealloc
    module procedure::recvDataItemScal
    module procedure::recvDataItemVect
    module procedure::recvDataItemVectRealloc
    module procedure::recvDataSetScal
    module procedure::recvDataSetVect
    module procedure::recvDataSetVectRealloc
    module procedure::recvNodeScal
    module procedure::recvNodeVect
    module procedure::recvNodeVectRealloc
    module procedure::recvPointScal
    module procedure::recvPointVect
    module procedure::recvPointVectRealloc
    module procedure::recvLineScal
    module procedure::recvLineVect
    module procedure::recvLineVectRealloc
    module procedure::recvFacetScal
    module procedure::recvFacetVect
    module procedure::recvFacetVectRealloc
    module procedure::recvBlockScal
    module procedure::recvBlockVect
    module procedure::recvBlockVectRealloc
  end interface
  public::recvData
  
  ! individual procedures
  public initMPI
  public finalMPI
  
contains
  
  !----------------------------
  ! initialize MPI environment
  !----------------------------
  subroutine initMPI()
    use mpi
    
    allocate(statMPI(MPI_status_size))
    call MPI_init(errMPI)
    call MPI_comm_rank(MPI_comm_world,pidMPI,errMPI)
    call MPI_comm_size(MPI_comm_world,sizeMPI,errMPI)
  end subroutine
  
  !--------------------------
  ! finalize MPI environment
  !--------------------------
  subroutine finalMPI()
    use mpi
    
    call MPI_finalize(errMPI)
  end subroutine
  
  !------------------------------------
  ! send scaler object of type integer
  !------------------------------------
  subroutine sendIntegerScal(obj,dest,tag)
    use mpi
    integer,intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
  end subroutine
  
  !---------------------------------------
  ! receive scaler object of type integer
  !---------------------------------------
  subroutine recvIntegerScal(obj,source,tag)
    use mpi
    integer,intent(in)::obj
    integer,intent(in)::source,tag
    
    call MPI_recv(obj,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
  end subroutine
  
  !------------------------------------
  ! send vector object of type integer
  !------------------------------------
  subroutine sendIntegerVect(obj,dest,tag)
    use mpi
    integer,intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(n>0)then
      call MPI_send(obj,n,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    end if
  end subroutine
  
  !---------------------------------------
  ! receive vector object of type integer
  !---------------------------------------
  subroutine recvIntegerVect(obj,source,tag)
    use mpi
    integer,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    if(n>0)then
      call MPI_recv(obj,n,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    end if
  end subroutine
  
  !------------------------------------------------------
  ! receive vector object of type integer and reallocate
  !------------------------------------------------------
  subroutine recvIntegerVectRealloc(obj,source,tag,realloc)
    use mpi
    integer,allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      if(n>0)then
        call MPI_recv(obj,n,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      end if
    else
      call recvIntegerVect(obj,source,tag)
    end if
  end subroutine
  
  !------------------------------------
  ! send matrix object of type integer
  !------------------------------------
  subroutine sendIntegerMat(obj,dest,tag)
    use mpi
    integer,intent(in)::obj(:,:)
    integer,intent(in)::dest,tag
    
    m=size(obj,1)
    n=size(obj,2)
    call MPI_send(m,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(m>0.and.n>0)then
      do i=1,m
        call MPI_send(obj(i,:),n,MPI_integer,dest,tag,MPI_comm_world,errMPI)
      end do
    end if
  end subroutine
  
  !---------------------------------------
  ! receive matrix object of type integer
  !---------------------------------------
  subroutine recvIntegerMat(obj,source,tag)
    use mpi
    integer,intent(inout)::obj(:,:)
    integer,intent(in)::source,tag
    integer,allocatable::buffInteger(:)
    
    call MPI_recv(m,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    if(m>0.and.n>0)then
      allocate(buffInteger(n))
      do i=1,m
        call MPI_recv(buffInteger,n,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
        obj(i,:)=buffInteger(:)
      end do
      deallocate(buffInteger)
    end if
  end subroutine
  
  !------------------------------------------------------
  ! receive matrix object of type integer and reallocate
  !------------------------------------------------------
  subroutine recvIntegerMatRealloc(obj,source,tag,realloc)
    use mpi
    integer,allocatable,intent(inout)::obj(:,:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    integer,allocatable::buffInteger(:)
    
    if(realloc)then
      call MPI_recv(m,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(m,n))
      if(m>0.and.n>0)then
        allocate(buffInteger(n))
        do i=1,m
          call MPI_recv(buffInteger,n,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
          obj(i,:)=buffInteger(:)
        end do
        deallocate(buffInteger)
      end if
    else
      call recvIntegerMat(obj,source,tag)
    end if
  end subroutine
  
  !---------------------------------------------
  ! send scaler object of type double precision
  !---------------------------------------------
  subroutine sendDoubleScal(obj,dest,tag)
    use mpi
    double precision,intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj,1,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
  end subroutine
  
  !------------------------------------------------
  ! receive scaler object of type double precision
  !------------------------------------------------
  subroutine recvDoubleScal(obj,source,tag)
    use mpi
    double precision,intent(in)::obj
    integer,intent(in)::source,tag
    
    call MPI_recv(obj,1,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
  end subroutine
  
  !---------------------------------------------
  ! send vector object of type double precision
  !---------------------------------------------
  subroutine sendDoubleVect(obj,dest,tag)
    use mpi
    double precision,intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(n>0)then
      call MPI_send(obj,n,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
    end if
  end subroutine
  
  !------------------------------------------------
  ! receive vector object of type double precision
  !------------------------------------------------
  subroutine recvDoubleVect(obj,source,tag)
    use mpi
    double precision,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    if(n>0)then
      call MPI_recv(obj,n,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
    end if
  end subroutine
  
  !---------------------------------------------------------------
  ! receive vector object of type double precision and reallocate
  !---------------------------------------------------------------
  subroutine recvDoubleVectRealloc(obj,source,tag,realloc)
    use mpi
    double precision,allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      if(n>0)then
        call MPI_recv(obj,n,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
      end if
    else
      call recvDoubleVect(obj,source,tag)
    end if
  end subroutine
  
  !---------------------------------------------
  ! send matrix object of type double precision
  !---------------------------------------------
  subroutine sendDoubleMat(obj,dest,tag)
    use mpi
    double precision,intent(in)::obj(:,:)
    integer,intent(in)::dest,tag
    
    m=size(obj,1)
    n=size(obj,2)
    call MPI_send(m,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(m>0.and.n>0)then
      do i=1,m
        call MPI_send(obj(i,:),n,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
      end do
    end if
  end subroutine
  
  !------------------------------------------------
  ! receive matrix object of type double precision
  !------------------------------------------------
  subroutine recvDoubleMat(obj,source,tag)
    use mpi
    double precision,intent(inout)::obj(:,:)
    integer,intent(in)::source,tag
    double precision,allocatable::buffDouble(:)
    
    call MPI_recv(m,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    if(m>0.and.n>0)then
      allocate(buffDouble(n))
      do i=1,m
        call MPI_recv(buffDouble,n,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
        obj(i,:)=buffDouble(:)
      end do
      deallocate(buffDouble)
    end if
  end subroutine
  
  !---------------------------------------------------------------
  ! receive matrix object of type double precision and reallocate
  !---------------------------------------------------------------
  subroutine recvDoubleMatRealloc(obj,source,tag,realloc)
    use mpi
    double precision,allocatable,intent(inout)::obj(:,:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    double precision,allocatable::buffDouble(:)
    
    if(realloc)then
      call MPI_recv(m,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(m,n))
      if(m>0.and.n>0)then
        allocate(buffDouble(n))
        do i=1,m
          call MPI_recv(buffDouble,n,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
          obj(i,:)=buffDouble(:)
        end do
        deallocate(buffDouble)
      end if
    else
      call recvDoubleMat(obj,source,tag)
    end if
  end subroutine
  
  !-----------------------------------------
  ! send scaler object of type typeDataItem
  !-----------------------------------------
  subroutine sendDataItemScal(obj,dest,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataItem),intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj%DataName,DATA_NAME_LENGTH,MPI_character,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%DataType,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    select case(obj%DataType)
      case(VAL_TYPE)
        call MPI_send(obj%Val,1,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
      case(TAB1D_TYPE)
        if(allocated(obj%Tab1d))then
          call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
          call sendDoubleMat(obj%Tab1d,dest,tag)
        else
          call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
        end if
      case default
    end select
  end subroutine
  
  !--------------------------------------------
  ! receive scaler object of type typeDataItem
  !--------------------------------------------
  subroutine recvDataItemScal(obj,source,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataItem),intent(inout)::obj
    integer,intent(in)::source,tag
    character(DATA_NAME_LENGTH) buffName
    integer buffInteger
    double precision buffDouble
    logical isAllocated
    
    call MPI_recv(buffName,DATA_NAME_LENGTH,MPI_character,&
    &             source,tag,MPI_comm_world,statMPI,errMPI)
    obj%DataName=buffName
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%DataType=buffInteger
    select case(obj%DataType)
      case(VAL_TYPE)
        call MPI_recv(buffDouble,1,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
        obj%Val=buffDouble
      case(TAB1D_TYPE)
        call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
        if(isAllocated)then
          call recvDoubleMatRealloc(obj%Tab1d,source,tag,realloc=.true.)
        end if
      case default
    end select
  end subroutine
  
  !-----------------------------------------
  ! send vector object of type typeDataItem
  !-----------------------------------------
  subroutine sendDataItemVect(obj,dest,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataItem),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendDataItemScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !--------------------------------------------
  ! receive vector object of type typeDataItem
  !--------------------------------------------
  subroutine recvDataItemVect(obj,source,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataItem),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvDataItemScal(obj(i),source,tag)
    end do
  end subroutine
  
  !-----------------------------------------------------------
  ! receive vector object of type typeDataItem and reallocate
  !-----------------------------------------------------------
  subroutine recvDataItemVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleMiscDataStruct
    type(typeDataItem),pointer,intent(inout)::obj(:) !TODO: what about allocatable arrays
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(associated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvDataItemScal(obj(i),source,tag)
      end do
    else
      call recvDataItemVect(obj,source,tag)
    end if
  end subroutine
  
  !----------------------------------------
  ! send scaler object of type typeDataSet
  !----------------------------------------
  subroutine sendDataSetScal(obj,dest,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataSet),intent(in)::obj
    integer,intent(in)::dest,tag
    
    if(associated(obj%DataItem))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendDataItemVect(obj%DataItem,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
  end subroutine
  
  !-------------------------------------------
  ! receive scaler object of type typeDataSet
  !-------------------------------------------
  subroutine recvDataSetScal(obj,source,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataSet),intent(inout)::obj
    integer,intent(in)::source,tag
    logical isAllocated
    
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvDataItemVectRealloc(obj%DataItem,source,tag,realloc=.true.)
      obj%ptrLast=>obj%DataItem(size(obj%DataItem))
    end if
  end subroutine
  
  !----------------------------------------
  ! send vector object of type typeDataSet
  !----------------------------------------
  subroutine sendDataSetVect(obj,dest,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataSet),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendDataSetScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !-------------------------------------------
  ! receive vector object of type typeDataSet
  !-------------------------------------------
  subroutine recvDataSetVect(obj,source,tag)
    use mpi
    use moduleMiscDataStruct
    type(typeDataSet),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvDataSetScal(obj(i),source,tag)
    end do
  end subroutine
  
  !----------------------------------------------------------
  ! receive vector object of type typeDataSet and reallocate
  !----------------------------------------------------------
  subroutine recvDataSetVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleMiscDataStruct
    type(typeDataSet),allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvDataSetScal(obj(i),source,tag)
      end do
    else
      call recvDataSetVect(obj,source,tag)
    end if
  end subroutine
  
  !-------------------------------------
  ! send scaler object of type typeNode
  !-------------------------------------
  subroutine sendNodeScal(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeNode),intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj%Ind,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call sendDoubleVect(obj%Pos,dest,tag)
    if(allocated(obj%FacetInd))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerVect(obj%FacetInd,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    if(allocated(obj%BlockInd))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerVect(obj%BlockInd,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
  end subroutine
  
  !----------------------------------------
  ! receive scaler object of type typeNode
  !----------------------------------------
  subroutine recvNodeScal(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeNode),intent(inout)::obj
    integer,intent(in)::source,tag
    integer buffInteger
    logical isAllocated
    
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Ind=buffInteger
    call recvDoubleVect(obj%Pos,source,tag)
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerVectRealloc(obj%FacetInd,source,tag,realloc=.true.)
    end if
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerVectRealloc(obj%BlockInd,source,tag,realloc=.true.)
    end if
  end subroutine
  
  !-------------------------------------
  ! send vector object of type typeNode
  !-------------------------------------
  subroutine sendNodeVect(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeNode),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendNodeScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !----------------------------------------
  ! receive vector object of type typeNode
  !----------------------------------------
  subroutine recvNodeVect(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeNode),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvNodeScal(obj(i),source,tag)
    end do
  end subroutine
  
  !-------------------------------------------------------
  ! receive vector object of type typeNode and reallocate
  !-------------------------------------------------------
  subroutine recvNodeVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleGrid
    type(typeNode),allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvNodeScal(obj(i),source,tag)
      end do
    else
      call recvNodeVect(obj,source,tag)
    end if
  end subroutine
  
  !--------------------------------------
  ! send scaler object of type typePoint
  !--------------------------------------
  subroutine sendPointScal(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typePoint),intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj%Ind,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%NodeInd,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%GeoEnti,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call sendDoubleVect(obj%Pos,dest,tag)
  end subroutine
  
  !-----------------------------------------
  ! receive scaler object of type typePoint
  !-----------------------------------------
  subroutine recvPointScal(obj,source,tag)
    use mpi
    use moduleGrid
    type(typePoint),intent(inout)::obj
    integer,intent(in)::source,tag
    integer buffInteger
    
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Ind=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%NodeInd=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%GeoEnti=buffInteger
    call recvDoubleVect(obj%Pos,source,tag)
  end subroutine
  
  !--------------------------------------
  ! send vector object of type typePoint
  !--------------------------------------
  subroutine sendPointVect(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typePoint),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendPointScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !-----------------------------------------
  ! receive vector object of type typePoint
  !-----------------------------------------
  subroutine recvPointVect(obj,source,tag)
    use mpi
    use moduleGrid
    type(typePoint),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvPointScal(obj(i),source,tag)
    end do
  end subroutine
  
  !--------------------------------------------------------
  ! receive vector object of type typePoint and reallocate
  !--------------------------------------------------------
  subroutine recvPointVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleGrid
    type(typePoint),allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvPointScal(obj(i),source,tag)
      end do
    else
      call recvPointVect(obj,source,tag)
    end if
  end subroutine
  
  !-------------------------------------
  ! send scaler object of type typeLine
  !-------------------------------------
  subroutine sendLineScal(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeLine),intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj%Ind,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%ShapeType,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%NodeNum,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call sendIntegerVect(obj%NodeInd,dest,tag)
    call MPI_send(obj%GeoEnti,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call sendDoubleVect(obj%PC,dest,tag)
    call MPI_send(obj%Length,1,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
  end subroutine
  
  !-----------------------------------------
  ! receive scaler object of type typeLine
  !-----------------------------------------
  subroutine recvLineScal(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeLine),intent(inout)::obj
    integer,intent(in)::source,tag
    integer buffInteger
    double precision buffDouble
    
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Ind=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%ShapeType=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%NodeNum=buffInteger
    call recvIntegerVectRealloc(obj%NodeInd,source,tag,realloc=.true.)
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%GeoEnti=buffInteger
    call recvDoubleVect(obj%PC,source,tag)
    call MPI_recv(buffDouble,1,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Length=buffDouble
  end subroutine
  
  !--------------------------------------
  ! send vector object of type typeLine
  !--------------------------------------
  subroutine sendLineVect(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeLine),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendLineScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !-----------------------------------------
  ! receive vector object of type typeLine
  !-----------------------------------------
  subroutine recvLineVect(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeLine),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvLineScal(obj(i),source,tag)
    end do
  end subroutine
  
  !-------------------------------------------------------
  ! receive vector object of type typeLine and reallocate
  !-------------------------------------------------------
  subroutine recvLineVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleGrid
    type(typeLine),allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvLineScal(obj(i),source,tag)
      end do
    else
      call recvLineVect(obj,source,tag)
    end if
  end subroutine
  
  !--------------------------------------
  ! send scaler object of type typeFacet
  !--------------------------------------
  subroutine sendFacetScal(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeFacet),intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj%Ind,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%ShapeType,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%NodeNum,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(allocated(obj%NodeInd))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerVect(obj%NodeInd,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    call MPI_send(obj%GeoEnti,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call sendIntegerVect(obj%NeibBlock,dest,tag)
    call sendDoubleVect(obj%PC,dest,tag)
    call MPI_send(obj%Area,1,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
    call sendDoubleVect(obj%Norm,dest,tag)
  end subroutine
  
  !-----------------------------------------
  ! receive scaler object of type typeFacet
  !-----------------------------------------
  subroutine recvFacetScal(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeFacet),intent(inout)::obj
    integer,intent(in)::source,tag
    integer buffInteger
    double precision buffDouble
    logical isAllocated
    
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Ind=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%ShapeType=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%NodeNum=buffInteger
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerVectRealloc(obj%NodeInd,source,tag,realloc=.true.)
    end if
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%GeoEnti=buffInteger
    call recvIntegerVect(obj%NeibBlock,source,tag)
    call recvDoubleVect(obj%PC,source,tag)
    call MPI_recv(buffDouble,1,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Area=buffDouble
    call recvDoubleVect(obj%Norm,source,tag)
  end subroutine
  
  !--------------------------------------
  ! send vector object of type typeFacet
  !--------------------------------------
  subroutine sendFacetVect(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeFacet),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendFacetScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !-----------------------------------------
  ! receive vector object of type typeFacet
  !-----------------------------------------
  subroutine recvFacetVect(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeFacet),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvFacetScal(obj(i),source,tag)
    end do
  end subroutine
  
  !--------------------------------------------------------
  ! receive vector object of type typeFacet and reallocate
  !--------------------------------------------------------
  subroutine recvFacetVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleGrid
    type(typeFacet),allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvFacetScal(obj(i),source,tag)
      end do
    else
      call recvFacetVect(obj,source,tag)
    end if
  end subroutine
  
  !--------------------------------------
  ! send scaler object of type typeBlock
  !--------------------------------------
  subroutine sendBlockScal(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeBlock),intent(in)::obj
    integer,intent(in)::dest,tag
    
    call MPI_send(obj%Ind,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%ShapeType,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%NodeNum,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    call MPI_send(obj%SurfNum,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(allocated(obj%NodeInd))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerVect(obj%NodeInd,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    call MPI_send(obj%GeoEnti,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    if(allocated(obj%Prt))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerVect(obj%Prt,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    if(allocated(obj%Neib))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerVect(obj%Neib,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    call sendDoubleVect(obj%PC,dest,tag)
    call MPI_send(obj%Vol,1,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
    if(allocated(obj%SurfNodeInd))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendIntegerMat(obj%SurfNodeInd,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    if(allocated(obj%SurfPC))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendDoubleMat(obj%SurfPC,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    if(allocated(obj%SurfArea))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendDoubleVect(obj%SurfArea,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
    if(allocated(obj%SurfNorm))then
      call MPI_send(.true.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
      call sendDoubleMat(obj%SurfNorm,dest,tag)
    else
      call MPI_send(.false.,1,MPI_logical,dest,tag,MPI_comm_world,errMPI)
    end if
  end subroutine
  
  !-----------------------------------------
  ! receive scaler object of type typeBlock
  !-----------------------------------------
  subroutine recvBlockScal(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeBlock),intent(inout)::obj
    integer,intent(in)::source,tag
    integer buffInteger
    double precision buffDouble
    logical isAllocated
    
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Ind=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%ShapeType=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%NodeNum=buffInteger
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%SurfNum=buffInteger
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerVectRealloc(obj%NodeInd,source,tag,realloc=.true.)
    end if
    call MPI_recv(buffInteger,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%GeoEnti=buffInteger
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerVectRealloc(obj%Prt,source,tag,realloc=.true.)
    end if
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerVectRealloc(obj%Neib,source,tag,realloc=.true.)
    end if
    call recvDoubleVect(obj%PC,source,tag)
    call MPI_recv(buffDouble,1,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
    obj%Vol=buffDouble
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvIntegerMatRealloc(obj%SurfNodeInd,source,tag,realloc=.true.)
    end if
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvDoubleMatRealloc(obj%SurfPC,source,tag,realloc=.true.)
    end if
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvDoubleVectRealloc(obj%SurfArea,source,tag,realloc=.true.)
    end if
    call MPI_recv(isAllocated,1,MPI_logical,source,tag,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvDoubleMatRealloc(obj%SurfNorm,source,tag,realloc=.true.)
    end if
  end subroutine
  
  !--------------------------------------
  ! send vector object of type typeBlock
  !--------------------------------------
  subroutine sendBlockVect(obj,dest,tag)
    use mpi
    use moduleGrid
    type(typeBlock),intent(in)::obj(:)
    integer,intent(in)::dest,tag
    
    n=size(obj)
    call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
    do i=1,n
      call sendBlockScal(obj(i),dest,tag)
    end do
  end subroutine
  
  !-----------------------------------------
  ! receive vector object of type typeBlock
  !-----------------------------------------
  subroutine recvBlockVect(obj,source,tag)
    use mpi
    use moduleGrid
    type(typeBlock),intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    do i=1,n
      call recvBlockScal(obj(i),source,tag)
    end do
  end subroutine
  
  !--------------------------------------------------------
  ! receive vector object of type typeBlock and reallocate
  !--------------------------------------------------------
  subroutine recvBlockVectRealloc(obj,source,tag,realloc)
    use mpi
    use moduleGrid
    type(typeBlock),allocatable,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    logical,intent(in)::realloc
    
    if(realloc)then
      call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
      if(allocated(obj))then
        deallocate(obj)
      end if
      allocate(obj(n))
      do i=1,n
        call recvBlockScal(obj(i),source,tag)
      end do
    else
      call recvBlockVect(obj,source,tag)
    end if
  end subroutine
  
end module
