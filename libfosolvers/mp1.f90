!----------------------------------------------------------------------------- best with 100 columns

!****************************
! MPI related data (level 1)
!****************************
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
  end interface
  public::sendData
  
  !----------------------
  ! receive generic data
  !----------------------
  interface recvData
    module procedure::recvIntegerScal
    module procedure::recvIntegerVect
    module procedure::recvIntegerVectRealloc
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
        n=size(obj%Tab1d,1)
        call MPI_send(n,1,MPI_integer,dest,tag,MPI_comm_world,errMPI)
        do i=1,n
          do j=1,2
            call MPI_send(obj%Tab1d(i,j),1,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
          end do
        end do
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
        call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
        if(allocated(obj%Tab1d))then
          deallocate(obj%Tab1d)
        end if
        allocate(obj%Tab1d(n,2))
        do i=1,n
          do j=1,2
            call MPI_recv(buffDouble,1,MPI_double_precision,&
            &             source,tag,MPI_comm_world,statMPI,errMPI)
            obj%Tab1d(i,j)=buffDouble
          end do
        end do
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
  
end module
