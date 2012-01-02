!----------------------------------------------------------------------------- best with 100 columns

!******************
! MPI related data
!******************
module moduleMP
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
    module procedure::sendDataItemScal
    module procedure::sendDataItemVect
  end interface
  public::sendData
  
  !---------------------
  ! recive generic data
  !---------------------
  interface recvData
    module procedure::recvIntegerScal
    module procedure::recvIntegerVect
    module procedure::recvIntegerVectRealloc
    module procedure::recvDoubleScal
    module procedure::recvDoubleVect
    module procedure::recvDoubleVectRealloc
    module procedure::recvDataItemScal
    module procedure::recvDataItemVect
    module procedure::recvDataItemVectRealloc
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
    call MPI_send(obj,n,MPI_integer,dest,tag,MPI_comm_world,errMPI)
  end subroutine
  
  !---------------------------------------
  ! receive vector object of type integer
  !---------------------------------------
  subroutine recvIntegerVect(obj,source,tag)
    use mpi
    integer,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    call MPI_recv(obj,n,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
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
      call MPI_recv(obj,n,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
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
    call MPI_send(obj,n,MPI_double_precision,dest,tag,MPI_comm_world,errMPI)
  end subroutine
  
  !------------------------------------------------
  ! receive vector object of type double precision
  !------------------------------------------------
  subroutine recvDoubleVect(obj,source,tag)
    use mpi
    double precision,intent(inout)::obj(:)
    integer,intent(in)::source,tag
    
    call MPI_recv(n,1,MPI_integer,source,tag,MPI_comm_world,statMPI,errMPI)
    call MPI_recv(obj,n,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
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
      call MPI_recv(obj,n,MPI_double_precision,source,tag,MPI_comm_world,statMPI,errMPI)
    else
      call recvDoubleVect(obj,source,tag)
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
  
end module
