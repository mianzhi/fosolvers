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
  end interface
  public::sendData
  
  !---------------------
  ! recive generic data
  !---------------------
  interface recvData
    module procedure::recvDoubleScal
    module procedure::recvDoubleVect
    module procedure::recvDoubleVectRealloc
    module procedure::recvIntegerScal
    module procedure::recvIntegerVect
    module procedure::recvIntegerVectRealloc
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
  
end module
