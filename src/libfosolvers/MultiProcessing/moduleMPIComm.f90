!----------------------------------------------------------------------------- best with 100 columns

!> MPI communication
module moduleMPIComm
  private
  
  ! variables
  integer,public,save::pidMPI !< MPI process id (from 0 to sizeMPI-1)
  integer,public,save::sizeMPI !< number of MPI processes
  
  ! individual procedures
  public initMPI
  public finalMPI
  
  !> send generic data
  interface sendDat
    module procedure::sendIScal
    module procedure::sendDScal
    module procedure::sendIVect
    module procedure::sendDVect
    module procedure::sendIMat
    module procedure::sendDMat
  end interface
  public::sendDat
  
  !> receive generic data
  interface recvDat
    module procedure::recvIScal
    module procedure::recvDScal
    module procedure::recvIVect
    module procedure::recvDVect
    module procedure::recvIMat
    module procedure::recvDMat
  end interface
  public::recvDat
  
contains
  
  !> initialize MPI environment
  subroutine initMPI()
    use mpi
    integer ierr
    
    call MPI_init(ierr)
    call MPI_comm_rank(MPI_comm_world,pidMPI,ierr)
    call MPI_comm_size(MPI_comm_world,sizeMPI,ierr)
  end subroutine
  
  !> finalize MPI environment
  subroutine finalMPI()
    use mpi
    integer ierr
    
    call MPI_finalize(ierr)
  end subroutine
  
  !> send integer scaler object
  subroutine sendIScal(obj,dest)
    use mpi
    integer,intent(in)::obj !< the integer scalar to be sent
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(obj,1,MPI_integer,dest,0,MPI_comm_world,ierr)
  end subroutine
  
  !> receive integer scaler object
  subroutine recvIScal(obj,source)
    use mpi
    integer,intent(out)::obj !< the integer scalar to be received
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(obj,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
  end subroutine
  
  !> send double scaler object
  subroutine sendDScal(obj,dest)
    use mpi
    double precision,intent(in)::obj !< the double scalar to be sent
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(obj,1,MPI_double_precision,dest,0,MPI_comm_world,ierr)
  end subroutine
  
  !> receive double scaler object
  subroutine recvDScal(obj,source)
    use mpi
    double precision,intent(out)::obj !< the double scalar to be received
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(obj,1,MPI_double_precision,source,0,MPI_comm_world,stat,ierr)
  end subroutine
  
  !> send arbitrary-rank integer array
  subroutine sendIAnyRankArr(arr,sizeArr,dest)
    use mpi
    integer,intent(in)::arr(*) !< the integer array to be sent
    integer,intent(in)::sizeArr(:) !< the size of all dimensions of arr
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(sizeArr(:),size(sizeArr),MPI_integer,dest,0,MPI_comm_world,ierr)
    n=product(sizeArr)
    call MPI_send(arr(1:n),n,MPI_integer,dest,0,MPI_comm_world,ierr)
  end subroutine
  
  !> receive arbitrary-rank integer array
  subroutine recvIAnyRankArr(arr,sizeArr,source)
    use mpi
    integer,intent(out)::arr(*) !< the integer array to be sent
    integer,intent(inout)::sizeArr(:) !< the size of all dimensions of arr
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(sizeArr(:),size(sizeArr),MPI_integer,source,0,MPI_comm_world,stat,ierr)
    n=product(sizeArr)
    call MPI_recv(arr(1:n),n,MPI_integer,source,0,MPI_comm_world,stat,ierr)
  end subroutine
  
  !> send arbitrary-rank double array
  subroutine sendDAnyRankArr(arr,sizeArr,dest)
    use mpi
    double precision,intent(in)::arr(*) !< the double array to be sent
    integer,intent(in)::sizeArr(:) !< the size of all dimensions of arr
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(sizeArr(:),size(sizeArr),MPI_integer,dest,0,MPI_comm_world,ierr)
    n=product(sizeArr)
    call MPI_send(arr(1:n),n,MPI_double_precision,dest,0,MPI_comm_world,ierr)
  end subroutine
  
  !> receive arbitrary-rank double array
  subroutine recvDAnyRankArr(arr,sizeArr,source)
    use mpi
    double precision,intent(out)::arr(*) !< the double array to be sent
    integer,intent(inout)::sizeArr(:) !< the size of all dimensions of arr
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(sizeArr(:),size(sizeArr),MPI_integer,source,0,MPI_comm_world,stat,ierr)
    n=product(sizeArr)
    call MPI_recv(arr(1:n),n,MPI_double_precision,source,0,MPI_comm_world,stat,ierr)
  end subroutine
  
  !> send integer vector object
  subroutine sendIVect(obj,dest)
    integer,intent(in)::obj(:) !< the integer vector to be sent
    integer,intent(in)::dest !< the destination process
    
    call sendIAnyRankArr(obj,[size(obj)],dest)
  end subroutine
  
  !> receive integer vector object
  subroutine recvIVect(obj,source)
    integer,intent(inout)::obj(:) !< the integer vector to be received
    integer,intent(in)::source !< the source process
    integer sz(1)
    
    call recvIAnyRankArr(obj,sz,source)
  end subroutine
  
  !> send double vector object
  subroutine sendDVect(obj,dest)
    double precision,intent(in)::obj(:) !< the double vector to be sent
    integer,intent(in)::dest !< the destination process
    
    call sendDAnyRankArr(obj,[size(obj)],dest)
  end subroutine
  
  !> receive double vector object
  subroutine recvDVect(obj,source)
    double precision,intent(inout)::obj(:) !< the double vector to be received
    integer,intent(in)::source !< the source process
    integer sz(1)
    
    call recvDAnyRankArr(obj,sz,source)
  end subroutine
  
  !> send integer matrix object
  subroutine sendIMat(obj,dest)
    integer,intent(in)::obj(:,:) !< the integer matrix to be sent
    integer,intent(in)::dest !< the destination process
    
    call sendIAnyRankArr(obj,[size(obj,1),size(obj,2)],dest)
  end subroutine
  
  !> receive integer matrix object
  subroutine recvIMat(obj,source)
    integer,intent(inout)::obj(:,:) !< the integer matrix to be received
    integer,intent(in)::source !< the source process
    integer sz(2)
    
    call recvIAnyRankArr(obj,sz,source)
  end subroutine
  
  !> send double matrix object
  subroutine sendDMat(obj,dest)
    double precision,intent(in)::obj(:,:) !< the double matrix to be sent
    integer,intent(in)::dest !< the destination process
    
    call sendDAnyRankArr(obj,[size(obj,1),size(obj,2)],dest)
  end subroutine
  
  !> receive double matrix object
  subroutine recvDMat(obj,source)
    double precision,intent(inout)::obj(:,:) !< the double matrix to be received
    integer,intent(in)::source !< the source process
    integer sz(2)
    
    call recvDAnyRankArr(obj,sz,source)
  end subroutine
  
end module
