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
    module procedure::recvIVectRealloc
    module procedure::recvDVect
    module procedure::recvDVectRealloc
    module procedure::recvIMat
    module procedure::recvIMatRealloc
    module procedure::recvDMat
    module procedure::recvDMatRealloc
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
  subroutine sendIArr(arr,n,dest)
    use mpi
    integer,intent(in)::arr(*) !< the integer array to be sent
    integer,intent(in)::n !< the size of arr
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(n,1,MPI_integer,dest,0,MPI_comm_world,ierr)
    call MPI_send(arr(1:n),n,MPI_integer,dest,0,MPI_comm_world,ierr)
  end subroutine
  
  !> receive arbitrary-rank integer array
  subroutine recvIArr(arr,n,source)
    use mpi
    integer,intent(out)::arr(*) !< the integer array to be sent
    integer,intent(out)::n !< the size of arr
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(n,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    call MPI_recv(arr(1:n),n,MPI_integer,source,0,MPI_comm_world,stat,ierr)
  end subroutine
  
  !> send arbitrary-rank double array
  subroutine sendDArr(arr,n,dest)
    use mpi
    double precision,intent(in)::arr(*) !< the double array to be sent
    integer,intent(in)::n !< the size of arr
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(n,1,MPI_integer,dest,0,MPI_comm_world,ierr)
    call MPI_send(arr(1:n),n,MPI_double_precision,dest,0,MPI_comm_world,ierr)
  end subroutine
  
  !> receive arbitrary-rank double array
  subroutine recvDArr(arr,n,source)
    use mpi
    double precision,intent(out)::arr(*) !< the double array to be sent
    integer,intent(out)::n !< the size of arr
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(n,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    call MPI_recv(arr(1:n),n,MPI_double_precision,source,0,MPI_comm_world,stat,ierr)
  end subroutine
  
  !> send integer vector object
  subroutine sendIVect(obj,dest)
    use mpi
    integer,intent(in)::obj(:) !< the integer vector to be sent
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(size(obj),1,MPI_integer,dest,0,MPI_comm_world,ierr)
    call sendIArr(obj,size(obj),dest)
  end subroutine
  
  !> receive integer vector object
  subroutine recvIVect(obj,source)
    use mpi
    integer,intent(out)::obj(:) !< the integer vector to be received
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(n,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    call recvIArr(obj,n,source)
  end subroutine
  
  !> receive integer vector object and reallocate if specified
  subroutine recvIVectRealloc(obj,source,realloc)
    use mpi
    use moduleBasicDataStruct
    integer,intent(inout),allocatable::obj(:) !< the integer vector to be received
    integer,intent(in)::source !< the source process
    logical,intent(in)::realloc !< whether the obj should be reallocated
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(n,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    if(realloc)then
      call reallocArr(obj,n)
    end if
    call recvIArr(obj,n,source)
  end subroutine
  
  !> send double vector object
  subroutine sendDVect(obj,dest)
    use mpi
    double precision,intent(in)::obj(:) !< the double vector to be sent
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send(size(obj),1,MPI_integer,dest,0,MPI_comm_world,ierr)
    call sendDArr(obj,size(obj),dest)
  end subroutine
  
  !> receive double vector object
  subroutine recvDVect(obj,source)
    use mpi
    double precision,intent(out)::obj(:) !< the double vector to be received
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(n,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    call recvDArr(obj,n,source)
  end subroutine
  
  !> receive double vector object and reallocate if specified
  subroutine recvDVectRealloc(obj,source,realloc)
    use mpi
    use moduleBasicDataStruct
    double precision,intent(inout),allocatable::obj(:) !< the double vector to be received
    integer,intent(in)::source !< the source process
    logical,intent(in)::realloc !< whether the obj should be reallocated
    integer ierr,stat(MPI_status_size)
    
    call MPI_recv(n,1,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    if(realloc)then
      call reallocArr(obj,n)
    end if
    call recvDArr(obj,n,source)
  end subroutine
  
  !> send integer matrix object
  subroutine sendIMat(obj,dest)
    use mpi
    integer,intent(in)::obj(:,:) !< the integer matrix to be sent
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send([size(obj,1),size(obj,2)],2,MPI_integer,dest,0,MPI_comm_world,ierr)
    call sendIArr(obj,size(obj),dest)
  end subroutine
  
  !> receive integer matrix object
  subroutine recvIMat(obj,source)
    use mpi
    integer,intent(out)::obj(:,:) !< the integer matrix to be received
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size),shp(2)
    
    call MPI_recv(shp,2,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    call recvIArr(obj,n,source)
  end subroutine
  
  !> receive integer matrix object and reallocate if specified
  subroutine recvIMatRealloc(obj,source,realloc)
    use mpi
    use moduleBasicDataStruct
    integer,intent(inout),allocatable::obj(:,:) !< the integer matrix to be received
    integer,intent(in)::source !< the source process
    logical,intent(in)::realloc !< whether the obj should be reallocated
    integer ierr,stat(MPI_status_size),shp(2)
    
    call MPI_recv(shp,2,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    if(realloc)then
      call reallocArr(obj,shp(1),shp(2))
    end if
    call recvIArr(obj,n,source)
  end subroutine
  
  !> send double matrix object
  subroutine sendDMat(obj,dest)
    use mpi
    double precision,intent(in)::obj(:,:) !< the double matrix to be sent
    integer,intent(in)::dest !< the destination process
    integer ierr
    
    call MPI_send([size(obj,1),size(obj,2)],2,MPI_integer,dest,0,MPI_comm_world,ierr)
    call sendDArr(obj,size(obj),dest)
  end subroutine
  
  !> receive double matrix object
  subroutine recvDMat(obj,source)
    use mpi
    double precision,intent(out)::obj(:,:) !< the double matrix to be received
    integer,intent(in)::source !< the source process
    integer ierr,stat(MPI_status_size),shp(2)
    
    call MPI_recv(shp,2,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    call recvDArr(obj,n,source)
  end subroutine
  
  !> receive double matrix object and reallocate if specified
  subroutine recvDMatRealloc(obj,source,realloc)
    use mpi
    use moduleBasicDataStruct
    double precision,intent(inout),allocatable::obj(:,:) !< the double matrix to be received
    integer,intent(in)::source !< the source process
    logical,intent(in)::realloc !< whether the obj should be reallocated
    integer ierr,stat(MPI_status_size),shp(2)
    
    call MPI_recv(shp,2,MPI_integer,source,0,MPI_comm_world,stat,ierr)
    if(realloc)then
      call reallocArr(obj,shp(1),shp(2))
    end if
    call recvDArr(obj,n,source)
  end subroutine
  
end module
