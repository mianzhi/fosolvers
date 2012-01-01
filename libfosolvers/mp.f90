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
  
  ! procedures
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
  
end module
