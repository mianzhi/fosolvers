!----------------------------------------------------------------------------- best with 100 columns

function basicFEM1() result(ierr)
  use modPolyFeGrid
  use modFileIO
  use modSparse
  use modBasicFEM
  type(polyFeGrid)::grid
  type(multiFront)::A
  double precision,allocatable::phi(:),rhs(:)
  logical,allocatable::isDirichlet(:)
  integer::ierr
  
  ierr=0
  open(10,file='data/tet10.msh',action='read')
    call readGMSH(10,grid)
  close(10)
  call grid%up()
  
  allocate(phi(grid%nN))
  allocate(rhs(grid%nN))
  allocate(isDirichlet(grid%nN))
  
  rhs(:)=0d0
  isDirichlet(:)=.false.
  do i=grid%nC+1,grid%nE
    if(grid%gid(i)==1)then
      isDirichlet(grid%iNE(1:grid%nNE(i),i))=.true.
      rhs((grid%iNE(1:grid%nNE(i),i)))=0d0
    end if
    if(grid%gid(i)==2)then
      isDirichlet(grid%iNE(1:grid%nNE(i),i))=.true.
      rhs((grid%iNE(1:grid%nNE(i),i)))=100d0
    end if
  end do
  call A%init(grid%nN,size(grid%iNE,1)**2*grid%nE)
  call findLaplacian(grid,A,isDirichlet)
  
  call A%fact()
  call A%solve(rhs,phi)
  
  do i=1,grid%nN
    if(abs(phi(i)-grid%pN(1,i)*100d0)>1d-6)then
      ierr=1
    end if
  end do
  
end function
