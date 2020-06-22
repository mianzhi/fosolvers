!----------------------------------------------------------------------------- best with 100 columns

function edgeFEM1() result(ierr)
  use modPolyEdgeFeGrid
  use modFileIO
  use modSparse
  use modEdgeFEM
  type(polyEdgeFeGrid)::grid
  type(multiFront)::curlCurl
  double precision,allocatable::A(:),rhs(:),J(:,:),H(:,:),full(:,:)
  logical,allocatable::isInactive(:)
  integer::ierr
  
  ierr=0
  open(10,file='data/wireInBox.msh',action='read')
    call readGMSH(10,grid)
  close(10)
  call grid%up()
  
  allocate(A(grid%nEdge))
  allocate(rhs(grid%nEdge))
  allocate(isInactive(grid%nEdge))
  allocate(J(3,grid%nC))
  allocate(H(3,grid%nC))
  allocate(full(3,grid%nE))
  
  J(:,:)=0d0
  rhs(:)=0d0
  isInactive(:)=.false.
  do i=1,grid%nC
    if(grid%gid(i)==2)then
      J(:,i)=[0d0,0d0,1d3]
    end if
  end do
  call findVolVectSrc(grid,J,rhs)
  do i=grid%nC+1,grid%nE
    if(grid%gid(i)==11.or.grid%gid(i)==12.or.grid%gid(i)==15.or.grid%gid(i)==17.or.&
    &  grid%gid(i)==13.or.grid%gid(i)==14.or.grid%gid(i)==18.or.grid%gid(i)==16)then
      isInactive(grid%iEdgeE(1:grid%nEdgeE(i),i))=.true.
      rhs(grid%iEdgeE(1:grid%nEdgeE(i),i))=0d0
    end if
  end do
  call curlCurl%init(grid%nEdge,size(grid%iEdgeE,1)**2*grid%nC)
  call findCurlCurl(grid,curlCurl,isInactive)
  
  call curlCurl%fact()
  call curlCurl%solve(rhs,A)
  
  call findCellCurl(grid,A,H)
  
  open(10,file='edgeFEM1.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  full(:,:)=0d0
  full(:,1:grid%nC)=H(:,:)
  call writeVTK(10,'H',full)
  full(:,:)=0d0
  full(:,1:grid%nC)=J(:,:)
  call writeVTK(10,'J',full)
  close(10)
  
end function
