!----------------------------------------------------------------------------- best with 100 columns

function writeVTK1() result(ierr)
  use modFileIO
  use modPolyMesh
  use modOtGrid
  integer ierr
  type(polyMesh)::mesh
  type(otGrid)::grid
  double precision::r(4096)
  
  ierr=0
  open(10,file='data/mesh_simple.gts',action='read')
  call readGTS(10,mesh)
  close(10)
  open(10,file='mesh_simple.vtk',action='write')
  call writeVTK(10,mesh)
  close(10)
  call grid%init([0d0,0d0,0d0],1d0,lvl=4)
  do i=1,grid%nC
    r(i)=norm2(grid%p(i))
  end do
  open(10,file='grid_octree.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,'test',r)
  close(10)
  call mesh%clear()
  call grid%clear()
end function
