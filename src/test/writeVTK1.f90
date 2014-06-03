!----------------------------------------------------------------------------- best with 100 columns

function writeVTK1() result(ierr)
  use modFileIO
  use modPolyMesh
  integer ierr
  type(polyMesh)::mesh
  
  ierr=0
  open(10,file='data/mesh_simple.gts',action='read')
  call readGTS(10,mesh)
  close(10)
  open(10,file='mesh_simple.vtk',action='write')
  call writeVTK(10,mesh)
  close(10)
  call mesh%clear()
end function
