!----------------------------------------------------------------------------- best with 100 columns

function readVTK1() result(ierr)
  use modFileIO
  use modPolyMesh
  integer ierr
  type(polyMesh)::mesh
  
  ierr=0
  open(10,file='data/mesh_simple.vtk',action='read')
  call readVTK(10,mesh)
  close(10)
  if(any(mesh%iNE/=reshape([1,3,2,1,4,3,1,2,4,2,3,4],shape(mesh%iNE))))then
    ierr=ierr+1
  end if
  if(any(mesh%sE/=TRI))then
    ierr=ierr+10
  end if
  if(any(mesh%nNE/=TRI_N))then
    ierr=ierr+100
  end if
  call mesh%clear()
end function
