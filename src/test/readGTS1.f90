!----------------------------------------------------------------------------- best with 100 columns

function readGTS1()
  use modFileIO
  use modPolyMesh
  integer readGTS1
  type(polyMesh)::mesh
  
  readGTS1=0
  open(10,file='data/mesh_simple.gts',action='read')
  call readGTS(10,mesh)
  if(any(mesh%iNE/=reshape([1,2,3,1,4,2,2,4,3,3,4,1],shape(mesh%iNE))))then
    readGTS1=readGTS1+1
  end if
  call mesh%clear()
  close(10)
end function
