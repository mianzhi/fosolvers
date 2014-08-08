!----------------------------------------------------------------------------- best with 100 columns

function polyMesh1()
  use modFileIO
  use modPolyMesh
  integer polyMesh1
  type(polyMesh)::mesh
  
  polyMesh1=0
  open(10,file='data/mesh_simple.gts',action='read')
  call readGTS(10,mesh)
  close(10)
  call mesh%up()
  if(any(abs(mesh%norm-reshape([0d0,0d0,-1d0,&
  &  -0.81649720746616528d0,0.47140365723169048d0,0.33333301988933012d0,&
  &  0.81649720746616528d0,0.47140365723169048d0,0.33333301988933012d0,&
  &  0.d0,-0.94280905929685954d0,0.33333328322831940d0],shape(mesh%iNE)))>1d-14))then
    polyMesh1=polyMesh1+1
  end if
  if(.not.norm2(matmul(mesh%norm,mesh%a))<=1d-14)then
    polyMesh1=polyMesh1+10
  end if
  call mesh%clear()
end function
