!----------------------------------------------------------------------------- best with 100 columns

!> the momentum residual function (not coupled with pressure)
function resMom(testU1d)
  use moduleGrid
  use moduleFVMGrad
  use miscNS
  double precision testU1d(:) !< the test velocity (unwrapped as 1d array)
  double precision resMom(size(testU1d)) !< the momentum residual function
  double precision testU(DIMS,size(testU1d)/DIMS),gradU(DIMS,DIMS,size(testU1d)/DIMS)
  
  testU=reshape(testU1d,[DIMS,grid%nNode])
  gradU=findGrad(testU,grid,BIND_NODE)
  
  resMom=reshape(testU-u,[DIMS*grid%nNode])
end function
