!----------------------------------------------------------------------------- best with 100 columns

!> foburn1d main program
program foburn1d
  use miscBurn1D
  use moduleGrid1D
  use moduleFVMConvect
  use moduleFVMDiffus
  
  call grid%genUniform(0d0,1d0,100)
  call setEnv()
  t=0d0
  dt=1d-5
  tFinal=1d0
  
  ! advance in time
  do while(t<tFinal)
    
    t=t+dt
  end do
  
  ! clean up
  call clearEnv()
  
end program
