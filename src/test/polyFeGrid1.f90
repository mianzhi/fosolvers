!----------------------------------------------------------------------------- best with 100 columns

function polyFeGrid1() result(ierr)
  use modFileIO
  use modPolyFeGrid
  integer ierr
  type(polyFeGrid)::grid
  
  ierr=0
  open(10,file='data/tet10.msh',action='read')
  call readGMSH(10,grid)
  call grid%up()
  do i=1,grid%nC
    if(any(abs(grid%detJ(1:4,i)-0.25d0)>1d-9)) ierr=ierr+1
  end do
  do i=grid%nC+1,grid%nE
    if(any(abs(grid%detJ(1:3,i)-0.5d0)>1d-9)) ierr=ierr+100
  end do
  call grid%up()
  close(10)
end function
