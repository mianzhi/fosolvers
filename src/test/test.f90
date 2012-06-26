!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleMPIComm
  type(typeGrid)::grid
  integer a,v(2,3)
  double precision b,d(3)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH1.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    a=2
    b=3d0
    v=reshape([1,2,3,4,5,6],[2,3])
    d=[1d0,2d0,3d0]
    call sendDat(a,1)
    call sendDat(b,1)
    call sendDat(v,1)
    call sendDat(d,1)
  else
    call recvDat(a,0)
    call recvDat(b,0)
    call recvDat(v,0)
    call recvDat(d,0)
    write(*,*),a,b
    write(*,*),v(:,1),v(:,2),v(:,3)
    write(*,*),d
    !open(13,file='rst.msh',status='replace')
    !call writeGMSH(13,grid)
    !close(13)
  end if
  call finalMPI()
end program
