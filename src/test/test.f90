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
  integer a,c(2),e(4,2)
  double precision b,d(3),f(3,2)
  integer,allocatable::cc(:),ee(:,:)
  double precision,allocatable::dd(:),ff(:,:)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH1.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    a=2
    b=2.4d0
    c=[4,2]
    d=[4d0,7d0,19d0]
    e=reshape([1,2,3,4,5,6,7,8],[4,2])
    f=reshape([1d0,2d0,3d0,4d0,5d0,6d0],[3,2])
    call sendDat(a,1)
    call sendDat(b,1)
    call sendDat(c,1)
    call sendDat(c,1)
    call sendDat(d,1)
    call sendDat(d,1)
    call sendDat(e,1)
    call sendDat(e,1)
    call sendDat(f,1)
    call sendDat(f,1)
  else
    call recvDat(a,0)
    call recvDat(b,0)
    call recvDat(c,0)
    call recvDat(cc,0,realloc=.true.)
    call recvDat(d,0)
    call recvDat(dd,0,realloc=.true.)
    call recvDat(e,0)
    call recvDat(ee,0,realloc=.true.)
    call recvDat(f,0)
    call recvDat(ff,0,realloc=.true.)
    write(*,*),a
    write(*,*),b
    write(*,*),c
    write(*,*),cc
    write(*,*),d
    write(*,*),dd
    do i=1,size(e,2)
      write(*,*),e(:,i)
      write(*,*),ee(:,i)
    end do
    do i=1,size(f,2)
      write(*,*),f(:,i)
      write(*,*),ff(:,i)
    end do
    !write(*,*),v(:,1),v(:,2),v(:,3)
    !open(13,file='rst.msh',status='replace')
    !call writeGMSH(13,grid)
    !close(13)
  end if
  call finalMPI()
end program
