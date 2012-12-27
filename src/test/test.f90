!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridInspection
  use moduleGridOperation
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleCondition
  use moduleNonlinearSolve
  use moduleMPIComm
  type(typeGrid)::grid
  double precision box(DIMS,2),u(20)
  type(typeCondition),allocatable::condition(:)
  external::testfun
  
  u(:)=10d0
  ProblemFunc=>testfun
  call solveNonlinear(u)
  do i=1,size(u)
    write(*,*),i,u(i)
  end do
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    box=findBoundBox(grid)
    write(*,*),box(:,1)
    write(*,*),box(:,2)
    
    open(13,file='bin/condition1',status='old')
    call readCondition(13,condition)
    write(*,*),condition(1)%Ent,condition(1)%dat%get('a_longer_name')
    write(*,*),condition(2)%Ent,condition(2)%dat%get('str')
  else
  end if
  call finalMPI()
end program

function testfun(r)
  double precision r(:)
  double precision testfun(size(r))
  n=size(r)
  forall(i=2:n-1)
    testfun(i)=-r(i-1)+2d0*r(i)-r(i+1)-1d0
  end forall
  testfun(1)=r(1)+2d0*r(1)-r(2)-1d0
  testfun(n)=-r(n-1)+2d0*r(n)+r(n)-1d0
end function
