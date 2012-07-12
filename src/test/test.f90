!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleMPIComm
  type(typeGrid)::grid
  type(typeGenDatLlist)::ll,temp
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH4.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    temp%key=['a']
    temp%datType=TAB1D_TYPE
    allocate(temp%dat(3,0:1),source=reshape([0d0,1d0,2d0,1d1,2d1,3d1],[3,2]))
    call ll%push(temp)
    temp%key=['b','c']
    call ll%push(temp)
    temp%key=['d','e','f']
    call ll%push(temp)
    write(*,*),ll%test('bc'),ll%test('ef')
    write(*,*),ll%get('a',[-0.3d0]),ll%get('def',[5.1d0]),ll%get('123')
  else
  end if
  call finalMPI()
end program
