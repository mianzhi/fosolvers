!----------------------------------------------------------------------------- best with 100 columns

!******************************
! declare the global variables
!******************************
module globalvar
  ! mesh related variables
  double precision,save,allocatable::Pos(:,:)
  integer,save,allocatable::lPoint(:,:),lLine(:,:),lTri(:,:),lQuad(:,:),lTet(:,:),lHex(:,:),&
  &                         iVol(:),neighbour(:,:)
  integer,save::nNode,nEle,nPoint,nLine,nTri,nQuad,nTet,nHex,faceTet(4,3)
  double precision,save::PosRange(3,2)
  ! condition related variables
  double precision,save::tFinal,tWrite,Tinit
  integer,save::ncodData
  double precision,save,allocatable::vRadSurf(:,:),vConvSurf(:,:),codData(:,:)
  integer,save,allocatable::vTempSurf(:),vFluxSurf(:),vSource(:)
  logical,save,allocatable::lTempSurf(:),lFluxSurf(:),lRadSurf(:),lConvSurf(:),lSource(:)
  ! material related variables
  integer,save,allocatable::nrhoTab(:),ncTab(:),nkTab(:),neTab(:)
  double precision,save,allocatable::rhoTab(:,:,:),cTab(:,:,:),kTab(:,:,:),eTab(:,:,:)
  character*3,save,allocatable::mtlid(:)
  ! FVM solving related variables & results
  double precision,save::t,dt,tolerance
  double precision,save,allocatable::Temp(:),e(:),rho(:),c(:),cond(:),&
  &                                  Vol(:),Paux(:,:,:),gradT(:,:)
  integer,save::countWrite,numIt(3)
  double precision::beta
  
end module

!**************
! main program
!**************
program heat3dFVM
  !----------------------------------------------
  ! declare and initialize variables & functions
  !----------------------------------------------
  use globalvar
  double precision::dtGen,volTet
  
  ! read mesh file
  call readmsh()
  if(nHex/=0)then
    write(*,'(a)'),'WARNING: this program does not support Hex. element'
  end if
  if(nTet==0)then
    write(*,'(a)'),'ERROR: no solid element is found in "grid.msh"'
    stop
  end if
  ! read simulation conditions
  call readcod()
  ! read material properties
  call readmtl()
  
  !--------------
  ! prepare data
  !--------------
  ! find neighbour cells (including ghost cells)
  faceTet(1,:)=[1,3,2] ! surface node list
  faceTet(2,:)=[1,2,4]
  faceTet(3,:)=[1,4,3]
  faceTet(4,:)=[2,3,4]
  call findNeighbour()
  ! build the table of specific energy e(T)[J/kg]
  call eTabGen()
  ! find cell volume
  allocate(Vol(1:nTet))
  do i=1,nTet
    Vol(i)=volTet(Pos(lTet(i,1),:),Pos(lTet(i,2),:),Pos(lTet(i,3),:),Pos(lTet(i,4),:))
  end do
  ! initialize variables and parameters
  call init()
  ! find auxiliary points' offset
  call auxGen()
  ! write initial data
  t=0d0
  countWrite=0
  call writerst()
  
  !-------
  ! solve
  !-------
  write(*,'(/,a,g10.4,a)'),'simulaiton time:',tFinal,'[sec]'
  
  numIt=[15,20,50] ! lower & upper limit of target iteration number, maximum iteration number
  beta=0.37d0 ! initial under-relaxation factor
  tolerance=1d-3 ! in [K]
  l=-1 ! l is the number of times of under-relaxation iteration; l<0 indicates the 1st run
  do while(t<tFinal)
    dt=dtGen(l) ! dt is adjusted according to l of the last time step
    call CrankNicolson(l) ! note: l is updated by CrankNicolson()
    ! write results
    if(t/tWrite>=countWrite)then
      call writerst()
    end if
    ! progress bar
    call progBar(t,tFinal,l,dt)
  end do
  
  write(*,'(/,a)'),'all done!'
  
  !----------
  ! clean up
  !----------
  deallocate(Pos)
  deallocate(lPoint)
  deallocate(lLine)
  deallocate(lTri)
  deallocate(lQuad)
  deallocate(lTet)
  deallocate(lHex)
  deallocate(iVol)
  deallocate(mtlid)
  deallocate(lTempSurf)
  deallocate(lFluxSurf)
  deallocate(lRadSurf)
  deallocate(lConvSurf)
  deallocate(lSource)
  deallocate(vTempSurf)
  deallocate(vFluxSurf)
  deallocate(vRadSurf)
  deallocate(vConvSurf)
  deallocate(vSource)
  deallocate(codData)
  deallocate(nrhoTab)
  deallocate(ncTab)
  deallocate(nkTab)
  deallocate(neTab)
  deallocate(rhoTab)
  deallocate(cTab)
  deallocate(kTab)
  deallocate(eTab)
  deallocate(neighbour)
  deallocate(Temp)
  deallocate(gradT)
  deallocate(rho)
  deallocate(c)
  deallocate(cond)
  deallocate(e)
  deallocate(Vol)
  deallocate(Paux)
  
end program
