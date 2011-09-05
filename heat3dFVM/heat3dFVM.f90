!----------------------------------------------------------------------------- best with 100 columns

!******************************
! declare the global variables
!******************************
module globalvar
  ! condition related variables
  double precision,save::tWrite,Tinit
  integer,save::ncodData
  double precision,save,allocatable::vRadSurf(:,:),vConvSurf(:,:),codData(:,:)
  integer,save,allocatable::vTempSurf(:),vFluxSurf(:),vSource(:)
  logical,save,allocatable::lTempSurf(:),lFluxSurf(:),lRadSurf(:),lConvSurf(:),lSource(:)
  ! material related variables
  integer,save,allocatable::nrhoTab(:),ncTab(:),nkTab(:),neTab(:),iVol(:)
  double precision,save,allocatable::rhoTab(:,:,:),cTab(:,:,:),kTab(:,:,:),eTab(:,:,:)
  character(3),save,allocatable::mtlid(:)
  ! FVM solving related variables & results
  double precision,save::dt,tolerance
  double precision,save,allocatable::Temp(:),e(:),rho(:),c(:),cond(:),&
  &                                  Vol(:),Paux(:,:,:),gradT(:,:)
  integer,save::numIt(3)
  double precision,save::beta
end module

!**************
! main program
!**************
program heat3dFVM
  use moduleGrid
  use moduleWrite
  use globalvar
  
  double precision::finddt ! functions will be called
  
  character(100) gridfname,codfname,mtlfname,rstfname ! for safer passing into subroutines
  write(gridfname,'(a)'),'grid.msh'
  write(codfname,'(a)'),'conditions.cod'
  write(mtlfname,'(a)'),'materials.mtl'
  write(rstfname,'(a)'),'rst.msh'
  
  call readmsh(gridfname,10)
  call sortEle()
  call updateFacetPara()
  call updateElePara()
  call initWriteEnv(0,0,0,1,0,0)
  
  if(nEle==0)then
    write(*,'(a)'),'ERROR: no solid element is found in "grid.msh"'
    stop
  end if
  
  ! read simulation conditions
  call readcod(codfname,11)
  ! read material properties
  call readmtl(mtlfname,12)
  
  !---------------------
  ! prepare for solving
  !---------------------
  ! build the table of specific energy e(T)[J/kg]
  call geneTab()
  ! initialize variables and parameters
  call initVar()
  ! find auxiliary points' offset
  call genPaux()
  ! write initial data
  t=0d0
  nWrite=0
  rstEleScal(1,:)=Temp(1:nEle)
  call writerst(rstfname,15,.true.)
  
  !-------
  ! solve
  !-------
  write(*,'(/,a,g10.4,a)'),'simulaiton time:',tFinal,'[sec]'
  numIt=[15,20,50] ! lower & upper limit of target iteration number, maximum iteration number
  beta=0.25d0 ! initial under-relaxation factor
  tolerance=1d-3 ! in [K]
  l=-1 ! l is the number of times of under-relaxation iteration; l<0 indicates the 1st run
  do while(t<tFinal)
    dt=finddt(l) ! dt is adjusted according to l of the last time step
    call CrankNicolson(l) ! note: l is updated by CrankNicolson()
    ! write results
    if(t/tWrite>=nWrite)then
      rstEleScal(1,:)=Temp(1:nEle)
      call writerst(rstfname,15,.true.)
    end if
    ! progress bar
    call showProgCN(t,tFinal,l,dt)
  end do
  
  write(*,'(/,a)'),'all done!'
end program
