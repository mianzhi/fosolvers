!----------------------------------------------------------------------------- best with 100 columns

!****************************
! read simulation conditions
!****************************
subroutine readcod(fname,codfile)
  use moduleGrid
  use moduleWrite
  use globalvar
  
  character(100),intent(in)::fname
  integer,intent(in)::codfile
  integer readerr
  character(200) temp_string
  double precision temp_double_vect(2)
  integer temp_int
  
  open(codfile,file=fname,status='old')
  readerr=0
  
  allocate(lTempSurf(nFacet))
  allocate(lFluxSurf(nFacet))
  allocate(lRadSurf(nFacet))
  allocate(lConvSurf(nFacet))
  allocate(lSource(nEle))
  allocate(vTempSurf(nFacet))
  allocate(vFluxSurf(nFacet))
  allocate(vRadSurf(nFacet,2))
  allocate(vConvSurf(nFacet,2))
  allocate(vSource(nEle))
  
  lTempSurf(:)=.false.
  lFluxSurf(:)=.false.
  lRadSurf(:)=.false.
  lConvSurf(:)=.false.
  lSource(:)=.false.
  vTempSurf(:)=0
  vFluxSurf(:)=0
  vRadSurf(:,:)=0d0
  vConvSurf(:,:)=0d0
  vSource(:)=0
  
  m=0 ! number of variables to be linked
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(codfile,'(a200)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    ! read simulation time
    if(temp_string(1:4)=='$Sim')then
      read(codfile,*,iostat=readerr),tFinal,tWrite
    end if
    ! read material label
    if(temp_string(1:4)=='$Mat')then
      read(codfile,*,iostat=readerr),n
      allocate(iVol(n))
      allocate(mtlid(n))
      do i=1,n
        read(codfile,*,iostat=readerr),iVol(i),mtlid(i)
      end do
    end if
    ! read initial temperature
    if(temp_string(1:4)=='$Ini')then
      read(codfile,*,iostat=readerr),Tinit
    end if
    ! read surface temperature
    if(temp_string(1:9)=='$TempSurf')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_int
        do i=1,nFacet
          if(Facet(i)%GeoEnti==l)then
            lTempSurf(i)=.true.
            vTempSurf(i)=temp_int
          end if
        end do
      end do
    end if
    ! read surface heat flux
    if(temp_string(1:9)=='$FluxSurf')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_int
        do i=1,nFacet
          if(Facet(i)%GeoEnti==l)then
            lFluxSurf(i)=.true.
            vFluxSurf(i)=temp_int
          end if
        end do
      end do
    end if
    ! read surface radiation
    if(temp_string(1:8)=='$RadSurf')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect
        do i=1,nFacet
          if(Facet(i)%GeoEnti==l)then
            lRadSurf(i)=.true.
            vRadSurf(i,:)=temp_double_vect
          end if
        end do
      end do
    end if
    ! read surface convection
    if(temp_string(1:9)=='$ConvSurf')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect
        do i=1,nFacet
          if(Facet(i)%GeoEnti==l)then
            lConvSurf(i)=.true.
            vConvSurf(i,:)=temp_double_vect
          end if
        end do
      end do
    end if
    ! read volumetric source
    if(temp_string(1:11)=='$SourceBody')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_int
        do i=1,nEle
          if(Ele(i)%GeoEnti==l)then
            if(lSource(i).eqv..true.)then
              write(*,'(a)'),'WARNING: heat source over defined.'
            end if
            lSource(i)=.true.
            vSource(i)=temp_int
          end if
        end do
      end do
    end if
    ! read source for element
    if(temp_string(1:10)=='$SourceEle')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_int
        do i=1,nEle
          if(i==l)then
            if(lSource(i).eqv..true.)then
              write(*,'(a)'),'WARNING: heat source over defined.'
            end if
            lSource(i)=.true.
            vSource(i)=temp_int
          end if
        end do
      end do
    end if
    ! read time-variant data
    if(temp_string(1:5)=='$Data')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),ncodData,m
        allocate(codData(0:m,ncodData))
        do i=1,ncodData
          read(codfile,*,iostat=readerr),codData(:,i)
        end do
      end do
    end if
  end do
  
  ! check if linked variables are valid
  k=max(maxval(vTempSurf),maxval(vFluxSurf),maxval(vSource))
  if(k>m)then
    write(*,'(a,i5)'),'ERROR: (in "conditions.cod") invalid variable index:',k
    stop
  end if
  
  close(codfile)
end subroutine

!**************************
! read material properties
!**************************
subroutine readmtl(fname,mtlfile)
  use globalvar
  
  character(100),intent(in)::fname
  integer,intent(in)::mtlfile
  integer readerr
  character(100) temp_string
  character(20) mtlname
  
  n=size(mtlid) ! number of volumes need to assign a material for
  allocate(nrhoTab(n))
  allocate(ncTab(n))
  allocate(nkTab(n))
  
  open(mtlfile,file=fname,status='old')
  readerr=0
  
  do i=1,n ! for all sub-regions
    do while(readerr==0)
      ! skip the irrelevant lines
      do while(readerr==0)
        read(mtlfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:1)=='$')then
          exit
        end if
      end do
      ! check if finished but have not found the material
      if(readerr/=0)then
        write(*,'(a,a,a,i4)'),'ERROR: can not find material "',mtlid(i),'" for sub-region ',iVol(i)
        stop
      end if
      ! find required data size
      if(temp_string(2:4)==mtlid(i))then ! if matched
        mtlname=temp_string(6:25)
        write(*,'(a,a,a,i4)'),'assign material "',trim(mtlname),'" to sub-region ',iVol(i)
        read(mtlfile,*,iostat=readerr),nrhoTab(i)
        do j=1,nrhoTab(i) ! skip data
          read(mtlfile,*,iostat=readerr),temp_string
        end do
        read(mtlfile,*,iostat=readerr),ncTab(i)
        do j=1,ncTab(i) ! skip data
          read(mtlfile,*,iostat=readerr),temp_string
        end do
        read(mtlfile,*,iostat=readerr),nkTab(i)
        do j=1,nkTab(i) ! skip data
          read(mtlfile,*,iostat=readerr),temp_string
        end do
        exit
      end if
    end do
    ! return to the biginning to read for the next sub-region
    rewind(mtlfile)
    readerr=0
  end do
  
  allocate(rhoTab(n,2,maxval(nrhoTab,1)))
  allocate(cTab(n,2,maxval(ncTab,1)))
  allocate(kTab(n,2,maxval(nkTab,1)))
    
  do i=1,n !for all sub-regions
    do while(readerr==0)
      ! skip the irrelevant lines
      do while(readerr==0)
        read(mtlfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:1)=='$')then
          exit
        end if
      end do
      ! read the property tables
      if(temp_string(2:4)==mtlid(i))then ! if matched
        read(mtlfile,*,iostat=readerr),nrhoTab(i)
        do j=1,nrhoTab(i)
          read(mtlfile,*,iostat=readerr),rhoTab(i,:,j)
        end do
        read(mtlfile,*,iostat=readerr),ncTab(i)
        do j=1,ncTab(i)
          read(mtlfile,*,iostat=readerr),cTab(i,:,j)
        end do
        read(mtlfile,*,iostat=readerr),nkTab(i)
        do j=1,nkTab(i)
          read(mtlfile,*,iostat=readerr),kTab(i,:,j)
        end do
        exit
      end if
    end do
    
    ! return to the biginning to read for the next sub-region
    rewind(mtlfile)
    readerr=0
  end do
  
  close(mtlfile)
end subroutine
