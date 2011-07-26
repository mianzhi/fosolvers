!----------------------------------------------------------------------------- best with 100 columns

!***************
! write results
!***************
subroutine writerst()
  use globalvar,only:nNode,nTet,lTet,Pos,countWrite,t,tWrite,tFinal,Temp
  integer ifile
  character*200 tempstring
  
  ifile=10
  if(countWrite==0)then ! if it is the 1st time writing results
    open(ifile,file='rst.msh',status='replace')
    ! write header
    write(ifile,'(a)'),'$MeshFormat'
    write(ifile,'(a)'),'2.2 0 8'
    write(ifile,'(a)'),'$EndMeshFormat'
    ! write nodes
    write(ifile,'(a)'),'$Nodes'
    write(tempstring,*),nNode
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
    do i=1,nNode
      write(tempstring,'(i8,1x,g14.7,1x,g14.7,1x,g14.7)'),i,Pos(i,1),Pos(i,2),Pos(i,3)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndNodes'
    ! write elements
    write(ifile,'(a)'),'$Elements'
    write(tempstring,*),nTet
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
    do i=1,nTet
      write(tempstring,'(i8,1x,i1,1x,i1,1x,i1,1x,i5,1x,i8,1x,i8,1x,i8,1x,i8)'),i,4,2,0,lTet(i,5),&
      &    lTet(i,1),lTet(i,2),lTet(i,3),lTet(i,4)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndElements'
  end if
  
  write(ifile,'(a)'),'$ElementData'
  write(ifile,'(i1)'),1 ! 1 string tag
  write(ifile,'(a)'),'"Temperature"' ! name of the data-set
  write(ifile,'(i1)'),1 ! 1 real tag
  write(tempstring,'(g14.7)'),t ! time
  write(ifile,'(a)'),trim(adjustl(tempstring))
  write(ifile,'(i1)'),3 ! 3 integer tags
  write(tempstring,'(i5)'),countWrite ! time step index
  write(ifile,'(a)'),trim(adjustl(tempstring))
  write(ifile,'(i1)'),1 ! 1 component scaler
  write(tempstring,*),nTet
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nTet
    write(tempstring,'(i8,1x,g14.7)'),i,Temp(i)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElementData'
  
  countWrite=countWrite+1
  
  if(t>tFinal)then ! if it is the last time writing results
    close(ifile)
  end if
end subroutine

!****************************
! read simulation conditions
!****************************
subroutine readcod()
  use globalvar,only:tFinal,tWrite,iVol,mtlid,Tinit,&
  &                  vTempSurf,vFluxSurf,vRadSurf,vConvSurf,vSource,codData,ncodData,&
  &                  lTempSurf,lFluxSurf,lRadSurf,lConvSurf,lSource,&
  &                  lPoint,lLine,lTri,lQuad,lTet,lHex,&
  &                  nPoint,nLine,nTri,nQuad,nTet,nHex,nNode,&
  &                  Pos
  
  integer codfile,readerr
  character*200 temp_string
  double precision temp_double_vect(2),temp_int
  
  codfile=13
  
  open(codfile,file='conditions.cod',status='old')
  readerr=0
  
  allocate(lTempSurf(nTri))
  allocate(lFluxSurf(nTri))
  allocate(lRadSurf(nTri))
  allocate(lConvSurf(nTri))
  allocate(lSource(nTet))
  allocate(vTempSurf(nTri))
  allocate(vFluxSurf(nTri))
  allocate(vRadSurf(nTri,2))
  allocate(vConvSurf(nTri,2))
  allocate(vSource(nTet))
  
  lTempSurf(:)=.false.
  lFluxSurf(:)=.false.
  lRadSurf(:)=.false.
  lConvSurf(:)=.false.
  lSource(:)=.false.
  vTempSurf(:)=0d0
  vFluxSurf(:)=0d0
  vRadSurf(:,:)=0d0
  vConvSurf(:,:)=0d0
  vSource(:)=0d0
  
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
        do i=1,nTri
          if(lTri(i,4)==l)then
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
        do i=1,nTri
          if(lTri(i,4)==l)then
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
        do i=1,nTri
          if(lTri(i,4)==l)then
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
        do i=1,nTri
          if(lTri(i,4)==l)then
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
        do i=1,nTet
          if(lTet(i,5)==l)then
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
        do i=1,nTet
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
subroutine readmtl()
  use globalvar,only:iVol,mtlid,nrhoTab,ncTab,nkTab,rhoTab,cTab,kTab
  
  integer mtlfile,readerr
  character*100 temp_string
  character*20 mtlname
  
  n=size(mtlid) ! number of volumes need to assign a material for
  allocate(nrhoTab(n))
  allocate(ncTab(n))
  allocate(nkTab(n))
  
  mtlfile=12
  
  open(mtlfile,file='materials.mtl',status='old')
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

!****************
! read mesh file
!****************
subroutine readmsh()
  use globalvar,only:Pos,&
  &                  lPoint,lLine,lTri,lQuad,lTet,lHex,&
  &                  nNode,nEle,nPoint,nLine,nTri,nQuad,nTet,nHex,&
  &                  PosRange
  
  integer gridfile,readerr,np,nt
  character*100 temp_string
  integer temp_int_vect1(20)
  integer,allocatable::temp_int_mat1(:,:)
  
  np=0
  
  gridfile=11
  
  open(gridfile,file='grid.msh',status='old')
  readerr=0
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(gridfile,'(a100)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    ! read node data
    if(temp_string(1:4)=='$NOD'.or.temp_string(1:4)=='$Nod')then
      read(gridfile,*,iostat=readerr),nNode
      allocate(Pos(nNode,3))
      read(gridfile,*,iostat=readerr),j,Pos(1,:)
      PosRange(:,1)=Pos(1,:)
      PosRange(:,2)=Pos(1,:)
      do i=2,nNode
        read(gridfile,*,iostat=readerr),j,Pos(i,:)
        PosRange(1,1)=min(PosRange(1,1),Pos(i,1))
        PosRange(1,2)=max(PosRange(1,2),Pos(i,1))
        PosRange(2,1)=min(PosRange(2,1),Pos(i,2))
        PosRange(2,2)=max(PosRange(2,2),Pos(i,2))
        PosRange(3,1)=min(PosRange(3,1),Pos(i,3))
        PosRange(3,2)=max(PosRange(3,2),Pos(i,3))
        if(i/=j)then
          write(*,'(a)'),'WARNING: node data may be not in sequence'
        end if
      end do
      cycle
    end if
    ! read element data
    if(temp_string(1:4)=='$ELM'.or.temp_string(1:4)=='$Ele')then
      if(temp_string(1:4)=='$ELM')then
        write(*,'(a)'),'WARNING: not capable with GMSH version 1 format'
        write(*,'(a)'),'  please try using the latest version of GMSH'
      end if
      read(gridfile,*,iostat=readerr),nEle
      allocate(lPoint(nEle,2))
      allocate(lLine(nEle,3))
      allocate(lTri(nEle,4))
      allocate(lQuad(nEle,5))
      allocate(lTet(nEle,5))
      allocate(lHex(nEle,9))
      !note the additional byte to storage geometric entity index
      nPoint=0
      nLine=0
      nTri=0
      nQuad=0
      nTet=0
      nHex=0
      do i=1,nEle
        read(gridfile,'(a)',iostat=readerr),temp_string
        read(temp_string,*),j,k,nt
        if(i/=j)then
          write(*,'(a)'),'WARNING: element data may be not in sequence'
        end if
        select case(k) !type of element
          case(15) !point
            np=1
          case(1) !line
            np=2
          case(2) !tri
            np=3
          case(3) !quad
            np=4
          case(4) !tet
            np=4
          case(5) !hex
            np=8
        end select
        read(temp_string,*),temp_int_vect1(1:3+nt+np)
        select case(temp_int_vect1(2)) !type of element
          case(15) !point
            nPoint=nPoint+1
            lPoint(nPoint,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lPoint(nPoint,np+1)=temp_int_vect1(3+nt)
          case(1) !line
            nLine=nLine+1
            lLine(nLine,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lLine(nLine,np+1)=temp_int_vect1(3+nt)
          case(2) !tri
            nTri=nTri+1
            lTri(nTri,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lTri(nTri,np+1)=temp_int_vect1(3+nt)
          case(3) !quad
            nQuad=nQuad+1
            lQuad(nQuad,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lQuad(nQuad,np+1)=temp_int_vect1(3+nt)
          case(4) !tet
            nTet=nTet+1
            lTet(nTet,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lTet(nTet,np+1)=temp_int_vect1(3+nt)
          case(5) !hex
            nHex=nHex+1
            lHex(nHex,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lHex(nHex,np+1)=temp_int_vect1(3+nt)
        end select
      end do
      !trim lPoint
      if(nPoint>0)then
        allocate(temp_int_mat1(nPoint,2))
        temp_int_mat1=lPoint(1:nPoint,:)
        deallocate(lPoint)
        allocate(lPoint(nPoint,2))
        lPoint=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lPoint)
        allocate(lPoint(1,1))
      end if
      !trim lLine
      if(nLine>0)then
        allocate(temp_int_mat1(nLine,3))
        temp_int_mat1=lLine(1:nLine,:)
        deallocate(lLine)
        allocate(lLine(nLine,3))
        lLine=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lLine)
        allocate(lLine(1,1))
      end if
      !trim lTri
      if(nTri>0)then
        allocate(temp_int_mat1(nTri,4))
        temp_int_mat1=lTri(1:nTri,:)
        deallocate(lTri)
        allocate(lTri(nTri,4))
        lTri=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lTri)
        allocate(lTri(1,1))
      end if
      !trim lQuad
      if(nQuad>0)then
        allocate(temp_int_mat1(nQuad,5))
        temp_int_mat1=lQuad(1:nQuad,:)
        deallocate(lQuad)
        allocate(lQuad(nQuad,5))
        lQuad=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lQuad)
        allocate(lQuad(1,1))
      end if
      !trim lTet
      if(nTet>0)then
        allocate(temp_int_mat1(nTet,5))
        temp_int_mat1=lTet(1:nTet,:)
        deallocate(lTet)
        allocate(lTet(nTet,5))
        lTet=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lTet)
        allocate(lTet(1,1))
      end if
      !trim lHex
      if(nHex>0)then
        allocate(temp_int_mat1(nHex,9))
        temp_int_mat1=lHex(1:nHex,:)
        deallocate(lHex)
        allocate(lHex(nHex,9))
        lHex=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lHex)
        allocate(lHex(1,1))
      end if
      cycle
    end if
  end do
  
  write(*,'(a)'),'grid summary:'
  write(*,'(a,i8)'),'  number of nodes: ',nNode
  write(*,'(a,i8)'),'  number of points:',nPoint
  write(*,'(a,i8)'),'  number of lines: ',nLine
  write(*,'(a,i8)'),'  number of tri.:  ',nTri
  write(*,'(a,i8)'),'  number of quad.: ',nQuad
  write(*,'(a,i8)'),'  number of tet.:  ',nTet
  write(*,'(a,i8,/)'),'  number of hex.:  ',nHex
  
  close(gridfile)
end subroutine
