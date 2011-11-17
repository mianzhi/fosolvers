!----------------------------------------------------------------------------- best with 100 columns

!*********************************
! result output related variables
!*********************************
module moduleWrite
  use moduleMiscDataStruct

  ! output control related variables
  double precision,save::t,tFinal
  integer,save::nWrite
  ! data to be output
  type(typePtrScalArray),allocatable,save::rstNodeScal(:),rstEleScal(:)
  type(typePtrVectArray),allocatable,save::rstNodeVect(:),rstEleVect(:)
  type(typePtrTensArray),allocatable,save::rstNodeTens(:),rstEleTens(:)
  integer,save::nrstNodeScal,nrstNodeVect,nrstNodeTens,nrstEleScal,nrstEleVect,nrstEleTens
end module

!****************
! read mesh file
!****************
subroutine readmsh(fname,gridfile)
  use moduleGrid
  
  integer gridfile,readerr,np,nt
  character(100) temp_string,fname
  integer temp_int_vect1(50)
  type(typePoint),allocatable::tempPoint(:)
  type(typeLine),allocatable::tempLine(:)
  type(typeTri),allocatable::tempTri(:)
  type(typeQuad),allocatable::tempQuad(:)
  type(typeTet),allocatable::tempTet(:)
  type(typeHex),allocatable::tempHex(:)
  
  np=0
  
  open(gridfile,file=fname,status='old')
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
      allocate(Node(nNode))
      do i=1,nNode
        read(gridfile,*,iostat=readerr),j,Node(i)%Pos(:)
        if(i/=j)then
          write(*,'(a)'),'WARNING: node data may be not in sequence'
        end if
      end do
      cycle
    end if
    ! read element (includes facet & point) data
    if(temp_string(1:4)=='$ELM'.or.temp_string(1:4)=='$Ele')then
      if(temp_string(1:4)=='$ELM')then
        write(*,'(a)'),'WARNING: not capable with GMSH version 1 format'
        write(*,'(a)'),'HINT: please try using the latest version of GMSH'
      end if
      ! NOTE: the value of nEle will be changed later
      read(gridfile,*,iostat=readerr),nEle
      ! NOTE: allocate sufficient space to each object array
      allocate(Point(nEle))
      allocate(Line(nEle))
      allocate(Tri(nEle))
      allocate(Quad(nEle))
      allocate(Tet(nEle))
      allocate(Hex(nEle))
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
        select case(k) ! type of element
          case(POINT_TYPE) ! point
            np=POINT_NODE_NUM
          case(LINE_TYPE) ! line
            np=LINE_NODE_NUM
          case(TRI_TYPE) ! tri
            np=TRI_NODE_NUM
          case(QUAD_TYPE) ! quad
            np=QUAD_NODE_NUM
          case(TET_TYPE) ! tet
            np=TET_NODE_NUM
          case(HEX_TYPE) ! hex
            np=HEX_NODE_NUM
        end select
        read(temp_string,*),temp_int_vect1(1:3+nt+np)
        select case(temp_int_vect1(2)) ! type of element
          case(POINT_TYPE) ! point
            nPoint=nPoint+1
            Point(nPoint)%NodeInd=temp_int_vect1(3+nt+1)
            if(nt>=2)then
              Point(nPoint)%GeoEnti=temp_int_vect1(5)
            else
              Point(nPoint)%GeoEnti=0
            end if
          case(LINE_TYPE) ! line
            nLine=nLine+1
            Line(nLine)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            if(nt>=2)then
              Line(nLine)%GeoEnti=temp_int_vect1(5)
            else
              Line(nLine)%GeoEnti=0
            end if
          case(TRI_TYPE) ! tri
            nTri=nTri+1
            Tri(nTri)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            if(nt>=2)then
              Tri(nTri)%GeoEnti=temp_int_vect1(5)
            else
              Tri(nTri)%GeoEnti=0
            end if
          case(QUAD_TYPE) ! quad
            nQuad=nQuad+1
            Quad(nQuad)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            if(nt>=2)then
              Quad(nQuad)%GeoEnti=temp_int_vect1(5)
            else
              Quad(nQuad)%GeoEnti=0
            end if
          case(TET_TYPE) ! tet
            nTet=nTet+1
            Tet(nTet)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            if(nt>=2)then
              Tet(nTet)%GeoEnti=temp_int_vect1(5)
              if(nt>=4)then
                Tet(nTet)%Prt=maxval(temp_int_vect1(7:6+temp_int_vect1(6)))
              else
                Tet(nTet)%Prt=1
              end if
            else
              Tet(nTet)%GeoEnti=0
            end if
          case(HEX_TYPE) ! hex
            nHex=nHex+1
            Hex(nHex)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            if(nt>=2)then
              Hex(nHex)%GeoEnti=temp_int_vect1(5)
              if(nt>=4)then
                Hex(nHex)%Prt=maxval(temp_int_vect1(7:6+temp_int_vect1(6)))
              else
                Hex(nHex)%Prt=1
              end if
            else
              Hex(nHex)%GeoEnti=0
            end if
        end select
      end do
      ! trim Point
      if(nPoint>0)then
        allocate(tempPoint(nPoint))
        tempPoint(:)=Point(1:nPoint)
        call move_alloc(tempPoint,Point)
      else
        deallocate(Point)
        allocate(Point(1))
      end if
      ! trim Line
      if(nLine>0)then
        allocate(tempLine(nLine))
        tempLine(:)=Line(1:nLine)
        call move_alloc(tempLine,Line)
      else
        deallocate(Line)
        allocate(Line(1))
      end if
      ! trim Tri
      if(nTri>0)then
        allocate(tempTri(nTri))
        tempTri(:)=Tri(1:nTri)
        call move_alloc(tempTri,Tri)
      else
        deallocate(Tri)
        allocate(Tri(1))
      end if
      ! trim Quad
      if(nQuad>0)then
        allocate(tempQuad(nQuad))
        tempQuad(:)=Quad(1:nQuad)
        call move_alloc(tempQuad,Quad)
      else
        deallocate(Quad)
        allocate(Quad(1))
      end if
      ! trim Tet
      if(nTet>0)then
        allocate(tempTet(nTet))
        tempTet(:)=Tet(1:nTet)
        call move_alloc(tempTet,Tet)
      else
        deallocate(Tet)
        allocate(Tet(1))
      end if
      ! trim Hex
      if(nHex>0)then
        allocate(tempHex(nHex))
        tempHex(:)=Hex(1:nHex)
        call move_alloc(tempHex,Hex)
      else
        deallocate(Hex)
        allocate(Hex(1))
      end if
      ! construct Facet and Ele
      nFacet=nTri+nQuad
      allocate(Facet(nFacet))
      forall(i=1:nTri)
        Facet(i)%ShapeType=TRI_TYPE
        Facet(i)%ShapeInd=i
        Facet(i)%NodeNum=TRI_NODE_NUM
      end forall
      forall(i=1:nQuad)
        Facet(nTri+i)%ShapeType=QUAD_TYPE
        Facet(nTri+i)%ShapeInd=i
        Facet(nTri+i)%NodeNum=QUAD_NODE_NUM
      end forall
      nEle=nTet+nHex
      allocate(Ele(nEle))
      forall(i=1:nTet)
        Ele(i)%ShapeType=TET_TYPE
        Ele(i)%ShapeInd=i
        Ele(i)%NodeNum=TET_NODE_NUM
        Ele(i)%SurfNum=TET_SURF_NUM
      end forall
      forall(i=1:nHex)
        Ele(nTet+i)%ShapeType=HEX_TYPE
        Ele(nTet+i)%ShapeInd=i
        Ele(nTet+i)%NodeNum=HEX_NODE_NUM
        Ele(nTet+i)%SurfNum=HEX_SURF_NUM
      end forall
      cycle
    end if
  end do
  
  ! auxiliary operations
  nPrt=1
  if(nTet>0)then
    nPrt=max(nPrt,maxval(Tet(:)%Prt))
  end if
  if(nHex>0)then
    nPrt=max(nPrt,maxval(Hex(:)%Prt))
  end if
  
  write(*,'(a)'),'grid summary:'
  write(*,'(a,i8)'),'  number of nodes:        ',nNode
  write(*,'(a,i8)'),'  number of points:       ',nPoint
  write(*,'(a,i8)'),'  number of lines:        ',nLine
  write(*,'(a,i8)'),'  number of tri.:         ',nTri
  write(*,'(a,i8)'),'  number of quad.:        ',nQuad
  write(*,'(a,i8)'),'  number of tet.:         ',nTet
  write(*,'(a,i8)'),'  number of hex.:         ',nHex
  write(*,'(a,i8,/)'),'  number of partitions.:  ',nPrt
  
  close(gridfile)
end subroutine

!****************
! read data file
!****************
subroutine readdata(fname,datafile)
  use moduleMiscDataStruct
  
  integer datafile,readerr
  character(100) temp_string,fname
  
  open(datafile,file=fname,status='old')
  readerr=0
  
  l=0
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(datafile,'(a100)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    
    ! Note: the data is counted but not recorded during the 1st reading
    ! read 1-dimensional table
    if(temp_string(1:6)=='$Tab1D')then
      do while(readerr==0)
        read(datafile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        l=l+1
      end do
    end if
  end do
  
  rewind(datafile)
  readerr=0
  allocate(dataTab1d(l/3))
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(datafile,'(a100)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    
    i=1
    ! read 1-dimensional table
    if(temp_string(1:6)=='$Tab1D')then
      do while(readerr==0)
        read(datafile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),dataTab1d(i)%length
        allocate(dataTab1d(i)%x(dataTab1d(i)%length))
        allocate(dataTab1d(i)%y(dataTab1d(i)%length))
        read(datafile,*,iostat=readerr),(dataTab1d(i)%x(k),k=1,dataTab1d(i)%length)
        read(datafile,*,iostat=readerr),(dataTab1d(i)%y(k),k=1,dataTab1d(i)%length)
        i=i+1
        
      end do
    end if
  end do
  
  close(datafile)
end subroutine

!****************************
! read simulation conditions
!****************************
subroutine readcod(fname,codfile)
  use moduleGrid
  use moduleWrite
  use moduleCond
  
  character(100),intent(in)::fname
  integer,intent(in)::codfile
  integer readerr,GeoEntity,nVal,nTab,iCond
  character(200) temp_string
  logical isNodeDone(nNode),isFacetDone(nFacet),isEleDone(nEle)
  type(typeCond)::tempCond
  
  allocate(CondNode(nNode))
  allocate(CondFacet(nFacet))
  allocate(CondEle(nEle))
  
  open(codfile,file=fname,status='old')
  readerr=0
  
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
    ! read node conditions
    if(temp_string(1:9)=='$NodeCond')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),tempCond%what,GeoEntity,nVal,nTab
        if(nVal>0)then
          allocate(tempCond%Val(nVal))
        end if
        if(nTab>0)then
          allocate(tempCond%Tab(nTab))
        end if
        read(temp_string,*),tempCond%what,GeoEntity,nVal,nTab,&
        &                   (tempCond%Val(i),i=1,nVal),(tempCond%Tab(i),i=1,nTab)
        isNodeDone(:)=.false.
        do i=1,nPoint
          if(Point(i)%GeoEnti==GeoEntity.and.(.not.isNodeDone(Point(i)%NodeInd)))then
            iCond=CondNode(Point(i)%NodeInd)%getSpace()
            CondNode(Point(i)%NodeInd)%Cond(iCond)=tempCond
            isNodeDone(Point(i)%NodeInd)=.true.
          end if
        end do
        do i=1,nFacet
          do j=1,Facet(i)%NodeNum
            if(Facet(i)%GeoEnti==GeoEntity.and.(.not.isNodeDone(Facet(i)%NodeInd(j))))then
              iCond=CondNode(Facet(i)%NodeInd(j))%getSpace()
              CondNode(Facet(i)%NodeInd(j))%Cond(iCond)=tempCond
              isNodeDone(Facet(i)%NodeInd(j))=.true.
            end if
          end do
        end do
        do i=1,nEle
          do j=1,Ele(i)%NodeNum
            if(Ele(i)%GeoEnti==GeoEntity.and.(.not.isNodeDone(Ele(i)%NodeInd(j))))then
              iCond=CondNode(Ele(i)%NodeInd(j))%getSpace()
              CondNode(Ele(i)%NodeInd(j))%Cond(iCond)=tempCond
              isNodeDone(Ele(i)%NodeInd(j))=.true.
            end if
          end do
        end do
        ! clean ups
        if(allocated(tempCond%Val))then
          deallocate(tempCond%Val)
        end if
        if(allocated(tempCond%Tab))then
          deallocate(tempCond%Tab)
        end if
      end do
    end if
    ! read facet conditions
    if(temp_string(1:10)=='$FacetCond')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),tempCond%what,GeoEntity,nVal,nTab
        if(nVal>0)then
          allocate(tempCond%Val(nVal))
        end if
        if(nTab>0)then
          allocate(tempCond%Tab(nTab))
        end if
        read(temp_string,*),tempCond%what,GeoEntity,nVal,nTab,&
        &                   (tempCond%Val(i),i=1,nVal),(tempCond%Tab(i),i=1,nTab)
        isFacetDone(:)=.false.
        do i=1,nFacet
          if(Facet(i)%GeoEnti==GeoEntity.and.(.not.isFacetDone(i)))then
            iCond=CondFacet(i)%getSpace()
            CondFacet(i)%Cond(iCond)=tempCond
            isFacetDone(i)=.true.
          end if
        end do
        ! clean ups
        if(allocated(tempCond%Val))then
          deallocate(tempCond%Val)
        end if
        if(allocated(tempCond%Tab))then
          deallocate(tempCond%Tab)
        end if
      end do
    end if
    ! read element conditions
    if(temp_string(1:8)=='$EleCond')then
      do while(readerr==0)
        read(codfile,'(a200)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),tempCond%what,GeoEntity,nVal,nTab
        if(nVal>0)then
          allocate(tempCond%Val(nVal))
        end if
        if(nTab>0)then
          allocate(tempCond%Tab(nTab))
        end if
        read(temp_string,*),tempCond%what,GeoEntity,nVal,nTab,&
        &                   (tempCond%Val(i),i=1,nVal),(tempCond%Tab(i),i=1,nTab)
        isEleDone(:)=.false.
        do i=1,nEle
          if(Ele(i)%GeoEnti==GeoEntity.and.(.not.isEleDone(i)))then
            iCond=CondEle(i)%getSpace()
            CondEle(i)%Cond(iCond)=tempCond
            isEleDone(i)=.true.
          end if
        end do
        ! clean ups
        if(allocated(tempCond%Val))then
          deallocate(tempCond%Val)
        end if
        if(allocated(tempCond%Tab))then
          deallocate(tempCond%Tab)
        end if
      end do
    end if
  end do
  
  close(codfile)
end subroutine

!************************************
! initialize data output environment
!************************************
! a: number of scaler data sets at nodes
! b: number of vector data sets at nodes
! c: number of tensor data sets at nodes
! d: number of scaler data sets at the center of elements
! e: number of vector data sets at the center of elements
! f: number of tensor data sets at the center of elements
subroutine initWriteEnv(a,b,c,d,e,f)
  use moduleGrid
  use moduleWrite
  integer,intent(in)::a,b,c,d,e,f
  ! allocate output data space
  nrstNodeScal=a
  nrstNodeVect=b
  nrstNodeTens=c
  nrstEleScal=d
  nrstEleVect=e
  nrstEleTens=f
  allocate(rstNodeScal(a))
  allocate(rstNodeVect(b))
  allocate(rstNodeTens(c))
  allocate(rstEleScal(d))
  allocate(rstEleVect(e))
  allocate(rstEleTens(f))
  ! set output control variables
  nWrite=0
end subroutine

!***************
! write results
!***************
! hold=.true.: hold the file when exiting subroutine; good for writing a time span
! hold=.false.: close the file when exiting subroutine; good for writing a snapshot
subroutine writerst(fname,ifile,hold)
  use moduleGrid
  use moduleWrite
  integer ifile
  character(100) fname
  character(400) tempstring
  logical hold
  
  if(nWrite==0.or..not.hold)then ! if it is the 1st time writing results
    open(ifile,file=fname,status='replace')
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
      write(tempstring,*),&
      &    i,Node(i)%Pos(:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndNodes'
    ! write elements
    write(ifile,'(a)'),'$Elements'
    write(tempstring,*),nEle+nPoint+nLine+nFacet
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
    ! Note: write solid elements first so that the element index would be consistant
    do i=1,nEle
      write(tempstring,*),&
      &    i,Ele(i)%ShapeType,2,0,Ele(i)%GeoEnti,(Ele(i)%NodeInd(j),j=1,Ele(i)%NodeNum)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    do i=1,nPoint
      write(tempstring,*),&
      &    nEle+i,POINT_TYPE,2,0,Point(i)%GeoEnti,Point(i)%NodeInd
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    do i=1,nLine
      write(tempstring,*),&
      &    nEle+nPoint+i,LINE_TYPE,2,0,Line(i)%GeoEnti,Line(i)%NodeInd(:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    do i=1,nFacet
      write(tempstring,*),&
      &    nEle+nPoint+nLine+i,Facet(i)%ShapeType,2,0,Facet(i)%GeoEnti,&
      &    (Facet(i)%NodeInd(j),j=1,Facet(i)%NodeNum)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndElements'
  end if
  
  ! write scalers at nodes
  do i=1,nrstNodeScal
    write(ifile,'(a)'),'$NodeData'
    write(ifile,'(i1)'),1 ! 1 string tag
    write(ifile,'(a,i2,a)'),'"scalNode',i,'"' ! name of the data-set
    write(ifile,'(i1)'),1 ! 1 real tag
    write(ifile,*),t ! time
    write(ifile,'(i1)'),3 ! 3 integer tags
    write(ifile,*),merge(nWrite,1,hold) ! time step index
    write(ifile,'(i1)'),1 ! 1 component scaler
    write(ifile,*),nNode
    do j=1,nNode
      write(tempstring,*),j,rstNodeScal(i)%ptr(j)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndNodeData'
  end do
  ! write vectors at nodes
  do i=1,nrstNodeVect
    write(ifile,'(a)'),'$NodeData'
    write(ifile,'(i1)'),1 ! 1 string tag
    write(ifile,'(a,i2,a)'),'"vectNode',i,'"' ! name of the data-set
    write(ifile,'(i1)'),1 ! 1 real tag
    write(ifile,*),t ! time
    write(ifile,'(i1)'),3 ! 3 integer tags
    write(ifile,*),merge(nWrite,1,hold) ! time step index
    write(ifile,'(i1)'),3 ! 3 components vector
    write(ifile,*),nNode
    do j=1,nNode
      write(tempstring,*),j,rstNodeVect(i)%ptr(j,:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndNodeData'
  end do
  ! write tensors at nodes
  do i=1,nrstNodeTens
    write(ifile,'(a)'),'$NodeData'
    write(ifile,'(i1)'),1 ! 1 string tag
    write(ifile,'(a,i2,a)'),'"tensNode',i,'"' ! name of the data-set
    write(ifile,'(i1)'),1 ! 1 real tag
    write(ifile,*),t ! time
    write(ifile,'(i1)'),3 ! 3 integer tags
    write(ifile,*),merge(nWrite,1,hold) ! time step index
    write(ifile,'(i1)'),9 ! 9 components tensor
    write(ifile,*),nNode
    do j=1,nNode
      write(tempstring,*),j,rstNodeTens(i)%ptr(j,1,:),rstNodeTens(i)%ptr(j,2,:),&
      &                     rstNodeTens(i)%ptr(j,3,:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndNodeData'
  end do
  ! write scalers at the center of elemnets
  do i=1,nrstEleScal
    write(ifile,'(a)'),'$ElementData'
    write(ifile,'(i1)'),1 ! 1 string tag
    write(ifile,'(a,i2,a)'),'"scalEle',i,'"' ! name of the data-set
    write(ifile,'(i1)'),1 ! 1 real tag
    write(ifile,*),t ! time
    write(ifile,'(i1)'),3 ! 3 integer tags
    write(ifile,*),merge(nWrite,1,hold) ! time step index
    write(ifile,'(i1)'),1 ! 1 component scaler
    write(ifile,*),nEle
    do j=1,nEle
      write(tempstring,*),j,rstEleScal(i)%ptr(j)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndElementData'
  end do
  ! write vectors at the center of elements
  do i=1,nrstEleVect
    write(ifile,'(a)'),'$ElementData'
    write(ifile,'(i1)'),1 ! 1 string tag
    write(ifile,'(a,i2,a)'),'"vectEle',i,'"' ! name of the data-set
    write(ifile,'(i1)'),1 ! 1 real tag
    write(ifile,*),t ! time
    write(ifile,'(i1)'),3 ! 3 integer tags
    write(ifile,*),merge(nWrite,1,hold) ! time step index
    write(ifile,'(i1)'),3 ! 3 components vector
    write(ifile,*),nEle
    do j=1,nEle
      write(tempstring,*),j,rstEleVect(i)%ptr(j,:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndElementData'
  end do
  ! write tensors at the center of elements
  do i=1,nrstEleTens
    write(ifile,'(a)'),'$ElementData'
    write(ifile,'(i1)'),1 ! 1 string tag
    write(ifile,'(a,i2,a)'),'"tensEle',i,'"' ! name of the data-set
    write(ifile,'(i1)'),1 ! 1 real tag
    write(ifile,*),t ! time
    write(ifile,'(i1)'),3 ! 3 integer tags
    write(ifile,*),merge(nWrite,1,hold) ! time step index
    write(ifile,'(i1)'),9 ! 9 components tensor
    write(ifile,*),nEle
    do j=1,nEle
      write(tempstring,*),j,rstEleTens(i)%ptr(j,1,:),rstEleTens(i)%ptr(j,2,:),&
      &                     rstEleTens(i)%ptr(j,3,:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndElementData'
  end do
  
  if(hold)then
    nWrite=nWrite+1
  end if
  
  if(t>tFinal.or..not.hold)then ! if it is the last time writing results
    close(ifile)
  end if
end subroutine
