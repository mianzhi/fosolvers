!----------------------------------------------------------------------------- best with 100 columns

!****************
! read mesh file
!****************
subroutine readmsh(fname,gridfile)
  use moduleGrid
  
  integer gridfile,readerr,np,nt
  character*100 temp_string,fname
  integer temp_int_vect1(20)
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
          case(15) ! point
            np=1
          case(1) ! line
            np=2
          case(2) ! tri
            np=3
          case(3) ! quad
            np=4
          case(4) ! tet
            np=4
          case(5) ! hex
            np=8
        end select
        read(temp_string,*),temp_int_vect1(1:3+nt+np)
        select case(temp_int_vect1(2)) ! type of element
          case(15) ! point
            nPoint=nPoint+1
            Point(nPoint)%NodeInd=temp_int_vect1(3+nt+1)
            Point(nPoint)%GeoEnti=temp_int_vect1(3+nt)
          case(1) ! line
            nLine=nLine+1
            Line(nLine)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            Line(nLine)%GeoEnti=temp_int_vect1(3+nt)
          case(2) ! tri
            nTri=nTri+1
            Tri(nTri)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            Tri(nTri)%GeoEnti=temp_int_vect1(3+nt)
          case(3) ! quad
            nQuad=nQuad+1
            Quad(nQuad)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            Quad(nQuad)%GeoEnti=temp_int_vect1(3+nt)
          case(4) ! tet
            nTet=nTet+1
            Tet(nTet)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            Tet(nTet)%GeoEnti=temp_int_vect1(3+nt)
          case(5) ! hex
            nHex=nHex+1
            Hex(nHex)%NodeInd(:)=temp_int_vect1(3+nt+1:3+nt+np)
            Hex(nHex)%GeoEnti=temp_int_vect1(3+nt)
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
        Facet(i)%ShapeType=2
        Facet(i)%ShapeInd=i
        Facet(i)%NodeNum=3
      end forall
      forall(i=1:nQuad)
        Facet(nTri+i)%ShapeType=3
        Facet(nTri+i)%ShapeInd=i
        Facet(nTri+i)%NodeNum=4
      end forall
      nEle=nTet+nHex
      allocate(Ele(nEle))
      forall(i=1:nTet)
        Ele(i)%ShapeType=4
        Ele(i)%ShapeInd=i
        Ele(i)%NodeNum=4
        Ele(i)%SurfNum=4
      end forall
      forall(i=1:nHex)
        Ele(nTet+i)%ShapeType=5
        Ele(nTet+i)%ShapeInd=i
        Ele(nTet+i)%NodeNum=8
        Ele(nTet+i)%SurfNum=6
      end forall
      cycle
    end if
  end do
  
  ! auxiliary operations
  SurfTabTet(1,:)=[1,3,2] ! table of surface nodes for all kinds of elements
  SurfTabTet(2,:)=[1,2,4]
  SurfTabTet(3,:)=[1,4,3]
  SurfTabTet(4,:)=[2,3,4]
  SurfTabHex(1,:)=[2,3,7,6]
  SurfTabHex(2,:)=[1,5,8,4]
  SurfTabHex(3,:)=[3,4,8,7]
  SurfTabHex(4,:)=[1,2,6,5]
  SurfTabHex(5,:)=[5,6,7,8]
  SurfTabHex(6,:)=[1,4,3,2]
  
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

!*********************************
! result output related variables
!*********************************
module moduleWrite
  ! output control related variables
  double precision,save::t,tFinal
  integer,save::nWrite
  ! data to be output
  double precision,allocatable,save::rstNodeScal(:,:),rstNodeVect(:,:,:),rstNodeTens(:,:,:),&
  &                                  rstEleScal(:,:),rstEleVect(:,:,:),rstEleTens(:,:,:)
  integer,save::nrstNodeScal,nrstNodeVect,nrstNodeTens,nrstEleScal,nrstEleVect,nrstEleTens
end module

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
  allocate(rstNodeScal(a,nNode))
  allocate(rstNodeVect(a,nNode,3))
  allocate(rstNodeTens(a,nNode,9))
  allocate(rstEleScal(a,nEle))
  allocate(rstEleVect(a,nEle,3))
  allocate(rstEleTens(a,nEle,9))
  ! set output control variables
  nWrite=0
end subroutine

!**********************************
! write results during a time span
!**********************************
subroutine writerstSpan(fname,ifile)
  use moduleGrid
  use moduleWrite
  integer ifile
  character*400 tempstring,fname
  
  if(nWrite==0)then ! if it is the 1st time writing results
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
      &    i,Node(i)%Pos(1),Node(i)%Pos(2),Node(i)%Pos(3)
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
      &    nEle+i,15,2,0,Point(i)%GeoEnti,Point(i)%NodeInd
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    do i=1,nLine
      write(tempstring,*),&
      &    nEle+nPoint+i,1,2,0,Line(i)%GeoEnti,Line(i)%NodeInd(:)
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
    write(ifile,*),nWrite ! time step index
    write(ifile,'(i1)'),1 ! 1 component scaler
    write(ifile,*),nNode
    do j=1,nNode
      write(tempstring,*),j,rstNodeScal(i,j)
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
    write(ifile,*),nWrite ! time step index
    write(ifile,'(i1)'),3 ! 3 components vector
    write(ifile,*),nNode
    do j=1,nNode
      write(tempstring,*),j,rstNodeVect(i,j,:)
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
    write(ifile,*),nWrite ! time step index
    write(ifile,'(i1)'),9 ! 9 components tensor
    write(ifile,*),nNode
    do j=1,nNode
      write(tempstring,*),j,rstNodeTens(i,j,:)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end do
    write(ifile,'(a)'),'$EndNodeData'
  end do
  
  
  nWrite=nWrite+1
  
  if(t>tFinal)then ! if it is the last time writing results
    close(ifile)
  end if
end subroutine
