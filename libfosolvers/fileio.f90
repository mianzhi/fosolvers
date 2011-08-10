!----------------------------------------------------------------------------- best with 100 columns

!****************
! read mesh file
!****************
subroutine readmsh()
  use moduleGrid
  
  integer gridfile,readerr,np,nt
  character*100 temp_string
  integer temp_int_vect1(20)
  type(typePoint),allocatable::tempPoint(:)
  type(typeLine),allocatable::tempLine(:)
  type(typeTri),allocatable::tempTri(:)
  type(typeQuad),allocatable::tempQuad(:)
  type(typeTet),allocatable::tempTet(:)
  type(typeHex),allocatable::tempHex(:)
  
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
