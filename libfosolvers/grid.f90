!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! data and elementary procedures for grid
!*****************************************
module moduleGrid
  private
  
  ! make some functions properly accessible (both privately & publicly)
  interface
    function find3PNorm(P1,P2,P3)
      double precision,intent(in)::P1(3),P2(3),P3(3)
      double precision find3PNorm(3)
    end function
  end interface
  public::find3PNorm
  
  ! table of surface nodes for all kinds of elements
  integer,public,save::SurfTabTet(4,3),SurfTabHex(6,4)
  
  ! bounding box of the geometry
  double precision,public,save::BoundBox(2,3)
  
  !----------
  ! typeNode
  !----------
  type,public::typeNode
    double precision Pos(3)
  end type
  type(typeNode),public,allocatable,save::Node(:)
  integer,public,save::nNode
  
  !-----------
  ! typePoint
  !-----------
  type,public::typePoint
    integer NodeInd
    integer GeoEnti
  contains
    procedure,public::findPos=>findPointPos
  end type
  type(typePoint),public,allocatable,save::Point(:)
  integer,public,save::nPoint
  
  !----------
  ! typeLine
  !----------
  type,public::typeLine
    integer NodeInd(2)
    integer GeoEnti
  contains
    procedure,public::findPC=>findLinePC
    procedure,public::findLength=>findLineLength
  end type
  type(typeLine),public,allocatable,save::Line(:)
  integer,public,save::nLine
  
  !---------
  ! typeTri
  !---------
  type,public::typeTri
    integer NodeInd(3)
    integer GeoEnti
  contains
    procedure,public::findPC=>findTriPC
    procedure,public::findArea=>findTriArea
    procedure,public::findNorm=>findTriNorm
    procedure,public::getNeibEle=>getTriNeibEle
  end type
  type(typeTri),public,allocatable,save::Tri(:)
  integer,public,save::nTri
  
  !----------
  ! typeQuad
  !----------
  type,public::typeQuad
    integer NodeInd(4)
    integer GeoEnti
  contains
    procedure,public::findPC=>findQuadPC
    procedure,public::findArea=>findQuadArea
    procedure,public::findNorm=>findQuadNorm
    procedure,public::getNeibEle=>getQuadNeibEle
  end type
  type(typeQuad),public,allocatable,save::Quad(:)
  integer,public,save::nQuad
  
  !---------
  ! typeTet
  !---------
  type,public::typeTet
    integer NodeInd(4)
    integer GeoEnti
    integer Prt
  contains
    procedure,public::findPC=>findTetPC
    procedure,public::findVol=>findTetVol
    procedure,public::getNeib=>getTetNeib
    procedure,public::findSurfPC=>findTetSurfPC
    procedure,public::findSurfArea=>findTetSurfArea
    procedure,public::findSurfNorm=>findTetSurfNorm
  end type
  type(typeTet),public,allocatable,save::Tet(:)
  integer,public,save::nTet
  
  !---------
  ! typeHex
  !---------
  type,public::typeHex
    integer NodeInd(8)
    integer GeoEnti
    integer Prt
  contains
    procedure,public::findPC=>findHexPC
    procedure,public::findVol=>findHexVol
    procedure,public::getNeib=>getHexNeib
    procedure,public::findSurfPC=>findHexSurfPC
    procedure,public::findSurfArea=>findHexSurfArea
    procedure,public::findSurfNorm=>findHexSurfNorm
  end type
  type(typeHex),public,allocatable,save::Hex(:)
  integer,public,save::nHex
  
  !-----------
  ! typeFacet
  !-----------
  type,public::typeFacet
    integer ShapeType
    ! type of shape:
    ! 2: 3-node triangle
    ! 3: 4-node quadrilateral
    integer ShapeInd
    integer NodeNum
    integer NodeInd(15) ! 15 is the maximum possible number of nodes a facet can have
    integer GeoEnti
    integer NeibEle(2) ! a facet may have 2 neighbour elements
    double precision PC(3)
    double precision Area
    double precision Norm(3)
  contains
    procedure,public::getNodeInd=>getFacetNodeInd
    procedure,public::findPC=>findFacetPC
    procedure,public::findArea=>findFacetArea
    procedure,public::findNorm=>findFacetNorm
    procedure,public::getGeoEnti=>getFacetGeoEnti
    procedure,public::getNeibEle=>getFacetNeibEle
  end type
  type(typeFacet),public,allocatable,save::Facet(:)
  integer,public,save::nFacet
  
  !---------
  ! typeEle
  !---------
  type,public::typeEle
    integer ShapeType
    ! type of shape:
    ! 4: 4-node tetrahedron
    ! 5: 8-node hexahedron
    integer ShapeInd
    integer NodeNum
    integer SurfNum
    integer NodeInd(27) ! 27 is the maximum possible number of nodes an element can have
    integer GeoEnti
    integer Prt
    integer Neib(6) ! 6 is the maximum possible number of neighbours
    double precision PC(3)
    double precision Vol
    double precision SurfPC(6,3)
    double precision SurfArea(6)
    double precision SurfNorm(6,3)
  contains
    procedure,public::getNodeInd=>getEleNodeInd
    procedure,public::findPC=>findElePC
    procedure,public::findVol=>findEleVol
    procedure,public::getGeoEnti=>getEleGeoEnti
    procedure,public::getPrt=>getElePrt
    procedure,public::getNeib=>getEleNeib
    procedure,public::findSurfPC=>findEleSurfPC
    procedure,public::findSurfArea=>findEleSurfArea
    procedure,public::findSurfNorm=>findEleSurfNorm
  end type
  type(typeEle),public,allocatable,save::Ele(:)
  integer,public,save::nEle
  
contains
  
  !---------------------------------
  ! find the position of this point
  !---------------------------------
  function findPointPos(this)
    class(typePoint),intent(in)::this
    double precision findPointPos(3)
    findPointPos(:)=Node(this%NodeInd)%Pos(:)
  end function
  
  !------------------------------
  ! find the center of this line
  !------------------------------
  function findLinePC(this)
    class(typeLine),intent(in)::this
    double precision findLinePC(3)
    findLinePC(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:))/2d0
  end function
  
  !------------------------------
  ! find the length of this line
  !------------------------------
  function findLineLength(this)
    class(typeLine),intent(in)::this
    double precision findLineLength,vect(3)
    vect(:)=Node(this%NodeInd(2))%Pos(:)-Node(this%NodeInd(1))%Pos(:)
    findLineLength=norm2(vect)
  end function
  
  !----------------------------------
  ! find the center of this triangle
  !----------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findTriPC(this)
    class(typeTri),intent(in)::this
    double precision findTriPC(3)
    findTriPC(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &            +Node(this%NodeInd(3))%Pos(:))/3d0
  end function
  
  !--------------------------------
  ! find the area of this triangle
  !--------------------------------
  function findTriArea(this)
    class(typeTri),intent(in)::this
    double precision findTriArea,find3PArea
    findTriArea=find3PArea(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                      Node(this%NodeInd(3))%Pos(:))
  end function
  
  !-----------------------------------------
  ! find the normal vector of this triangle
  !-----------------------------------------
  function findTriNorm(this)
    class(typeTri),intent(in)::this
    double precision findTriNorm(3)
    findTriNorm(:)=find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                         Node(this%NodeInd(3))%Pos(:))
  end function
  
  !--------------------------------------------
  ! get the neighbour element of this triangle
  !--------------------------------------------
  ! Note: one triangle may have 2 neighbour elements
  function getTriNeibEle(this)
    class(typeTri),intent(in)::this
    integer getTriNeibEle(2),lPTet(4)
    logical maskTri(3)
    getTriNeibEle(:)=0
    k=0
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(4) ! the neigbour element can be a tetrahedron
          maskTri(:)=.false.
          lPTet(:)=Tet(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,3
            if(any(lPTet(:)==this%NodeInd(j)))then
              maskTri(j)=.true.
            end if
          end do
          if(all(maskTri(:)))then
            k=k+1
            getTriNeibEle(k)=i
            if(k==2)then
              exit
            end if
          end if
        case default
      end select
    end do
  end function
  
  !---------------------------------------
  ! find the center of this quadrilateral
  !---------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findQuadPC(this)
    class(typeQuad),intent(in)::this
    double precision findQuadPC(3)
    findQuadPC(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &             +Node(this%NodeInd(3))%Pos(:)+Node(this%NodeInd(4))%Pos(:))/4d0
  end function
  
  !-------------------------------------
  ! find the area of this quadrilateral
  !-------------------------------------
  function findQuadArea(this)
    class(typeQuad),intent(in)::this
    double precision findQuadArea,find3PArea
    findQuadArea=find3PArea(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                       Node(this%NodeInd(3))%Pos(:))&
    &           +find3PArea(Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:),&
    &                       Node(this%NodeInd(1))%Pos(:))
  end function
  
  !----------------------------------------------
  ! find the normal vector of this quadrilateral
  !----------------------------------------------
  function findQuadNorm(this)
    class(typeQuad),intent(in)::this
    double precision findQuadNorm(3)
    findQuadNorm(:)=find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                          Node(this%NodeInd(3))%Pos(:))
  end function
  
  !-------------------------------------------------
  ! get the neighbour element of this quadrilateral
  !-------------------------------------------------
  ! Note: one quadirlateral may have 2 neighbour elements
  function getQuadNeibEle(this)
    class(typeQuad),intent(in)::this
    integer getQuadNeibEle(2),lPHex(8)
    logical maskQuad(4)
    getQuadNeibEle(:)=0
    k=0
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(5) ! the neigbour element can be a hexahedron
          maskQuad(:)=.false.
          lPHex(:)=Hex(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,4
            if(any(lPHex(:)==this%NodeInd(j)))then
              maskQuad(j)=.true.
            end if
          end do
          if(all(maskQuad(:)))then
            k=k+1
            getQuadNeibEle(k)=i
            if(k==2)then
              exit
            end if
          end if
        case default
      end select
    end do
  end function
  
  !-------------------------------------
  ! find the center of this tetrahedron
  !-------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findTetPC(this)
    class(typeTet),intent(in)::this
    double precision findTetPC(3)
    findTetPC(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &            +Node(this%NodeInd(3))%Pos(:)+Node(this%NodeInd(4))%Pos(:))/4d0
  end function
  
  !-------------------------------------
  ! find the volume of this tetrahedron
  !-------------------------------------
  function findTetVol(this)
    class(typeTet),intent(in)::this
    double precision findTetVol,find4PVol
    findTetVol=find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                    Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:))
  end function
  
  !-------------------------------------------------------------
  ! get the k_th neighbour element or facet of this tetrahedron
  !-------------------------------------------------------------
  ! Note: if facet are found, the function will return the negative value of the facet index.
  ! Note: we are getting neighbour element or facet, not neighbour hexahedron or quadrilateral
  !       index, though the neighbour element (facet) is likely to be a tetrahedron (triangle).
  function getTetNeib(this,k)
    class(typeTet),intent(in)::this
    integer,intent(in)::k
    integer getTetNeib,lPTet(4),lPTri(3)
    logical maskTet(4)
    getTetNeib=0
    if(k<1.or.k>4)then
      write(*,'(a,i2,a)'),'ERROR: a tetrahedron can not have ',k,'th surface'
      stop
    end if
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(4) ! the neigbour element can be a tetrahedron
          maskTet(:)=.false.
          lPTet(:)=Tet(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,4
            if(any(lPTet(:)==this%NodeInd(j)))then
              maskTet(j)=.true.
            end if
          end do
          if(all(maskTet(SurfTabTet(k,:))).and.any(maskTet(:).eqv..false.))then
            getTetNeib=i
            exit
          end if
        case default
      end select
    end do
    do i=1,nFacet
      select case(Facet(i)%ShapeType)
        case(2) ! the neigbour facet can be a triangle
          maskTet(:)=.false.
          lPTri(:)=Tri(Facet(i)%ShapeInd)%NodeInd(:)
          do j=1,4
            if(any(lPTri(:)==this%NodeInd(j)))then
              maskTet(j)=.true.
            end if
          end do
          if(all(maskTet(SurfTabTet(k,:))))then
            if(getTetNeib==0.or.Tri(Facet(i)%ShapeInd)%GeoEnti<0)then
              getTetNeib=-i
            end if
            ! Note: different partitions are splited, while different geometric entities are not
            exit
          end if
        case default
      end select
    end do
  end function
  
  !----------------------------------------------------
  ! find the center of this tetrahedron's k_th surface
  !----------------------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findTetSurfPC(this,k)
    class(typeTet),intent(in)::this
    integer,intent(in)::k
    double precision findTetSurfPC(3)
    if(k<1.or.k>4)then
      write(*,'(a,i2,a)'),'ERROR: a tetrahedron can not have ',k,'th surface'
      stop
    end if
    findTetSurfPC(:)=(Node(this%NodeInd(SurfTabTet(k,1)))%Pos(:)&
    &                +Node(this%NodeInd(SurfTabTet(k,2)))%Pos(:)&
    &                +Node(this%NodeInd(SurfTabTet(k,3)))%Pos(:))/3d0
  end function
  
  !--------------------------------------------------
  ! find the area of this tetrahedron's k_th surface
  !--------------------------------------------------
  function findTetSurfArea(this,k)
    class(typeTet),intent(in)::this
    integer,intent(in)::k
    double precision findTetSurfArea
    if(k<1.or.k>4)then
      write(*,'(a,i2,a)'),'ERROR: a tetrahedron can not have ',k,'th surface'
      stop
    end if
    findTetSurfArea=find3PArea(Node(this%NodeInd(SurfTabTet(k,1)))%Pos(:),&
    &                          Node(this%NodeInd(SurfTabTet(k,2)))%Pos(:),&
    &                          Node(this%NodeInd(SurfTabTet(k,3)))%Pos(:))
  end function
  
  !-----------------------------------------------------------
  ! find the normal vector of this tetrahedron's k_th surface
  !-----------------------------------------------------------
  function findTetSurfNorm(this,k)
    class(typeTet),intent(in)::this
    integer,intent(in)::k
    double precision findTetSurfNorm(3)
    if(k<1.or.k>4)then
      write(*,'(a,i2,a)'),'ERROR: a tetrahedron can not have ',k,'th surface'
      stop
    end if
    findTetSurfNorm(:)=find3PNorm(Node(this%NodeInd(SurfTabTet(k,1)))%Pos(:),&
    &                             Node(this%NodeInd(SurfTabTet(k,2)))%Pos(:),&
    &                             Node(this%NodeInd(SurfTabTet(k,3)))%Pos(:))
  end function
  
  !------------------------------------
  ! find the center of this hexahedron
  !------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findHexPC(this)
    class(typeHex),intent(in)::this
    double precision findHexPC(3)
    findHexPC(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &            +Node(this%NodeInd(3))%Pos(:)+Node(this%NodeInd(4))%Pos(:)&
    &            +Node(this%NodeInd(5))%Pos(:)+Node(this%NodeInd(6))%Pos(:)&
    &            +Node(this%NodeInd(7))%Pos(:)+Node(this%NodeInd(8))%Pos(:))/8d0
  end function
  
  !------------------------------------
  ! find the volume of this hexahedron
  !------------------------------------
  function findHexVol(this)
    class(typeHex),intent(in)::this
    double precision findHexVol,find4PVol
    findHexVol=find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                    Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(6))%Pos(:))&
    &         +find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(3))%Pos(:),&
    &                    Node(this%NodeInd(4))%Pos(:),Node(this%NodeInd(6))%Pos(:))&
    &         +find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(4))%Pos(:),&
    &                    Node(this%NodeInd(5))%Pos(:),Node(this%NodeInd(6))%Pos(:))&
    &         +find4PVol(Node(this%NodeInd(4))%Pos(:),Node(this%NodeInd(5))%Pos(:),&
    &                    Node(this%NodeInd(6))%Pos(:),Node(this%NodeInd(7))%Pos(:))&
    &         +find4PVol(Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:),&
    &                    Node(this%NodeInd(6))%Pos(:),Node(this%NodeInd(7))%Pos(:))&
    &         +find4PVol(Node(this%NodeInd(4))%Pos(:),Node(this%NodeInd(5))%Pos(:),&
    &                    Node(this%NodeInd(7))%Pos(:),Node(this%NodeInd(8))%Pos(:))
  end function
  
  !------------------------------------------------------------
  ! get the k_th neighbour element or facet of this hexahedron
  !------------------------------------------------------------
  ! Note: if facet are found, the function will return the negative value of the facet index.
  ! Note: we are getting neighbour element or facet, not neighbour hexahedron or quadrilateral
  !       index, though the neighbour element (facet) is likely to be a hexahedron (quadrilateral).
  function getHexNeib(this,k)
    class(typeHex),intent(in)::this
    integer,intent(in)::k
    integer getHexNeib,lPHex(8),lPQuad(4)
    logical maskHex(8)
    getHexNeib=0
    if(k<1.or.k>6)then
      write(*,'(a,i2,a)'),'ERROR: a hexahedron can not have ',k,'th surface'
      stop
    end if
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(5) ! the neigbour element can be a hexahedron
          maskHex(:)=.false.
          lPHex(:)=Hex(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,8
            if(any(lPHex(:)==this%NodeInd(j)))then
              maskHex(j)=.true.
            end if
          end do
          if(all(maskHex(SurfTabHex(k,:))).and.any(maskHex(:).eqv..false.))then
            getHexNeib=i
            exit
          end if
        case default
      end select
    end do
    do i=1,nFacet
      select case(Facet(i)%ShapeType)
        case(3) ! the neigbour facet can be a quadrilateral
          maskHex(:)=.false.
          lPQuad(:)=Quad(Facet(i)%ShapeInd)%NodeInd(:)
          do j=1,8
            if(any(lPQuad(:)==this%NodeInd(j)))then
              maskHex(j)=.true.
            end if
          end do
          if(all(maskHex(SurfTabHex(k,:))))then
            if(getHexNeib==0.or.Quad(Facet(i)%ShapeInd)%GeoEnti<0)then
              getHexNeib=-i
            end if
            ! Note: different partitions are splited, while different geometric entities are not
            exit
          end if
        case default
      end select
    end do
  end function
  
  !---------------------------------------------------
  ! find the center of this hexahedron's k_th surface
  !---------------------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findHexSurfPC(this,k)
    class(typeHex),intent(in)::this
    integer,intent(in)::k
    double precision findHexSurfPC(3)
    if(k<1.or.k>6)then
      write(*,'(a,i2,a)'),'ERROR: a hexahedron can not have ',k,'th surface'
      stop
    end if
    findHexSurfPC(:)=(Node(this%NodeInd(SurfTabHex(k,1)))%Pos(:)&
    &                +Node(this%NodeInd(SurfTabHex(k,2)))%Pos(:)&
    &                +Node(this%NodeInd(SurfTabHex(k,3)))%Pos(:)&
    &                +Node(this%NodeInd(SurfTabHex(k,4)))%Pos(:))/4d0
  end function
  
  !-------------------------------------------------
  ! find the area of this hexahedron's k_th surface
  !-------------------------------------------------
  function findHexSurfArea(this,k)
    class(typeHex),intent(in)::this
    integer,intent(in)::k
    double precision findHexSurfArea
    if(k<1.or.k>6)then
      write(*,'(a,i2,a)'),'ERROR: a hexahedron can not have ',k,'th surface'
      stop
    end if
    findHexSurfArea=find3PArea(Node(this%NodeInd(SurfTabHex(k,1)))%Pos(:),&
    &                          Node(this%NodeInd(SurfTabHex(k,2)))%Pos(:),&
    &                          Node(this%NodeInd(SurfTabHex(k,3)))%Pos(:))&
    &              +find3PArea(Node(this%NodeInd(SurfTabHex(k,3)))%Pos(:),&
    &                          Node(this%NodeInd(SurfTabHex(k,4)))%Pos(:),&
    &                          Node(this%NodeInd(SurfTabHex(k,1)))%Pos(:))
  end function
  
  !----------------------------------------------------------
  ! find the normal vector of this hexahedron's k_th surface
  !----------------------------------------------------------
  function findHexSurfNorm(this,k)
    class(typeHex),intent(in)::this
    integer,intent(in)::k
    double precision findHexSurfNorm(3)
    if(k<1.or.k>6)then
      write(*,'(a,i2,a)'),'ERROR: a hexahedron can not have ',k,'th surface'
      stop
    end if
    findHexSurfNorm(:)=find3PNorm(Node(this%NodeInd(SurfTabHex(k,1)))%Pos(:),&
    &                             Node(this%NodeInd(SurfTabHex(k,2)))%Pos(:),&
    &                             Node(this%NodeInd(SurfTabHex(k,3)))%Pos(:))
  end function
    
  !-------------------------------------------
  ! get the list of node indics of this facet
  !-------------------------------------------
  function getFacetNodeInd(this)
    class(typeFacet),intent(in)::this
    integer getFacetNodeInd(15) ! 15 is the maximum possible number of nodes an facet can have
    getFacetNodeInd(:)=0
    select case(this%ShapeType)
      case(2)
        getFacetNodeInd(1:this%NodeNum)=Tri(this%ShapeInd)%NodeInd(:)
      case(3)
        getFacetNodeInd(1:this%NodeNum)=Quad(this%ShapeInd)%NodeInd(:)
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !-------------------------------
  ! find the center of this facet
  !-------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findFacetPC(this)
    class(typeFacet),intent(in)::this
    double precision findFacetPC(3)
    select case(this%ShapeType)
      case(2)
        findFacetPC(:)=Tri(this%ShapeInd)%findPC()
      case(3)
        findFacetPC(:)=Quad(this%ShapeInd)%findPC()
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !-----------------------------
  ! find the area of this facet
  !-----------------------------
  function findFacetArea(this)
    class(typeFacet),intent(in)::this
    double precision findFacetArea
    select case(this%ShapeType)
      case(2)
        findFacetArea=Tri(this%ShapeInd)%findArea()
      case(3)
        findFacetArea=Quad(this%ShapeInd)%findArea()
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !--------------------------------------
  ! find the normal vector of this facet
  !--------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findFacetNorm(this)
    class(typeFacet),intent(in)::this
    double precision findFacetNorm(3)
    select case(this%ShapeType)
      case(2)
        findFacetNorm(:)=Tri(this%ShapeInd)%findNorm()
      case(3)
        findFacetNorm(:)=Quad(this%ShapeInd)%findNorm()
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !------------------------------------------
  ! get the geometrical entity of this facet
  !------------------------------------------
  function getFacetGeoEnti(this)
    class(typeFacet),intent(in)::this
    integer getFacetGeoEnti
    select case(this%ShapeType)
      case(2)
        getFacetGeoEnti=Tri(this%ShapeInd)%GeoEnti
      case(3)
        getFacetGeoEnti=Quad(this%ShapeInd)%GeoEnti
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !----------------------------------------
  ! get the neighbour element of this facet
  !-----------------------------------------
  ! Note: a facet may have 2 neighbour elements
  function getFacetNeibEle(this)
    class(typeFacet),intent(in)::this
    integer getFacetNeibEle(2)
    select case(this%ShapeType)
      case(2)
        getFacetNeibEle(:)=Tri(this%ShapeInd)%getNeibEle()
      case(3)
        getFacetNeibEle(:)=Quad(this%ShapeInd)%getNeibEle()
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !---------------------------------------------
  ! get the list of node indics of this element
  !---------------------------------------------
  function getEleNodeInd(this)
    class(typeEle),intent(in)::this
    integer getEleNodeInd(27) ! 27 is the maximum possible number of nodes an element can have
    getEleNodeInd(:)=0
    select case(this%ShapeType)
      case(4)
        getEleNodeInd(1:this%NodeNum)=Tet(this%ShapeInd)%NodeInd(:)
      case(5)
        getEleNodeInd(1:this%NodeNum)=Hex(this%ShapeInd)%NodeInd(:)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !---------------------------------
  ! find the center of this element
  !---------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findElePC(this)
    class(typeEle),intent(in)::this
    double precision findElePC(3)
    select case(this%ShapeType)
      case(4)
        findElePC(:)=Tet(this%ShapeInd)%findPC()
      case(5)
        findElePC(:)=Hex(this%ShapeInd)%findPC()
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !---------------------------------
  ! find the volume of this element
  !---------------------------------
  function findEleVol(this)
    class(typeEle),intent(in)::this
    double precision findEleVol
    select case(this%ShapeType)
      case(4)
        findEleVol=Tet(this%ShapeInd)%findVol()
      case(5)
        findEleVol=Hex(this%ShapeInd)%findVol()
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !--------------------------------------------
  ! get the geometrical entity of this element
  !--------------------------------------------
  function getEleGeoEnti(this)
    class(typeEle),intent(in)::this
    integer getEleGeoEnti
    select case(this%ShapeType)
      case(4)
        getEleGeoEnti=Tet(this%ShapeInd)%GeoEnti
      case(5)
        getEleGeoEnti=Hex(this%ShapeInd)%GeoEnti
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !-----------------------------------------
  ! get the partition index of this element
  !-----------------------------------------
  function getElePrt(this)
    class(typeEle),intent(in)::this
    integer getElePrt
    select case(this%ShapeType)
      case(4)
        getElePrt=Tet(this%ShapeInd)%Prt
      case(5)
        getElePrt=Hex(this%ShapeInd)%Prt
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !------------------------------------------------------
  ! get the k_th neighbour element index of this element
  !------------------------------------------------------
  function getEleNeib(this,k)
    class(typeEle),intent(in)::this
    integer,intent(in)::k
    integer getEleNeib
    getEleNeib=0
    select case(this%ShapeType)
      case(4)
        getEleNeib=Tet(this%ShapeInd)%getNeib(k)
      case(5)
        getEleNeib=Hex(this%ShapeInd)%getNeib(k)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !------------------------------------------------
  ! find the center of this element's k_th surface
  !------------------------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findEleSurfPC(this,k)
    class(typeEle),intent(in)::this
    integer,intent(in)::k
    double precision findEleSurfPC(3)
    select case(this%ShapeType)
      case(4)
        findEleSurfPC(:)=Tet(this%ShapeInd)%findSurfPC(k)
      case(5)
        findEleSurfPC(:)=Hex(this%ShapeInd)%findSurfPC(k)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !----------------------------------------------
  ! find the area of this element's k_th surface
  !----------------------------------------------
  function findEleSurfArea(this,k)
    class(typeEle),intent(in)::this
    integer,intent(in)::k
    double precision findEleSurfArea
    select case(this%ShapeType)
      case(4)
        findEleSurfArea=Tet(this%ShapeInd)%findSurfArea(k)
      case(5)
        findEleSurfArea=Hex(this%ShapeInd)%findSurfArea(k)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !-------------------------------------------------------
  ! find the normal vector of this element's k_th surface
  !-------------------------------------------------------
  function findEleSurfNorm(this,k)
    class(typeEle),intent(in)::this
    integer,intent(in)::k
    double precision findEleSurfNorm(3)
    select case(this%ShapeType)
      case(4)
        findEleSurfNorm(:)=Tet(this%ShapeInd)%findSurfNorm(k)
      case(5)
        findEleSurfNorm(:)=Hex(this%ShapeInd)%findSurfNorm(k)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
end module
