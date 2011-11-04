!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! data and elementary procedures for grid
!*****************************************
module moduleGrid
  private
  
  ! some constants
  integer,parameter,public::DIMS=3
  
  integer,parameter,public::POINT_TYPE=15
  integer,parameter,public::LINE_TYPE=1
  integer,parameter,public::TRI_TYPE=2
  integer,parameter,public::QUAD_TYPE=3
  integer,parameter,public::TET_TYPE=4
  integer,parameter,public::HEX_TYPE=5
  
  integer,parameter,public::POINT_NODE_NUM=1
  integer,parameter,public::LINE_NODE_NUM=2
  integer,parameter,public::TRI_NODE_NUM=3
  integer,parameter,public::QUAD_NODE_NUM=4
  integer,parameter,public::TET_NODE_NUM=4
  integer,parameter,public::HEX_NODE_NUM=8
  
  integer,parameter,public::TET_SURF_NUM=4
  integer,parameter,public::HEX_SURF_NUM=6
  
  integer,parameter,public::FACET_MAX_NODE_NUM=15
  integer,parameter,public::ELE_MAX_NODE_NUM=27
  
  integer,parameter,public::FACET_NEIB_ELE_NUM=2
  
  integer,parameter,public::ELE_MAX_SURF_NUM=6
  
  ! table of surface nodes for all kinds of elements
  integer,public,save::SurfTabTet(TET_SURF_NUM,TRI_NODE_NUM),SurfTabHex(HEX_SURF_NUM,QUAD_NODE_NUM)
  
  ! bounding box of the geometry
  double precision,public,save::BoundBox(2,DIMS)
  
  ! make some functions properly accessible (both privately & publicly)
  interface
    function find3PNorm(P1,P2,P3)
      double precision,intent(in)::P1(3),P2(3),P3(3)
      double precision find3PNorm(3)
    end function
  end interface
  public::find3PNorm
  
  !----------
  ! typeNode
  !----------
  type,public::typeNode
    double precision Pos(DIMS)
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
    integer NodeInd(LINE_NODE_NUM)
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
    integer NodeInd(TRI_NODE_NUM)
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
    integer NodeInd(QUAD_NODE_NUM)
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
    integer NodeInd(TET_NODE_NUM)
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
    integer NodeInd(HEX_NODE_NUM)
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
    integer NodeInd(FACET_MAX_NODE_NUM)
    integer GeoEnti
    integer NeibEle(FACET_NEIB_ELE_NUM) ! a facet may have 2 neighbour elements
    double precision PC(DIMS)
    double precision Area
    double precision Norm(DIMS)
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
    integer NodeInd(ELE_MAX_NODE_NUM)
    integer GeoEnti
    integer Prt
    integer Neib(ELE_MAX_SURF_NUM)
    double precision PC(DIMS)
    double precision Vol
    double precision SurfPC(ELE_MAX_SURF_NUM,DIMS)
    double precision SurfArea(ELE_MAX_SURF_NUM)
    double precision SurfNorm(ELE_MAX_SURF_NUM,DIMS)
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
    double precision findPointPos(DIMS)
    findPointPos(:)=Node(this%NodeInd)%Pos(:)
  end function
  
  !------------------------------
  ! find the center of this line
  !------------------------------
  function findLinePC(this)
    class(typeLine),intent(in)::this
    double precision findLinePC(DIMS)
    findLinePC(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:))/2d0
  end function
  
  !------------------------------
  ! find the length of this line
  !------------------------------
  function findLineLength(this)
    class(typeLine),intent(in)::this
    double precision findLineLength,vect(DIMS)
    vect(:)=Node(this%NodeInd(2))%Pos(:)-Node(this%NodeInd(1))%Pos(:)
    findLineLength=norm2(vect)
  end function
  
  !----------------------------------
  ! find the center of this triangle
  !----------------------------------
  ! Note: actually we are finding the average position of the vertices
  function findTriPC(this)
    class(typeTri),intent(in)::this
    double precision findTriPC(DIMS)
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
    double precision findTriNorm(DIMS)
    findTriNorm(:)=find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                         Node(this%NodeInd(3))%Pos(:))
  end function
  
  !--------------------------------------------
  ! get the neighbour element of this triangle
  !--------------------------------------------
  ! Note: one triangle may have 2 neighbour elements
  function getTriNeibEle(this)
    class(typeTri),intent(in)::this
    integer getTriNeibEle(FACET_NEIB_ELE_NUM),lPTet(TET_NODE_NUM)
    logical maskTri(TRI_NODE_NUM)
    getTriNeibEle(:)=0
    k=0
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(TET_TYPE) ! the neigbour element can be a tetrahedron
          maskTri(:)=.false.
          lPTet(:)=Tet(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,TRI_NODE_NUM
            if(any(lPTet(:)==this%NodeInd(j)))then
              maskTri(j)=.true.
            end if
          end do
          if(all(maskTri(:)))then
            k=k+1
            getTriNeibEle(k)=i
            if(k==FACET_NEIB_ELE_NUM)then
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
    double precision findQuadPC(DIMS)
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
    double precision findQuadNorm(DIMS)
    findQuadNorm(:)=find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                          Node(this%NodeInd(3))%Pos(:))
  end function
  
  !-------------------------------------------------
  ! get the neighbour element of this quadrilateral
  !-------------------------------------------------
  ! Note: one quadirlateral may have 2 neighbour elements
  function getQuadNeibEle(this)
    class(typeQuad),intent(in)::this
    integer getQuadNeibEle(FACET_NEIB_ELE_NUM),lPHex(HEX_NODE_NUM)
    logical maskQuad(QUAD_NODE_NUM)
    getQuadNeibEle(:)=0
    k=0
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(HEX_TYPE) ! the neigbour element can be a hexahedron
          maskQuad(:)=.false.
          lPHex(:)=Hex(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,QUAD_NODE_NUM
            if(any(lPHex(:)==this%NodeInd(j)))then
              maskQuad(j)=.true.
            end if
          end do
          if(all(maskQuad(:)))then
            k=k+1
            getQuadNeibEle(k)=i
            if(k==FACET_NEIB_ELE_NUM)then
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
    double precision findTetPC(DIMS)
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
    integer getTetNeib,lPTet(TET_NODE_NUM),lPTri(TRI_NODE_NUM)
    logical maskTet(TET_NODE_NUM)
    getTetNeib=0
    if(k<1.or.k>TET_SURF_NUM)then
      write(*,'(a,i2,a)'),'ERROR: a tetrahedron can not have ',k,'th surface'
      stop
    end if
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(TET_TYPE) ! the neigbour element can be a tetrahedron
          maskTet(:)=.false.
          lPTet(:)=Tet(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,TET_NODE_NUM
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
        case(TRI_TYPE) ! the neigbour facet can be a triangle
          maskTet(:)=.false.
          lPTri(:)=Tri(Facet(i)%ShapeInd)%NodeInd(:)
          do j=1,TET_NODE_NUM
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
    double precision findTetSurfPC(DIMS)
    if(k<1.or.k>TET_SURF_NUM)then
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
    if(k<1.or.k>TET_SURF_NUM)then
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
    double precision findTetSurfNorm(DIMS)
    if(k<1.or.k>TET_SURF_NUM)then
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
    double precision findHexPC(DIMS)
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
    integer getHexNeib,lPHex(HEX_NODE_NUM),lPQuad(QUAD_NODE_NUM)
    logical maskHex(HEX_NODE_NUM)
    getHexNeib=0
    if(k<1.or.k>HEX_SURF_NUM)then
      write(*,'(a,i2,a)'),'ERROR: a hexahedron can not have ',k,'th surface'
      stop
    end if
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(HEX_TYPE) ! the neigbour element can be a hexahedron
          maskHex(:)=.false.
          lPHex(:)=Hex(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,HEX_NODE_NUM
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
        case(QUAD_TYPE) ! the neigbour facet can be a quadrilateral
          maskHex(:)=.false.
          lPQuad(:)=Quad(Facet(i)%ShapeInd)%NodeInd(:)
          do j=1,HEX_NODE_NUM
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
    if(k<1.or.k>HEX_SURF_NUM)then
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
    if(k<1.or.k>HEX_SURF_NUM)then
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
    double precision findHexSurfNorm(DIMS)
    if(k<1.or.k>HEX_SURF_NUM)then
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
    integer getFacetNodeInd(FACET_MAX_NODE_NUM)
    getFacetNodeInd(:)=0
    select case(this%ShapeType)
      case(TRI_TYPE)
        getFacetNodeInd(1:this%NodeNum)=Tri(this%ShapeInd)%NodeInd(:)
      case(QUAD_TYPE)
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
    double precision findFacetPC(DIMS)
    select case(this%ShapeType)
      case(TRI_TYPE)
        findFacetPC(:)=Tri(this%ShapeInd)%findPC()
      case(QUAD_TYPE)
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
      case(TRI_TYPE)
        findFacetArea=Tri(this%ShapeInd)%findArea()
      case(QUAD_TYPE)
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
    double precision findFacetNorm(DIMS)
    select case(this%ShapeType)
      case(TRI_TYPE)
        findFacetNorm(:)=Tri(this%ShapeInd)%findNorm()
      case(QUAD_TYPE)
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
      case(TRI_TYPE)
        getFacetGeoEnti=Tri(this%ShapeInd)%GeoEnti
      case(QUAD_TYPE)
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
    integer getFacetNeibEle(FACET_NEIB_ELE_NUM)
    select case(this%ShapeType)
      case(TRI_TYPE)
        getFacetNeibEle(:)=Tri(this%ShapeInd)%getNeibEle()
      case(QUAD_TYPE)
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
    integer getEleNodeInd(ELE_MAX_NODE_NUM)
    getEleNodeInd(:)=0
    select case(this%ShapeType)
      case(TET_TYPE)
        getEleNodeInd(1:this%NodeNum)=Tet(this%ShapeInd)%NodeInd(:)
      case(HEX_TYPE)
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
    double precision findElePC(DIMS)
    select case(this%ShapeType)
      case(TET_TYPE)
        findElePC(:)=Tet(this%ShapeInd)%findPC()
      case(HEX_TYPE)
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
      case(TET_TYPE)
        findEleVol=Tet(this%ShapeInd)%findVol()
      case(HEX_TYPE)
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
      case(TET_TYPE)
        getEleGeoEnti=Tet(this%ShapeInd)%GeoEnti
      case(HEX_TYPE)
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
      case(TET_TYPE)
        getElePrt=Tet(this%ShapeInd)%Prt
      case(HEX_TYPE)
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
      case(TET_TYPE)
        getEleNeib=Tet(this%ShapeInd)%getNeib(k)
      case(HEX_TYPE)
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
    double precision findEleSurfPC(DIMS)
    select case(this%ShapeType)
      case(TET_TYPE)
        findEleSurfPC(:)=Tet(this%ShapeInd)%findSurfPC(k)
      case(HEX_TYPE)
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
      case(TET_TYPE)
        findEleSurfArea=Tet(this%ShapeInd)%findSurfArea(k)
      case(HEX_TYPE)
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
    double precision findEleSurfNorm(DIMS)
    select case(this%ShapeType)
      case(TET_TYPE)
        findEleSurfNorm(:)=Tet(this%ShapeInd)%findSurfNorm(k)
      case(HEX_TYPE)
        findEleSurfNorm(:)=Hex(this%ShapeInd)%findSurfNorm(k)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end function
  
end module
