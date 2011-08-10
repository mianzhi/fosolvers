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
  end type
  type(typeQuad),public,allocatable,save::Quad(:)
  integer,public,save::nQuad
  
  !---------
  ! typeTet
  !---------
  type,public::typeTet
    integer NodeInd(4)
    integer GeoEnti
  contains
    procedure,public::findPC=>findTetPC
    procedure,public::findVol=>findTetVol
    procedure,public::getNeib=>getTetNeib
  end type
  type(typeTet),public,allocatable,save::Tet(:)
  integer,public,save::nTet
  
  !---------
  ! typeHex
  !---------
  type,public::typeHex
    integer NodeInd(8)
    integer GeoEnti
  contains
    procedure,public::findPC=>findHexPC
    procedure,public::findVol=>findHexVol
    procedure,public::getNeib=>getHexNeib
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
  contains
    procedure,public::findPC=>findFacetPC
    procedure,public::findArea=>findFacetArea
    procedure,public::findNorm=>findFacetNorm
    procedure,public::getGeoEnti=>getFacetGeoEnti
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
    integer Neib(6) ! 6 is the maximum possible number of neighbours
  contains
    procedure,public::findPC=>findElePC
    procedure,public::findVol=>findEleVol
    procedure,public::getGeoEnti=>getEleGeoEnti
    procedure,public::getNeib=>getEleNeib
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
    findLineLength=sqrt(dot_product(vect,vect))
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
  
  !----------------------------------------------------------
  ! get the k_th neighbour element index of this tetrahedron
  !----------------------------------------------------------
  ! Note: we are getting neighbour element index, not neighbour tetrahedron index, though the
  !       neighbour element is likely to be a tetrahedron
  function getTetNeib(this,k)
    class(typeTet),intent(in)::this
    integer,intent(in)::k
    integer getTetNeib,lPTet(4)
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
            if(this%NodeInd(j)==lPTet(1).or.this%NodeInd(j)==lPTet(2).or.&
            &  this%NodeInd(j)==lPTet(3).or.this%NodeInd(j)==lPTet(4))then
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
  
  !---------------------------------------------------------
  ! get the k_th neighbour element index of this hexahedron
  !---------------------------------------------------------
  ! Note: we are getting neighbour element index, not neighbour hexahedron index, though the
  !       neighbour element is likely to be a hexahedron
  function getHexNeib(this,k)
    class(typeHex),intent(in)::this
    integer,intent(in)::k
    integer getHexNeib,lPHex(8)
    logical maskHex(8)
    getHexNeib=0
    if(k<1.or.k>8)then
      write(*,'(a,i2,a)'),'ERROR: a hexahedron can not have ',k,'th surface'
      stop
    end if
    do i=1,nEle
      select case(Ele(i)%ShapeType)
        case(5) ! the neigbour element can be a hexahedron
          maskHex(:)=.false.
          lPHex(:)=Hex(Ele(i)%ShapeInd)%NodeInd(:)
          do j=1,8
            if(this%NodeInd(j)==lPHex(1).or.this%NodeInd(j)==lPHex(2).or.&
            &  this%NodeInd(j)==lPHex(3).or.this%NodeInd(j)==lPHex(4).or.&
            &  this%NodeInd(j)==lPHex(5).or.this%NodeInd(j)==lPHex(6).or.&
            &  this%NodeInd(j)==lPHex(7).or.this%NodeInd(j)==lPHex(8))then
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
  
end module

!******************************************************
! find the area of a triangle having P1~P3 as vertices
!******************************************************
function find3PArea(P1,P2,P3)
  double precision,intent(in)::P1(3),P2(3),P3(3)
  double precision find3PArea,a(3),b(3)
  a=P2-P1
  b=P3-P1
  find3PArea=((a(2)*b(3)-a(3)*b(2))**2d0&
  &          +(a(3)*b(1)-a(1)*b(3))**2d0&
  &          +(a(1)*b(2)-a(2)*b(1))**2d0)**0.5d0/2d0
end function

!***************************************************************
! find the normal vector of a triangle having P1~P3 as vertices
!***************************************************************
function find3PNorm(P1,P2,P3)
  double precision,intent(in)::P1(3),P2(3),P3(3)
  double precision find3PNorm(3),a(3),b(3),c(3)
  a(:)=P2(:)-P1(:)
  b(:)=P3(:)-P2(:)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  find3PNorm(:)=c(:)/sqrt(dot_product(c(:),c(:)))
end function

!***********************************************************
! find the volume of a tetrahedron having P1~P4 as vertices
!***********************************************************
function find4PVol(P1,P2,P3,P4)
  double precision,intent(in)::P1(3),P2(3),P3(3),P4(3)
  double precision find4PVol,y1z2,y1z3,y1z4,y2z1,y2z3,y2z4,y3z1,y3z2,y3z4,y4z1,y4z2,y4z3
  y1z2=P1(2)*P2(3)
  y1z3=P1(2)*P3(3)
  y1z4=P1(2)*P4(3)
  y2z1=P2(2)*P1(3)
  y2z3=P2(2)*P3(3)
  y2z4=P2(2)*P4(3)
  y3z1=P3(2)*P1(3)
  y3z2=P3(2)*P2(3)
  y3z4=P3(2)*P4(3)
  y4z1=P4(2)*P1(3)
  y4z2=P4(2)*P2(3)
  y4z3=P4(2)*P3(3)
  find4PVol=((P2(1)*(y3z4-y4z3)-P3(1)*(y2z4-y4z2)+P4(1)*(y2z3-y3z2))&
  &         -(P1(1)*(y3z4-y4z3)-P3(1)*(y1z4-y4z1)+P4(1)*(y1z3-y3z1))&
  &         +(P1(1)*(y2z4-y4z2)-P2(1)*(y1z4-y4z1)+P4(1)*(y1z2-y2z1))&
  &         -(P1(1)*(y2z3-y3z2)-P2(1)*(y1z3-y3z1)+P3(1)*(y1z2-y2z1)))/6d0
  find4PVol=abs(find4PVol)
end function
