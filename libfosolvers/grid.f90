!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! data and elementary procedures for grid
!*****************************************
module moduleGrid
  private
  
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
  contains
    procedure,public::findPC=>findFacetPC
    procedure,public::findArea=>findFacetArea
    procedure,public::findNorm=>findFacetNorm
    procedure,public::findGeoEnti=>findFacetGeoEnti
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
  contains
    procedure,public::findPC=>findElePC
    procedure,public::findVol=>findEleVol
    procedure,public::findGeoEnti=>findEleGeoEnti
  end type
  type(typeEle),public,allocatable,save::Ele(:)
  integer,public,save::nEle
  
contains
  
  !---------------------------------
  ! find the position of this point
  !---------------------------------
  subroutine findPointPos(this,rst)
    class(typePoint),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=Node(this%NodeInd)%Pos(:)
  end subroutine
  
  !------------------------------
  ! find the center of this line
  !------------------------------
  subroutine findLinePC(this,rst)
    class(typeLine),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:))/2d0
  end subroutine
  
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
  subroutine findTriPC(this,rst)
    class(typeTri),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &      +Node(this%NodeInd(3))%Pos(:))/3d0
  end subroutine
  
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
  subroutine findTriNorm(this,rst)
    class(typeTri),intent(in)::this
    double precision,intent(out)::rst(3)
    call find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &               Node(this%NodeInd(3))%Pos(:),rst)
  end subroutine
  
  !---------------------------------------
  ! find the center of this quadrilateral
  !---------------------------------------
  ! Note: actually we are finding the average position of the vertices
  subroutine findQuadPC(this,rst)
    class(typeQuad),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &      +Node(this%NodeInd(3))%Pos(:)+Node(this%NodeInd(4))%Pos(:))/4d0
  end subroutine
  
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
  subroutine findQuadNorm(this,rst)
    class(typeQuad),intent(in)::this
    double precision,intent(out)::rst(3)
    call find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &               Node(this%NodeInd(3))%Pos(:),rst)
  end subroutine
  
  !-------------------------------------
  ! find the center of this tetrahedron
  !-------------------------------------
  ! Note: actually we are finding the average position of the vertices
  subroutine findTetPC(this,rst)
    class(typeTet),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &      +Node(this%NodeInd(3))%Pos(:)+Node(this%NodeInd(4))%Pos(:))/4d0
  end subroutine
  
  !-------------------------------------
  ! find the volume of this tetrahedron
  !-------------------------------------
  function findTetVol(this)
    class(typeTet),intent(in)::this
    double precision findTetVol,find4PVol
    findTetVol=find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                    Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:))
  end function
  
  !------------------------------------
  ! find the center of this hexahedron
  !------------------------------------
  ! Note: actually we are finding the average position of the vertices
  subroutine findHexPC(this,rst)
    class(typeHex),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=(Node(this%NodeInd(1))%Pos(:)+Node(this%NodeInd(2))%Pos(:)&
    &      +Node(this%NodeInd(3))%Pos(:)+Node(this%NodeInd(4))%Pos(:)&
    &      +Node(this%NodeInd(5))%Pos(:)+Node(this%NodeInd(6))%Pos(:)&
    &      +Node(this%NodeInd(7))%Pos(:)+Node(this%NodeInd(8))%Pos(:))/8d0
  end subroutine
  
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
  
  !-------------------------------
  ! find the center of this facet
  !-------------------------------
  ! Note: actually we are finding the average position of the vertices
  subroutine findFacetPC(this,rst)
    class(typeFacet),intent(in)::this
    double precision,intent(out)::rst(3)
    select case(this%ShapeType)
      case(2)
        call Tri(this%ShapeInd)%findPC(rst(:))
      case(3)
        call Quad(this%ShapeInd)%findPC(rst(:))
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end subroutine
  
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
  subroutine findFacetNorm(this,rst)
    class(typeFacet),intent(in)::this
    double precision,intent(out)::rst(3)
    select case(this%ShapeType)
      case(2)
        call Tri(this%ShapeInd)%findNorm(rst(:))
      case(3)
        call Quad(this%ShapeInd)%findNorm(rst(:))
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end subroutine
  
  !-------------------------------------------
  ! find the geometrical entity of this facet
  !-------------------------------------------
  function findFacetGeoEnti(this)
    class(typeFacet),intent(in)::this
    integer findFacetGeoEnti
    select case(this%ShapeType)
      case(2)
        findFacetGeoEnti=Tri(this%ShapeInd)%GeoEnti
      case(3)
        findFacetGeoEnti=Quad(this%ShapeInd)%GeoEnti
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet shapeType: ',this%shapeType
        stop
    end select
  end function
  
  !---------------------------------
  ! find the center of this element
  !---------------------------------
  ! Note: actually we are finding the average position of the vertices
  subroutine findElePC(this,rst)
    class(typeEle),intent(in)::this
    double precision,intent(out)::rst(3)
    select case(this%ShapeType)
      case(4)
        call Tet(this%ShapeInd)%findPC(rst(:))
      case(5)
        call Hex(this%ShapeInd)%findPC(rst(:))
      case default
        write(*,'(a,i2)'),'ERROR: unknown element shapeType: ',this%shapeType
        stop
    end select
  end subroutine
  
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
  
  !---------------------------------------------
  ! find the geometrical entity of this element
  !---------------------------------------------
  function findEleGeoEnti(this)
    class(typeEle),intent(in)::this
    integer findEleGeoEnti
    select case(this%ShapeType)
      case(4)
        findEleGeoEnti=Tet(this%ShapeInd)%GeoEnti
      case(5)
        findEleGeoEnti=Hex(this%ShapeInd)%GeoEnti
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
subroutine find3PNorm(P1,P2,P3,rst)
  double precision,intent(in)::P1(3),P2(3),P3(3)
  double precision,intent(out)::rst(3)
  double precision a(3),b(3)
  a(:)=P2(:)-P1(:)
  b(:)=P3(:)-P2(:)
  rst(1)=a(2)*b(3)-a(3)*b(2)
  rst(2)=a(3)*b(1)-a(1)*b(3)
  rst(3)=a(1)*b(2)-a(2)*b(1)
  rst(:)=rst(:)/sqrt(dot_product(rst(:),rst(:)))
end subroutine

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
