!----------------------------------------------------------------------------- best with 100 columns

!> simple computational geometry procedures
module moduleSimpleGeometry
  private
  
  ! private constants
  integer,parameter::DIMS=3 !< dimensions
  integer,parameter::NV_TRI=3 !< number of vertices per triangle
  integer,parameter::NV_QUAD=4 !< number of vertices per quadrilateral
  integer,parameter::NV_TET=4 !< number of vertices per tetrahedron
  integer,parameter::NV_HEX=8 !< number of vertices per hexahedron
  
  ! public procedures
  public find3PArea
  public find4PArea
  public find3PNorm
  public find4PVol
  public find8PVol
  
contains
  
  !> find the area of a triangle
  pure function find3PArea(P)
    double precision,intent(in)::P(DIMS,NV_TRI) !< positions of the 3 vertices
    double precision find3PArea !< area of the triangle
    double precision a(DIMS),b(DIMS)
    a(:)=P(:,2)-P(:,1)
    b(:)=P(:,3)-P(:,1)
    find3PArea=((a(2)*b(3)-a(3)*b(2))**2d0&
    &          +(a(3)*b(1)-a(1)*b(3))**2d0&
    &          +(a(1)*b(2)-a(2)*b(1))**2d0)**0.5d0/2d0
  end function
  
  !> find the area of a quadrilateral
  pure function find4PArea(P)
    double precision,intent(in)::P(DIMS,NV_QUAD) !< positions of the 4 vertices
    double precision find4PArea !< area of the triangle
    find4PArea=find3PArea(P(:,[1,2,3]))+find3PArea(P(:,[3,4,1]))
  end function
  
  !> find the normal vector of a triangle
  pure function find3PNorm(P)
    double precision,intent(in)::P(DIMS,NV_TRI) !< positions of the 3 vertices
    double precision find3PNorm(DIMS) !< normal vector of the triangle
    double precision a(DIMS),b(DIMS),c(DIMS)
    a(:)=P(:,2)-P(:,1)
    b(:)=P(:,3)-P(:,2)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    find3PNorm(:)=c(:)/norm2(c)
  end function
  
  !> find the volume of a tetrahedron
  pure function find4PVol(P)
    double precision,intent(in)::P(DIMS,NV_TET) !< positions of the 4 vertices
    double precision find4PVol !< volume of the tetrahedron
    double precision y1z2,y1z3,y1z4,y2z1,y2z3,y2z4,y3z1,y3z2,y3z4,y4z1,y4z2,y4z3
    y1z2=P(2,1)*P(3,2)
    y1z3=P(2,1)*P(3,3)
    y1z4=P(2,1)*P(3,4)
    y2z1=P(2,2)*P(3,1)
    y2z3=P(2,2)*P(3,3)
    y2z4=P(2,2)*P(3,4)
    y3z1=P(2,3)*P(3,1)
    y3z2=P(2,3)*P(3,2)
    y3z4=P(2,3)*P(3,4)
    y4z1=P(2,4)*P(3,1)
    y4z2=P(2,4)*P(3,2)
    y4z3=P(2,4)*P(3,3)
    find4PVol=abs((P(1,2)*(y3z4-y4z3)-P(1,3)*(y2z4-y4z2)+P(1,4)*(y2z3-y3z2))&
    &            -(P(1,1)*(y3z4-y4z3)-P(1,3)*(y1z4-y4z1)+P(1,4)*(y1z3-y3z1))&
    &            +(P(1,1)*(y2z4-y4z2)-P(1,2)*(y1z4-y4z1)+P(1,4)*(y1z2-y2z1))&
    &            -(P(1,1)*(y2z3-y3z2)-P(1,2)*(y1z3-y3z1)+P(1,3)*(y1z2-y2z1)))/6d0
  end function
  
  !> find the volume of a hexahedron
  pure function find8PVol(P)
    double precision,intent(in)::P(DIMS,NV_HEX) !< positions of the 8 vertices
    double precision find8PVol !< volume of the tetrahedron
    
    find8PVol=find4PVol(P(:,[1,2,3,6]))+find4PVol(P(:,[1,3,4,6]))+find4PVol(P(:,[1,4,5,6]))&
    &        +find4PVol(P(:,[4,5,6,7]))+find4PVol(P(:,[3,4,6,7]))+find4PVol(P(:,[4,5,7,8]))
  end function
  
end module
