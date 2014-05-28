!----------------------------------------------------------------------------- best with 100 columns

!> geometry module
module modGeometry
  public
  
  ! constants
  integer,private,parameter::DIMS=3 !< dimensions
  
  ! long operations are interfaced
  interface
  end interface
  
contains
  
  !> find the area of a triangle
  pure function a3p(p)
    double precision,intent(in)::p(DIMS,3) !< positions of the 3 vertices
    double precision a3p !< area of the triangle
    double precision a(DIMS),b(DIMS)
    
    a(:)=p(:,2)-p(:,1)
    b(:)=p(:,3)-p(:,1)
    a3p=((a(2)*b(3)-a(3)*b(2))**2d0&
    &   +(a(3)*b(1)-a(1)*b(3))**2d0&
    &   +(a(1)*b(2)-a(2)*b(1))**2d0)**0.5d0/2d0
  end function
  
  !> find the area of a quadrilateral
  pure function a4p(p)
    double precision,intent(in)::p(DIMS,4) !< positions of the 4 vertices
    double precision a4p !< area of the triangle
    
    a4p=a3p(p(:,[1,2,3]))+a3p(p(:,[3,4,1]))
  end function
  
  !> find the normal vector of a triangle
  pure function n3p(p)
    double precision,intent(in)::p(DIMS,3) !< positions of the 3 vertices
    double precision n3p(DIMS) !< normal vector of the triangle
    double precision a(DIMS),b(DIMS),c(DIMS)
    
    a(:)=p(:,2)-p(:,1)
    b(:)=p(:,3)-p(:,2)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    n3p(:)=c(:)/norm2(c)
  end function
  
  !> find the volume of a tetrahedron
  pure function v4p(p)
    double precision,intent(in)::p(DIMS,4) !< positions of the 4 vertices
    double precision v4p !< volume of the tetrahedron
    double precision y1z2,y1z3,y1z4,y2z1,y2z3,y2z4,y3z1,y3z2,y3z4,y4z1,y4z2,y4z3
    
    y1z2=p(2,1)*p(3,2)
    y1z3=p(2,1)*p(3,3)
    y1z4=p(2,1)*p(3,4)
    y2z1=p(2,2)*p(3,1)
    y2z3=p(2,2)*p(3,3)
    y2z4=p(2,2)*p(3,4)
    y3z1=p(2,3)*p(3,1)
    y3z2=p(2,3)*p(3,2)
    y3z4=p(2,3)*p(3,4)
    y4z1=p(2,4)*p(3,1)
    y4z2=p(2,4)*p(3,2)
    y4z3=p(2,4)*p(3,3)
    v4p=abs((p(1,2)*(y3z4-y4z3)-p(1,3)*(y2z4-y4z2)+p(1,4)*(y2z3-y3z2))&
    &      -(p(1,1)*(y3z4-y4z3)-p(1,3)*(y1z4-y4z1)+p(1,4)*(y1z3-y3z1))&
    &      +(p(1,1)*(y2z4-y4z2)-p(1,2)*(y1z4-y4z1)+p(1,4)*(y1z2-y2z1))&
    &      -(p(1,1)*(y2z3-y3z2)-p(1,2)*(y1z3-y3z1)+p(1,3)*(y1z2-y2z1)))/6d0
  end function
  
  !> find the volume of a hexahedron
  pure function v8p(p)
    double precision,intent(in)::p(DIMS,8) !< positions of the 8 vertices
    double precision v8p !< volume of the tetrahedron
    
    v8p=v4p(p(:,[1,2,3,6]))+v4p(p(:,[1,3,4,6]))+v4p(p(:,[1,4,5,6]))&
    &  +v4p(p(:,[4,5,6,7]))+v4p(p(:,[3,4,6,7]))+v4p(p(:,[4,5,7,8]))
  end function
  
end module
