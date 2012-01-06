!----------------------------------------------------------------------------- best with 100 columns

!*****************
! simple geometry
!*****************
module moduleSimpleGeometry
  private
  
  ! procedures
  public find3PArea
  public find3PNorm
  public find4PVol
  
contains

  !------------------------------------------------------
  ! find the area of a triangle having P1~P3 as vertices
  !------------------------------------------------------
  pure function find3PArea(P1,P2,P3)
    double precision,intent(in)::P1(3),P2(3),P3(3)
    double precision find3PArea,a(3),b(3)
    a=P2-P1
    b=P3-P1
    find3PArea=((a(2)*b(3)-a(3)*b(2))**2d0&
    &          +(a(3)*b(1)-a(1)*b(3))**2d0&
    &          +(a(1)*b(2)-a(2)*b(1))**2d0)**0.5d0/2d0
  end function
  
  !---------------------------------------------------------------
  ! find the normal vector of a triangle having P1~P3 as vertices
  !---------------------------------------------------------------
  pure function find3PNorm(P1,P2,P3)
    double precision,intent(in)::P1(3),P2(3),P3(3)
    double precision find3PNorm(3),a(3),b(3),c(3)
    a(:)=P2(:)-P1(:)
    b(:)=P3(:)-P2(:)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    find3PNorm(:)=c(:)/norm2(c(:))
  end function
  
  !-----------------------------------------------------------
  ! find the volume of a tetrahedron having P1~P4 as vertices
  !-----------------------------------------------------------
  pure function find4PVol(P1,P2,P3,P4)
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
  
end module
