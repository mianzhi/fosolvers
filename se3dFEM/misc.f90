!----------------------------------------------------------------------------- best with 100 columns

!**************************************
! find the stress tensor for element k
!**************************************
subroutine findstress(k)
  use globalvar,only:strain,stress,lambda,mu
  double precision C(6,6)
  C(:,:)=0d0
  C(1:3,1:3)=lambda
  do i=1,3
    C(i,i)=C(i,i)+2d0*mu
    C(i+3,i+3)=mu
  end do
  stress(k,:)=matmul(C,strain(k,:))
end subroutine

!***********************************
! find the strain tensor for Tet. k
!***********************************
subroutine strainTet(strain,k)
  use globalvar,only:Pos,lTet,disp
  double precision strain(6)
  double precision gradPhi(4,3),invP(4,4),R(6,12),U(12)
  integer k
  U(1:3)=disp(lTet(k,1),:)
  U(4:6)=disp(lTet(k,2),:)
  U(7:9)=disp(lTet(k,3),:)
  U(10:12)=disp(lTet(k,4),:)
  invP(1,:)=[1d0,1d0,1d0,1d0]
  invP(2:4,1)=Pos(lTet(k,1),:)
  invP(2:4,2)=Pos(lTet(k,2),:)
  invP(2:4,3)=Pos(lTet(k,3),:)
  invP(2:4,4)=Pos(lTet(k,4),:)
  gradPhi(1,:)=[0d0,0d0,0d0]
  gradPhi(2,:)=[1d0,0d0,0d0]
  gradPhi(3,:)=[0d0,1d0,0d0]
  gradPhi(4,:)=[0d0,0d0,1d0]
  call inv(invP,4)
  gradPhi=matmul(invP,gradPhi)
  R(:,:)=0d0
  R([1,4,5],[1,4,7,10])=transpose(gradPhi)
  R([4,2,6],[2,5,8,11])=transpose(gradPhi)
  R([5,6,3],[3,6,9,12])=transpose(gradPhi)
  strain=matmul(R,U)
end subroutine

!***********************************
! find the strain tensor for Hex. m
!***********************************
subroutine strainHex(strain,m)
  use globalvar,only:Pos,lHex,disp
  double precision strain(6)
  double precision gp(5),gw(5),xr(8),yr(8),zr(8),x(8),y(8),z(8)
  double precision gradPhi_r(3,8),Jac(3,3),gradPhi(3,8),R(6,24),U(24)
  integer m
  U(1:3)=disp(lHex(m,1),:)
  U(4:6)=disp(lHex(m,2),:)
  U(7:9)=disp(lHex(m,3),:)
  U(10:12)=disp(lHex(m,4),:)
  U(13:15)=disp(lHex(m,5),:)
  U(16:18)=disp(lHex(m,6),:)
  U(19:21)=disp(lHex(m,7),:)
  U(22:24)=disp(lHex(m,8),:)
  gp(1)=-1d0/3d0*sqrt(5d0+2d0*sqrt(10d0/7d0))
  gp(2)=-1d0/3d0*sqrt(5d0-2d0*sqrt(10d0/7d0))
  gp(3)=0d0
  gp(4)=-gp(2)
  gp(5)=-gp(1)
  gw([1,5])=(322d0-13d0*sqrt(70d0))/900d0
  gw([2,4])=(322d0+13d0*sqrt(70d0))/900d0
  gw(3)=128d0/225d0
  xr=[-1d0,1d0,1d0,-1d0,-1d0,1d0,1d0,-1d0]
  yr=[-1d0,-1d0,1d0,1d0,-1d0,-1d0,1d0,1d0]
  zr=[-1d0,-1d0,-1d0,-1d0,1d0,1d0,1d0,1d0]
  x(:)=Pos(lHex(m,1:8),1)
  y(:)=Pos(lHex(m,1:8),2)
  z(:)=Pos(lHex(m,1:8),3)
  strain(:)=0d0
  do i=1,5
    do j=1,5
      do k=1,5
        do l=1,8
          gradPhi_r(1,l)=1d0/8d0*xr(l)*(1+yr(l)*gp(j))*(1+zr(l)*gp(k))
          gradPhi_r(2,l)=1d0/8d0*(1+xr(l)*gp(i))*yr(l)*(1+zr(l)*gp(k))
          gradPhi_r(3,l)=1d0/8d0*(1+xr(l)*gp(i))*(1+yr(l)*gp(j))*zr(l)
        end do
        Jac(:,1)=matmul(gradPhi_r,x)
        Jac(:,2)=matmul(gradPhi_r,y)
        Jac(:,3)=matmul(gradPhi_r,z)
        call inv(Jac,3)
        gradPhi=matmul(Jac,gradPhi_r)
        R(:,:)=0d0
        R([1,4,5],[1,4,7,10,13,16,19,22])=gradPhi
        R([4,2,6],[2,5,8,11,14,17,20,23])=gradPhi
        R([5,6,3],[3,6,9,12,15,18,21,24])=gradPhi
        ! note that here we are finding weighted average instead of integral
        strain=strain+gw(i)*gw(j)*gw(k)/8d0*matmul(R,U)
      end do
    end do
  end do
end subroutine

!************************************************
! sort the elements (to minimize the band-width)
!************************************************
subroutine sortEle(m)
  use globalvar,only:lEle,lTet,lHex,nEle,nTet,nHex,Pos
  integer m,NodeList(8),temp_int(2)
  double precision sortData(nEle),temp_double
  ! initialize
  do i=1,nEle
    n=lEle(i,1)
    if(n==4)then
      NodeList(1:n)=lTet(lEle(i,2),1:n)
    else
    if(n==8)then
      NodeList(1:n)=lHex(lEle(i,2),1:n)
    end if
    end if
    sortData(i)=sum([(Pos(NodeList(i),m),i=1,n)])/dble(n)
  end do
  ! insertion sort
  do i=2,nEle
    do j=i,2,-1
      if(sortData(j-1)<=sortData(j))then
        exit
      end if
      temp_int(:)=lEle(j,:)
      lEle(j,:)=lEle(j-1,:)
      lEle(j-1,:)=temp_int(:)
      temp_double=sortData(j)
      sortData(j)=sortData(j-1)
      sortData(j-1)=temp_double
    end do
  end do
end subroutine

!*****************************************
! find a free space in the frontal matrix
!*****************************************
function findSpace(F2N,width)
  integer width,findSpace
  integer F2N(width)
  findSpace=0
  do i=1,width
    if(F2N(i)==0)then
      findSpace=i
      exit
    end if
  end do
  if(findSpace==0)then
    write(*,'(a)'),'ERROR: frontal matrix over flow.'
    stop
  end if
end function

!************************
! volume of a hexahedron
!************************
function volHex(P1,P2,P3,P4,P5,P6,P7,P8)
  double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3),volHex
  double precision volTet
  volHex=volTet(P1,P2,P3,P6)+volTet(P1,P3,P4,P6)+volTet(P1,P4,P5,P6)&
  &     +volTet(P4,P5,P6,P7)+volTet(P3,P4,P6,P7)+volTet(P4,P5,P7,P8)
end function

!*************************
! volume of a tetrahedron
!*************************
function volTet(P1,P2,P3,P4)
  double precision P1(3),P2(3),P3(3),P4(3),volTet
  double precision y1z2,y1z3,y1z4,y2z1,y2z3,y2z4,y3z1,y3z2,y3z4,y4z1,y4z2,y4z3
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
  volTet=((P2(1)*(y3z4-y4z3)-P3(1)*(y2z4-y4z2)+P4(1)*(y2z3-y3z2))&
  &      -(P1(1)*(y3z4-y4z3)-P3(1)*(y1z4-y4z1)+P4(1)*(y1z3-y3z1))&
  &      +(P1(1)*(y2z4-y4z2)-P2(1)*(y1z4-y4z1)+P4(1)*(y1z2-y2z1))&
  &      -(P1(1)*(y2z3-y3z2)-P2(1)*(y1z3-y3z1)+P3(1)*(y1z2-y2z1)))/6d0
  volTet=abs(volTet)
end function

!****************************
! normal vector of a polygon
!****************************
subroutine normal(rst,P1,P2,P3)
  double precision rst(3),P1(3),P2(3),P3(3)
  double precision a(3),b(3)
  a(:)=P2(:)-P1(:)
  b(:)=P3(:)-P2(:)
  rst(1)=a(2)*b(3)-a(3)*b(2)
  rst(2)=a(3)*b(1)-a(1)*b(3)
  rst(3)=a(1)*b(2)-a(2)*b(1)
  rst(:)=rst(:)/sqrt(dot_product(rst(:),rst(:)))
end subroutine

!****************************
! surface area of a triangle
!****************************
function areaTri(P1,P2,P3)
  double precision areaTri,P1(3),P2(3),P3(3)
  double precision a(3),b(3)
  a=P2-P1
  b=P3-P1
  areaTri=((a(2)*b(3)-a(3)*b(2))**2d0&
  &       +(a(3)*b(1)-a(1)*b(3))**2d0&
  &       +(a(1)*b(2)-a(2)*b(1))**2d0)**0.5d0/2d0
end function

!*********************
! inverse of a matrix
!*********************
subroutine inv(A,n)
  integer n,im
  double precision A(n,n)
  double precision M(n,2*n),temp_vect(2*n)
  M(:,:)=0d0
  M(1:n,1:n)=A
  do i=1,n
    M(i,n+i)=1d0
  end do
  ! forward elimination
  do i=1,n-1
    im=i
    ! find pivot
    do j=i+1,n
      if(abs(M(j,i))>abs(M(im,i)))then
        im=j
      end if
    end do
    if(abs(M(im,i))+1d0<=1d0)then
      write(*,'(a)'),'ERROR: matrix can not be inversed.'
      stop
    end if
    ! switch
    if(im/=i)then
      temp_vect(i:2*n)=M(i,i:2*n)
      M(i,i:2*n)=M(im,i:2*n)
      M(im,i:2*n)=temp_vect(i:2*n)
    end if
    ! eliminate
    do j=i+1,n
      M(j,i+1:2*n)=M(j,i+1:2*n)-M(j,i)/M(i,i)*M(i,i+1:2*n)
    end do
    M(i+1:n,i)=0d0
  end do
  ! backward substitution
  do i=n,2,-1
    do j=1,i-1
      M(j,i+1:2*n)=M(j,i+1:2*n)-M(j,i)/M(i,i)*M(i,i+1:2*n)
    end do
    M(1:i-1,i)=0d0
  end do
  do i=1,n
    M(i,n+1:2*n)=M(i,n+1:2*n)/M(i,i)
  end do
  A=M(1:n,n+1:2*n)
end subroutine

!**************
! progress bar
!**************
subroutine progBar(v,vtot)
  double precision v,vtot
  real prog
  prog=v/vtot*100d0
  write(*,'(6a,f5.1,a,$)'),(char(8),i=1,6),prog,'%'
end subroutine
