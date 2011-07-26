!----------------------------------------------------------------------------- best with 100 columns

!***********************************
! stiffness matrix for Tet. element
!***********************************
subroutine stiTet(ATet,P1,P2,P3,P4,lambda,mu)
  double precision ATet(12,12),P1(3),P2(3),P3(3),P4(3),lambda,mu
  double precision gradPhi(4,3),invP(4,4),C(6,6),R(6,12)
  double precision volTet
  
  invP(1,:)=[1d0,1d0,1d0,1d0]
  invP(2,:)=[P1(1),P2(1),P3(1),P4(1)]
  invP(3,:)=[P1(2),P2(2),P3(2),P4(2)]
  invP(4,:)=[P1(3),P2(3),P3(3),P4(3)]
  gradPhi(1,:)=[0d0,0d0,0d0]
  gradPhi(2,:)=[1d0,0d0,0d0]
  gradPhi(3,:)=[0d0,1d0,0d0]
  gradPhi(4,:)=[0d0,0d0,1d0]
  
  call inv(invP,4)
  gradPhi=matmul(invP,gradPhi)
  
  C(:,:)=0d0
  C(1:3,1:3)=lambda
  do i=1,3
    C(i,i)=C(i,i)+2*mu
    C(i+3,i+3)=mu
  end do
  
  R(:,:)=0d0
  
  R([1,4,5],[1,4,7,10])=transpose(gradPhi)
  R([4,2,6],[2,5,8,11])=transpose(gradPhi)
  R([5,6,3],[3,6,9,12])=transpose(gradPhi)
  
  ATet=volTet(P1,P2,P3,P4)*matmul(matmul(transpose(R),C),R)
end subroutine

!***********************************
! stiffness matrix for Hex. element
!***********************************
subroutine stiHex(AHex,P1,P2,P3,P4,P5,P6,P7,P8,lambda,mu)
  double precision AHex(24,24),P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3),lambda,mu
  double precision xr(8),yr(8),zr(8),x(8),y(8),z(8)
  double precision gp(5),gw(5)
  double precision C(6,6),gradPhi_r(3,8),Jac(3,3),gradPhi(3,8),R(6,24),detJ
  
  AHex(:,:)=0d0
  
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
  x=[P1(1),P2(1),P3(1),P4(1),P5(1),P6(1),P7(1),P8(1)]
  y=[P1(2),P2(2),P3(2),P4(2),P5(2),P6(2),P7(2),P8(2)]
  z=[P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3)]
  
  C(:,:)=0d0
  C(1:3,1:3)=lambda
  do i=1,3
    C(i,i)=C(i,i)+2*mu
    C(i+3,i+3)=mu
  end do
  
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
        detJ=abs(Jac(1,1)*(Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2))&
        &       -Jac(1,2)*(Jac(2,1)*Jac(3,3)-Jac(2,3)*Jac(3,1))&
        &       +Jac(1,3)*(Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1)))
        call inv(Jac,3)
        gradPhi=matmul(Jac,gradPhi_r)
        R(:,:)=0d0
        R([1,4,5],[1,4,7,10,13,16,19,22])=gradPhi
        R([4,2,6],[2,5,8,11,14,17,20,23])=gradPhi
        R([5,6,3],[3,6,9,12,15,18,21,24])=gradPhi
        AHex=AHex+gw(i)*gw(j)*gw(k)*detJ*matmul(matmul(transpose(R),C),R)
      end do
    end do
  end do
end subroutine
