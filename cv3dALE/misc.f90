!author: Wang, Mianzhi
!miscellaneous functions:
!	voltrans(), volume_hex(), volume_tetra(), area_quad_2d(), dist_3d(), dist_2d()
!miscellaneous subroutines:
!	shearstress(), projarea_subcell(), subcell(), area_quad_3d()

subroutine shearstress(P1,P2,P3,P4,P5,P6,P7,P8,U1,U2,U3,U4,U5,U6,U7,U8,rst)
!calculate the 3 upper-triangular elements of stress tensor
!input:
!	P1~8: [1*3] position of vertices
!	U1~8: [1*3] velocity at vertices
use globalvar,only:visc
double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3),&
				&U1(3),U2(3),U3(3),U4(3),U5(3),U6(3),U7(3),U8(3),&
				&rst(3),temp1(3,3),temp2(3,3),temp,dist_3d,temp3(3,3)
temp1(1,:)=[U2(1)+U3(1)+U6(1)+U7(1)-U1(1)-U4(1)-U5(1)-U8(1),&
		   &U5(1)+U6(1)+U7(1)+U8(1)-U1(1)-U2(1)-U3(1)-U4(1),&
		   &U3(1)+U4(1)+U7(1)+U8(1)-U1(1)-U2(1)-U5(1)-U6(1)]
temp1(2,:)=[U2(2)+U3(2)+U6(2)+U7(2)-U1(2)-U4(2)-U5(2)-U8(2),&
		   &U5(2)+U6(2)+U7(2)+U8(2)-U1(2)-U2(2)-U3(2)-U4(2),&
		   &U3(2)+U4(2)+U7(2)+U8(2)-U1(2)-U2(2)-U5(2)-U6(2)]
temp1(3,:)=[U2(3)+U3(3)+U6(3)+U7(3)-U1(3)-U4(3)-U5(3)-U8(3),&
		   &U5(3)+U6(3)+U7(3)+U8(3)-U1(3)-U2(3)-U3(3)-U4(3),&
		   &U3(3)+U4(3)+U7(3)+U8(3)-U1(3)-U2(3)-U5(3)-U6(3)]
temp1=temp1/4
temp=dist_3d((P2+P3+P6+P7)/4,(P1+P4+P5+P8)/4)
temp2(1,:)=((P2+P3+P6+P7)-(P1+P4+P5+P8))/4/temp/temp
temp=dist_3d((P5+P6+P7+P8)/4,(P1+P2+P3+P4)/4)
temp2(2,:)=((P5+P6+P7+P8)-(P1+P2+P3+P4))/4/temp/temp
temp=dist_3d((P3+P4+P7+P8)/4,(P1+P2+P5+P6)/4)
temp2(3,:)=((P3+P4+P7+P8)-(P1+P2+P5+P6))/4/temp/temp
temp3=matmul(temp1,temp2) ! velocity gradient tensor
rst=[temp3(1,2)+temp3(2,1),temp3(1,3)+temp3(3,1),temp3(2,3)+temp3(3,2)]
rst=visc*rst/2
end subroutine

function voltrans(P1,P2,P3,P4,pp1,pp2,pp3,pp4)
!calculate volume transfer of a cell surface during rezoning
!input:
!	P1~4: [1*3] position of vertices before rezoning
!	pp1~4: [1*3] position of vertices after resoning
!note: sequence of input must be CCW (see from inside) to gurantee the correct sign
!called function & subroutines
!	area_quad_3d
	double precision P1(3),P2(3),P3(3),P4(3),pp1(3),pp2(3),pp3(3),pp4(3),voltrans,temp(3)
	double precision area1(3),area2(3) !normal area of surface before & after rezoning
	call area_quad_3d(P1,P2,P3,P4,area1)
	call area_quad_3d(pp1,pp2,pp3,pp4,area2)
	temp=(pp1+pp2+pp3+pp4)-(P1+P2+P3+P4)
	voltrans=-dot_product((area1+area2)/2,temp/4)
end function

subroutine projarea_subcell(P1,P2,P3,P4,P5,P6,P7,P8,n,rst)
!calculate vertice cell's subcells' projection area on which stress exerts
!input:
!	P1~8: [1*3] position of vertices
!	n: the index of the subcell (1~8)
!output:
!	rst: [1*3] result
!called functions & subroutines:
!	subcell(), area_quad_2d()
	double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3)
	double precision rst(3),area_quad_2d
	double precision ptab(8,3) !store the position of candidate points
	integer n
	integer itab(8,12) !point table for selecting the 3 surfaces
	ptab(1,:)=P1
	ptab(2,:)=P2
	ptab(3,:)=P3
	ptab(4,:)=P4
	ptab(5,:)=P5
	ptab(6,:)=P6
	ptab(7,:)=P7
	ptab(8,:)=P8
	call subcell(ptab(1,:),ptab(2,:),ptab(3,:),ptab(4,:),ptab(5,:),ptab(6,:),ptab(7,:),ptab(8,:),n)
	!for different n:
	itab(1,:)=[1,4,3,2,1,5,8,4,1,2,6,5]
	itab(2,:)=[1,4,3,2,2,3,7,6,1,2,6,5]
	itab(3,:)=[1,4,3,2,2,3,7,6,3,4,8,7]
	itab(4,:)=[1,4,3,2,1,5,8,4,3,4,8,7]
	itab(5,:)=[5,6,7,8,1,5,8,4,1,2,6,5]
	itab(6,:)=[5,6,7,8,2,3,7,6,1,2,6,5]
	itab(7,:)=[5,6,7,8,2,3,7,6,3,4,8,7]
	itab(8,:)=[5,6,7,8,1,5,8,4,3,4,8,7]
	rst(1)=area_quad_2d(ptab(itab(n,1),[2,3]),ptab(itab(n,2),[2,3]),&
					   &ptab(itab(n,3),[2,3]),ptab(itab(n,4),[2,3]))&
		 &+area_quad_2d(ptab(itab(n,5),[2,3]),ptab(itab(n,6),[2,3]),&
					   &ptab(itab(n,7),[2,3]),ptab(itab(n,8),[2,3]))&
		 &+area_quad_2d(ptab(itab(n,9),[2,3]),ptab(itab(n,10),[2,3]),&
					   &ptab(itab(n,11),[2,3]),ptab(itab(n,12),[2,3]))
	rst(2)=area_quad_2d(ptab(itab(n,1),[3,1]),ptab(itab(n,2),[3,1]),&
					   &ptab(itab(n,3),[3,1]),ptab(itab(n,4),[3,1]))&
		 &+area_quad_2d(ptab(itab(n,5),[3,1]),ptab(itab(n,6),[3,1]),&
					   &ptab(itab(n,7),[3,1]),ptab(itab(n,8),[3,1]))&
		 &+area_quad_2d(ptab(itab(n,9),[3,1]),ptab(itab(n,10),[3,1]),&
					   &ptab(itab(n,11),[3,1]),ptab(itab(n,12),[3,1]))
	rst(3)=area_quad_2d(ptab(itab(n,1),[1,2]),ptab(itab(n,2),[1,2]),&
					   &ptab(itab(n,3),[1,2]),ptab(itab(n,4),[1,2]))&
		 &+area_quad_2d(ptab(itab(n,5),[1,2]),ptab(itab(n,6),[1,2]),&
					   &ptab(itab(n,7),[1,2]),ptab(itab(n,8),[1,2]))&
		 &+area_quad_2d(ptab(itab(n,9),[1,2]),ptab(itab(n,10),[1,2]),&
					   &ptab(itab(n,11),[1,2]),ptab(itab(n,12),[1,2]))
end subroutine

subroutine subcell(P1,P2,P3,P4,P5,P6,P7,P8,n)
!transfer into the n-th sub-cell of a vertice cell
!input:
!	P1~8: [1*3] position of vertices
!	n: the index of the subcell (1~8)
	double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3)
	double precision ptab(27,3) !store the position of candidate points
	integer n
	integer itab(8,8) !table for selecting candidte points for different n
	ptab(1,:)=P1
	ptab(2,:)=P2
	ptab(3,:)=P3
	ptab(4,:)=P4
	ptab(5,:)=P5
	ptab(6,:)=P6
	ptab(7,:)=P7
	ptab(8,:)=P8
	ptab(9,:)=(P1+P2)/2
	ptab(10,:)=(P2+P3)/2
	ptab(11,:)=(P3+P4)/2
	ptab(12,:)=(P4+P1)/2
	ptab(13,:)=(P1+P5)/2
	ptab(14,:)=(P2+P6)/2
	ptab(15,:)=(P3+P7)/2
	ptab(16,:)=(P4+P8)/2
	ptab(17,:)=(P5+P6)/2
	ptab(18,:)=(P6+P7)/2
	ptab(19,:)=(P7+P8)/2
	ptab(20,:)=(P8+P5)/2
	ptab(21,:)=(P2+P3+P7+P6)/4
	ptab(22,:)=(P1+P5+P8+P4)/4
	ptab(23,:)=(P5+P6+P7+P8)/4
	ptab(24,:)=(P1+P4+P3+P2)/4
	ptab(25,:)=(P3+P4+P8+P7)/4
	ptab(26,:)=(P1+P2+P6+P5)/4
	ptab(27,:)=(P1+P2+P3+P4+P5+P6+P7+P8)/8
	itab(1,:)=[27,21,15,25,23,18, 7,19]
	itab(2,:)=[22,27,25,16,20,23,19, 8]
	itab(3,:)=[13,26,27,22, 5,17,23,20]
	itab(4,:)=[26,14,21,27,17, 6,18,23]
	itab(5,:)=[24,10, 3,11,27,21,15,25]
	itab(6,:)=[12,24,11, 4,22,27,25,16]
	itab(7,:)=[ 1, 9,24,12,13,26,27,22]
	itab(8,:)=[ 9, 2,10,24,26,14,21,27]
	P1=ptab(itab(n,1),:)
	P2=ptab(itab(n,2),:)
	P3=ptab(itab(n,3),:)
	P4=ptab(itab(n,4),:)
	P5=ptab(itab(n,5),:)
	P6=ptab(itab(n,6),:)
	P7=ptab(itab(n,7),:)
	P8=ptab(itab(n,8),:)
end subroutine

function volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
!calculate the volume of a hexahedron
!input:
!	P1~8: [1*3] position of vertices
!called functions, subrutines:
!	volume_tetra()
	double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3),volume_hex
	double precision volume_tetra
	volume_hex=volume_tetra(P1,P2,P3,P6)+volume_tetra(P1,P3,P4,P6)+volume_tetra(P1,P4,P5,P6)&
			 &+volume_tetra(P4,P5,P6,P7)+volume_tetra(P3,P4,P6,P7)+volume_tetra(P4,P5,P7,P8)
end function

function volume_tetra(P1,P2,P3,P4)
!calculate the volume of a tetrahedron
!input:
!	P1~4: [1*3] position of vertices
	double precision P1(3),P2(3),P3(3),P4(3),volume_tetra
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
	volume_tetra=((P2(1)*(y3z4-y4z3)-P3(1)*(y2z1-y4z2)+P4(1)*(y2z3-y3z2))&
				&-(P1(1)*(y3z4-y4z3)-P3(1)*(y1z4-y4z1)+P4(1)*(y1z3-y3z1))&
				&+(P1(1)*(y2z1-y4z2)-P2(1)*(y1z4-y4z1)+P4(1)*(y1z2-y2z1))&
				&-(P1(1)*(y2z3-y3z2)-P2(1)*(y1z3-y3z1)+P3(1)*(y1z2-y2z1)))/6.0
	volume_tetra=abs(volume_tetra)
end function

subroutine area_quad_3d(P1,P2,P3,P4,rst)
!calculate the 3d normal area of a quadrangle
!input:
!	P1~4: [1*3] position of points
!output:
!	rst: [1*3] result
!called functions:
!	dist_3d()
!note: sequence of input must be CCW (see from inside) to gurantee the correct sign
	double precision P1(3),P2(3),P3(3),P4(3),rst(3),area,temp,dist_3d
	area=dist_3d((P1+P2)/2,(P2+P3)/2)*&
		&dist_3d((P1+P2+P2+P3)/4,(P1+P2+P3+P4)/4)*4 !the normal area
	rst(1)=(P2(2)-P1(2))*(P4(3)-P1(3))-(P4(2)-P1(2))*(P2(3)-P1(3))
	rst(2)=(P2(3)-P1(3))*(P4(1)-P1(1))-(P4(3)-P1(3))*(P2(1)-P1(1))
	rst(3)=(P2(1)-P1(1))*(P4(2)-P1(2))-(P4(1)-P1(1))*(P2(2)-P1(2))
	temp=sqrt(rst(1)*rst(1)&
			&+rst(2)*rst(2)&
			&+rst(3)*rst(3))
	if(temp==0.0)then
		rst(1)=(P4(2)-P3(2))*(P2(3)-P3(3))-(P2(2)-P3(2))*(P4(3)-P3(3))
		rst(2)=(P4(3)-P3(3))*(P2(1)-P3(1))-(P2(3)-P3(3))*(P4(1)-P3(1))
		rst(3)=(P4(1)-P3(1))*(P2(2)-P3(2))-(P2(1)-P3(1))*(P4(2)-P3(2))
		temp=sqrt(rst(1)*rst(1)&
			&+rst(2)*rst(2)&
			&+rst(3)*rst(3))
	end if
	rst=rst/temp*area
end subroutine

function area_quad_2d(P1,P2,P3,P4)
!calculate the 2d area of a quadrangle
!input:
!	P1~4: [1*2] position of points
!called functions:
!	dist_2d()
!note: sequence of input must be CCW (see from inside) to gurantee the correct sign
	double precision P1(2),P2(2),P3(2),P4(2),area_quad_2d,dist_2d,a,b,c,p
	a=dist_2d((P1+P2)/2,(P2+P3)/2)
	b=dist_2d((P3+P4)/2,(P2+P3)/2)
	c=dist_2d((P1+P2)/2,(P3+P4)/2)
	p=(a+b+c)/2
	area_quad_2d=4*sqrt(p*(p-a)*(p-b)*(p-c))
	if((P2(1)-P1(1))*(P4(2)-P1(2))-(P4(1)-P1(1))*(P2(2)-P1(2))<0)then
		area_quad_2d=-area_quad_2d
	else
		if((P2(1)-P1(1))*(P4(2)-P1(2))-(P4(1)-P1(1))*(P2(2)-P1(2))==0)then
			if((P4(1)-P3(1))*(P2(2)-P3(2))-(P2(1)-P3(1))*(P4(2)-P3(2))<0)then
				area_quad_2d=-area_quad_2d
			end if
		end if
	end if
end function

function dist_3d(P1,P2)
!calculate the 3d distance between 2 points
!input:
!	P1,2: [1*3] position of points
	double precision P1(3),P2(3),dist_3d
	double precision temp1,temp2,temp3
	temp1=P1(1)-P2(1)
	temp2=P1(2)-P2(2)
	temp3=P1(3)-P2(3)
	dist_3d=sqrt(temp1*temp1+temp2*temp2+temp3*temp3)
end function

function dist_2d(P1,P2)
!calculate the 2d distance between 2 points
!input:
!	P1,2: [1*2] position of points
	double precision P1(2),P2(2),dist_2d
	double precision temp1,temp2
	temp1=P1(1)-P2(1)
	temp2=P1(2)-P2(2)
	dist_2d=sqrt(temp1*temp1+temp2*temp2)
end function

subroutine try(a)
	integer a(2,2)
	a(1,:)=[5,6]
end subroutine

