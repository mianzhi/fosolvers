!author: Wang, Mianzhi
!ALE solver for compressible viscous 3d flow

subroutine cv3dALE(dt,step,MaxIt,tolerance)
!input:
!	dt: lenghth of time step
!	step: current step number
!	MaxIt: maximum times of iteration for implicit step
!	tolerance: torlerance of implicit computation
!
!output:
!	no ouputs; the global data will be updated
!
!vertices denotaiton:	face/neighbor-cell denotation:	neighbor-cell of vertices:
!                                       5
!     4------8                     @------@                   4        8
!     |\     |\                    |\  2  |\
!     | 3------7                   | @------@                   3  \|    7
!     1 | ---5 |                 4 @ | ---@ | 3                  ---@---
!      \|     \|                    \|   1 \|                 1     |\ 5
!       2------6                     @------@
!                                       6                       2        6
!neighbor-vertices of vertices:
!              5
!            2
!             \|
!        4  ---@---  3
!              |\
!                1
!              6
!
!data structure:
!structural grid:
!	   cell list: element(cell#,14)=integer[vertice1~8,neighbor-cell1~6]
!	vertice list: vertice(vertice#,14)=integer[neighbor-vertice1~6,neighbor-cell1~8]
!
!primative variables:
!	 position of vertices: px,py,pz(vertice#)=double
!	 velocity at vertices: vu,vv,vw(vertice#)=double
!	    pressure of cells: Pr(cell#)=double
!	total energy of cells: TE(cell#)=double
!	        mass of cells: Mass(cell#)=double
!	        I.E. of cells: IE(cell#)=double
!        density of cells: rr(cell#)=double
!	      volume of cells: Vol(cell#)=double

!declare variables
!global data
use globalvar
!input
real dt,to
integer step,MaxIt
!functions
double precision volume_hex,area_quad_2d,voltrans
!auxiliary
integer ncell,nvertice !number of cells & vertices
double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3) !position of 8 auxiliary points
double precision,allocatable :: Massvertice(:),Volvertice(:),Momvertice(:,:)
							!mass/volume/momentum of vertices
double precision,allocatable :: rvu(:),rvv(:),rvw(:) !velocity of vertices
													 !(from vertices' defined position)
double precision temp_double1
double precision,allocatable :: temp_double_vect1(:),temp_double_vect2(:),temp_double_vect3(:),&
							   &temp_double_vect4(:),temp_double_vect5(:)
double precision,allocatable :: temp_double_mat1(:,:),temp_double_mat2(:,:),temp_double_mat3(:,:),&
							   &temp_double_mat4(:,:),temp_double_mat5(:,:),temp_double_mat6(:,:)
double precision,allocatable :: temp_double_3d1(:,:,:)
integer,allocatable :: temp_int_mat1(:,:),temp_int_mat2(:,:)
ncell=size(element,1) !number of cells
nvertice=size(vertice,1) !number of vertices
allocate(Massvertice(nvertice))
allocate(Volvertice(nvertice))
allocate(Momvertice(nvertice,3))
allocate(rvu(nvertice))
allocate(rvv(nvertice))
allocate(rvw(nvertice))

write(*,"(A,I4,A)"),'> STEP ',step,' STARTS'

!!!!!!!!!!!!!!!!!!!!!!!
! initial calculation !
!!!!!!!!!!!!!!!!!!!!!!!
!volume & mass of CELLs from points' position
!-------------------------------------
do i=1,ncell
	j=element(i,1) !index of a vertice
	P1=[px(j),py(j),pz(j)]
	j=element(i,2)
	P2=[px(j),py(j),pz(j)]
	j=element(i,3)
	P3=[px(j),py(j),pz(j)]
	j=element(i,4)
	P4=[px(j),py(j),pz(j)]
	j=element(i,5)
	P5=[px(j),py(j),pz(j)]
	j=element(i,6)
	P6=[px(j),py(j),pz(j)]
	j=element(i,7)
	P7=[px(j),py(j),pz(j)]
	j=element(i,8)
	P8=[px(j),py(j),pz(j)]
	Vol(i)=volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
	Mass(i)=Vol(i)*rr(i)
end do
!mass of VERTICEs from volume & density of the 8 neighbor cells
!--------------------------------------------------------------
do i=1,nvertice
	Massvertice(i)=0.0
	Volvertice(i)=0.0
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		do l=1,8
			k=vertice(i,l+6) !index of a cell
			j=element(k,1)
			P1=[px(j),py(j),pz(j)]
			j=element(k,2)
			P2=[px(j),py(j),pz(j)]
			j=element(k,3)
			P3=[px(j),py(j),pz(j)]
			j=element(k,4)
			P4=[px(j),py(j),pz(j)]
			j=element(k,5)
			P5=[px(j),py(j),pz(j)]
			j=element(k,6)
			P6=[px(j),py(j),pz(j)]
			j=element(k,7)
			P7=[px(j),py(j),pz(j)]
			j=element(k,8)
			P8=[px(j),py(j),pz(j)]
			call subcell(P1,P2,P3,P4,P5,P6,P7,P8,l)
			temp_double1=volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
			Volvertice(i)=Volvertice(i)+temp_double1
			Massvertice(i)=Massvertice(i)+temp_double1*rr(k)
		end do
	end if
end do
!total energy of CELLs from internal engergy of cells & velocity of vertices (KE)
!--------------------------------------------------------------------------------
do i=1,ncell
	temp_double1=0
	do j=1,8
		k=element(i,j)
		temp_double1=vu(k)*vu(k)+vv(k)*vv(k)+vw(k)*vw(k)+temp_double1
	end do
	TE(i)=IE(i)+temp_double1/16
end do
write(*,*),'initial computation done'

!!!!!!!!!!!!!!
! Lagrangian !
!!!!!!!!!!!!!!
!phase 1, part a
!===============
!apply BC
!--------
 call BC_npns(step,dt)
 call BC_velo(refv)
!TODO: create other BC
!velocity of VERTICEs from body force & stress (normal & shear)
!--------------------------------------------------------------
allocate(temp_double_vect1(3)) !projection area of sub-cell in 3 directions
allocate(temp_double_mat1(8,3)) !velocity at 8 vertices
allocate(temp_double_mat2(ncell,3)) !the upper-triangular elements of stress tensor
do i=1,ncell
	j=element(i,1)
	P1=[px(j),py(j),pz(j)]
	temp_double_mat1(1,:)=[vu(j),vv(j),vw(j)]
	j=element(i,2)
	P2=[px(j),py(j),pz(j)]
	temp_double_mat1(2,:)=[vu(j),vv(j),vw(j)]
	j=element(i,3)
	P3=[px(j),py(j),pz(j)]
	temp_double_mat1(3,:)=[vu(j),vv(j),vw(j)]
	j=element(i,4)
	P4=[px(j),py(j),pz(j)]
	temp_double_mat1(4,:)=[vu(j),vv(j),vw(j)]
	j=element(i,5)
	P5=[px(j),py(j),pz(j)]
	temp_double_mat1(5,:)=[vu(j),vv(j),vw(j)]
	j=element(i,6)
	P6=[px(j),py(j),pz(j)]
	temp_double_mat1(6,:)=[vu(j),vv(j),vw(j)]
	j=element(i,7)
	P7=[px(j),py(j),pz(j)]
	temp_double_mat1(7,:)=[vu(j),vv(j),vw(j)]
	j=element(i,8)
	P8=[px(j),py(j),pz(j)]
	temp_double_mat1(8,:)=[vu(j),vv(j),vw(j)]
	call shearstress(P1,P2,P3,P4,P5,P6,P7,P8,&
					&temp_double_mat1(1,:),temp_double_mat1(2,:),&
					&temp_double_mat1(3,:),temp_double_mat1(4,:),&
					&temp_double_mat1(5,:),temp_double_mat1(6,:),&
					&temp_double_mat1(7,:),temp_double_mat1(8,:),temp_double_mat2(i,:))
end do
do i=1,nvertice
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
	 	do l=1,8
	 		temp_double_vect1(:)=0 !store the projection area of the subcell
	 		k=vertice(i,l+6)	!index of a neighbor cell
			j=element(k,1)
			P1=[px(j),py(j),pz(j)]
			j=element(k,2)
			P2=[px(j),py(j),pz(j)]
			j=element(k,3)
			P3=[px(j),py(j),pz(j)]
			j=element(k,4)
			P4=[px(j),py(j),pz(j)]
			j=element(k,5)
			P5=[px(j),py(j),pz(j)]
			j=element(k,6)
			P6=[px(j),py(j),pz(j)]
			j=element(k,7)
			P7=[px(j),py(j),pz(j)]
			j=element(k,8)
			P8=[px(j),py(j),pz(j)]
			call projarea_subcell(P1,P2,P3,P4,P5,P6,P7,P8,l,temp_double_vect1) !projection area
			vu(i)=vu(i)+dt/Massvertice(i)*(Pr(k)*temp_double_vect1(1)&
						 &-temp_double_mat2(k,1)*temp_double_vect1(2)&
						 &-temp_double_mat2(k,2)*temp_double_vect1(3))
			vv(i)=vv(i)+dt/Massvertice(i)*(Pr(k)*temp_double_vect1(2)&
						 &-temp_double_mat2(k,1)*temp_double_vect1(1)&
						 &-temp_double_mat2(k,3)*temp_double_vect1(3))
			vw(i)=vw(i)+dt/Massvertice(i)*(Pr(k)*temp_double_vect1(3)&
						 &-temp_double_mat2(k,2)*temp_double_vect1(1)&
						 &-temp_double_mat2(k,3)*temp_double_vect1(2))
	 	end do
	end if
end do
!total energy of CELLs from work done by shear stress & body forces (normal stress NOT inculeded)
!------------------------------------------------------------------------------------------------
allocate(temp_double_mat3(8,3)) !the position of 8 vertices
allocate(temp_int_mat1(6,4)) !table of suface's vertices
temp_int_mat1(1,:)=[2,3,7,6]
temp_int_mat1(2,:)=[1,5,8,4]
temp_int_mat1(3,:)=[5,6,7,8]
temp_int_mat1(4,:)=[1,4,3,2]
temp_int_mat1(5,:)=[3,4,8,7]
temp_int_mat1(6,:)=[1,2,6,5]
do i=1,ncell
	do j=1,8
		temp_double_mat3(j,:)=[px(element(i,j)),py(element(i,j)),pz(element(i,j))]
		temp_double_mat1(j,:)=[vu(element(i,j)),vv(element(i,j)),vw(element(i,j))]
	end do
	do j=1,6
		call area_quad_3d(temp_double_mat3(temp_int_mat1(j,1),:),&
						 &temp_double_mat3(temp_int_mat1(j,2),:),&
						 &temp_double_mat3(temp_int_mat1(j,3),:),&
						 &temp_double_mat3(temp_int_mat1(j,4),:),temp_double_vect1)
		if(element(i,j+8)/=0)then
			TE(i)=TE(i)-dt/Mass(i)*((temp_double_mat1(temp_int_mat1(j,1),1)& !in x
								   &+temp_double_mat1(temp_int_mat1(j,2),1)&
								   &+temp_double_mat1(temp_int_mat1(j,3),1)&
								   &+temp_double_mat1(temp_int_mat1(j,4),1))/4*&
								   &(temp_double_vect1(2)*(temp_double_mat2(i,1)&
														 &+temp_double_mat2(element(i,j+8),1))/2&
								   &+temp_double_vect1(3)*(temp_double_mat2(i,2)&
														 &+temp_double_mat2(element(i,j+8),2))/2)&
								  &+(temp_double_mat1(temp_int_mat1(j,1),2)& !in y
								   &+temp_double_mat1(temp_int_mat1(j,2),2)&
								   &+temp_double_mat1(temp_int_mat1(j,3),2)&
								   &+temp_double_mat1(temp_int_mat1(j,4),2))/4*&
								   &(temp_double_vect1(1)*(temp_double_mat2(i,1)&
														 &+temp_double_mat2(element(i,j+8),1))/2&
								   &+temp_double_vect1(3)*(temp_double_mat2(i,3)&
														 &+temp_double_mat2(element(i,j+8),3))/2)&
								  &+(temp_double_mat1(temp_int_mat1(j,1),3)& !in z
								   &+temp_double_mat1(temp_int_mat1(j,2),3)&
								   &+temp_double_mat1(temp_int_mat1(j,3),3)&
								   &+temp_double_mat1(temp_int_mat1(j,4),3))/4*&
								   &(temp_double_vect1(1)*(temp_double_mat2(i,2)&
														 &+temp_double_mat2(element(i,j+8),2))/2&
								   &+temp_double_vect1(2)*(temp_double_mat2(i,3)&
														 &+temp_double_mat2(element(i,j+8),3))/2))
		else
			TE(i)=TE(i)-dt/Mass(i)*((temp_double_mat1(temp_int_mat1(j,1),1)&
								   &+temp_double_mat1(temp_int_mat1(j,2),1)&
								   &+temp_double_mat1(temp_int_mat1(j,3),1)&
								   &+temp_double_mat1(temp_int_mat1(j,4),1))/4*&
								   &(temp_double_vect1(2)*temp_double_mat2(i,1)&
								   &+temp_double_vect1(3)*temp_double_mat2(i,2))& !in x
								  &+(temp_double_mat1(temp_int_mat1(j,1),2)&
								   &+temp_double_mat1(temp_int_mat1(j,2),2)&
								   &+temp_double_mat1(temp_int_mat1(j,3),2)&
								   &+temp_double_mat1(temp_int_mat1(j,4),2))/4*&
								   &(temp_double_vect1(1)*temp_double_mat2(i,1)&
								   &+temp_double_vect1(3)*temp_double_mat2(i,3))& !in y
								  &+(temp_double_mat1(temp_int_mat1(j,1),3)&
								   &+temp_double_mat1(temp_int_mat1(j,2),3)&
								   &+temp_double_mat1(temp_int_mat1(j,3),3)&
								   &+temp_double_mat1(temp_int_mat1(j,4),3))/4*&
								   &(temp_double_vect1(1)*temp_double_mat2(i,2)&
								   &+temp_double_vect1(2)*temp_double_mat2(i,3))) !in z
		end if
	end do
end do
deallocate(temp_double_vect1)
deallocate(temp_double_mat1)
deallocate(temp_double_mat2)
deallocate(temp_double_mat3)
deallocate(temp_int_mat1)
write(*,*),'Lagrangian phase 1, part a done'

!phase 2, implicit
!=================
write(*,'(a)'),'starting iteration for implicit computation'
allocate(temp_double_vect1(ncell)) !store the temporary volume
allocate(temp_double_vect2(ncell)) !store the temporary density
allocate(temp_double_vect3(ncell)) !store the temporary internal energy
allocate(temp_double_vect4(ncell)) !store the temporary delta P
allocate(temp_double_vect5(3)) !projection area in 3 directions
allocate(temp_double_mat1(nvertice,3)) !store the temporary position of vertices
do l=1,MaxIt
	write(*,"(A,$)"),'='
	!position of VERTICEs form new velocity of vertices
	!--------------------------------------------------
	do i=1,nvertice
		temp_double_mat1(i,1)=px(i)+vu(i)*dt
		temp_double_mat1(i,2)=py(i)+vv(i)*dt
		temp_double_mat1(i,3)=pz(i)+vw(i)*dt
	end do
	!volume of CELLs from new position of vertices
	!---------------------------------------------
	do i=1,ncell
		j=element(i,1) !index of a vertice
		P1=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,2)
		P2=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,3)
		P3=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,4)
		P4=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,5)
		P5=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,6)
		P6=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,7)
		P7=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		j=element(i,8)
		P8=[temp_double_mat1(j,1),temp_double_mat1(j,2),temp_double_mat1(j,3)]
		temp_double_vect1(i)=volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
	end do
	!density of CELLs from new volume of cells
	!-----------------------------------------
	temp_double_vect2=rr*Vol/temp_double_vect1
	!internal energy of CELLs from multiple cell parameters
	!------------------------------------------------------
	temp_double_vect3=IE+Pr/rr*(1.0-temp_double_vect1/Vol)
	!pressure corrector of CELLs from previous results
	!-------------------------------------------------
	KT=temp_double_vect3/700
	!TODO:generate a non-uniform, non-linear relationship between internal energy & temperature
	temp_double_vect4=rr*idealgas/molmass*KT-Pr
	Pr=Pr+temp_double_vect4
	!velocity of VERTICEs due to change of normal stress (pressure)
	!--------------------------------------------------------------
	do i=1,nvertice
		if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
		 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		 	temp_double_vect5(:)=0
		 	do m=1,8
		 		temp_double_vect5(:)=0 !store the projection area of the subcell
		 		k=vertice(i,m+6)	!index of a neighbor cell
				j=element(k,1)
				P1=[px(j),py(j),pz(j)]
				j=element(k,2)
				P2=[px(j),py(j),pz(j)]
				j=element(k,3)
				P3=[px(j),py(j),pz(j)]
				j=element(k,4)
				P4=[px(j),py(j),pz(j)]
				j=element(k,5)
				P5=[px(j),py(j),pz(j)]
				j=element(k,6)
				P6=[px(j),py(j),pz(j)]
				j=element(k,7)
				P7=[px(j),py(j),pz(j)]
				j=element(k,8)
				P8=[px(j),py(j),pz(j)]
				call projarea_subcell(P1,P2,P3,P4,P5,P6,P7,P8,m,temp_double_vect5)
				vu(i)=vu(i)+dt/Massvertice(i)*temp_double_vect4(k)*temp_double_vect5(1)
				vv(i)=vv(i)+dt/Massvertice(i)*temp_double_vect4(k)*temp_double_vect5(2)
				vw(i)=vw(i)+dt/Massvertice(i)*temp_double_vect4(k)*temp_double_vect5(3)
		 	end do
		end if
	end do
	!apply BC
	!--------
	call BC_npns(step,dt)
	call BC_velo(refv)
	!check error
	!-----------
	if(sum(abs(temp_double_vect4/Pr))/size(Pr)<=tolerance)then
	  write(*,'(a)'),''
		write(*,"(A,I2,A)"),'implicit computation (phase 2) done after ',l,' times of iteration'
		exit
	end if
end do
deallocate(temp_double_vect1)
deallocate(temp_double_vect2)
deallocate(temp_double_vect3)
deallocate(temp_double_vect4)
deallocate(temp_double_vect5)
deallocate(temp_double_mat1)

!phase 1, part b
!===============
!total energy of CELLs from work done by normal stress only
!(calculate the surface (of cells) pressure first)
!----------------------------------------------------------
allocate(temp_double_mat1(2,3)) !store the projection area & surface velocity in x,y,z
allocate(temp_double_mat2(8,3)) !store the position of points
allocate(temp_double_mat3(8,3)) !store the velocity at points
allocate(temp_int_mat1(6,4)) !store the table of surface point
temp_int_mat1(1,:)=[2,3,7,6]
temp_int_mat1(2,:)=[1,5,8,4]
temp_int_mat1(3,:)=[5,6,7,8]
temp_int_mat1(4,:)=[1,4,3,2]
temp_int_mat1(5,:)=[3,4,8,7]
temp_int_mat1(6,:)=[1,2,6,5]
do i=1,ncell
	temp_double1=0 !to store the pressure at cell surface
	do j=1,8
		k=element(i,j)
		temp_double_mat2(j,:)=[px(k),py(k),pz(k)]
		temp_double_mat3(j,:)=[vu(k),vv(k),vw(k)]
	end do
	do j=1,6
		!surface pressure
		if(element(i,j+8)/=0)then
			temp_double1=(Pr(i)*Mass(element(i,j+8))+Pr(element(i,j+8))*Mass(i))&
								&/(Mass(element(i,j+8))+Mass(i))
		else
			temp_double1=Pr(i)
		end if
		!surface prejection area
		temp_double_mat1(1,1)=area_quad_2d(temp_double_mat2(temp_int_mat1(j,1),[2,3]),&
										  &temp_double_mat2(temp_int_mat1(j,2),[2,3]),&
										  &temp_double_mat2(temp_int_mat1(j,3),[2,3]),&
										  &temp_double_mat2(temp_int_mat1(j,4),[2,3]))
		temp_double_mat1(1,2)=area_quad_2d(temp_double_mat2(temp_int_mat1(j,1),[3,1]),&
										  &temp_double_mat2(temp_int_mat1(j,2),[3,1]),&
										  &temp_double_mat2(temp_int_mat1(j,3),[3,1]),&
										  &temp_double_mat2(temp_int_mat1(j,4),[3,1]))
		temp_double_mat1(1,3)=area_quad_2d(temp_double_mat2(temp_int_mat1(j,1),[1,2]),&
										  &temp_double_mat2(temp_int_mat1(j,2),[1,2]),&
										  &temp_double_mat2(temp_int_mat1(j,3),[1,2]),&
										  &temp_double_mat2(temp_int_mat1(j,4),[1,2]))
		!surface velocity
		temp_double_mat1(2,1)=(temp_double_mat3(temp_int_mat1(j,1),1)&
							 &+temp_double_mat3(temp_int_mat1(j,2),1)&
							 &+temp_double_mat3(temp_int_mat1(j,3),1)&
							 &+temp_double_mat3(temp_int_mat1(j,4),1))/4
		temp_double_mat1(2,2)=(temp_double_mat3(temp_int_mat1(j,1),2)&
							 &+temp_double_mat3(temp_int_mat1(j,2),2)&
							 &+temp_double_mat3(temp_int_mat1(j,3),2)&
							 &+temp_double_mat3(temp_int_mat1(j,4),2))/4
		temp_double_mat1(2,3)=(temp_double_mat3(temp_int_mat1(j,1),3)&
							 &+temp_double_mat3(temp_int_mat1(j,2),3)&
							 &+temp_double_mat3(temp_int_mat1(j,3),3)&
							 &+temp_double_mat3(temp_int_mat1(j,4),3))/4
		TE(i)=TE(i)+dt/Mass(i)*(temp_double1*sum(product(temp_double_mat1,dim=1)))
	end do
end do
deallocate(temp_double_mat1)
deallocate(temp_double_mat2)
deallocate(temp_double_mat3)
deallocate(temp_int_mat1)
!position of VERTICEs from velocity
!----------------------------------
px=px+vu*dt
py=py+vv*dt
pz=pz+vw*dt
!volume of CELLs from points' position
!-------------------------------------
do i=1,ncell
	j=element(i,1) !index of a vertice
	P1=[px(j),py(j),pz(j)]
	j=element(i,2)
	P2=[px(j),py(j),pz(j)]
	j=element(i,3)
	P3=[px(j),py(j),pz(j)]
	j=element(i,4)
	P4=[px(j),py(j),pz(j)]
	j=element(i,5)
	P5=[px(j),py(j),pz(j)]
	j=element(i,6)
	P6=[px(j),py(j),pz(j)]
	j=element(i,7)
	P7=[px(j),py(j),pz(j)]
	j=element(i,8)
	P8=[px(j),py(j),pz(j)]
	Vol(i)=volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
end do
!density of the CELLs from mass of the cell and its new volume
!-------------------------------------------------------------
rr=Mass/Vol
!apply BC
!--------
 call BC_velo(refv)
 call BC_npns(step,dt)
write(*,*),'Lagrangian phase 1, part b done'

!!!!!!!!!!!!
! Eulerian !
!!!!!!!!!!!!
!rezone
!======
!reference velocity of the VERTICEs
!----------------------------------
rvu=(ppx(step,:)-px)/dt
rvv=(ppy(step,:)-py)/dt
rvw=(ppz(step,:)-pz)/dt
!volume exchange between CELLs due to rezoning
!---------------------------------------------
allocate(temp_double_mat1(ncell,6)) !to store the delta V
allocate(temp_double_mat2(8,3)) !the table of vertices position before rezoning
allocate(temp_double_mat3(8,3)) !the table of vertices position after rezoning 
allocate(temp_int_mat1(6,4)) !the table of vertices of different surfaces
temp_double_mat1(:,:)=0.0
temp_int_mat1(1,:)=[2,3,7,6]
temp_int_mat1(2,:)=[1,5,8,4]
temp_int_mat1(3,:)=[5,6,7,8]
temp_int_mat1(4,:)=[1,4,3,2]
temp_int_mat1(5,:)=[3,4,8,7]
temp_int_mat1(6,:)=[1,2,6,5]
do i=1,ncell
	do j=1,8
		k=element(i,j)
		temp_double_mat2(j,:)=[px(k),py(k),pz(k)]
		temp_double_mat3(j,:)=[ppx(step,k),ppy(step,k),ppz(step,k)]
	end do
	do j=1,6
		if(temp_double_mat1(i,j)==0.0)then
			temp_double_mat1(i,j)=voltrans(temp_double_mat2(temp_int_mat1(j,1),:),&
										  &temp_double_mat2(temp_int_mat1(j,2),:),&
										  &temp_double_mat2(temp_int_mat1(j,3),:),&
										  &temp_double_mat2(temp_int_mat1(j,4),:),&
										  &temp_double_mat3(temp_int_mat1(j,1),:),&
										  &temp_double_mat3(temp_int_mat1(j,2),:),&
										  &temp_double_mat3(temp_int_mat1(j,3),:),&
										  &temp_double_mat3(temp_int_mat1(j,4),:))
			if(element(i,j+8)/=0)then
				if(mod(j,2)/=0)then
					temp_double_mat1(element(i,j+8),j+1)=-temp_double_mat1(i,j)
				else
					temp_double_mat1(element(i,j+8),j-1)=-temp_double_mat1(i,j)
				end if
			end if
		end if
	end do
end do
deallocate(temp_double_mat2)
deallocate(temp_double_mat3)
deallocate(temp_int_mat1)
!mass exchange & total energy (M*E) exchange between CELLs due to rezoning
!-------------------------------------------------------------------------
allocate(temp_double_mat2(ncell,6)) !to store delta M
allocate(temp_double_mat3(ncell,6)) !to store delta M*E
temp_double_mat2(:,:)=0.0
temp_double_mat3(:,:)=0.0
do i=1,ncell
	do j=1,6
		if(temp_double_mat2(i,j)==0.0)then
			temp_double1=temp_double_mat1(i,j) !read the delta V
			if(temp_double1<=0)then
				temp_double_mat2(i,j)=temp_double_mat1(i,j)*rr(i)
				temp_double_mat3(i,j)=temp_double_mat2(i,j)*TE(i)
			else
				if(element(i,j+8)/=0)then
					temp_double_mat2(i,j)=temp_double_mat1(i,j)*rr(element(i,j+8))
					temp_double_mat3(i,j)=temp_double_mat2(i,j)*TE(element(i,j+8))
				else
					temp_double_mat2(i,j)=temp_double_mat1(i,j)*rr(i)
					temp_double_mat3(i,j)=temp_double_mat2(i,j)*TE(i)
				end if
			end if
			if(element(i,j+8)/=0)then
				if(mod(j,2)/=0)then
					temp_double_mat2(element(i,j+8),j+1)=-temp_double_mat2(i,j)
					temp_double_mat3(element(i,j+8),j+1)=-temp_double_mat3(i,j)
				else
					temp_double_mat2(element(i,j+8),j-1)=-temp_double_mat2(i,j)
					temp_double_mat3(element(i,j+8),j-1)=-temp_double_mat3(i,j)
				end if
			end if
		end if
	end do
end do
!volume exchange between VERTICEs due to rezoning
!------------------------------------------------
allocate(temp_double_mat4(nvertice,6)) !store the delta Vvertice
allocate(temp_double_mat5(8,3)) !the table of vertice's position before rezoning
allocate(temp_double_mat6(8,3)) !the table of vertice's position after rezoning
allocate(temp_int_mat1(6,4)) !store the table of neighbor cells
							 !	is also the table of vertices of a surface
temp_int_mat1(1,:)=[2,3,7,6]
temp_int_mat1(2,:)=[1,5,8,4]
temp_int_mat1(3,:)=[5,6,7,8]
temp_int_mat1(4,:)=[1,4,3,2]
temp_int_mat1(5,:)=[3,4,8,7]
temp_int_mat1(6,:)=[1,2,6,5]
temp_double_mat4(:,:)=0.0
do i=1,nvertice
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		do j=1,6
			if(temp_double_mat4(i,j)==0.0)then
				do k=1,4
					l=vertice(i,temp_int_mat1(j,k)+6) !index of a cell
					do m=1,8
						temp_double_mat5(m,:)=[px(element(l,m)),&
											  &py(element(l,m)),&
											  &pz(element(l,m))]
						temp_double_mat6(m,:)=[ppx(step,element(l,m)),&
											  &ppy(step,element(l,m)),&
											  &ppz(step,element(l,m))]
					end do
					call subcell(temp_double_mat5(1,:),&
								&temp_double_mat5(2,:),&
								&temp_double_mat5(3,:),&
								&temp_double_mat5(4,:),&
								&temp_double_mat5(5,:),&
								&temp_double_mat5(6,:),&
								&temp_double_mat5(7,:),&
								&temp_double_mat5(8,:),&
								&temp_int_mat1(j,k))
					call subcell(temp_double_mat6(1,:),&
								&temp_double_mat6(2,:),&
								&temp_double_mat6(3,:),&
								&temp_double_mat6(4,:),&
								&temp_double_mat6(5,:),&
								&temp_double_mat6(6,:),&
								&temp_double_mat6(7,:),&
								&temp_double_mat6(8,:),&
								&temp_int_mat1(j,k))
					temp_double_mat4(i,j)=temp_double_mat4(i,j)&
										 &+voltrans(temp_double_mat5(temp_int_mat1(j,1),:),&
										 		   &temp_double_mat5(temp_int_mat1(j,2),:),&
										 		   &temp_double_mat5(temp_int_mat1(j,3),:),&
										 		   &temp_double_mat5(temp_int_mat1(j,4),:),&
										 		   &temp_double_mat6(temp_int_mat1(j,1),:),&
										 		   &temp_double_mat6(temp_int_mat1(j,2),:),&
										 		   &temp_double_mat6(temp_int_mat1(j,3),:),&
										 		   &temp_double_mat6(temp_int_mat1(j,4),:))
				end do
				if(mod(j,2)/=0)then
					temp_double_mat4(vertice(i,j),j+1)=-temp_double_mat4(i,j)
				else
					temp_double_mat4(vertice(i,j),j-1)=-temp_double_mat4(i,j)
				end if
			end if
		end do
	end if
end do
deallocate(temp_double_mat5)
deallocate(temp_double_mat6)
deallocate(temp_int_mat1)
!momentum exchange between VERTICEs due to rezoning
!--------------------------------------------------
allocate(temp_double_3d1(nvertice,6,3)) !store the momentum exchange between vertices
temp_double_3d1(:,:,:)=0
do i=1,nvertice
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		do j=1,6
			if(temp_double_3d1(i,j,1)==0.0)then
				temp_double1=temp_double_mat4(i,j) !read the delta V of vertice
				if(temp_double1<=0)then
					temp_double_3d1(i,j,:)=temp_double1*Massvertice(i)/Volvertice(i)
					temp_double_3d1(i,j,1)=temp_double_3d1(i,j,1)*vu(i)
					temp_double_3d1(i,j,2)=temp_double_3d1(i,j,2)*vv(i)
					temp_double_3d1(i,j,3)=temp_double_3d1(i,j,3)*vw(i)
				else
					temp_double_3d1(i,j,:)=temp_double1*Massvertice(vertice(i,j))/Volvertice(i)
					temp_double_3d1(i,j,1)=temp_double_3d1(i,j,1)*vu(vertice(i,j))
					temp_double_3d1(i,j,2)=temp_double_3d1(i,j,2)*vv(vertice(i,j))
					temp_double_3d1(i,j,3)=temp_double_3d1(i,j,3)*vw(vertice(i,j))
				end if
				if(mod(j,2)/=0)then
					temp_double_3d1(vertice(i,j),j+1,:)=-temp_double_3d1(i,j,:)
				else
					temp_double_3d1(vertice(i,j),j-1,:)=-temp_double_3d1(i,j,:)
				end if
			end if
		end do
	end if
end do
!new positions of the VERTICEs
!-----------------------------
px=ppx(step,:)
py=ppy(step,:)
pz=ppz(step,:)
!new volume of CELLs
!-------------------
do i=1,ncell
	j=element(i,1) !index of a vertice
	P1=[px(j),py(j),pz(j)]
	j=element(i,2)
	P2=[px(j),py(j),pz(j)]
	j=element(i,3)
	P3=[px(j),py(j),pz(j)]
	j=element(i,4)
	P4=[px(j),py(j),pz(j)]
	j=element(i,5)
	P5=[px(j),py(j),pz(j)]
	j=element(i,6)
	P6=[px(j),py(j),pz(j)]
	j=element(i,7)
	P7=[px(j),py(j),pz(j)]
	j=element(i,8)
	P8=[px(j),py(j),pz(j)]
	Vol(i)=volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
end do
!new mass and total energy of CELLs
!-----------------
do i=1,ncell
	TE(i)=TE(i)*Mass(i)
	Mass(i)=Mass(i)+sum(temp_double_mat2(i,:))
	TE(i)=(TE(i)+sum(temp_double_mat3(i,:)))/Mass(i)
end do
!new density of CELLs from new volume & new mass of cells
!--------------------------------------------------------
rr=Mass/Vol
!new momentum of VERTICEs
!---------------------------------
Momvertice(:,:)=0.0
do i=1,nvertice
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		Momvertice(i,:)=Massvertice(i)*[vu(i),vv(i),vw(i)]
		Momvertice(i,1)=Momvertice(i,1)+sum(temp_double_3d1(i,:,1))
		Momvertice(i,2)=Momvertice(i,2)+sum(temp_double_3d1(i,:,2))
		Momvertice(i,3)=Momvertice(i,3)+sum(temp_double_3d1(i,:,3))
	end if
end do
!new mass of VERTICEs
!--------------------
do i=1,nvertice
	Massvertice(i)=0.0
	Volvertice(i)=0.0
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		do l=1,8
			k=vertice(i,l+6) !index of a cell
			j=element(k,1)
			P1=[px(j),py(j),pz(j)]
			j=element(k,2)
			P2=[px(j),py(j),pz(j)]
			j=element(k,3)
			P3=[px(j),py(j),pz(j)]
			j=element(k,4)
			P4=[px(j),py(j),pz(j)]
			j=element(k,5)
			P5=[px(j),py(j),pz(j)]
			j=element(k,6)
			P6=[px(j),py(j),pz(j)]
			j=element(k,7)
			P7=[px(j),py(j),pz(j)]
			j=element(k,8)
			P8=[px(j),py(j),pz(j)]
			call subcell(P1,P2,P3,P4,P5,P6,P7,P8,l)
			temp_double1=volume_hex(P1,P2,P3,P4,P5,P6,P7,P8)
			Volvertice(i)=Volvertice(i)+temp_double1
			Massvertice(i)=Massvertice(i)+temp_double1*rr(k)
		end do
	end if
end do
!new momentum of VERTICEs
!------------------------
do i=1,nvertice
	if(vertice(i,7)/=0.and.vertice(i,8)/=0.and.vertice(i,9)/=0.and.vertice(i,10)/=0.and.&
	 &vertice(i,11)/=0.and.vertice(i,12)/=0.and.vertice(i,13)/=0.and.vertice(i,14)/=0)then
		vu(i)=Momvertice(i,1)/Massvertice(i)
		vv(i)=Momvertice(i,2)/Massvertice(i)
		vw(i)=Momvertice(i,3)/Massvertice(i)
	end if
end do
!deallocate temporary data
!-------------------------
deallocate(temp_double_mat1)
deallocate(temp_double_mat2)
deallocate(temp_double_mat3)
deallocate(temp_double_mat4)
deallocate(temp_double_3d1)
write(*,*),'Eulerian computation done'

!!!!!!!!!!!!!
! Auxiliary !
!!!!!!!!!!!!!
!internal energy of CELLs from total energy of cells & velocity at vertices (KE)
!-------------------------------------------------------------------------------
do i=1,ncell
	temp_double1=0
	do j=1,8
		k=element(i,j)
		temp_double1=vu(k)*vu(k)+vv(k)*vv(k)+vw(k)*vw(k)+temp_double1
	end do
	IE(i)=TE(i)-temp_double1/16
end do
!pressure of CELLs from density & internal energy of cells
!---------------------------------------------------------
KT=IE/700
!TODO:generate a non-uniform, non-linear relationship between internal energy & temperature
Pr=rr*idealgas/molmass*KT
write(*,*),'Auxiliary computaion done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate Temporary Variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(Massvertice)
deallocate(Volvertice)
deallocate(Momvertice)
deallocate(rvu)
deallocate(rvv)
deallocate(rvw)
end subroutine
