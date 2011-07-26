!author: Wang, Mianzhi
!ALE solver for compressible viscous 3d flow
!main program (FOR DEBUG ONLY, NEED TO BE IMPROVED)

!vertices denotaiton:	face/neighbor-cell denotation:	neighbor-cell of vertices:
!                                       5
!     4------8                     @------@                   4        8
!     |\     |\                    |\  2  |\
!     | 3------7                   | @------@                   3  \|    7
!     1 | ---5 |                 4 @ | ---@ | 3                  ---@---
!      \|     \|                    \|   1 \|                 1     |\ 5
!       2------6                     @------@
!                                       6                       2        6
!
!data structure:
!structural grid:
!	   cell list: element(cell#,16)=integer[vertice1~8,neighbor-cell1~8]
!	vertice list: vertice(vertice#,16)=integer[neighbor-vertice1~8,neighbor-cell1~8]
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

module globalvar
!declare the global variables
	integer,save,allocatable :: element(:,:)
	integer,save,allocatable :: vertice(:,:)
	integer,save,allocatable :: BClist_npns(:),BClist_velo(:)
	double precision,save,allocatable :: px(:),py(:),pz(:),vu(:),vv(:),vw(:)
	double precision,save,allocatable :: Pr(:),TE(:),Mass(:),IE(:),rr(:),Vol(:),KT(:)
	double precision,save,allocatable :: ppx(:,:),ppy(:,:),ppz(:,:)
	double precision,save :: idealgas,molmass,visc
	double precision,save :: refv(3)
end module

program main
!main program

!declare & initialize variables
use globalvar
real dt,tolerance
integer ncell,nvertice
integer MaxStep,MaxIt
real hx,hy,hz
integer counter1,counter2,step
integer,allocatable :: temp_int_vect1(:),temp_int_vect2(:)
l=13
m=13
n=13
ncell=l*m*n
nvertice=(l+1)*(m+1)*(n+1)
MaxStep=10
MaxIt=500
tolerance=1e-7
allocate(element(ncell,14))
allocate(vertice(nvertice,14))
allocate(px(nvertice))
allocate(py(nvertice))
allocate(pz(nvertice))
allocate(vu(nvertice))
allocate(vv(nvertice))
allocate(vw(nvertice))
allocate(Pr(ncell))
allocate(TE(ncell))
allocate(Mass(ncell))
allocate(IE(ncell))
allocate(rr(ncell))
allocate(Vol(ncell))
allocate(KT(ncell))
allocate(ppx(MaxStep,nvertice))
allocate(ppy(MaxStep,nvertice))
allocate(ppz(MaxStep,nvertice))
allocate(temp_int_vect1(nvertice))
allocate(temp_int_vect2(nvertice))

idealgas=8.314472
molmass=29e-3
visc=10 !proper viscosity to debug the viscous calculations
!visc=17.9e-6 !the viscosity of air

hx=0.1
hy=0.1
hz=0.1
dt=2e-4
refv=[0.5,1.0,0.0]

 call Grid_cube(l,m,n,hx,hy,hz)
write(*,*),'>> grid generated'

!set moveable grid (in this case the grid is stationary)
do step=1,MaxStep
	do i=1,nvertice
		ppx(step,i)=px(i)
		ppy(step,i)=py(i)
		ppz(step,i)=pz(i)
	end do
end do

!set IC
vu(:)=0.0
vv(:)=0.0
vw(:)=0.0
rr(:)=1.293
IE(:)=273*700
TE=IE
Pr(:)=101204
write(*,*),'>> IC set'

!find boundary vertices
counter1=1
counter2=1
do i=1,l+1
	do j=1,m+1
		do k=1,n+1
			if(i==l+1)then
				temp_int_vect1(counter1)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
				counter1=counter1+1
			else
				if(i==1)then
					temp_int_vect1(counter1)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
					counter1=counter1+1
				else
					if(j==m+1)then
						temp_int_vect1(counter1)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
						counter1=counter1+1
					else
						if(j==1)then
							temp_int_vect1(counter1)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
							counter1=counter1+1
						else
							if(k==n+1)then
								temp_int_vect2(counter2)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
								counter2=counter2+1
							else
								if(k==1)then
									temp_int_vect1(counter1)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
									counter1=counter1+1
								end if
							end if
						end if
					end if
				end if
			end if
		end do
	end do
end do
allocate(BClist_npns(counter1-1))
allocate(BClist_velo(counter2-1))
BClist_npns=temp_int_vect1(1:counter1-1)
BClist_velo=temp_int_vect2(1:counter2-1)
write(*,"(A,I8,A)"),'>>',counter1+counter2-2,' boundary vertices found'

do step=1,MaxStep
	call cv3dALE(dt,step,MaxIt,tolerance)
end do

write(*,*),'>> calculation done'

 call writedata()

!deallocate the allocatable variables
deallocate(element)
deallocate(vertice)
deallocate(px)
deallocate(py)
deallocate(pz)
deallocate(vu)
deallocate(vv)
deallocate(vw)
deallocate(Pr)
deallocate(TE)
deallocate(Mass)
deallocate(IE)
deallocate(rr)
deallocate(Vol)
deallocate(KT)
deallocate(ppx)
deallocate(ppy)
deallocate(ppz)
deallocate(BClist_npns)
deallocate(BClist_velo)
deallocate(temp_int_vect1)
deallocate(temp_int_vect2)

end program
