!author: Wang, Mianzhi
!subroutines:
!	Grid_cube(),BC_npns(),BC_velo()

subroutine Grid_cube(l,m,n,hx,hy,hz)
!generate grid for cubic region
!input:
!	l,m,n: number of steps in 3 directions
!	hx,hy,hz: step size
use globalvar,only:element,vertice,px,py,pz
integer l,m,n
real hx,hy,hz
element(:,:)=0
vertice(:,:)=0
px(:)=0
py(:)=0
pz(:)=0
do i=1,l+1 !generate points
	do j=1,m+1
		do k=1,n+1
			px((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k)=hx*(i-1)
			py((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k)=hy*(j-1)
			pz((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k)=hz*(k-1)
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,7)=(i-2)*m*n+(j-2)*n+k-1
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,8)=(i-1)*m*n+(j-2)*n+k-1
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,9)=(i-1)*m*n+(j-2)*n+k
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,10)=(i-2)*m*n+(j-2)*n+k
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,11)=(i-2)*m*n+(j-1)*n+k-1
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,12)=(i-1)*m*n+(j-1)*n+k-1
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,13)=(i-1)*m*n+(j-1)*n+k
			vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,14)=(i-2)*m*n+(j-1)*n+k
			if(i/=l+1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,1)=(i)*(m+1)*(n+1)+(j-1)*(n+1)+k
			end if
			if(i/=1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,2)=(i-2)*(m+1)*(n+1)+(j-1)*(n+1)+k
			end if
			if(j/=m+1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,3)=(i-1)*(m+1)*(n+1)+(j)*(n+1)+k
			end if
			if(j/=1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,4)=(i-1)*(m+1)*(n+1)+(j-2)*(n+1)+k
			end if
			if(k/=n+1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,5)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k+1
			end if
			if(k/=1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,6)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k-1
			end if
		end do
	end do
end do
do i=1,l+1
	do j=1,m+1
		do k=1,n+1
			if(i==l+1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,8)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,9)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,12)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,13)=0
			end if
			if(i==1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,7)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,10)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,11)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,14)=0
			end if
			if(j==m+1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,11)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,12)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,13)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,14)=0
			end if
			if(j==1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,7)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,8)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,9)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,10)=0
			end if
			if(k==n+1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,9)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,10)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,13)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,14)=0
			end if
			if(k==1)then
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,7)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,8)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,11)=0
				vertice((i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k,12)=0
			end if
		end do
	end do
end do
do i=1,l !generate cells
	do j=1,m
		do k=1,n
			element((i-1)*m*n+(j-1)*n+k,1)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k
			element((i-1)*m*n+(j-1)*n+k,2)=(i)*(m+1)*(n+1)+(j-1)*(n+1)+k
			element((i-1)*m*n+(j-1)*n+k,3)=(i)*(m+1)*(n+1)+(j-1)*(n+1)+k+1
			element((i-1)*m*n+(j-1)*n+k,4)=(i-1)*(m+1)*(n+1)+(j-1)*(n+1)+k+1
			element((i-1)*m*n+(j-1)*n+k,5)=(i-1)*(m+1)*(n+1)+(j)*(n+1)+k
			element((i-1)*m*n+(j-1)*n+k,6)=(i)*(m+1)*(n+1)+(j)*(n+1)+k
			element((i-1)*m*n+(j-1)*n+k,7)=(i)*(m+1)*(n+1)+(j)*(n+1)+k+1
			element((i-1)*m*n+(j-1)*n+k,8)=(i-1)*(m+1)*(n+1)+(j)*(n+1)+k+1
			if(i/=l)then
				element((i-1)*m*n+(j-1)*n+k,9)=(i)*m*n+(j-1)*n+k
			end if
			if(i/=1)then
				element((i-1)*m*n+(j-1)*n+k,10)=(i-2)*m*n+(j-1)*n+k
			end if
			if(j/=m)then
				element((i-1)*m*n+(j-1)*n+k,11)=(i-1)*m*n+(j)*n+k
			end if
			if(j/=1)then
				element((i-1)*m*n+(j-1)*n+k,12)=(i-1)*m*n+(j-2)*n+k
			end if
			if(k/=n)then
				element((i-1)*m*n+(j-1)*n+k,13)=(i-1)*m*n+(j-1)*n+k+1
			end if
			if(k/=1)then
				element((i-1)*m*n+(j-1)*n+k,14)=(i-1)*m*n+(j-1)*n+k-1
			end if
		end do
	end do
end do
end subroutine

subroutine BC_npns(step,dt)
!apply non-penetration & non-slip BC
!input:
!	step: # of current step
!	dt: step size
!output:
!	no ouputs: the velocity of vertices in the list will be updated according to the grid velocity
use globalvar,only:ppx,ppy,ppz,vu,vv,vw,BClist_npns
integer step
real dt
l=size(ppx,1)
n=size(BClist_npns,1)
do i=1,n
vu(BClist_npns(i))=0
vv(BClist_npns(i))=0
vw(BClist_npns(i))=0
!	if(step<l)then
!		vu(BClist_npns(i))=(ppx(step+1,BClist_npns(i))-ppx(step,BClist_npns(i)))/dt
!		vv(BClist_npns(i))=(ppy(step+1,BClist_npns(i))-ppy(step,BClist_npns(i)))/dt
!		vw(BClist_npns(i))=(ppz(step+1,BClist_npns(i))-ppz(step,BClist_npns(i)))/dt
!	else
!		vu(BClist_npns(i))=(ppx(step,BClist_npns(i))-ppx(step-1,BClist_npns(i)))/dt
!		vv(BClist_npns(i))=(ppy(step,BClist_npns(i))-ppy(step-1,BClist_npns(i)))/dt
!		vw(BClist_npns(i))=(ppz(step,BClist_npns(i))-ppz(step-1,BClist_npns(i)))/dt
!	end if
end do
end subroutine

subroutine BC_velo(refv)
!apply velosity BC
!input:
!	refv: [1*3] velocity to be set
!output:
!	no ouputs: the velocity of vertices in the list will be updated according to the grid velocity
use globalvar,only:vu,vv,vw,BClist_velo
double precision refv(3)
n=size(BClist_velo,1)
do i=1,n
	vu(BClist_velo(i))=refv(1)
	vv(BClist_velo(i))=refv(2)
	vw(BClist_velo(i))=refv(3)
end do
end subroutine
