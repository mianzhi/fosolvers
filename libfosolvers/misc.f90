!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! update the bounding box of the geometry
!*****************************************
subroutine updateBoundBox()
  use moduleGrid
  ! BoundBox:
  !   xmin | ymin | zmin
  !   xmax | ymax | zmax
  BoundBox(1,:)=Node(1)%Pos(:)
  BoundBox(2,:)=Node(1)%Pos(:)
  do i=2,nNode
    BoundBox(1,1)=min(BoundBox(1,1),Node(i)%Pos(1))
    BoundBox(1,2)=min(BoundBox(1,2),Node(i)%Pos(2))
    BoundBox(1,3)=min(BoundBox(1,3),Node(i)%Pos(3))
    BoundBox(2,1)=max(BoundBox(2,1),Node(i)%Pos(1))
    BoundBox(2,2)=max(BoundBox(2,2),Node(i)%Pos(2))
    BoundBox(2,3)=max(BoundBox(2,3),Node(i)%Pos(3))
  end do
end subroutine

!*************************************
! update the parameters of all facets
!*************************************
subroutine updateFacetPara()
  use moduleGrid
  !$omp parallel do
  do i=1,nFacet
    ! position of the center
    Facet(i)%PC(:)=Facet(i)%findPC()
    ! area
    Facet(i)%Area=Facet(i)%findArea()
    ! normal vector
    Facet(i)%Norm=Facet(i)%findNorm()
  end do
  !$omp end parallel do
end subroutine

!***************************************
! update the parameters of all elements
!***************************************
subroutine updateElePara()
  use moduleGrid
  !$omp parallel do
  do i=1,nEle
    ! neighbour
    Ele(i)%Neib(:)=0
    do j=1,Ele(i)%SurfNum
      Ele(i)%Neib(j)=Ele(i)%getNeib(j)
    end do
    ! position of the center
    Ele(i)%PC(:)=Ele(i)%findPC()
    ! volume
    Ele(i)%Vol=Ele(i)%findVol()
    ! position of the center of all surfaces
    Ele(i)%SurfPC(:,:)=0d0
    do j=1,Ele(i)%SurfNum
      Ele(i)%SurfPC(j,:)=Ele(i)%findSurfPC(j)
    end do
    ! surface area of all surfaces
    Ele(i)%SurfArea(:)=0d0
    do j=1,Ele(i)%SurfNum
      Ele(i)%SurfArea(j)=Ele(i)%findSurfArea(j)
    end do
    ! normal vectors of all surfaces
    Ele(i)%SurfNorm(:,:)=0d0
    do j=1,Ele(i)%SurfNum
      Ele(i)%SurfNorm(j,:)=Ele(i)%findSurfNorm(j)
    end do
  end do
  !$omp end parallel do
end subroutine

!*******************
! sort the elements
!*******************
subroutine sortEle()
  use moduleGrid
  double precision sortData(nEle),tempData,tempPos(3)
  integer m
  type(typeEle)::tempEle
  ! decide which direction we will sort along
  ! Note: updateBoundBox() is embeded because of its small consumption of computation time
  call updateBoundBox()
  m=maxloc(BoundBox(:,2)-BoundBox(:,1),1)
  ! find the elements' center (according to which they will be sorted)
  do i=1,nEle
    tempPos(:)=Ele(i)%findPC()
    sortData(i)=tempPos(m)
  end do
  ! insertion sort
  do i=2,nEle
    do j=i,2,-1
      if(sortData(j-1)<=sortData(j))then
        exit
      end if
      tempEle=Ele(j)
      Ele(j)=Ele(j-1)
      Ele(j-1)=tempEle
      tempData=sortData(j)
      sortData(j)=sortData(j-1)
      sortData(j-1)=tempData
    end do
  end do
end subroutine

!*****************************************************************************
! find the gradient of vector v of dimension m at the center of the element k
!*****************************************************************************
subroutine findEleGradAny(k,v,m,rst)
  use moduleGrid
  integer,intent(in)::k,m
  double precision,intent(in)::v(nEle,m)
  double precision,intent(out)::rst(m,3)
  double precision,allocatable::dx(:,:),dv(:,:),s(:)
  integer,parameter::lwork=900 ! according to returned work(1) with lwork=-1
  double precision distsq,work(lwork)
  integer nEq,rank,iwork(50),er
  nEq=count(Ele(k)%Neib(:)>0)
  allocate(dx(nEq,3))
  allocate(dv(max(nEq,3),m))
  allocate(s(min(nEq,3)))
  ! weighted least-squares gradient evaluation
  j=0
  do i=1,Ele(k)%SurfNum
    if(Ele(k)%Neib(i)>0)then
      j=j+1
      dx(j,:)=Ele(Ele(k)%Neib(i))%PC(:)-Ele(k)%PC(:)
      dv(j,:)=v(Ele(k)%Neib(i),:)-v(k,:)
      distsq=dot_product(dx(j,:),dx(j,:))
      dx(j,:)=dx(j,:)/distsq
      dv(j,:)=dv(j,:)/distsq
    end if
  end do
  call DGELSD(nEq,3,m,dx,nEq,dv,max(nEq,3),s,-1d0,rank,work,lwork,iwork,er) ! calling lapack
  if(er/=0)then
    write(*,'(a,i2)'),'Error: Lapack returning error flag: ',er
    stop
  end if
  rst(:,:)=transpose(dv(1:3,1:m))
end subroutine

!**************************************************************
! find the gradient of scaler v at the center of the element k
!**************************************************************
subroutine findEleGradScal(k,v,rst)
  use moduleGrid
  integer,intent(in)::k
  double precision,intent(in)::v(nEle)
  double precision,intent(out)::rst(3)
  call findEleGradAny(k,v,1,rst)
end subroutine

!*****************************************************************************
! find the gradient of vector v of dimension 3 at the center of the element k
!*****************************************************************************
subroutine findEleGradVect(k,v,rst)
  use moduleGrid
  integer,intent(in)::k
  double precision,intent(in)::v(nEle,3)
  double precision,intent(out)::rst(3,3)
  call findEleGradAny(k,v,3,rst)
end subroutine

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
function find3PNorm(P1,P2,P3)
  double precision,intent(in)::P1(3),P2(3),P3(3)
  double precision find3PNorm(3),a(3),b(3),c(3)
  a(:)=P2(:)-P1(:)
  b(:)=P3(:)-P2(:)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  find3PNorm(:)=c(:)/sqrt(dot_product(c(:),c(:)))
end function

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

!***************
! show progress
!***************
subroutine showProg(v,vtot)
  double precision,intent(in)::v,vtot
  double precision prog
  prog=v/vtot*100d0
  write(*,'(15a,a,f5.1,a,$)'),(char(8),i=1,15),'progress:',prog,'%'
end subroutine
