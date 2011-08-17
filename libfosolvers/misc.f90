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

!***************
! show progress
!***************
subroutine showProg(v,vtot)
  double precision,intent(in)::v,vtot
  double precision prog
  prog=v/vtot*100d0
  write(*,'(15a,a,f5.1,a,$)'),(char(8),i=1,15),'progress:',prog,'%'
end subroutine
