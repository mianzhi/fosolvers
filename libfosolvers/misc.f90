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

!**************************************
! update the neighbour of all elements
!**************************************
subroutine updateNeib()
  use moduleGrid
  !$omp parallel do
  do i=1,nEle
    Ele(i)%Neib(:)=0
    do j=1,Ele(i)%SurfNum
      Ele(i)%Neib(j)=Ele(i)%getNeib(j)
    end do
  end do
  !$omp end parallel do
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
