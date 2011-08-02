!----------------------------------------------------------------------------- best with 100 columns

!***************************************
! find the bounding box of the geometry
!***************************************
subroutine findBound(rst)
  use moduleGrid
  double precision,intent(out)::rst(2,3)
  ! rst:
  !   xmin | ymin | zmin
  !   xmax | ymax | zmax
  rst(1,:)=Node(1)%Pos(:)
  rst(2,:)=Node(1)%Pos(:)
  do i=2,nNode
    rst(1,1)=min(rst(1,1),Node(i)%Pos(1))
    rst(1,2)=min(rst(1,2),Node(i)%Pos(2))
    rst(1,3)=min(rst(1,3),Node(i)%Pos(3))
    rst(2,1)=max(rst(2,1),Node(i)%Pos(1))
    rst(2,2)=max(rst(2,2),Node(i)%Pos(2))
    rst(2,3)=max(rst(2,3),Node(i)%Pos(3))
  end do
end subroutine

!***************
! show progress
!***************
subroutine showProg(v,vtot)
  double precision,intent(in)::v,vtot
  real prog
  prog=v/vtot*100d0
  write(*,'(15a,a,f5.1,a,$)'),(char(8),i=1,15),'progress:',prog,'%'
end subroutine
