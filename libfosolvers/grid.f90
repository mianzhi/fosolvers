!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! data and elementary procedures for grid
!*****************************************
module moduleGrid
  private
  
  !----------
  ! typeNode
  !----------
  type,public::typeNode
    double precision Pos(3)
  end type
  type(typeNode),public,allocatable,save::Node(:)
  integer,public,save::nNode
  
  !-----------
  ! typePoint
  !-----------
  type,public::typePoint
    integer NodeInd
    integer GeoEnti
  contains
    procedure,public::findPos=>findPointPos
  end type
  type(typePoint),public,allocatable,save::Point(:)
  integer,public,save::nPoint
  
contains
  
  !---------------------------------
  ! find the position of this point
  !---------------------------------
  subroutine findPointPos(this,rst)
    class(typePoint),intent(in)::this
    double precision,intent(out)::rst(3)
    rst(:)=Node(this%NodeInd)%Pos(:)
  end subroutine
  
end module

program test
  use moduleGrid
  double precision rst(3)
  allocate(Node(2))
  Node(1)%Pos=[1,2,3]
  Node(2)%Pos=[4,5,6]
  allocate(Point(3))
  Point(1)%NodeInd=1
  Point(2)%NodeInd=2
  Point(3)%NodeInd=1
  do i=1,3
    call Point(i)%findPos(rst)
    write(*,*),rst
  end do
  write(*,*),nNode,nPoint
end program
