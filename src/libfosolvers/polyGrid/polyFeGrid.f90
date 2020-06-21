!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron and polygon finite element grid module
module modPolyFeGrid
  use modPolyGrid
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> polyhedron and polygon finite volume grid type
  type,extends(polyGrid),public::polyFeGrid
    double precision,allocatable::invJ(:,:,:,:) !< inverse mapping Jacobian at quadrature points
    double precision,allocatable::detJ(:,:) !< determinant of mapping Jacobian at quadrature points
  contains
    procedure,public::clear=>clearPolyFeGrid
    procedure,public::up=>upPolyFeGrid
    final::purgePolyFeGrid
  end type
  
contains
  
  !> clear this polyFeGrid
  elemental subroutine clearPolyFeGrid(this)
    class(polyFeGrid),intent(inout)::this !< this polyFeGrid
    
    call this%polyGrid%clear()
    if(allocated(this%invJ)) deallocate(this%invJ)
    if(allocated(this%detJ)) deallocate(this%detJ)
  end subroutine
  
  !> update this polyFeGrid
  subroutine upPolyFeGrid(this)
    use modSMQ
    class(polyFeGrid),intent(inout)::this !< this polyFeGrid
    double precision,allocatable::invJ(:,:,:,:)
    double precision,allocatable::detJ(:,:)
    double precision::Jacobian(DIMS,DIMS)
    integer,parameter::MAX_QP_PER_E=11 !< maximum 11 quadrature points per element
    integer::nQp
    
    if(.not.this%isUp)then
      call initSMQ()
      call this%polyGrid%up()
      allocate(invJ(DIMS,DIMS,MAX_QP_PER_E,this%nE))
      allocate(detJ(MAX_QP_PER_E,this%nE))
      nQp=0
      do i=1,this%nE
        select case(this%sE(i))
        case(TET)
          n=size(TET_QW)
          nQp=max(n,nQp)
          do j=1,n
            Jacobian=mapJTet(this%pN(:,this%iNE(1:this%nNE(i),i)),TET_QP(:,j))
            call findInvDet3by3(Jacobian,invJ(:,:,j,i),detJ(j,i))
          end do
        case(TET10)
          n=size(TET10_QW)
          nQp=max(n,nQp)
          do j=1,n
            Jacobian=mapJTet10(this%pN(:,this%iNE(1:this%nNE(i),i)),TET10_QP(:,j))
            call findInvDet3by3(Jacobian,invJ(:,:,j,i),detJ(j,i))
          end do
        case(TRI)
          n=size(TRI_QW)
          nQp=max(n,nQp)
          do j=1,n
            Jacobian(1:2,:)=mapJTri(this%pN(:,this%iNE(1:this%nNE(i),i)),TRI_QP(:,j))
            ! save d(xi)/d(xx1) cross d(xi)/d(yy1) at the 3rd row
            Jacobian(3,1)=Jacobian(1,2)*Jacobian(2,3)-Jacobian(1,3)*Jacobian(2,2)
            Jacobian(3,2)=Jacobian(1,3)*Jacobian(2,1)-Jacobian(1,1)*Jacobian(2,3)
            Jacobian(3,3)=Jacobian(1,1)*Jacobian(2,2)-Jacobian(1,2)*Jacobian(2,1)
            invJ(:,:,j,i)=0d0 ! Jacobian is non-square
            detJ(j,i)=norm2(Jacobian(3,:)) ! magnitude of d(xi)/d(xx1) cross d(xi)/d(yy1)
          end do
        case(TRI6)
          n=size(TRI6_QW)
          nQp=max(n,nQp)
          do j=1,n
            Jacobian(1:2,:)=mapJTri6(this%pN(:,this%iNE(1:this%nNE(i),i)),TRI6_QP(:,j))
            ! save d(xi)/d(xx1) cross d(xi)/d(yy1) at the 3rd row
            Jacobian(3,1)=Jacobian(1,2)*Jacobian(2,3)-Jacobian(1,3)*Jacobian(2,2)
            Jacobian(3,2)=Jacobian(1,3)*Jacobian(2,1)-Jacobian(1,1)*Jacobian(2,3)
            Jacobian(3,3)=Jacobian(1,1)*Jacobian(2,2)-Jacobian(1,2)*Jacobian(2,1)
            invJ(:,:,j,i)=0d0 ! Jacobian is non-square
            detJ(j,i)=norm2(Jacobian(3,:)) ! magnitude of d(xi)/d(xx1) cross d(xi)/d(yy1)
          end do
        case default
        end select
      end do
      if(allocated(this%invJ)) deallocate(this%invJ)
      if(allocated(this%detJ)) deallocate(this%detJ)
      allocate(this%invJ,source=invJ(:,:,1:nQp,:))
      allocate(this%detJ,source=detJ(1:nQp,:))
      deallocate(invJ)
      deallocate(detJ)
    end if
  contains
    
    !> find inverse and determinant of 3 by 3 matrix
    subroutine findInvDet3by3(A,invA,detA)
      double precision,intent(in)::A(3,3) !< 3 by 3 matrix
      double precision,intent(inout)::invA(3,3) !< inverse of A
      double precision,intent(inout)::detA !< determinant of A

      detA=A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)&
      &    +A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
      invA(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      invA(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      invA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      invA(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      invA(2,2)=(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      invA(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      invA(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      invA(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      invA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
      invA(:,:)=invA(:,:)/detA
    end subroutine
    
  end subroutine
  
  !> destructor of polyFeGrid
  elemental subroutine purgePolyFeGrid(this)
    type(polyFeGrid),intent(inout)::this !< this polyFeGrid
    
    call this%clear()
  end subroutine
  
end module
