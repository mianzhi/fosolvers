!----------------------------------------------------------------------------- best with 100 columns

!> shape, mapping, and quadrature
!> NOTE: initSMQ() must be called to have full functionalists
module modSMQ
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  double precision,parameter::TET10_C(10,10)=reshape(& !< TET10 shape function factors
  &  [ 2d0, 2d0, 2d0, 4d0, 4d0, 4d0,-3d0,-3d0,-3d0, 1d0,&
  &    2d0, 0d0, 0d0, 0d0, 0d0, 0d0,-1d0, 0d0, 0d0, 0d0,&
  &    0d0, 2d0, 0d0, 0d0, 0d0, 0d0, 0d0,-1d0, 0d0, 0d0,&
  &    0d0, 0d0, 2d0, 0d0, 0d0, 0d0, 0d0, 0d0,-1d0, 0d0,&
  &   -4d0, 0d0, 0d0,-4d0,-4d0, 0d0, 4d0, 0d0, 0d0, 0d0,&
  &    0d0, 0d0, 0d0, 4d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
  &    0d0,-4d0, 0d0,-4d0, 0d0,-4d0, 0d0, 4d0, 0d0, 0d0,&
  &    0d0, 0d0,-4d0, 0d0,-4d0,-4d0, 0d0, 0d0, 4d0, 0d0,&
  &    0d0, 0d0, 0d0, 0d0, 0d0, 4d0, 0d0, 0d0, 0d0, 0d0,&
  &    0d0, 0d0, 0d0, 0d0, 4d0, 0d0, 0d0, 0d0, 0d0, 0d0],[10,10])
  double precision,parameter,public::TET10_QP(DIMS,4)=reshape(& !< TET10 quadrature points
  &  [5.85410196624968515d-1,1.38196601125010504d-1,1.38196601125010504d-1,&
  &   1.38196601125010504d-1,5.85410196624968515d-1,1.38196601125010504d-1,&
  &   1.38196601125010504d-1,1.38196601125010504d-1,5.85410196624968515d-1,&
  &   1.38196601125010504d-1,1.38196601125010504d-1,1.38196601125010504d-1],[DIMS,4])
  double precision,parameter,public::TET10_QW(4)=& !< TET10 quadrature weights
  &  [4.16666666666666667d-2,4.16666666666666667d-2,4.16666666666666667d-2,4.16666666666666667d-2]
  
  double precision,parameter::TRI6_C(6,6)=reshape(& !< TRI6 shape function factors
  &  [ 2d0, 2d0, 4d0,-3d0,-3d0, 1d0,&
  &    2d0, 0d0, 0d0,-1d0, 0d0, 0d0,&
  &    0d0, 2d0, 0d0, 0d0,-1d0, 0d0,&
  &   -4d0, 0d0,-4d0, 4d0, 0d0, 0d0,&
  &    0d0, 0d0, 4d0, 0d0, 0d0, 0d0,&
  &    0d0,-4d0,-4d0, 0d0, 4d0, 0d0],[6,6])
  double precision,parameter,public::TRI6_QP(2,3)=reshape(& !< TRI6 quadrature points
  &  [1d0/6d0,1d0/6d0,2d0/3d0,1d0/6d0,1d0/6d0,2d0/3d0],[2,3])
  double precision,parameter,public::TRI6_QW(3)=& !<TRI6 quadrature weights
  &  [1d0/6d0,1d0/6d0,1d0/6d0]
  
  ! to be used as constant once initSMQ() is called
  double precision,public::TET10_SHAPE_QP(10,4) !< TET10 shape function at QP
  double precision,public::TET10_GRAD_QP(DIMS,10,4) !< TET10 shape grad at QP
  double precision,public::TRI6_SHAPE_QP(6,3) !< TRI6 shape function at QP
  double precision,public::TRI6_GRAD_QP(2,6,3) !< TRI6 shape grad at QP
  
  public::initSMQ
  
contains
  
  !> initialize the SMQ module
  subroutine initSMQ
    do i=1,4
      TET10_SHAPE_QP(:,i)=shapeTet10(TET10_QP(:,i))
      TET10_GRAD_QP(:,:,i)=gradShapeTet10(TET10_QP(:,i))
    end do
    do i=1,3
      TRI6_SHAPE_QP(:,i)=shapeTri6(TRI6_QP(:,i))
      TRI6_GRAD_QP(:,:,i)=gradShapeTri6(TRI6_QP(:,i))
    end do
  end subroutine
  
  !> shape function of reference TET10
  pure function shapeTet10(xx)
    double precision,intent(in)::xx(DIMS) !< test location in the reference coordinate
    double precision::shapeTet10(10) !< value of the 10 shape functions
    double precision::terms(10)
    
    terms(:)=[xx(1)**2,xx(2)**2,xx(3)**2,xx(1)*xx(2),xx(1)*xx(3),xx(2)*xx(3),xx(1),xx(2),xx(3),1d0]
    shapeTet10=matmul(terms,TET10_C)
  end function
  
  !> shape function gradient of reference TET10
  pure function gradShapeTet10(xx)
    double precision,intent(in)::xx(DIMS) !< test location in the reference coordinate
    double precision::gradShapeTet10(DIMS,10) !< gradient of the 10 shape functions
    double precision::terms(DIMS,10)
    
    terms(1,:)=[2*xx(1),0d0,0d0,xx(2),xx(3),0d0,1d0,0d0,0d0,0d0]
    terms(2,:)=[0d0,2*xx(2),0d0,xx(1),0d0,xx(3),0d0,1d0,0d0,0d0]
    terms(3,:)=[0d0,0d0,2*xx(3),0d0,xx(1),xx(2),0d0,0d0,1d0,0d0]
    gradShapeTet10=matmul(terms,TET10_C)
  end function
  
  !> Jacobian matrix of the iso-parametric map for TET10
  pure function mapJTet10(xN,xx)
    double precision,intent(in)::xN(DIMS,10) !< the location of the 10 nodes
    double precision,intent(in)::xx(DIMS) !< location in reference coordinate where J is evaluated
    double precision::mapJTet10(DIMS,DIMS) !< the Jacobian matrix
    
    mapJTet10=matmul(gradShapeTet10(xx),transpose(xN))
  end function
  
  !> shape function of reference TRI6
  pure function shapeTri6(xx)
    double precision,intent(in)::xx(2) !< test location in the reference coordinate
    double precision::shapeTri6(6) !< value of the 6 shape functions
    double precision::terms(6)
    
    terms(:)=[xx(1)**2,xx(2)**2,xx(1)*xx(2),xx(1),xx(2),1d0]
    shapeTri6=matmul(terms,TRI6_C)
  end function
  
  !> shape function gradient of reference TRI6
  pure function gradShapeTri6(xx)
    double precision,intent(in)::xx(2) !< test location in the reference coordinate
    double precision::gradShapeTri6(2,6) !< gradient of the 6 shape functions
    double precision::terms(2,6)
    
    terms(1,:)=[2*xx(1),0d0,xx(2),1d0,0d0,0d0]
    terms(2,:)=[0d0,2*xx(2),xx(1),0d0,1d0,0d0]
    gradShapeTri6=matmul(terms,TRI6_C)
  end function
  
  !> Jacobian matrix of the iso-parametric map for TRI6 (2 xx components by 3 x component)
  pure function mapJTri6(xN,xx)
    double precision,intent(in)::xN(DIMS,6) !< the location of the 6 nodes
    double precision,intent(in)::xx(2) !< location in reference coordinate where J is evaluated
    double precision::mapJTri6(2,DIMS) !< the Jacobian matrix
    
    mapJTri6=matmul(gradShapeTri6(xx),transpose(xN))
  end function
  
end module
