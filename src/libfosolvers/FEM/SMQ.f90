!----------------------------------------------------------------------------- best with 100 columns

!> shape, mapping, and quadrature
!> NOTE: initSMQ() must be called to have full functionalists
module modSMQ
  public
  
  ! constants
  integer,private,parameter::DIMS=3 !< dimensions
  
  double precision,parameter::TET_C(4,4)=reshape(& !< TET shape function factors
  &  [-1d0,-1d0,-1d0, 1d0,&
  &    1d0, 0d0, 0d0, 0d0,&
  &    0d0, 1d0, 0d0, 0d0,&
  &    0d0, 0d0, 1d0, 0d0],[4,4])
  double precision,parameter,public::TET_QP(DIMS,1)=reshape(& !< TET quadrature points
  &  [0.25d0,0.25d0,0.25d0],[DIMS,1])
  double precision,parameter,public::TET_QW(1)=& !< TET quadrature weights
  &  [1.66666666666666666d-1]
  
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
  
  double precision,parameter::TRI_C(3,3)=reshape(& !< TRI shape function factors
  &  [-1d0,-1d0, 1d0,&
  &    1d0, 0d0, 0d0,&
  &    0d0, 1d0, 0d0],[3,3])
  double precision,parameter,public::TRI_QP(2,1)=reshape(& !< TRI quadrature points
  &  [1d0/3d0,1d0/3d0],[2,1])
  double precision,parameter,public::TRI_QW(1)=& !<TRI quadrature weights
  &  [0.5d0]
  
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
  double precision,public::TET_SHAPE_QP(4,1) !< TET shape function at QP
  double precision,public::TET_GRAD_QP(DIMS,4,1) !< TET shape grad at QP
  double precision,public::TET10_SHAPE_QP(10,4) !< TET10 shape function at QP
  double precision,public::TET10_GRAD_QP(DIMS,10,4) !< TET10 shape grad at QP
  double precision,public::TRI_SHAPE_QP(3,1) !< TRI shape function at QP
  double precision,public::TRI_GRAD_QP(2,3,1) !< TRI shape grad at QP
  double precision,public::TRI6_SHAPE_QP(6,3) !< TRI6 shape function at QP
  double precision,public::TRI6_GRAD_QP(2,6,3) !< TRI6 shape grad at QP
  
contains
  
  !> initialize the SMQ module
  subroutine initSMQ
    do i=1,1
      TET_SHAPE_QP(:,i)=shapeTet(TET_QP(:,i))
      TET_GRAD_QP(:,:,i)=gradShapeTet(TET_QP(:,i))
    end do
    do i=1,4
      TET10_SHAPE_QP(:,i)=shapeTet10(TET10_QP(:,i))
      TET10_GRAD_QP(:,:,i)=gradShapeTet10(TET10_QP(:,i))
    end do
    do i=1,1
      TRI_SHAPE_QP(:,i)=shapeTri(TRI_QP(:,i))
      TRI_GRAD_QP(:,:,i)=gradShapeTri(TRI_QP(:,i))
    end do
    do i=1,3
      TRI6_SHAPE_QP(:,i)=shapeTri6(TRI6_QP(:,i))
      TRI6_GRAD_QP(:,:,i)=gradShapeTri6(TRI6_QP(:,i))
    end do
  end subroutine
  
  !> shape function of reference TET
  pure function shapeTet(xx)
    double precision,intent(in)::xx(DIMS) !< test location in the reference coordinate
    double precision::shapeTet(4) !< value of the 4 shape functions
    double precision::terms(4)
    
    terms(:)=[xx(1),xx(2),xx(3),1d0]
    shapeTet=matmul(terms,TET_C)
  end function
  
  !> shape function gradient of reference TET
  pure function gradShapeTet(xx)
    double precision,intent(in),optional::xx(DIMS) !< test location in the reference coordinate
    double precision::gradShapeTet(DIMS,4) !< gradient of the 4 shape functions
    double precision::terms(DIMS,4)
    
    terms(1,:)=[1d0,0d0,0d0,0d0]
    terms(2,:)=[0d0,1d0,0d0,0d0]
    terms(3,:)=[0d0,0d0,1d0,0d0]
    gradShapeTet=matmul(terms,TET_C)
    if(present(xx))then
    end if
  end function
  
  !> Jacobian matrix of the iso-parametric map for TET
  pure function mapJTet(xN,xx)
    double precision,intent(in)::xN(DIMS,4) !< the location of the 4 nodes
    double precision,intent(in)::xx(DIMS) !< location in reference coordinate where J is evaluated
    double precision::mapJTet(DIMS,DIMS) !< the Jacobian matrix
    double precision::grad(DIMS,4)
    
    grad=gradShapeTet(xx)
    mapJTet=matmul(grad,transpose(xN))
  end function
  
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
    double precision::grad(DIMS,10)
    
    grad=gradShapeTet10(xx)
    mapJTet10=matmul(grad,transpose(xN))
  end function
  
  !> shape function of reference TRI
  pure function shapeTri(xx)
    double precision,intent(in)::xx(2) !< test location in the reference coordinate
    double precision::shapeTri(3) !< value of the 3 shape functions
    double precision::terms(3)
    
    terms(:)=[xx(1),xx(2),1d0]
    shapeTri=matmul(terms,TRI_C)
  end function
  
  !> shape function gradient of reference TRI
  pure function gradShapeTri(xx)
    double precision,intent(in),optional::xx(2) !< test location in the reference coordinate
    double precision::gradShapeTri(2,3) !< gradient of the 3 shape functions
    double precision::terms(2,3)
    
    terms(1,:)=[1d0,0d0,0d0]
    terms(2,:)=[0d0,1d0,0d0]
    gradShapeTri=matmul(terms,TRI_C)
    if(present(xx))then
    end if
  end function
  
  !> Jacobian matrix of the iso-parametric map for TRI (2 xx components by 3 x components)
  pure function mapJTri(xN,xx)
    double precision,intent(in)::xN(DIMS,3) !< the location of the 3 nodes
    double precision,intent(in)::xx(2) !< location in reference coordinate where J is evaluated
    double precision::mapJTri(2,DIMS) !< the Jacobian matrix
    double precision::grad(2,3)
    
    grad=gradShapeTri(xx)
    mapJTri=matmul(grad,transpose(xN))
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
  
  !> Jacobian matrix of the iso-parametric map for TRI6 (2 xx components by 3 x components)
  pure function mapJTri6(xN,xx)
    double precision,intent(in)::xN(DIMS,6) !< the location of the 6 nodes
    double precision,intent(in)::xx(2) !< location in reference coordinate where J is evaluated
    double precision::mapJTri6(2,DIMS) !< the Jacobian matrix
    double precision::grad(2,6)
    
    grad=gradShapeTri6(xx)
    mapJTri6=matmul(grad,transpose(xN))
  end function
  
end module
