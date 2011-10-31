!----------------------------------------------------------------------------- best with 100 columns

!****************************
! global data related to MPI
!****************************
module moduleMPIvar
  use mpi
  integer,save::errMPI,pidMPI,sizeMPI,nblkMPI
  integer,allocatable,save::statMPI(:),llenMPI(:),ltypeMPI(:)
  integer(MPI_address_kind),allocatable,save::ldispMPI(:)
  integer(MPI_address_kind),save::baseMPI
  integer,save::typeNodeMPI,typePointMPI,typeLineMPI,typeTriMPI,typeQuadMPI,typeTetMPI,typeHexMPI,&
  &             typeFacetMPI,typeEleMPI
end module

!****************************
! initialize MPI environment
!****************************
subroutine initMPI()
  use moduleGrid
  use moduleMPIvar
  
  type(typeNode)::sampNode
  type(typePoint)::sampPoint
  type(typeLine)::sampLine
  type(typeTri)::sampTri
  type(typeQuad)::sampQuad
  type(typeTet)::sampTet
  type(typeHex)::sampHex
  type(typeFacet)::sampFacet
  type(typeEle)::sampEle
  
  allocate(statMPI(MPI_status_size))
  call MPI_init(errMPI)
  call MPI_comm_rank(MPI_comm_world,pidMPI,errMPI)
  call MPI_comm_size(MPI_comm_world,sizeMPI,errMPI)
  
  ! Note: all the types are splited into scalers to ensure correctness
  ! node struct
  nblkMPI=3
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,3
    call MPI_get_address(sampNode%Pos(i),ldispMPI(i),errMPI)
  end do
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_double_precision
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeNodeMPI,errMPI)
  call MPI_type_commit(typeNodeMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! point struct
  nblkMPI=2
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  call MPI_get_address(sampPoint%NodeInd,ldispMPI(1),errMPI)
  call MPI_get_address(sampPoint%GeoEnti,ldispMPI(2),errMPI)
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_integer
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typePointMPI,errMPI)
  call MPI_type_commit(typePointMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! line struct
  nblkMPI=3
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,2
    call MPI_get_address(sampLine%NodeInd(i),ldispMPI(i),errMPI)
  end do
  call MPI_get_address(sampLine%GeoEnti,ldispMPI(3),errMPI)
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_integer
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeLineMPI,errMPI)
  call MPI_type_commit(typeLineMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! triangle struct
  nblkMPI=4
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,3
    call MPI_get_address(sampTri%NodeInd(i),ldispMPI(i),errMPI)
  end do
  call MPI_get_address(sampTri%GeoEnti,ldispMPI(4),errMPI)
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_integer
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeTriMPI,errMPI)
  call MPI_type_commit(typeTriMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! quadrilateral struct
  nblkMPI=5
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,4
    call MPI_get_address(sampQuad%NodeInd(i),ldispMPI(i),errMPI)
  end do
  call MPI_get_address(sampQuad%GeoEnti,ldispMPI(5),errMPI)
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_integer
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeQuadMPI,errMPI)
  call MPI_type_commit(typeQuadMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! tetrahedron struct
  nblkMPI=6
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,4
    call MPI_get_address(sampTet%NodeInd(i),ldispMPI(i),errMPI)
  end do
  call MPI_get_address(sampTet%GeoEnti,ldispMPI(5),errMPI)
  call MPI_get_address(sampTet%Prt,ldispMPI(6),errMPI)
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_integer
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeTetMPI,errMPI)
  call MPI_type_commit(typeTetMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! hexahedron struct
  nblkMPI=10
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,8
    call MPI_get_address(sampHex%NodeInd(i),ldispMPI(i),errMPI)
  end do
  call MPI_get_address(sampHex%GeoEnti,ldispMPI(9),errMPI)
  call MPI_get_address(sampHex%Prt,ldispMPI(10),errMPI)
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_integer
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeHexMPI,errMPI)
  call MPI_type_commit(typeHexMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! facet struct
  nblkMPI=28
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  call MPI_get_address(sampFacet%ShapeType,ldispMPI(1),errMPI)
  call MPI_get_address(sampFacet%ShapeInd,ldispMPI(2),errMPI)
  call MPI_get_address(sampFacet%NodeNum,ldispMPI(3),errMPI)
  do i=1,15
    call MPI_get_address(sampFacet%NodeInd(i),ldispMPI(i+3),errMPI)
  end do
  call MPI_get_address(sampFacet%GeoEnti,ldispMPI(19),errMPI)
  do i=1,2
    call MPI_get_address(sampFacet%NeibEle(i),ldispMPI(i+19),errMPI)
  end do
  do i=1,3
    call MPI_get_address(sampFacet%PC(i),ldispMPI(i+21),errMPI)
  end do
  call MPI_get_address(sampFacet%Area,ldispMPI(25),errMPI)
  do i=1,3
    call MPI_get_address(sampFacet%Norm(i),ldispMPI(i+25),errMPI)
  end do
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(1:21)=MPI_integer
  ltypeMPI(22:28)=MPI_double_precision
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeFacetMPI,errMPI)
  call MPI_type_commit(typeFacetMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! element struct
  nblkMPI=85
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  call MPI_get_address(sampEle%ShapeType,ldispMPI(1),errMPI)
  call MPI_get_address(sampEle%ShapeInd,ldispMPI(2),errMPI)
  call MPI_get_address(sampEle%NodeNum,ldispMPI(3),errMPI)
  call MPI_get_address(sampEle%SurfNum,ldispMPI(4),errMPI)
  do i=1,27
    call MPI_get_address(sampEle%NodeInd(i),ldispMPI(i+4),errMPI)
  end do
  call MPI_get_address(sampEle%GeoEnti,ldispMPI(32),errMPI)
  call MPI_get_address(sampEle%Prt,ldispMPI(33),errMPI)
  do i=1,6
    call MPI_get_address(sampEle%Neib(i),ldispMPI(i+33),errMPI)
  end do
  do i=1,3
    call MPI_get_address(sampEle%PC(i),ldispMPI(i+39),errMPI)
  end do
  call MPI_get_address(sampEle%Vol,ldispMPI(43),errMPI)
  do i=1,6
    do j=1,3
      call MPI_get_address(sampEle%SurfPC(i,j),ldispMPI((i-1)*3+j+43),errMPI)
    end do
  end do
  do i=1,6
    call MPI_get_address(sampEle%SurfArea(i),ldispMPI(i+61),errMPI)
  end do
  do i=1,6
    do j=1,3
      call MPI_get_address(sampEle%SurfNorm(i,j),ldispMPI((i-1)*3+j+67),errMPI)
    end do
  end do
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(1:39)=MPI_integer
  ltypeMPI(40:85)=MPI_double_precision
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeEleMPI,errMPI)
  call MPI_type_commit(typeEleMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
end subroutine

!**************************
! finalize MPI environment
!**************************
subroutine finalMPI()
  use moduleMPIvar
  call MPI_type_free(typeNodeMPI,errMPI)
  call MPI_type_free(typePointMPI,errMPI)
  call MPI_type_free(typeLineMPI,errMPI)
  call MPI_type_free(typeTriMPI,errMPI)
  call MPI_type_free(typeQuadMPI,errMPI)
  call MPI_type_free(typeTetMPI,errMPI)
  call MPI_type_free(typeHexMPI,errMPI)
  call MPI_type_free(typeFacetMPI,errMPI)
  call MPI_type_free(typeEleMPI,errMPI)
  call MPI_finalize(errMPI)
end subroutine
