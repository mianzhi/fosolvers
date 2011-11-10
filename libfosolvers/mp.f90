!----------------------------------------------------------------------------- best with 100 columns

!****************************
! global data related to MPI
!****************************
module moduleMPIvar
  use mpi
  use moduleMiscDataStruct
  
  integer,parameter::ROOT_PID=0
  
  integer,save::errMPI,pidMPI,sizeMPI
  integer,allocatable,save::statMPI(:)
  integer,save::typeNodeMPI,typePointMPI,typeLineMPI,typeTriMPI,typeQuadMPI,typeTetMPI,typeHexMPI,&
  &             typeFacetMPI,typeEleMPI,typeCondMPI
  integer,allocatable,save::mapNode(:),mapPoint(:),mapLine(:),mapTri(:),mapQuad(:),mapTet(:),&
  &                         mapHex(:),mapFacet(:),mapEle(:)
  type(typePtrScalArray),allocatable,save::transNodeScal(:),transEleScal(:)
  type(typePtrVectArray),allocatable,save::transNodeVect(:),transEleVect(:)
  type(typePtrTensArray),allocatable,save::transNodeTens(:),transEleTens(:)
  integer,save::ntransNodeScal,ntransNodeVect,ntransNodeTens,ntransEleScal,ntransEleVect,&
  &             ntransEleTens
end module

!****************************
! initialize MPI environment
!****************************
subroutine initMPI()
  use moduleGrid
  use moduleCond
  use moduleMPIvar
  
  integer nblkMPI
  integer,allocatable::llenMPI(:),ltypeMPI(:)
  integer(MPI_address_kind),allocatable::ldispMPI(:)
  integer(MPI_address_kind)::baseMPI
  
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
  
  ! Note: all the types are split into scalers to ensure correctness
  ! node struct
  nblkMPI=DIMS
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,DIMS
    call MPI_get_address(sampNode%Pos(i),ldispMPI(i),errMPI)
  end do
  baseMPI=ldispMPI(1)
  ldispMPI(:)=ldispMPI(:)-baseMPI
  ltypeMPI(:)=MPI_double_precision
  call MPI_type_create_struct(nblkMPI,llenMPI,ldispMPI,ltypeMPI,typeNodeMPI,errMPI)
  call MPI_type_commit(typeNodeMPI,errMPI)
  deallocate(llenMPI,ldispMPI,ltypeMPI)
  ! point struct
  nblkMPI=POINT_NODE_NUM+1
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
  nblkMPI=LINE_NODE_NUM+1
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,LINE_NODE_NUM
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
  nblkMPI=TRI_NODE_NUM+1
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,TRI_NODE_NUM
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
  nblkMPI=QUAD_NODE_NUM+1
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,QUAD_NODE_NUM
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
  nblkMPI=TET_NODE_NUM+2
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,TET_NODE_NUM
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
  nblkMPI=HEX_NODE_NUM+2
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  do i=1,HEX_NODE_NUM
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
  nblkMPI=FACET_MAX_NODE_NUM+FACET_NEIB_ELE_NUM+DIMS+DIMS+5
  allocate(llenMPI(nblkMPI))
  allocate(ldispMPI(nblkMPI))
  allocate(ltypeMPI(nblkMPI))
  llenMPI(:)=1
  call MPI_get_address(sampFacet%ShapeType,ldispMPI(1),errMPI)
  call MPI_get_address(sampFacet%ShapeInd,ldispMPI(2),errMPI)
  call MPI_get_address(sampFacet%NodeNum,ldispMPI(3),errMPI)
  do i=1,FACET_MAX_NODE_NUM
    call MPI_get_address(sampFacet%NodeInd(i),ldispMPI(i+3),errMPI)
  end do
  call MPI_get_address(sampFacet%GeoEnti,ldispMPI(19),errMPI)
  do i=1,FACET_NEIB_ELE_NUM
    call MPI_get_address(sampFacet%NeibEle(i),ldispMPI(i+19),errMPI)
  end do
  do i=1,DIMS
    call MPI_get_address(sampFacet%PC(i),ldispMPI(i+21),errMPI)
  end do
  call MPI_get_address(sampFacet%Area,ldispMPI(25),errMPI)
  do i=1,DIMS
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
  do i=1,ELE_MAX_NODE_NUM
    call MPI_get_address(sampEle%NodeInd(i),ldispMPI(i+4),errMPI)
  end do
  call MPI_get_address(sampEle%GeoEnti,ldispMPI(32),errMPI)
  call MPI_get_address(sampEle%Prt,ldispMPI(33),errMPI)
  do i=1,ELE_MAX_SURF_NUM
    call MPI_get_address(sampEle%Neib(i),ldispMPI(i+33),errMPI)
  end do
  do i=1,DIMS
    call MPI_get_address(sampEle%PC(i),ldispMPI(i+39),errMPI)
  end do
  call MPI_get_address(sampEle%Vol,ldispMPI(43),errMPI)
  do i=1,ELE_MAX_SURF_NUM
    do j=1,DIMS
      call MPI_get_address(sampEle%SurfPC(i,j),ldispMPI((i-1)*3+j+43),errMPI)
    end do
  end do
  do i=1,ELE_MAX_SURF_NUM
    call MPI_get_address(sampEle%SurfArea(i),ldispMPI(i+61),errMPI)
  end do
  do i=1,ELE_MAX_SURF_NUM
    do j=1,DIMS
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

!********************************************
! broadcast static information (data tables)
!********************************************
subroutine bcastStatic()
  use moduleCond
  use moduleMPIvar
  
  integer n
  
  ! broadcast dataTab
  if(allocated(dataTab1d))then
    n=size(dataTab1d)
  else
    n=0
  end if
  call MPI_bcast(n,1,MPI_integer,ROOT_PID,MPI_comm_world,errMPI)
  do i=1,n
    call MPI_bcast(dataTab1d(i)%length,1,MPI_integer,ROOT_PID,MPI_comm_world,errMPI)
    call MPI_bcast(dataTab1d(i)%x,dataTab1d(i)%length,MPI_double_precision,&
    &              ROOT_PID,MPI_comm_world,errMPI)
    call MPI_bcast(dataTab1d(i)%y,dataTab1d(i)%length,MPI_double_precision,&
    &              ROOT_PID,MPI_comm_world,errMPI)
  end do
  
end subroutine

!*************************************************
! receive static information (Conditions,dataTab)
!*************************************************
subroutine recvStatic()
  use moduleCond
  use moduleMPIvar
  
  integer n
  
  ! receive dataTab
  call MPI_bcast(n,1,MPI_integer,ROOT_PID,MPI_comm_world,errMPI)
  if(n>0)then
    allocate(dataTab1d(n))
  end if
  do i=1,n
    call MPI_bcast(dataTab1d(i)%length,1,MPI_integer,ROOT_PID,MPI_comm_world,errMPI)
    allocate(dataTab1d(i)%x(dataTab1d(i)%length))
    allocate(dataTab1d(i)%y(dataTab1d(i)%length))
    call MPI_bcast(dataTab1d(i)%x,dataTab1d(i)%length,MPI_double_precision,&
    &              ROOT_PID,MPI_comm_world,errMPI)
    call MPI_bcast(dataTab1d(i)%y,dataTab1d(i)%length,MPI_double_precision,&
    &              ROOT_PID,MPI_comm_world,errMPI)
  end do
  
end subroutine

!****************************************************
! generate table of data to be wraped with partition
!****************************************************
subroutine genPrtDataTab(a,b,c,d,e,f)
  use moduleMPIvar
  integer,intent(in)::a,b,c,d,e,f
  ! allocate output data space
  ntransNodeScal=a
  ntransNodeVect=b
  ntransNodeTens=c
  ntransEleScal=d
  ntransEleVect=e
  ntransEleTens=f
  allocate(transNodeScal(a))
  allocate(transNodeVect(b))
  allocate(transNodeTens(c))
  allocate(transEleScal(d))
  allocate(transEleVect(e))
  allocate(transEleTens(f))
end subroutine

!*************************************
! distribute partition k to process p
!*************************************
! e.g.: to distribute element scaler data v(nEle) together with partition
!   call genPrtDataTab(0,0,0,1,0,0)
!   ...
!   transEleScal(1)%ptr=>v
!   ...
!   call distriPrt(k,p)
subroutine distriPrt(k,p)
  use moduleGrid
  use moduleCond
  use moduleMPIvar
  
  integer,intent(in)::k,p
  integer nkNode,nkPoint,nkLine,nkTri,nkQuad,nkTet,nkHex,nkFacet,nkEle
  type(typeNode),allocatable::buffNode(:)
  type(typePoint),allocatable::buffPoint(:)
  type(typeLine),allocatable::buffLine(:)
  type(typeTri),allocatable::buffTri(:)
  type(typeQuad),allocatable::buffQuad(:)
  type(typeTet),allocatable::buffTet(:)
  type(typeHex),allocatable::buffHex(:)
  type(typeFacet),allocatable::buffFacet(:)
  type(typeEle),allocatable::buffEle(:)
  integer gmapNode(nNode),gmapPoint(nPoint),gmapLine(nLine),gmapTri(nTri),gmapQuad(nQuad),&
  &       gmapTet(nTet),gmapHex(nHex),gmapFacet(nFacet),gmapEle(nEle)
  double precision,allocatable::buffDataScal(:),buffDataVect(:,:),buffDataTens(:,:,:)
  
  ! copy elements
  nkEle=count(Ele(:)%Prt==k)
  allocate(buffEle(nkEle))
  allocate(mapEle(nkEle))
  j=1
  do i=1,nEle
    if(Ele(i)%Prt==k)then
      mapEle(j)=i
      gmapEle(i)=j
      j=j+1
    end if
  end do
  buffEle(:)=Ele(mapEle(:))
  ! copy facets
  nkFacet=sum([(count(buffEle(i)%Neib(:)<0),i=1,nkEle)])
  allocate(buffFacet(nkFacet))
  allocate(mapFacet(nkFacet))
  l=1
  do i=1,nkEle
    do j=1,size(buffEle(i)%Neib)
      if(buffEle(i)%Neib(j)<0)then
        mapFacet(l)=-buffEle(i)%Neib(j)
        gmapFacet(-buffEle(i)%Neib(j))=l
        l=l+1
      end if
    end do
  end do
  buffFacet(:)=Facet(mapFacet)
  ! copy tetrahedrons
  nkTet=count(buffEle(:)%ShapeType==TET_TYPE)
  allocate(buffTet(nkTet))
  allocate(mapTet(nkTet))
  j=1
  do i=1,nkEle
    if(buffEle(i)%ShapeType==TET_TYPE)then
      mapTet(j)=buffEle(i)%ShapeInd
      gmapTet(buffEle(i)%ShapeInd)=j
      j=j+1
    end if
  end do
  buffTet(:)=Tet(mapTet)
  ! copy hexahedrons
  nkHex=count(buffEle(:)%ShapeType==HEX_TYPE)
  allocate(buffHex(nkHex))
  allocate(mapHex(nkHex))
  j=1
  do i=1,nkEle
    if(buffEle(i)%ShapeType==HEX_TYPE)then
      mapHex(j)=buffEle(i)%ShapeInd
      gmapHex(buffEle(i)%ShapeInd)=j
      j=j+1
    end if
  end do
  buffHex(:)=Hex(mapHex)
  ! copy triangles
  nkTri=count(buffFacet(:)%ShapeType==TRI_TYPE)
  allocate(buffTri(nkTri))
  allocate(mapTri(nkTri))
  j=1
  do i=1,nkFacet
    if(buffFacet(i)%ShapeType==TRI_TYPE)then
      mapTri(j)=buffFacet(i)%ShapeInd
      gmapTri(buffFacet(i)%ShapeInd)=j
      j=j+1
    end if
  end do
  buffTri(:)=Tri(mapTri)
  ! copy quadrilaterals
  nkQuad=count(buffFacet(:)%ShapeType==QUAD_TYPE)
  allocate(buffQuad(nkQuad))
  allocate(mapQuad(nkQuad))
  j=1
  do i=1,nkFacet
    if(buffFacet(i)%ShapeType==QUAD_TYPE)then
      mapQuad(j)=buffFacet(i)%ShapeInd
      gmapQuad(buffFacet(i)%ShapeInd)=j
      j=j+1
    end if
  end do
  buffQuad(:)=Quad(mapQuad)
  ! copy nodes
  gmapNode(:)=0
  l=0
  do i=1,nkEle
    do j=1,buffEle(i)%NodeNum
      if(gmapNode(buffEle(i)%NodeInd(j))==0)then
        l=l+1
        gmapNode(buffEle(i)%NodeInd(j))=l
      end if
    end do
  end do
  nkNode=count(gmapNode(:)>0)
  allocate(buffNode(nkNode))
  allocate(mapNode(nkNode))
  do i=1,nNode
    if(gmapNode(i)>0)then
      mapNode(gmapNode(i))=i
    end if
  end do
  buffNode(:)=Node(mapNode)
  ! copy points
  gmapPoint(:)=0
  j=0
  do i=1,nPoint
    if(gmapNode(Point(i)%NodeInd)>0)then
      j=j+1
      gmapPoint(i)=j
    end if
  end do
  nkPoint=count(gmapPoint(:)>0)
  allocate(buffPoint(nkPoint))
  allocate(mapPoint(nkPoint))
  do i=1,nPoint
    if(gmapPoint(i)>0)then
      mapPoint(gmapPoint(i))=i
    end if
  end do
  buffPoint(:)=Point(mapPoint)
  ! copy lines
  gmapLine(:)=0
  j=0
  do i=1,nLine
    if(all(gmapNode(Line(i)%NodeInd(:))>0))then
      j=j+1
      gmapLine(i)=j
    end if
  end do
  nkLine=count(gmapLine(:)>0)
  allocate(buffLine(nkLine))
  allocate(mapLine(nkLine))
  do i=1,nLine
    if(gmapLine(i)>0)then
      mapLine(gmapLine(i))=i
    end if
  end do
  buffLine(:)=Line(mapLine)
  
  ! correct points
  forall(i=1:nkPoint)
    buffPoint(i)%NodeInd=gmapNode(buffPoint(i)%NodeInd)
  end forall
  ! correct lines
  forall(i=1:nkLine)
    buffLine(i)%NodeInd(:)=gmapNode(buffLine(i)%NodeInd(:))
  end forall
  ! correct triangles
  forall(i=1:nkTri)
    buffTri(i)%NodeInd(:)=gmapNode(buffTri(i)%NodeInd(:))
  end forall
  ! correct quadrilaterals
  forall(i=1:nkQuad)
    buffQuad(i)%NodeInd(:)=gmapNode(buffQuad(i)%NodeInd(:))
  end forall
  ! correct tetrahedrons
  forall(i=1:nkTet)
    buffTet(i)%NodeInd(:)=gmapNode(buffTet(i)%NodeInd(:))
  end forall
  ! correct hexahedrons
  forall(i=1:nkHex)
    buffHex(i)%NodeInd(:)=gmapNode(buffHex(i)%NodeInd(:))
  end forall
  ! correct facets
  do i=1,nkFacet
    select case(buffFacet(i)%ShapeType)
      case(TRI_TYPE)
        buffFacet(i)%ShapeInd=gmapTri(buffFacet(i)%ShapeInd)
      case(QUAD_TYPE)
        buffFacet(i)%ShapeInd=gmapQuad(buffFacet(i)%ShapeInd)
      case default
        write(*,'(a,i2)'),'ERROR: unknown facet ShapeType: ',buffFacet(i)%ShapeType
    end select
    buffFacet(i)%NodeInd(1:buffFacet(i)%NodeNum)=&
    &           gmapNode(buffFacet(i)%NodeInd(1:buffFacet(i)%NodeNum))
    do j=1,FACET_NEIB_ELE_NUM
      if(buffFacet(i)%NeibEle(j)>0)then
        buffFacet(i)%NeibEle(j)=gmapEle(buffFacet(i)%NeibEle(j))
      else
        buffFacet(i)%NeibEle(j)=0
      end if
    end do
  end do
  ! correct elements
  do i=1,nkEle
    select case(buffEle(i)%ShapeType)
      case(TET_TYPE)
        buffEle(i)%ShapeInd=gmapTet(buffEle(i)%ShapeInd)
      case(HEX_TYPE)
        buffEle(i)%ShapeInd=gmapHex(buffEle(i)%ShapeInd)
      case default
        write(*,'(a,i2)'),'ERROR: unknown element ShapeType: ',buffEle(i)%ShapeType
    end select
    buffEle(i)%NodeInd(1:buffEle(i)%NodeNum)=gmapNode(buffEle(i)%NodeInd(1:buffEle(i)%NodeNum))
    do j=1,buffEle(i)%SurfNum
      if(buffEle(i)%Neib(j)>0)then
        buffEle(i)%Neib(j)=gmapEle(buffEle(i)%Neib(j))
      else
        buffEle(i)%Neib(j)=-gmapFacet(-buffEle(i)%Neib(j))
      end if
    end do
  end do
  
  ! send nodes
  call MPI_send(nkNode,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapNode,nkNode,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffNode,nkNode,typeNodeMPI,p,p,MPI_comm_world,errMPI)
  ! send points
  call MPI_send(nkPoint,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapPoint,nkPoint,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffPoint,nkPoint,typePointMPI,p,p,MPI_comm_world,errMPI)
  ! send lines
  call MPI_send(nkLine,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapLine,nkLine,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffLine,nkLine,typeLineMPI,p,p,MPI_comm_world,errMPI)
  ! send triangles
  call MPI_send(nkTri,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapTri,nkTri,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffTri,nkTri,typeTriMPI,p,p,MPI_comm_world,errMPI)
  ! send quadrilaterals
  call MPI_send(nkQuad,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapQuad,nkQuad,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffQuad,nkQuad,typeQuadMPI,p,p,MPI_comm_world,errMPI)
  ! send tetrahedrons
  call MPI_send(nkTet,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapTet,nkTet,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffTet,nkTet,typeTetMPI,p,p,MPI_comm_world,errMPI)
  ! send hexahedrons
  call MPI_send(nkHex,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapHex,nkHex,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffHex,nkHex,typeHexMPI,p,p,MPI_comm_world,errMPI)
  ! send facets
  call MPI_send(nkFacet,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapFacet,nkFacet,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffFacet,nkFacet,typeFacetMPI,p,p,MPI_comm_world,errMPI)
  ! send elements
  call MPI_send(nkEle,1,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(mapEle,nkEle,MPI_integer,p,p,MPI_comm_world,errMPI)
  call MPI_send(buffEle,nkEle,typeEleMPI,p,p,MPI_comm_world,errMPI)
  
  ! send conditions at nodes
  if(allocated(CondNode))then
    call MPI_send(.true.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
    do i=1,nkNode
      if(allocated(CondNode(mapNode(i))%Cond))then
        n=size(CondNode(mapNode(i))%Cond) ! number of conditions at i_th node of partition k
        call MPI_send(n,1,MPI_integer,p,p,MPI_comm_world,errMPI)
        do j=1,n
          call MPI_send(CondNode(mapNode(i))%Cond(j)%what,COND_NAME_LEN,MPI_character,&
          &             p,p,MPI_comm_world,errMPI)
          if(allocated(CondNode(mapNode(i))%Cond(j)%Val))then
            m=size(CondNode(mapNode(i))%Cond(j)%Val) ! number of data cells
            call MPI_send(m,1,MPI_integer,p,p,MPI_comm_world,errMPI)
            call MPI_send(CondNode(mapNode(i))%Cond(j)%Val,m,MPI_double_precision,&
            &             p,p,MPI_comm_world,errMPI)
          else
            call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
          end if
          if(allocated(CondNode(mapNode(i))%Cond(j)%Tab))then
            m=size(CondNode(mapNode(i))%Cond(j)%Tab) ! number of associated tables
            call MPI_send(m,1,MPI_integer,p,p,MPI_comm_world,errMPI)
            call MPI_send(CondNode(mapNode(i))%Cond(j)%Tab,m,MPI_integer,p,p,MPI_comm_world,errMPI)
          else
            call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
          end if
        end do
      else
        call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
      end if
    end do
  else
    call MPI_send(.false.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
  end if
  ! send conditions on facets
  if(allocated(CondFacet))then
    call MPI_send(.true.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
    do i=1,nkFacet
      if(allocated(CondFacet(mapFacet(i))%Cond))then
        n=size(CondFacet(mapFacet(i))%Cond) ! number of conditions at i_th facet of partition k
        call MPI_send(n,1,MPI_integer,p,p,MPI_comm_world,errMPI)
        do j=1,n
          call MPI_send(CondFacet(mapFacet(i))%Cond(j)%what,COND_NAME_LEN,MPI_character,&
          &             p,p,MPI_comm_world,errMPI)
          if(allocated(CondFacet(mapFacet(i))%Cond(j)%Val))then
            m=size(CondFacet(mapFacet(i))%Cond(j)%Val) ! number of data cells
            call MPI_send(m,1,MPI_integer,p,p,MPI_comm_world,errMPI)
            call MPI_send(CondFacet(mapFacet(i))%Cond(j)%Val,m,MPI_double_precision,&
            &             p,p,MPI_comm_world,errMPI)
          else
            call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
          end if
          if(allocated(CondFacet(mapFacet(i))%Cond(j)%Tab))then
            m=size(CondFacet(mapFacet(i))%Cond(j)%Tab) ! number of associated tables
            call MPI_send(m,1,MPI_integer,p,p,MPI_comm_world,errMPI)
            call MPI_send(CondFacet(mapFacet(i))%Cond(j)%Tab,m,MPI_integer,&
            &             p,p,MPI_comm_world,errMPI)
          else
            call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
          end if
        end do
      else
        call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
      end if
    end do
  else
    call MPI_send(.false.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
  end if
  ! send conditions in elements
  if(allocated(CondEle))then
    call MPI_send(.true.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
    do i=1,nkEle
      if(allocated(CondEle(mapEle(i))%Cond))then
        n=size(CondEle(mapEle(i))%Cond) ! number of conditions at i_th element of partition k
        call MPI_send(n,1,MPI_integer,p,p,MPI_comm_world,errMPI)
        do j=1,n
          call MPI_send(CondEle(mapEle(i))%Cond(j)%what,COND_NAME_LEN,MPI_character,&
          &             p,p,MPI_comm_world,errMPI)
          if(allocated(CondEle(mapEle(i))%Cond(j)%Val))then
            m=size(CondEle(mapEle(i))%Cond(j)%Val) ! number of data cells
            call MPI_send(m,1,MPI_integer,p,p,MPI_comm_world,errMPI)
            call MPI_send(CondEle(mapEle(i))%Cond(j)%Val,m,MPI_double_precision,&
            &             p,p,MPI_comm_world,errMPI)
          else
            call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
          end if
          if(allocated(CondEle(mapEle(i))%Cond(j)%Tab))then
            m=size(CondEle(mapEle(i))%Cond(j)%Tab) ! number of associated tables
            call MPI_send(m,1,MPI_integer,p,p,MPI_comm_world,errMPI)
            call MPI_send(CondEle(mapEle(i))%Cond(j)%Tab,m,MPI_integer,p,p,MPI_comm_world,errMPI)
          else
            call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
          end if
        end do
      else
        call MPI_send(0,1,MPI_integer,p,p,MPI_comm_world,errMPI)
      end if
    end do
  else
    call MPI_send(.false.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
  end if
  
  ! send scaler data on nodes
  allocate(buffDataScal(nkNode))
  do i=1,ntransNodeScal
    buffDataScal(:)=transNodeScal(i)%ptr(mapNode)
    call MPI_send(buffDataScal,nkNode,MPI_double_precision,p,p,MPI_comm_world,errMPI)
  end do
  deallocate(buffDataScal)
  ! send vector data on nodes
  allocate(buffDataVect(nkNode,DIMS))
  do i=1,ntransNodeVect
    buffDataVect(:,:)=transNodeVect(i)%ptr(mapNode,:)
    do j=1,DIMS
      call MPI_send(buffDataVect(:,j),nkNode,MPI_double_precision,p,p,MPI_comm_world,errMPI)
    end do
  end do
  deallocate(buffDataVect)
  ! send tensor data on nodes
  allocate(buffDataTens(nkNode,DIMS,DIMS))
  do i=1,ntransNodeTens
    buffDataTens(:,:,:)=transNodeTens(i)%ptr(mapNode,:,:)
    do j=1,DIMS
      do l=1,DIMS
        call MPI_send(buffDataTens(:,j,l),nkNode,MPI_double_precision,p,p,MPI_comm_world,errMPI)
      end do
    end do
  end do
  deallocate(buffDataTens)
  ! send scaler data in elements
  allocate(buffDataScal(nkEle))
  do i=1,ntransEleScal
    buffDataScal(:)=transEleScal(i)%ptr(mapEle)
    call MPI_send(buffDataScal,nkEle,MPI_double_precision,p,p,MPI_comm_world,errMPI)
  end do
  deallocate(buffDataScal)
  ! send vector data in elements
  allocate(buffDataVect(nkEle,DIMS))
  do i=1,ntransEleVect
    buffDataVect(:,:)=transEleVect(i)%ptr(mapEle,:)
    do j=1,DIMS
      call MPI_send(buffDataVect(:,j),nkEle,MPI_double_precision,p,p,MPI_comm_world,errMPI)
    end do
  end do
  deallocate(buffDataVect)
  ! send tensor data in elements
  allocate(buffDataTens(nkEle,DIMS,DIMS))
  do i=1,ntransEleTens
    buffDataTens(:,:,:)=transEleTens(i)%ptr(mapEle,:,:)
    do j=1,DIMS
      do l=1,DIMS
        call MPI_send(buffDataTens(:,j,l),nkEle,MPI_double_precision,p,p,MPI_comm_world,errMPI)
      end do
    end do
  end do
  deallocate(buffDataTens)
  
  ! clean up
  deallocate(buffNode,buffPoint,buffLine,buffTri,buffQuad,buffTet,buffHex,buffFacet,buffEle)
  deallocate(mapNode,mapPoint,mapLine,mapTri,mapQuad,mapTet,mapHex,mapFacet,mapEle)
end subroutine

!*******************
! receive partition
!*******************
! e.g.: to receive also the element data comes with the partition and save it to v
!   call recvPrt()
!   allocate(v(nEle))
!   v=transEleScal(1)%ptr
!   deallocate(transEleScal(1)%ptr)
subroutine recvPrt()
  use moduleGrid
  use moduleCond
  use moduleMPIvar
  
  logical flag
  
  ! receive nodes
  call MPI_recv(nNode,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Node(nNode))
  allocate(mapNode(nNode))
  call MPI_recv(mapNode,nNode,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Node,nNode,typeNodeMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive points
  call MPI_recv(nPoint,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Point(nPoint))
  allocate(mapPoint(nPoint))
  call MPI_recv(mapPoint,nPoint,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Point,nPoint,typePointMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive lines
  call MPI_recv(nLine,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Line(nLine))
  allocate(mapLine(nLine))
  call MPI_recv(mapLine,nLine,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Line,nLine,typeLineMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive triangles
  call MPI_recv(nTri,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Tri(nTri))
  allocate(mapTri(nTri))
  call MPI_recv(mapTri,nTri,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Tri,nTri,typeTriMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive quadrilaterals
  call MPI_recv(nQuad,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Quad(nQuad))
  allocate(mapQuad(nQuad))
  call MPI_recv(mapQuad,nQuad,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Quad,nQuad,typeQuadMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive tetrahedrons
  call MPI_recv(nTet,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Tet(nTet))
  allocate(mapTet(nTet))
  call MPI_recv(mapTet,nTet,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Tet,nTet,typeTetMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive hexahedrons
  call MPI_recv(nHex,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Hex(nHex))
  allocate(mapHex(nHex))
  call MPI_recv(mapHex,nHex,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Hex,nHex,typeHexMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive facets
  call MPI_recv(nFacet,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Facet(nFacet))
  allocate(mapFacet(nFacet))
  call MPI_recv(mapFacet,nFacet,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Facet,nFacet,typeFacetMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  ! receive elements
  call MPI_recv(nEle,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  allocate(Ele(nEle))
  allocate(mapEle(nEle))
  call MPI_recv(mapEle,nEle,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(Ele,nEle,typeEleMPI,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  
  ! receive conditions at nodes
  call MPI_recv(flag,1,MPI_logical,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  if(flag)then
    allocate(CondNode(nNode))
    do i=1,nNode
      call MPI_recv(n,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
      if(n>0)then
        allocate(CondNode(i)%Cond(n))
        do j=1,n
          call MPI_recv(CondNode(i)%Cond(j)%what,COND_NAME_LEN,MPI_character,&
          &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          call MPI_recv(m,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          if(m>0)then
            allocate(CondNode(i)%Cond(j)%Val(m))
            call MPI_recv(CondNode(i)%Cond(j)%Val,m,MPI_double_precision,&
            &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          end if
          call MPI_recv(m,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          if(m>0)then
            allocate(CondNode(i)%Cond(j)%Tab(m))
            call MPI_recv(CondNode(i)%Cond(j)%Tab,m,MPI_integer,&
            &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          end if
        end do
      end if
    end do
  end if
  ! receive conditions on facets
  call MPI_recv(flag,1,MPI_logical,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  if(flag)then
    allocate(CondFacet(nFacet))
    do i=1,nFacet
      call MPI_recv(n,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
      if(n>0)then
        allocate(CondFacet(i)%Cond(n))
        do j=1,n
          call MPI_recv(CondFacet(i)%Cond(j)%what,COND_NAME_LEN,MPI_character,&
          &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          call MPI_recv(m,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          if(m>0)then
            allocate(CondFacet(i)%Cond(j)%Val(m))
            call MPI_recv(CondFacet(i)%Cond(j)%Val,m,MPI_double_precision,&
            &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          end if
          call MPI_recv(m,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          if(m>0)then
            allocate(CondFacet(i)%Cond(j)%Tab(m))
            call MPI_recv(CondFacet(i)%Cond(j)%Tab,m,MPI_integer,&
            &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          end if
        end do
      end if
    end do
  end if
  ! receive conditions in elements
  call MPI_recv(flag,1,MPI_logical,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  if(flag)then
    allocate(CondEle(nEle))
    do i=1,nEle
      call MPI_recv(n,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
      if(n>0)then
        allocate(CondEle(i)%Cond(n))
        do j=1,n
          call MPI_recv(CondEle(i)%Cond(j)%what,COND_NAME_LEN,MPI_character,&
          &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          call MPI_recv(m,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          if(m>0)then
            allocate(CondEle(i)%Cond(j)%Val(m))
            call MPI_recv(CondEle(i)%Cond(j)%Val,m,MPI_double_precision,&
            &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          end if
          call MPI_recv(m,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          if(m>0)then
            allocate(CondEle(i)%Cond(j)%Tab(m))
            call MPI_recv(CondEle(i)%Cond(j)%Tab,m,MPI_integer,&
            &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
          end if
        end do
      end if
    end do
  end if
  
  ! receive scaler data at nodes
  do i=1,ntransNodeScal
    allocate(transNodeScal(i)%ptr(nNode))
    call MPI_recv(transNodeScal(i)%ptr,nNode,MPI_double_precision,&
    &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  end do
  ! receive vector data at nodes
  do i=1,ntransNodeVect
    allocate(transNodeVect(i)%ptr(nNode,DIMS))
    do j=1,DIMS
      call MPI_recv(transNodeVect(i)%ptr(:,j),nNode,MPI_double_precision,&
      &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
    end do
  end do
  ! receive tensor data at nodes
  do i=1,ntransNodeTens
    allocate(transNodeTens(i)%ptr(nNode,DIMS,DIMS))
    do j=1,DIMS
      do k=1,DIMS
        call MPI_recv(transNodeTens(i)%ptr(:,j,k),nNode,MPI_double_precision,&
        &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
      end do
    end do
  end do
  ! receive scaler data in elements
  do i=1,ntransEleScal
    allocate(transEleScal(i)%ptr(nEle))
    call MPI_recv(transEleScal(i)%ptr,nEle,MPI_double_precision,&
    &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
  end do
  ! receive vector data in elements
  do i=1,ntransEleVect
    allocate(transEleVect(i)%ptr(nEle,DIMS))
    do j=1,DIMS
      call MPI_recv(transEleVect(i)%ptr(:,j),nEle,MPI_double_precision,&
      &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
    end do
  end do
  ! receive tensor data in elements
  do i=1,ntransEleTens
    allocate(transEleTens(i)%ptr(nEle,DIMS,DIMS))
    do j=1,DIMS
      do k=1,DIMS
        call MPI_recv(transEleTens(i)%ptr(:,j,k),nEle,MPI_double_precision,&
        &             ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
      end do
    end do
  end do
end subroutine

!*************
! return data
!*************
subroutine retnData()
  use moduleMPIvar
  use moduleGrid
  
  ! send auxiliary data
  call MPI_send(nNode,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,errMPI)
  call MPI_send(nEle,1,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,errMPI)
  call MPI_send(mapNode,nNode,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,errMPI)
  call MPI_send(mapEle,nEle,MPI_integer,ROOT_PID,pidMPI,MPI_comm_world,errMPI)
  
  ! send scaler data on nodes
  do i=1,ntransNodeScal
    call MPI_send(transNodeScal(i)%ptr,nNode,MPI_double_precision,&
    &             ROOT_PID,pidMPI,MPI_comm_world,errMPI)
  end do
  ! send vector data on nodes
  do i=1,ntransNodeVect
    do j=1,DIMS
      call MPI_send(transNodeVect(i)%ptr(:,j),nNode,MPI_double_precision,&
      &             ROOT_PID,pidMPI,MPI_comm_world,errMPI)
    end do
  end do
  ! send tensor data on nodes
  do i=1,ntransNodeTens
    do j=1,DIMS
      do k=1,DIMS
        call MPI_send(transNodeTens(i)%ptr(:,j,k),nNode,MPI_double_precision,&
        &             ROOT_PID,pidMPI,MPI_comm_world,errMPI)
      end do
    end do
  end do
  ! send scaler data in elements
  do i=1,ntransEleScal
    call MPI_send(transEleScal(i)%ptr,nEle,MPI_double_precision,&
    &             ROOT_PID,pidMPI,MPI_comm_world,errMPI)
  end do
  ! send vector data in elements
  do i=1,ntransEleVect
    do j=1,DIMS
      call MPI_send(transEleVect(i)%ptr(:,j),nEle,MPI_double_precision,&
      &             ROOT_PID,pidMPI,MPI_comm_world,errMPI)
    end do
  end do
  ! send tensor data in elements
  do i=1,ntransEleTens
    do j=1,DIMS
      do k=1,DIMS
        call MPI_send(transEleTens(i)%ptr(:,j,k),nEle,MPI_double_precision,&
        &             ROOT_PID,pidMPI,MPI_comm_world,errMPI)
      end do
    end do
  end do
  
end subroutine

!**************************************************************************
! try to gather data from any process, received partition t from process p
!**************************************************************************
subroutine gathData(t,p)
  use moduleMPIvar
  use moduleGrid
  
  integer,intent(out)::t,p
  integer nkNode,nkEle
  double precision,allocatable::buffDataScal(:),buffDataVect(:,:),buffDataTens(:,:,:)
  
  ! receive auxiliary data
  call MPI_recv(nkNode,1,MPI_integer,MPI_any_source,MPI_any_tag,MPI_comm_world,statMPI,errMPI)
  p=statMPI(MPI_source)
  call MPI_recv(nkEle,1,MPI_integer,p,p,MPI_comm_world,statMPI,errMPI)
  allocate(mapNode(nkNode))
  allocate(mapEle(nkEle))
  call MPI_recv(mapNode,nkNode,MPI_integer,p,p,MPI_comm_world,statMPI,errMPI)
  call MPI_recv(mapEle,nkEle,MPI_integer,p,p,MPI_comm_world,statMPI,errMPI)
  t=Ele(mapEle(1))%Prt
    
  allocate(buffDataScal(nNode))
  allocate(buffDataVect(nNode,DIMS))
  allocate(buffDataTens(nNode,DIMS,DIMS))
  ! receive scaler data at nodes
  do i=1,ntransNodeScal
    call MPI_recv(buffDataScal,nkNode,MPI_double_precision,p,p,MPI_comm_world,statMPI,errMPI)
    transNodeScal(i)%ptr(mapNode)=buffDataScal
  end do
  ! receive vector data at nodes
  do i=1,ntransNodeVect
    do j=1,DIMS
      call MPI_recv(buffDataVect(:,j),nkNode,MPI_double_precision,p,p,MPI_comm_world,statMPI,errMPI)
    end do
    transNodeVect(i)%ptr(mapNode,:)=buffDataVect(:,:)
  end do
  ! receive tensor data at nodes
  do i=1,ntransNodeTens
    do j=1,DIMS
      do k=1,DIMS
        call MPI_recv(buffDataTens(:,j,k),nkNode,MPI_double_precision,&
        &             p,p,MPI_comm_world,statMPI,errMPI)
      end do
    end do
    transNodeTens(i)%ptr(mapNode,:,:)=buffDataTens(:,:,:)
  end do
  deallocate(buffDataScal,buffDataVect,buffDataTens)
  
  allocate(buffDataScal(nEle))
  allocate(buffDataVect(nEle,DIMS))
  allocate(buffDataTens(nEle,DIMS,DIMS))
  ! receive scaler data in elements
  do i=1,ntransEleScal
    call MPI_recv(buffDataScal,nkEle,MPI_double_precision,p,p,MPI_comm_world,statMPI,errMPI)
    transEleScal(i)%ptr(mapEle)=buffDataScal
  end do
  ! receive vector data in elements
  do i=1,ntransEleVect
    do j=1,DIMS
      call MPI_recv(buffDataVect(:,j),nkEle,MPI_double_precision,p,p,MPI_comm_world,statMPI,errMPI)
    end do
    transEleVect(i)%ptr(mapEle,:)=buffDataVect(:,:)
  end do
  ! receive tensor data in elements
  do i=1,ntransEleTens
    do j=1,DIMS
      do k=1,DIMS
        call MPI_recv(buffDataTens(:,j,k),nkEle,MPI_double_precision,&
        &             p,p,MPI_comm_world,statMPI,errMPI)
      end do
    end do
    transEleTens(i)%ptr(mapEle,:,:)=buffDataTens(:,:,:)
  end do
  deallocate(buffDataScal,buffDataVect,buffDataTens)
  
  ! clean up
  deallocate(mapNode,mapEle)
end subroutine