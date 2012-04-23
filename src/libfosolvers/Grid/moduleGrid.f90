!----------------------------------------------------------------------------- best with 100 columns

!> grid
module moduleGrid
  use moduleBasicDataStruct
  private
  
  ! constants
  integer,parameter,public::DIMS=3 !< dimensions
  
  integer,parameter,public::POINT_TYPE=15 !< point type code
  integer,parameter,public::LINE_TYPE=1  !< line type code
  integer,parameter,public::TRI_TYPE=2 !< tri type code
  integer,parameter,public::QUAD_TYPE=3 !< quad type code
  integer,parameter,public::TET_TYPE=4 !< tet type code
  integer,parameter,public::HEX_TYPE=5 !< hex type code
  
  integer,parameter,public::POINT_NODE_NUM=1 !< number of nodes per point
  integer,parameter,public::LINE_NODE_NUM=2 !< number of nodes per line
  integer,parameter,public::TRI_NODE_NUM=3 !< number of nodes per tri
  integer,parameter,public::QUAD_NODE_NUM=4 !< number of nodes per quad
  integer,parameter,public::TET_NODE_NUM=4 !< number of nodes per tet
  integer,parameter,public::HEX_NODE_NUM=8 !< number of nodes per hex
  
  !> grid data and procedures
  type,public::typeGrid
    ! basic grid data
    integer nNode,nPoint,nLine,nFacet,nBlock
    type(typeListIScal)::PhyEntiList !< list of physical entities
    type(typeListIScal)::PrtList !< list of partitions
    double precision,allocatable::NodePos(:,:) !< node position
    integer,allocatable::PointNodeInd(:) !< node index of the points
    type(typeListIScal),allocatable::LineNodeInd(:) !< node index of the lines
    type(typeListIScal),allocatable::FacetNodeInd(:) !< node index of the facets
    type(typeListIScal),allocatable::BlockNodeInd(:) !< node index of the blocks
    integer,allocatable::PointGeoEnti(:) !< geometric entity of the points
    integer,allocatable::LineGeoEnti(:) !< geometric entity of the lines
    integer,allocatable::FacetGeoEnti(:) !< geometric entity of the facets
    integer,allocatable::BlockGeoEnti(:) !< geometric entity of the blocks
    type(typeListIScal),allocatable::PointPhyEnti(:) !< physical entity of the points
    type(typeListIScal),allocatable::LinePhyEnti(:) !< physical entity of the lines
    type(typeListIScal),allocatable::FacetPhyEnti(:) !< physical entity of the facets
    type(typeListIScal),allocatable::BlockPhyEnti(:) !< physical entity of the blocks
    type(typeListIScal),allocatable::PointPrt(:) !< partition of the points
    type(typeListIScal),allocatable::LinePrt(:) !< partition of the lines
    type(typeListIScal),allocatable::FacetPrt(:) !< partition of the facets
    type(typeListIScal),allocatable::BlockPrt(:) !< partition of the blocks
  contains
    procedure,public::clear=>clearGrid
    ! grid input/output
!    procedure,public::readGMSH=>readGridGMSH
!    procedure,public::writeGMSH=>writeGridGMSH
    ! procedures using basic data
!    procedure,public::getNodePos
!    procedure,public::getPointNodeInd
!    procedure,public::getPointGeoEnti
!    procedure,public::getPointPhyEnti
!    procedure,public::getPointPrt
!    procedure,public::getLineNodeInd
!    procedure,public::getLineGeoEnti
!    procedure,public::getLinePhyEnti
!    procedure,public::getLinePrt
!    procedure,public::getFacetNodeInd
!    procedure,public::getFacetGeoEnti
!    procedure,public::getFacetPhyEnti
!    procedure,public::getFacetPrt
!    procedure,public::getBlockNodeInd
!    procedure,public::getBlockGeoEnti
!    procedure,public::getBlockPhyEnti
!    procedure,public::getBlockPrt
    !FIXME:final::purgeGrid
  end type
  
contains
  
  !> clear this Grid
  subroutine clearGrid(this)
    class(typeGrid),intent(inout)::this !< this Grid
    
    if(allocated(this%NodePos)) deallocate(this%NodePos)
    if(allocated(this%PointNodeInd)) deallocate(this%PointNodeInd)
    if(allocated(this%LineNodeInd)) deallocate(this%LineNodeInd)
    if(allocated(this%FacetNodeInd)) deallocate(this%FacetNodeInd)
    if(allocated(this%BlockNodeInd)) deallocate(this%BlockNodeInd)
    if(allocated(this%PointGeoEnti)) deallocate(this%PointGeoEnti)
    if(allocated(this%LineGeoEnti)) deallocate(this%LineGeoEnti)
    if(allocated(this%FacetGeoEnti)) deallocate(this%FacetGeoEnti)
    if(allocated(this%BlockGeoEnti)) deallocate(this%BlockGeoEnti)
    if(allocated(this%PointPhyEnti)) deallocate(this%PointPhyEnti)
    if(allocated(this%LinePhyEnti)) deallocate(this%LinePhyEnti)
    if(allocated(this%FacetPhyEnti)) deallocate(this%FacetPhyEnti)
    if(allocated(this%BlockPhyEnti)) deallocate(this%BlockPhyEnti)
    if(allocated(this%PointPrt)) deallocate(this%PointPrt)
    if(allocated(this%LinePrt)) deallocate(this%LinePrt)
    if(allocated(this%FacetPrt)) deallocate(this%FacetPrt)
    if(allocated(this%BlockPrt)) deallocate(this%BlockPrt)
  end subroutine
  
  !> read grid data from opened file id into this Grid
  subroutine readGridGMSH(this,id)
    use moduleBasicDataStruct
    class(typeGrid),intent(inout)::this !< this Grid
    integer,intent(in)::id !< the file id
    integer ierr
    integer,parameter::DFLT_STR_LEN=400
    character(DFLT_STR_LEN)::tempStr
    integer nEle,shapetype,nt,np,tempIVect(DFLT_STR_LEN/2)
    integer,allocatable::tempPointNodeInd(:)
    type(typeListIScal),allocatable::tempLineNodeInd(:)
    type(typeListIScal),allocatable::tempFacetNodeInd(:)
    type(typeListIScal),allocatable::tempBlockNodeInd(:)
    integer,allocatable::tempPointGeoEnti(:)
    integer,allocatable::tempLineGeoEnti(:)
    integer,allocatable::tempFacetGeoEnti(:)
    integer,allocatable::tempBlockGeoEnti(:)
    type(typeListIScal),allocatable::tempPointPhyEnti(:)
    type(typeListIScal),allocatable::tempLinePhyEnti(:)
    type(typeListIScal),allocatable::tempFacetPhyEnti(:)
    type(typeListIScal),allocatable::tempBlockPhyEnti(:)
    type(typeListIScal),allocatable::tempPointPrt(:)
    type(typeListIScal),allocatable::tempLinePrt(:)
    type(typeListIScal),allocatable::tempFacetPrt(:)
    type(typeListIScal),allocatable::tempBlockPrt(:)
    
    call this%clear()
    ierr=0
    rewind(id,iostat=ierr)
    
    do while(ierr==0)
      ! skip the irrelevant lines
      do while(ierr==0)
        read(id,*,iostat=ierr),tempStr
        if(tempStr(1:1)=='$') exit
      end do
      ! check if finished
      if(ierr/=0) exit
      ! read node data
      if(tempStr(1:6)=='$Nodes')then
        read(id,*,iostat=ierr),this%nNode
        allocate(this%NodePos(DIMS,this%nNode))
        do i=1,this%nNode
          read(id,*,iostat=ierr),j,this%NodePos(:,i)
        end do
        cycle
      end if
      ! read element (includes facet, line & point) data
      if(tempStr(1:9)=='$Elements')then
        read(id,*,iostat=ierr),nEle
        allocate(tempPointNodeInd(nEle))
        allocate(tempLineNodeInd(nEle))
        allocate(tempFacetNodeInd(nEle))
        allocate(tempBlockNodeInd(nEle))
        allocate(tempPointGeoEnti(nEle))
        allocate(tempLineGeoEnti(nEle))
        allocate(tempFacetGeoEnti(nEle))
        allocate(tempBlockGeoEnti(nEle))
        allocate(tempPointPhyEnti(nEle))
        allocate(tempLinePhyEnti(nEle))
        allocate(tempFacetPhyEnti(nEle))
        allocate(tempBlockPhyEnti(nEle))
        allocate(tempPointPrt(nEle))
        allocate(tempLinePrt(nEle))
        allocate(tempFacetPrt(nEle))
        allocate(tempBlockPrt(nEle))
        do i=1,nEle
          read(id,*,iostat=ierr),tempStr
          read(tempStr,*),j,shapetype,nt ! read the first 3 integers
          np=0
          select case(shapetype)
            case(POINT_TYPE)
              np=POINT_NODE_NUM
            case(LINE_TYPE)
              np=LINE_NODE_NUM
            case(TRI_TYPE)
              np=TRI_NODE_NUM
            case(QUAD_TYPE)
              np=QUAD_NODE_NUM
            case(TET_TYPE)
              np=TET_NODE_NUM
            case(HEX_TYPE)
              np=HEX_NODE_NUM
          end select
          read(tempStr,*),tempIVect(1:3+nt+np)
          if(tempIVect(2)==POINT_TYPE)then ! save point
            this%nPoint=this%nPoint+1
            tempPointNodeInd(this%nPoint)=tempIVect(3+nt+np)
            if(nt>=2)then
              tempPointGeoEnti(this%nPoint)=tempIVect(5)
              call tempPointPhyEnti(this%nPoint)%push(tempIVect(4))
              if(nt>=4)then
                call tempPointPrt(this%nPoint)%push(tempIVect(7:6+tempIVect(6)))
              else
                call tempPointPrt(this%nPoint)%push(0)
              end if
            else
              tempPointGeoEnti(this%nPoint)=0
              call tempPointPhyEnti(this%nPoint)%push(0)
              call tempPointPrt(this%nPoint)%push(0)
            end if
          else if(tempIVect(2)==LINE_TYPE)then ! save line
            this%nLine=this%nLine+1
            call tempLineNodeInd(this%nLine)%push(tempIVect(3+nt+1:3+nt+np))
            if(nt>=2)then
              tempLineGeoEnti(this%nLine)=tempIVect(5)
              call tempLinePhyEnti(this%nLine)%push(tempIVect(4))
              if(nt>=4)then
                call tempLinePrt(this%nLine)%push(tempIVect(7:6+tempIVect(6)))
              else
                call tempLinePrt(this%nLine)%push(0)
              end if
            else
              tempLineGeoEnti(this%nLine)=0
              call tempLinePhyEnti(this%nLine)%push(0)
              call tempLinePrt(this%nLine)%push(0)
            end if
          else if(tempIVect(2)==TRI_TYPE.or.tempIVect(2)==QUAD_TYPE)then ! save facet
            this%nFacet=this%nFacet+1
            call tempFacetNodeInd(this%nFacet)%push(tempIVect(3+nt+1:3+nt+np))
            if(nt>=2)then
              tempFacetGeoEnti(this%nFacet)=tempIVect(5)
              call tempFacetPhyEnti(this%nFacet)%push(tempIVect(4))
              if(nt>=4)then
                call tempFacetPrt(this%nFacet)%push(tempIVect(7:6+tempIVect(6)))
              else
                call tempFacetPrt(this%nFacet)%push(0)
              end if
            else
              tempFacetGeoEnti(this%nFacet)=0
              call tempFacetPhyEnti(this%nFacet)%push(0)
              call tempFacetPrt(this%nFacet)%push(0)
            end if
          else if(tempIVect(2)==TET_TYPE.or.tempIVect(2)==HEX_TYPE)then ! save block
            this%nBlock=this%nBlock+1
            call tempBlockNodeInd(this%nBlock)%push(tempIVect(3+nt+1:3+nt+np))
            if(nt>=2)then
              tempBlockGeoEnti(this%nBlock)=tempIVect(5)
              call tempBlockPhyEnti(this%nBlock)%push(tempIVect(4))
              if(nt>=4)then
                call tempBlockPrt(this%nBlock)%push(tempIVect(7:6+tempIVect(6)))
              else
                call tempBlockPrt(this%nBlock)%push(0)
              end if
            else
              tempBlockGeoEnti(this%nBlock)=0
              call tempBlockPhyEnti(this%nBlock)%push(0)
              call tempBlockPrt(this%nBlock)%push(0)
            end if
          end if
        end do
!        ! trim Point
!        if(nPoint>0)then
!          allocate(tempPoint(nPoint))
!          tempPoint(:)=Point(1:nPoint)
!          call move_alloc(tempPoint,Point)
!        else
!          deallocate(Point)
!        end if
!        ! trim Line
!        if(nLine>0)then
!          allocate(tempLine(nLine))
!          tempLine(:)=Line(1:nLine)
!          call move_alloc(tempLine,Line)
!        else
!          deallocate(Line)
!        end if
!        ! trim Facet
!        if(nFacet>0)then
!          allocate(tempFacet(nFacet))
!          tempFacet(:)=Facet(1:nFacet)
!          call move_alloc(tempFacet,Facet)
!        else
!          deallocate(Facet)
!        end if
!        ! trim Block
!        if(nBlock>0)then
!          allocate(tempBlock(nBlock))
!          tempBlock(:)=Block(1:nBlock)
!          call move_alloc(tempBlock,Block)
!        else
!          deallocate(Block)
!        end if
        exit
      end if
    end do
  end subroutine
  
end module
