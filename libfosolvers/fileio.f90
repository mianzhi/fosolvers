!----------------------------------------------------------------------------- best with 100 columns

!***********************
! file input and output
!***********************
module moduleFileIO
  use moduleMiscDataStruct
  private
  
  ! some constants
  integer,parameter,public::DEFAULT_STRING_LEN=200
  
  integer,parameter,public::RST_BIND_NODE=1
  integer,parameter,public::RST_BIND_FACET=2
  integer,parameter,public::RST_BIND_BLOCK=3
  
  ! output control related variables
  integer,public,save::nWrite=0
  ! data to be output
  type(typePtrScalArray),allocatable,public,save::rstNodeScal(:),rstFacetScal(:),rstBlockScal(:)
  type(typePtrVectArray),allocatable,public,save::rstNodeVect(:),rstFacetVect(:),rstBlockVect(:)
  type(typePtrTensArray),allocatable,public,save::rstNodeTens(:),rstFacetTens(:),rstBlockTens(:)
  
  ! generic add write list
  interface addWrite
    module procedure::addWriteScal
    module procedure::addWriteVect
    module procedure::addWriteTens
  end interface
  public addWrite
  
  ! other procedures
  public readMsh
  public writeRst
  public readMtl
  
contains
  
  !----------------
  ! read mesh file
  !----------------
  subroutine readMsh(fname,fid,verbose)
    use moduleGrid
    use moduleUtility
    character(*),intent(in)::fname
    integer,intent(in)::fid
    logical,optional,intent(in)::verbose
    integer readerr,nEle,shapetype,np,nt,tempIntVect1(50)
    character(DEFAULT_STRING_LEN) tempString
    logical info
    type(typePoint),allocatable::tempPoint(:)
    type(typeLine),allocatable::tempLine(:)
    type(typeFacet),allocatable::tempFacet(:)
    type(typeBlock),allocatable::tempBlock(:)
    
    if(.not.present(verbose))then
      info=.false.
    else
      info=verbose
    end if
    
    if(info)then
      call showNoAdv('reading msh file...')
    end if
    
    np=0
    
    open(fid,file=fname,status='old')
    readerr=0
    
    do while(readerr==0)
      ! skip the irrelevant lines
      do while(readerr==0)
        read(fid,*,iostat=readerr),tempString
        if(tempString(1:1)=='$')then
          exit
        end if
      end do
      ! check if finished
      if(readerr/=0)then
        exit
      end if
      ! read node data
      if(tempString(1:4)=='$NOD'.or.tempString(1:4)=='$Nod')then
        read(fid,*,iostat=readerr),nNode
        allocate(Node(nNode))
        do i=1,nNode
          read(fid,*,iostat=readerr),j,Node(i)%Pos(:)
          if(i/=j)then
            call showWarning('node data may be not in sequence.')
          end if
        end do
        cycle
      end if
      ! read element (includes facet, line & point) data
      if(tempString(1:4)=='$ELM'.or.tempString(1:4)=='$Ele')then
        if(tempString(1:4)=='$ELM')then
          call showWarning('not capable with GMSH version 1 format.')
        end if
        read(fid,*,iostat=readerr),nEle
        ! Note: allocate sufficient space to each object array
        allocate(Point(nEle))
        allocate(Line(nEle))
        allocate(Facet(nEle))
        allocate(Block(nEle))
        nPoint=0
        nLine=0
        nFacet=0
        nBlock=0
        do i=1,nEle
          read(fid,'(a)',iostat=readerr),tempString
          read(tempString,*),j,shapetype,nt ! read the first 3 integers
          if(i/=j)then
            call showWarning('element data may be not in sequence.')
          end if
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
          read(tempString,*),tempIntVect1(1:3+nt+np)
          select case(tempIntVect1(2)) ! type of element
            case(POINT_TYPE) ! point
              nPoint=nPoint+1
              Point(nPoint)%NodeInd=tempIntVect1(3+nt+1)
              if(nt>=2)then
                Point(nPoint)%GeoEnti=tempIntVect1(5)
              else
                Point(nPoint)%GeoEnti=0
              end if
            case(LINE_TYPE) ! line
              nLine=nLine+1
              Line(nLine)%NodeInd(:)=tempIntVect1(3+nt+1:3+nt+np)
              if(nt>=2)then
                Line(nLine)%GeoEnti=tempIntVect1(5)
              else
                Line(nLine)%GeoEnti=0
              end if
            case(TRI_TYPE) ! tri
              nFacet=nFacet+1
              call Facet(nFacet)%specify(TRI_TYPE)
              Facet(nFacet)%NodeInd(:)=tempIntVect1(3+nt+1:3+nt+np)
              if(nt>=2)then
                Facet(nFacet)%GeoEnti=tempIntVect1(5)
              else
                Facet(nFacet)%GeoEnti=0
              end if
            case(QUAD_TYPE) ! quad
              nFacet=nFacet+1
              call Facet(nFacet)%specify(QUAD_TYPE)
              Facet(nFacet)%NodeInd(:)=tempIntVect1(3+nt+1:3+nt+np)
              if(nt>=2)then
                Facet(nFacet)%GeoEnti=tempIntVect1(5)
              else
                Facet(nFacet)%GeoEnti=0
              end if
            case(TET_TYPE) ! tet
              nBlock=nBlock+1
              call Block(nBlock)%specify(TET_TYPE)
              Block(nBlock)%NodeInd(:)=tempIntVect1(3+nt+1:3+nt+np)
              if(nt>=2)then
                Block(nBlock)%GeoEnti=tempIntVect1(5)
                if(nt>=4)then
                  allocate(Block(nBlock)%Prt(tempIntVect1(6)))
                  Block(nBlock)%Prt(:)=tempIntVect1(7:6+tempIntVect1(6))
                else
                  allocate(Block(nBlock)%Prt(1))
                  Block(nBlock)%Prt=1
                end if
              else
                Block(nBlock)%GeoEnti=0
                allocate(Block(nBlock)%Prt(1))
                Block(nBlock)%Prt=1
              end if
            case(HEX_TYPE) ! hex
              nBlock=nBlock+1
              call Block(nBlock)%specify(HEX_TYPE)
              Block(nBlock)%NodeInd(:)=tempIntVect1(3+nt+1:3+nt+np)
              if(nt>=2)then
                Block(nBlock)%GeoEnti=tempIntVect1(5)
                if(nt>=4)then
                  allocate(Block(nBlock)%Prt(tempIntVect1(6)))
                  Block(nBlock)%Prt(:)=tempIntVect1(7:6+tempIntVect1(6))
                else
                  allocate(Block(nBlock)%Prt(1))
                  Block(nBlock)%Prt=1
                end if
              else
                Block(nBlock)%GeoEnti=0
                allocate(Block(nBlock)%Prt(1))
                Block(nBlock)%Prt=1
              end if
          end select
        end do
        ! trim Point
        if(nPoint>0)then
          allocate(tempPoint(nPoint))
          tempPoint(:)=Point(1:nPoint)
          call move_alloc(tempPoint,Point)
        else
          deallocate(Point)
        end if
        ! trim Line
        if(nLine>0)then
          allocate(tempLine(nLine))
          tempLine(:)=Line(1:nLine)
          call move_alloc(tempLine,Line)
        else
          deallocate(Line)
        end if
        ! trim Facet
        if(nFacet>0)then
          allocate(tempFacet(nFacet))
          tempFacet(:)=Facet(1:nFacet)
          call move_alloc(tempFacet,Facet)
        else
          deallocate(Facet)
        end if
        ! trim Block
        if(nBlock>0)then
          allocate(tempBlock(nBlock))
          tempBlock(:)=Block(1:nBlock)
          call move_alloc(tempBlock,Block)
        else
          deallocate(Block)
        end if
        exit
      end if
    end do
    
    ! auxiliary operations
    if(info)then
      call showNoAdv('processing grid...')
    end if
    call updateGrid()
    
    if(info)then
      call showNoAdv('grid summary:')
      write(*,'(a)'),''
      write(*,'(a,i8)'),'  number of nodes:        ',nNode
      write(*,'(a,i8)'),'  number of points:       ',nPoint
      write(*,'(a,i8)'),'  number of lines:        ',nLine
      write(*,'(a,i8)'),'  number of facets:       ',nFacet
      write(*,'(a,i8)'),'  number of blocks:       ',nBlock
      write(*,'(a,i8,/)'),'  number of partitions.:  ',nPrt
    end if
    
    close(fid)
  end subroutine
  
  !-------------------------------
  ! add scaler data to write list
  !-------------------------------
  subroutine addWriteScal(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:)
    integer,intent(in),optional::binding
    logical bindnode,bindfacet,bindblock
    type(typePtrScalArray),allocatable::temp(:)
    
    bindnode=.false.
    bindfacet=.false.
    bindblock=.false.
    if(present(binding))then
      select case(binding)
        case(RST_BIND_NODE)
          bindnode=.true.
        case(RST_BIND_FACET)
          bindfacet=.true.
        case(RST_BIND_BLOCK)
          bindblock=.true.
        case default
      end select
    else
      bindnode=(size(v,1)==nNode)
      bindfacet=(size(v,1)==nFacet)
      bindblock=(size(v,1)==nBlock)
    end if
    
    if(bindnode)then
      if(allocated(rstNodeScal))then
        allocate(temp(size(rstNodeScal)+1))
        temp(1:size(rstNodeScal))=rstNodeScal(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstNodeScal)
      else
        allocate(rstNodeScal(1))
        rstNodeScal(1)%ptr=>v
      end if
    end if
    if(bindfacet)then
      if(allocated(rstFacetScal))then
        allocate(temp(size(rstFacetScal)+1))
        temp(1:size(rstFacetScal))=rstFacetScal(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstFacetScal)
      else
        allocate(rstFacetScal(1))
        rstFacetScal(1)%ptr=>v
      end if
    end if
    if(bindblock)then
      if(allocated(rstBlockScal))then
        allocate(temp(size(rstBlockScal)+1))
        temp(1:size(rstBlockScal))=rstBlockScal(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstBlockScal)
      else
        allocate(rstBlockScal(1))
        rstBlockScal(1)%ptr=>v
      end if
    end if
    
    nWrite=0
  end subroutine
  
  !-------------------------------
  ! add vector data to write list
  !-------------------------------
  subroutine addWriteVect(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:,:)
    integer,intent(in),optional::binding
    logical bindnode,bindfacet,bindblock
    type(typePtrVectArray),allocatable::temp(:)
    
    bindnode=.false.
    bindfacet=.false.
    bindblock=.false.
    if(present(binding))then
      select case(binding)
        case(RST_BIND_NODE)
          bindnode=.true.
        case(RST_BIND_FACET)
          bindfacet=.true.
        case(RST_BIND_BLOCK)
          bindblock=.true.
        case default
      end select
    else
      bindnode=(size(v,1)==nNode)
      bindfacet=(size(v,1)==nFacet)
      bindblock=(size(v,1)==nBlock)
    end if
    
    if(bindnode)then
      if(allocated(rstNodeVect))then
        allocate(temp(size(rstNodeVect)+1))
        temp(1:size(rstNodeVect))=rstNodeVect(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstNodeVect)
      else
        allocate(rstNodeVect(1))
        rstNodeVect(1)%ptr=>v
      end if
    end if
    if(bindfacet)then
      if(allocated(rstFacetVect))then
        allocate(temp(size(rstFacetVect)+1))
        temp(1:size(rstFacetVect))=rstFacetVect(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstFacetVect)
      else
        allocate(rstFacetVect(1))
        rstFacetVect(1)%ptr=>v
      end if
    end if
    if(bindblock)then
      if(allocated(rstBlockVect))then
        allocate(temp(size(rstBlockVect)+1))
        temp(1:size(rstBlockVect))=rstBlockVect(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstBlockVect)
      else
        allocate(rstBlockVect(1))
        rstBlockVect(1)%ptr=>v
      end if
    end if
    
    nWrite=0
  end subroutine
  
  !-------------------------------
  ! add tensor data to write list
  !-------------------------------
  subroutine addWriteTens(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:,:,:)
    integer,intent(in),optional::binding
    logical bindnode,bindfacet,bindblock
    type(typePtrTensArray),allocatable::temp(:)
    
    bindnode=.false.
    bindfacet=.false.
    bindblock=.false.
    if(present(binding))then
      select case(binding)
        case(RST_BIND_NODE)
          bindnode=.true.
        case(RST_BIND_FACET)
          bindfacet=.true.
        case(RST_BIND_BLOCK)
          bindblock=.true.
        case default
      end select
    else
      bindnode=(size(v,1)==nNode)
      bindfacet=(size(v,1)==nFacet)
      bindblock=(size(v,1)==nBlock)
    end if
    
    if(bindnode)then
      if(allocated(rstNodeTens))then
        allocate(temp(size(rstNodeTens)+1))
        temp(1:size(rstNodeTens))=rstNodeTens(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstNodeTens)
      else
        allocate(rstNodeTens(1))
        rstNodeTens(1)%ptr=>v
      end if
    end if
    if(bindfacet)then
      if(allocated(rstFacetTens))then
        allocate(temp(size(rstFacetTens)+1))
        temp(1:size(rstFacetTens))=rstFacetTens(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstFacetTens)
      else
        allocate(rstFacetTens(1))
        rstFacetTens(1)%ptr=>v
      end if
    end if
    if(bindblock)then
      if(allocated(rstBlockTens))then
        allocate(temp(size(rstBlockTens)+1))
        temp(1:size(rstBlockTens))=rstBlockTens(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,rstBlockTens)
      else
        allocate(rstBlockTens(1))
        rstBlockTens(1)%ptr=>v
      end if
    end if
    
    nWrite=0
  end subroutine
  
  !---------------
  ! write results
  !---------------
  ! hold=.true.: hold the file when exiting subroutine; good for writing a time span
  ! hold=.false.: close the file when exiting subroutine; good for writing a snapshot
  subroutine writeRst(fname,fid,span)
    use moduleGrid
    use moduleTime
    character(*),intent(in)::fname
    integer,intent(in)::fid
    logical,optional,intent(in)::span
    character(DEFAULT_STRING_LEN) tempString
    logical hold
    
    if(.not.present(span))then
      hold=.false.
    else
      hold=span
    end if
    
    if(nWrite==0.or..not.hold)then ! if it is the 1_st time writing this file
      open(fid,file=fname,status='replace')
      ! write header
      write(fid,'(a)'),'$MeshFormat'
      write(fid,'(a)'),'2.2 0 8'
      write(fid,'(a)'),'$EndMeshFormat'
      ! write nodes
      write(fid,'(a)'),'$Nodes'
      write(tempString,*),nNode
      tempString=adjustl(tempString)
      write(fid,'(a)'),trim(tempString)
      do i=1,nNode
        write(tempString,*),i,Node(i)%Pos(:)
        tempString=adjustl(tempString)
        write(fid,'(a)'),trim(tempString)
      end do
      write(fid,'(a)'),'$EndNodes'
      ! write elements
      write(fid,'(a)'),'$Elements'
      write(tempString,*),nBlock+nPoint+nLine+nFacet
      tempString=adjustl(tempString)
      write(fid,'(a)'),trim(tempString)
      ! Note: write blocks first so that the block index would be consistent
      do i=1,nBlock
        write(tempString,*),i,Block(i)%ShapeType,2,0,Block(i)%GeoEnti,Block(i)%NodeInd(:)
        tempString=adjustl(tempString)
        write(fid,'(a)'),trim(tempString)
      end do
      do i=1,nPoint
        write(tempString,*),nBlock+i,POINT_TYPE,2,0,Point(i)%GeoEnti,Point(i)%NodeInd
        tempString=adjustl(tempString)
        write(fid,'(a)'),trim(tempString)
      end do
      do i=1,nLine
        write(tempString,*),nBlock+nPoint+i,LINE_TYPE,2,0,Line(i)%GeoEnti,Line(i)%NodeInd(:)
        tempString=adjustl(tempString)
        write(fid,'(a)'),trim(tempString)
      end do
      do i=1,nFacet
        write(tempString,*),nBlock+nPoint+nLine+i,Facet(i)%ShapeType,2,0,Facet(i)%GeoEnti,&
        &                   Facet(i)%NodeInd(:)
        tempString=adjustl(tempString)
        write(fid,'(a)'),trim(tempString)
      end do
      write(fid,'(a)'),'$EndElements'
    end if
    
    ! write scalers at nodes
    if(allocated(rstNodeScal))then
      do i=1,size(rstNodeScal)
        write(fid,'(a)'),'$NodeData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"scalNode',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),1 ! 1 component scaler
        write(fid,*),nNode
        do j=1,nNode
          write(tempString,*),j,rstNodeScal(i)%ptr(j)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndNodeData'
      end do
    end if
    ! write vectors at nodes
    if(allocated(rstNodeVect))then
      do i=1,size(rstNodeVect)
        write(fid,'(a)'),'$NodeData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"vectNode',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),3 ! 3 components vector
        write(fid,*),nNode
        do j=1,nNode
          write(tempString,*),j,rstNodeVect(i)%ptr(j,:)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndNodeData'
      end do
    end if
    ! write tensors at nodes
    if(allocated(rstNodeTens))then
      do i=1,size(rstNodeTens)
        write(fid,'(a)'),'$NodeData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"tensNode',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),9 ! 9 components tensor
        write(fid,*),nNode
        do j=1,nNode
          write(tempString,*),j,rstNodeTens(i)%ptr(j,1,:),rstNodeTens(i)%ptr(j,2,:),&
          &                     rstNodeTens(i)%ptr(j,3,:)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndNodeData'
      end do
    end if
    
    ! write scalers at the centre of blocks
    if(allocated(rstBlockScal))then
      do i=1,size(rstBlockScal)
        write(fid,'(a)'),'$ElementData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"scalBlock',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),1 ! 1 component scaler
        write(fid,*),nBlock
        do j=1,nBlock
          write(tempString,*),j,rstBlockScal(i)%ptr(j)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndElementData'
      end do
    end if
    ! write vectors at the centre of blocks
    if(allocated(rstBlockVect))then
      do i=1,size(rstBlockVect)
        write(fid,'(a)'),'$ElementData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"vectBlock',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),3 ! 3 components vector
        write(fid,*),nBlock
        do j=1,nBlock
          write(tempString,*),j,rstBlockVect(i)%ptr(j,:)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndElementData'
      end do
    end if
    ! write tensors at the centre of blocks
    if(allocated(rstBlockTens))then
      do i=1,size(rstBlockTens)
        write(fid,'(a)'),'$ElementData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"tensBlock',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),9 ! 9 components tensor
        write(fid,*),nBlock
        do j=1,nBlock
          write(tempString,*),j,rstBlockTens(i)%ptr(j,1,:),rstBlockTens(i)%ptr(j,2,:),&
          &                     rstBlockTens(i)%ptr(j,3,:)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndElementData'
      end do
    end if
    
    ! write scalers at the centre of facets
    if(allocated(rstFacetScal))then
      do i=1,size(rstFacetScal)
        write(fid,'(a)'),'$ElementData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"scalFacet',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),1 ! 1 component scaler
        write(fid,*),nFacet
        do j=1,nFacet
          write(tempString,*),nBlock+nPoint+nLine+j,rstFacetScal(i)%ptr(j)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndElementData'
      end do
    end if
    ! write vectors at the centre of blocks
    if(allocated(rstFacetVect))then
      do i=1,size(rstFacetVect)
        write(fid,'(a)'),'$ElementData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"vectFacet',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),3 ! 3 components vector
        write(fid,*),nFacet
        do j=1,nFacet
          write(tempString,*),nBlock+nPoint+nLine+j,rstFacetVect(i)%ptr(j,:)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndElementData'
      end do
    end if
    ! write tensors at the centre of blocks
    if(allocated(rstFacetTens))then
      do i=1,size(rstFacetTens)
        write(fid,'(a)'),'$ElementData'
        write(fid,'(i1)'),1 ! 1 string tag
        write(fid,'(a,i2,a)'),'"tensFacet',i,'"' ! name of the data-set
        write(fid,'(i1)'),1 ! 1 real tag
        write(fid,*),t ! time
        write(fid,'(i1)'),3 ! 3 integer tags
        write(fid,*),merge(nWrite,1,hold) ! time step index
        write(fid,'(i1)'),9 ! 9 components tensor
        write(fid,*),nFacet
        do j=1,nFacet
          write(tempString,*),nBlock+nPoint+nLine+j,rstFacetTens(i)%ptr(j,1,:),&
          &                                         rstFacetTens(i)%ptr(j,2,:),&
          &                                         rstFacetTens(i)%ptr(j,3,:)
          tempString=adjustl(tempString)
          write(fid,'(a)'),trim(tempString)
        end do
        write(fid,'(a)'),'$EndElementData'
      end do
    end if
    
    if(hold)then
      nWrite=nWrite+1
    end if
    
    if(t>tFinal.or..not.hold)then ! if it is the last time writing results
      close(fid)
    end if
  end subroutine
  
  !--------------------
  ! read material file
  !--------------------
  subroutine readMtl(fname,fid)
    use moduleMtl
    character(*),intent(in)::fname
    integer,intent(in)::fid
    integer readerr
    character(DEFAULT_STRING_LEN) tempString
    
    open(fid,file=fname,status='old')
    readerr=0
    
    read(fid,*,iostat=readerr),n
    allocate(Mtl(n))
    
    do while(readerr==0)
      ! skip the irrelevant lines
      do while(readerr==0)
        read(fid,*,iostat=readerr),tempString
        if(tempString(1:1)=='$')then
          exit
        end if
      end do
      ! check if finished
      if(readerr/=0)then
        exit
      end if
      ! read node data
      if(tempString(1:4)=='$Mtl'.or.tempString(1:4)=='$Mtl')then
        !TODO: read
        cycle
      end if
    end do
  end subroutine
  
end module

