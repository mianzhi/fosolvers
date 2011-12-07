!----------------------------------------------------------------------------- best with 100 columns

!***********************
! file input and output
!***********************
module moduleFileIO
  use moduleMiscDataStruct
  private
  
  ! some constants
  integer,parameter,public::DEFAULT_STRING_LEN=200
  
  ! output control related variables
  integer,public,save::nWrite
  ! data to be output
  type(typePtrScalArray),allocatable,public,save::rstNodeScal(:),rstFacetScal(:),rstBlockScal(:)
  type(typePtrVectArray),allocatable,public,save::rstNodeVect(:),rstFacetVect(:),rstBlockVect(:)
  type(typePtrTensArray),allocatable,public,save::rstNodeTens(:),rstFacetTens(:),rstBlockTens(:)
  
  ! procedures
  public readmsh
  public initWriteEnv
  public writerst
  
contains
  
  !----------------
  ! read mesh file
  !----------------
  subroutine readmsh(fname,fid,verbose)
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
  
  !------------------------------------
  ! initialize data output environment
  !------------------------------------
  ! a: number of scaler data sets at nodes
  ! b: number of vector data sets at nodes
  ! c: number of tensor data sets at nodes
  ! d: number of scaler data sets at the centre of facets
  ! e: number of vector data sets at the centre of facets
  ! f: number of tensor data sets at the centre of facets
  ! g: number of scaler data sets at the centre of blocks
  ! h: number of vector data sets at the centre of blocks
  ! i: number of tensor data sets at the centre of blocks
  subroutine initWriteEnv(a,b,c,d,e,f,g,h,i)
    integer,intent(in)::a,b,c,d,e,f,g,h,i
    
    ! allocate output data space
    allocate(rstNodeScal(a))
    allocate(rstNodeVect(b))
    allocate(rstNodeTens(c))
    allocate(rstFacetScal(d))
    allocate(rstFacetVect(e))
    allocate(rstFacetTens(f))
    allocate(rstBlockScal(g))
    allocate(rstBlockVect(h))
    allocate(rstBlockTens(i))
    ! set output control variables
    nWrite=0
  end subroutine
  
  !---------------
  ! write results
  !---------------
  ! hold=.true.: hold the file when exiting subroutine; good for writing a time span
  ! hold=.false.: close the file when exiting subroutine; good for writing a snapshot
  subroutine writerst(fname,fid,span)
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
    ! write vectors at nodes
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
    ! write tensors at nodes
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
    
    ! write scalers at the centre of blocks
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
    ! write vectors at the centre of blocks
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
    ! write tensors at the centre of blocks
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
    
    ! write scalers at the centre of facets
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
    ! write vectors at the centre of blocks
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
    ! write tensors at the centre of blocks
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
    
    if(hold)then
      nWrite=nWrite+1
    end if
    
    if(t>tFinal.or..not.hold)then ! if it is the last time writing results
      close(fid)
    end if
  end subroutine
  
end module

