!----------------------------------------------------------------------------- best with 100 columns

!***********************
! file input and output
!***********************
module moduleFileIO
  private
  
  ! some constants
  integer,parameter,public::DEFAULT_STRING_LEN=200
  
  ! procedures
  public readmsh
  
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
end module

