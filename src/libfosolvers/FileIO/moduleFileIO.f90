!----------------------------------------------------------------------------- best with 100 columns

!> file input/output
module moduleFileIO
  private
  
  !> read GMSH file
  interface readGMSH
    module procedure::readGMSHGrid
  end interface
  public::readGMSH
  
  !> write GMSH file
  interface writeGMSH
    module procedure::writeGMSHGrid
  end interface
  public::writeGMSH
  
contains
  
  !> read grid data from opened file id into grid
  subroutine readGMSHGrid(id,grid)
    use moduleBasicDataStruct
    use moduleSimpleSetLogic
    use moduleGrid
    integer,intent(in)::id !< the file id
    class(typeGrid),intent(inout)::grid !< resulting grid
    integer ierr
    integer,parameter::DFLT_STR_LEN=400
    character(DFLT_STR_LEN)::tempStr
    integer nEle,shapetype,nt,np,tempIVect(DFLT_STR_LEN/2),jP,jL,jF,jB
    type(typeEle),allocatable::tempEle(:)
    type(typeHtr1DIArr),allocatable::NodePrt(:)
    
    call grid%init()
    ierr=0
    rewind(id,iostat=ierr)
    
    do while(ierr==0)
      ! skip the irrelevant lines
      do while(ierr==0)
        read(id,*,iostat=ierr),tempStr
        if(tempStr(1:1)=='$')then
          exit
        end if
      end do
      ! check if finished
      if(ierr/=0)then
        exit
      end if
      ! read node position
      if(tempStr(1:6)=='$Nodes')then
        read(id,*,iostat=ierr),grid%nNode
        allocate(grid%NodePos(DIMS,grid%nNode))
        do i=1,grid%nNode
          read(id,*,iostat=ierr),j,grid%NodePos(:,i)
        end do
        cycle
      end if
      ! read elements
      if(tempStr(1:9)=='$Elements')then
        ! construct tempEle first
        read(id,*,iostat=ierr),nEle
        allocate(tempEle(nEle))
        do i=1,nEle
          read(id,'(a)',iostat=ierr),tempStr
          read(tempStr,*),j,shapetype,nt
          tempEle(i)%Shp=shapetype ! save Shp
          select case(shapetype)
          case(POINT_TYPE)
            np=POINT_NODE_NUM
            grid%nPoint=grid%nPoint+1
          case(LINE_TYPE)
            np=LINE_NODE_NUM
            grid%nLine=grid%nLine+1
          case(TRI_TYPE)
            np=TRI_NODE_NUM
            grid%nFacet=grid%nFacet+1
          case(QUAD_TYPE)
            np=QUAD_NODE_NUM
            grid%nFacet=grid%nFacet+1
          case(TET_TYPE)
            np=TET_NODE_NUM
            grid%nBlock=grid%nBlock+1
          case(HEX_TYPE)
            np=HEX_NODE_NUM
            grid%nBlock=grid%nBlock+1
          case default
            np=0
          end select
          tempEle(i)%nNode=np ! save nNode
          read(tempStr,*),tempIVect(1:3+nt+np)
          call pushArr(tempEle(i)%iNode,tempIVect(3+nt+1:3+nt+np)) ! save iNode
          if(nt>=2)then
            tempEle(i)%Ent=tempIVect(5) ! save Ent
            call pushArr(tempEle(i)%Dmn,tempIVect(4)) ! save Dmn
            if(nt>=4)then
              call pushArr(tempEle(i)%Prt,tempIVect(7:6+tempIVect(6))) ! save Prt
            else
              call pushArr(tempEle(i)%Prt,0)
            end if
          else
            tempEle(i)%Ent=0
            call pushArr(tempEle(i)%Dmn,0)
          end if
        end do
        ! copy data from tempEle to grid
        allocate(grid%Point(grid%nPoint))
        allocate(grid%Line(grid%nLine))
        allocate(grid%Facet(grid%nFacet))
        allocate(grid%Block(grid%nBlock))
        jP=0
        jL=0
        jF=0
        jB=0
        do i=1,nEle
          select case(tempEle(i)%Shp)
          case(POINT_TYPE)
            jP=jP+1
            grid%Point(jP)=tempEle(i)
          case(LINE_TYPE)
            jL=jL+1
            grid%Line(jL)=tempEle(i)
          case(TRI_TYPE,QUAD_TYPE)
            jF=jF+1
            grid%Facet(jF)=tempEle(i)
          case(TET_TYPE,HEX_TYPE)
            jB=jB+1
            grid%Block(jB)=tempEle(i)
          end select
        end do
        deallocate(tempEle)
        ! process partition
        if(maxval([(maxval(grid%Block(l)%Prt(:)),l=1,grid%nBlock)])>0)then
          allocate(NodePrt(grid%nNode))
          call grid%updateBlockNeib()
          do i=1,grid%nBlock
            do j=1,size(grid%Block(i)%Prt)
              if(grid%Block(i)%Prt(j)>0)then
                do k=1,grid%Block(i)%nNode
                  if(allocated(NodePrt(grid%Block(i)%iNode(k))%dat))then
                    call applUnion(NodePrt(grid%Block(i)%iNode(k))%dat,[grid%Block(i)%Prt(j)])
                  else
                    call pushArr(NodePrt(grid%Block(i)%iNode(k))%dat,grid%Block(i)%Prt(j))
                  end if
                end do
                do k=1,getBlockSurfNum(grid%Block(i))
                  if(grid%BlockNeibFacet(i)%dat(k)/=0)then
                    if(grid%Facet(grid%BlockNeibFacet(i)%dat(k))%Prt(1)==0)then
                      grid%Facet(grid%BlockNeibFacet(i)%dat(k))%Prt(1)=grid%Block(i)%Prt(j)
                    else
                      call applUnion(grid%Facet(grid%BlockNeibFacet(i)%dat(k))%Prt,&
                      &             [grid%Block(i)%Prt(j)])
                    end if
                  end if
                end do
              end if
            end do
          end do
          do i=1,grid%nPoint
            call reallocArr(grid%Point(i)%Prt,size(NodePrt(grid%Point(i)%iNode(1))%dat))
            grid%Point(i)%Prt(:)=NodePrt(grid%Point(i)%iNode(1))%dat(:)
          end do
          do i=1,grid%nLine
            call reallocArr(grid%Line(i)%Prt,size(NodePrt(grid%Line(i)%iNode(1))%dat))
            grid%Line(i)%Prt(:)=NodePrt(grid%Line(i)%iNode(1))%dat(:)
            do j=2,grid%Line(i)%nNode
              call applIntersection(grid%Line(i)%Prt,NodePrt(grid%Line(i)%iNode(j))%dat)
            end do
          end do
          deallocate(NodePrt)
        end if
        ! process physical domain TODO
        cycle
      end if
    end do
  end subroutine
  
  !> write grid data into opened file id
  subroutine writeGMSHGrid(id,grid)
    use moduleBasicDataStruct
    use moduleGrid
    integer,intent(in)::id !< the file id
    class(typeGrid),intent(in)::grid !< grid to be written
    integer,parameter::DFLT_STR_LEN=400
    character(DFLT_STR_LEN)::tempStr
    
    write(id,'(a)'),'$MeshFormat'
    write(id,'(a)'),'2.2 0 8'
    write(id,'(a)'),'$EndMeshFormat'
    write(id,'(a)'),'$Nodes'
    write(tempStr,*),grid%nNode
    write(id,'(a)'),trim(adjustl(tempStr))
    do i=1,grid%nNode
      write(tempStr,*),i,grid%NodePos(:,i)
      write(id,'(a)'),trim(adjustl(tempStr))
    end do
    write(id,'(a)'),'$EndNodes'
    write(id,'(a)'),'$Elements'
    write(tempStr,*),grid%nBlock+grid%nFacet+grid%nLine+grid%nPoint
    write(id,'(a)'),trim(adjustl(tempStr))
    do i=1,grid%nBlock
      write(tempStr,*),i,grid%Block(i)%Shp,2,0,grid%Block(i)%Ent,grid%Block(i)%iNode(:)
      write(id,'(a)'),trim(adjustl(tempStr))
    end do
    do i=1,grid%nFacet
      write(tempStr,*),grid%nBlock+i,grid%Facet(i)%Shp,2,0,grid%Facet(i)%Ent,grid%Facet(i)%iNode(:)
      write(id,'(a)'),trim(adjustl(tempStr))
    end do
    do i=1,grid%nLine
      write(tempStr,*),grid%nBlock+grid%nFacet+i,grid%Line(i)%Shp,2,0,grid%Line(i)%Ent,&
      &                grid%Line(i)%iNode(:)
      write(id,'(a)'),trim(adjustl(tempStr))
    end do
    do i=1,grid%nPoint
      write(tempStr,*),grid%nBlock+grid%nFacet+grid%nLine+i,grid%Point(i)%Shp,2,0,&
      &                grid%Point(i)%Ent,grid%Point(i)%iNode(:)
      write(id,'(a)'),trim(adjustl(tempStr))
    end do
    write(id,'(a)'),'$EndElements'
  end subroutine
  
end module
