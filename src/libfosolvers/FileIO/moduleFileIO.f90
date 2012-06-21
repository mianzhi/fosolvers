!----------------------------------------------------------------------------- best with 100 columns

!> file input/output
module moduleFileIO
  private
  
  ! constants
  integer,parameter::RANK_SCAL=1 !< rank of scalar
  integer,parameter::RANK_VECT=2 !< rank of vector
  integer,parameter::RANK_TENS=3 !< rank of tensor
  
  !> read GMSH file
  interface readGMSH
    module procedure::readGMSHGrid
  end interface
  public::readGMSH
  
  !> write GMSH file
  interface writeGMSH
    module procedure::writeGMSHGrid
    module procedure::writeGMSHScal
    module procedure::writeGMSHVect
    module procedure::writeGMSHTens
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
        grid%nPrt=maxval([(maxval(grid%Block(l)%Prt(:)),l=1,grid%nBlock)])
        if(grid%nPrt>0)then
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
        else
          grid%nPrt=1
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
  
  !> write scalar/vector/tensor v into opened file id
  subroutine writeGMSHScalVectTens(id,rank,v,grid,bind,vName,n,t)
    use moduleGrid
    integer,intent(in)::id !< the file id
    integer,intent(in)::rank !< rank of v
    double precision,intent(in)::v(*) !< data to be written
    type(typeGrid),intent(in)::grid !< grid on which v is defined
    integer,intent(in)::bind !< bind with Node/Facet/Block
    character(*),intent(in)::vName !< name of v
    integer,intent(in),optional::n !< time step index
    double precision,intent(in),optional::t !< time
    integer,parameter::DFLT_STR_LEN=400
    character(DFLT_STR_LEN)::tempStr
    integer nEntry,nItem
    
    select case(bind)
    case(BIND_NODE)
      write(id,'(a)'),'$NodeData'
    case(BIND_FACET,BIND_BLOCK)
      write(id,'(a)'),'$ElementData'
    case default
    end select
    write(id,'(i1)'),1 ! 1 string tag
    write(id,'(a,a,a)'),'"',vName,'"'
    write(id,'(i1)'),1 ! 1 real tag
    if(present(t))then
      write(id,*),t
    else
      write(id,*),0d0
    end if
    write(id,'(i1)'),3 ! 3 integer tags
    if(present(n))then
      write(id,*),n
    else
      write(id,*),0
    end if
    select case(rank)
    case(RANK_SCAL)
      nEntry=1
    case(RANK_VECT)
      nEntry=3
    case(RANK_TENS)
      nEntry=9
    case default
      nEntry=0
    end select
    write(id,'(i1)'),nEntry
    select case(bind)
    case(BIND_NODE)
      nItem=grid%nNode
    case(BIND_FACET)
      nItem=grid%nFacet
    case(BIND_BLOCK)
      nItem=grid%nBlock
    case default
      nItem=0
    end select
    write(id,*),nItem
    do i=1,nItem
      select case(bind)
      case(BIND_NODE,BIND_BLOCK)
        write(tempStr,*),i,v((i-1)*nEntry+1:(i-1)*nEntry+nEntry)
      case(BIND_FACET)
        write(tempStr,*),i+grid%nBlock,v((i-1)*nEntry+1:(i-1)*nEntry+nEntry)
      case default
      end select
      write(id,'(a)'),trim(adjustl(tempStr))
    end do
    select case(bind)
    case(BIND_NODE)
      write(id,'(a)'),'$EndNodeData'
    case(BIND_FACET,BIND_BLOCK)
      write(id,'(a)'),'$EndElementData'
    case default
    end select
  end subroutine
  
  !> write scalar v into opened file id
  subroutine writeGMSHScal(id,v,grid,bind,vName,n,t)
    use moduleGrid
    integer,intent(in)::id !< the file id
    double precision,intent(in)::v(:) !< scalar to be written
    type(typeGrid),intent(in)::grid !< grid on which v is defined
    integer,intent(in)::bind !< bind with Node/Facet/Block
    character(*),intent(in)::vName !< name of v
    integer,intent(in),optional::n !< time step index
    double precision,intent(in),optional::t !< time
    
    if(present(n).and.present(t))then
      call writeGMSHScalVectTens(id,RANK_SCAL,v,grid,bind,vName,n,t)
    else
      call writeGMSHScalVectTens(id,RANK_SCAL,v,grid,bind,vName)
    end if
  end subroutine
  
  !> write vector v into opened file id
  subroutine writeGMSHVect(id,v,grid,bind,vName,n,t)
    use moduleGrid
    integer,intent(in)::id !< the file id
    double precision,intent(in)::v(:,:) !< vector to be written
    type(typeGrid),intent(in)::grid !< grid on which v is defined
    integer,intent(in)::bind !< bind with Node/Facet/Block
    character(*),intent(in)::vName !< name of v
    integer,intent(in),optional::n !< time step index
    double precision,intent(in),optional::t !< time
    
    if(present(n).and.present(t))then
      call writeGMSHScalVectTens(id,RANK_VECT,v,grid,bind,vName,n,t)
    else
      call writeGMSHScalVectTens(id,RANK_VECT,v,grid,bind,vName)
    end if
  end subroutine
  
  !> write tensor v into opened file id
  subroutine writeGMSHTens(id,v,grid,bind,vName,n,t)
    use moduleGrid
    integer,intent(in)::id !< the file id
    double precision,intent(in)::v(:,:,:) !< tensor to be written
    type(typeGrid),intent(in)::grid !< grid on which v is defined
    integer,intent(in)::bind !< bind with Node/Facet/Block
    character(*),intent(in)::vName !< name of v
    integer,intent(in),optional::n !< time step index
    double precision,intent(in),optional::t !< time
    
    if(present(n).and.present(t))then
      call writeGMSHScalVectTens(id,RANK_TENS,v,grid,bind,vName,n,t)
    else
      call writeGMSHScalVectTens(id,RANK_TENS,v,grid,bind,vName)
    end if
  end subroutine
  
end module
