!----------------------------------------------------------------------------- best with 100 columns

!> file input/output
module moduleFileIO
  private
  
  !> read GMSH file
  interface readGMSH
    module procedure::readGMSHGrid
  end interface
  public::readGMSH
  
contains
  
  !> read grid data from opened file id into grid
  subroutine readGMSHGrid(id,grid)
    use moduleBasicDataStruct
    use moduleGrid
    integer,intent(in)::id !< the file id
    class(typeGrid),intent(inout)::grid !< resulting grid
    integer ierr
    integer,parameter::DFLT_STR_LEN=400
    character(DFLT_STR_LEN)::tempStr
    integer nEle,shapetype,nt,np,tempIVect(DFLT_STR_LEN/2),jP,jL,jF,jB
    type(typeEle),allocatable::tempEle(:)
    
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
        ! process physical domain and partition TODO
        cycle
      end if
    end do
  end subroutine
  
end module
