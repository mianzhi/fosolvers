!----------------------------------------------------------------------------- best with 100 columns

!> grid operations
module moduleGridOperation
  use moduleGrid
  private
  
  ! constants
  integer,parameter::NPLFB=5 !< node, point, line, facet & block, 5 things need to be mapped
  integer,parameter,public::MAP_NODE=1 !< map index of node
  integer,parameter,public::MAP_POINT=2 !< map index of Point
  integer,parameter,public::MAP_LINE=3 !< map index of Line
  integer,parameter,public::MAP_FACET=4 !< map index of Facet
  integer,parameter,public::MAP_BLOCK=5 !< map index of Block
  
  ! public procedures
  public::splitGridPrt
  
contains
  
  !> split the original grid into resulting grids with respect to partition
  subroutine splitGridPrt(org,rst,r2o,isOverlap)
    use moduleBasicDataStruct
    type(typeGrid),intent(in)::org !< the original grid
    type(typeGrid),intent(inout),allocatable::rst(:) !< the resulting grids
    type(typeHtr1DIArr),intent(inout),allocatable::r2o(:,:) !< mapping from rst to org
    logical,intent(in),optional::isOverlap
    logical overlap
    
    if(present(isOverlap))then
      overlap=isOverlap
    else
      overlap=.false.
    end if
    if(allocated(rst)) deallocate(rst)
    if(allocated(r2o)) deallocate(r2o)
    allocate(rst(org%nPrt))
    allocate(r2o(NPLFB,org%nPrt))
    do i=1,org%nPrt
      call rst(i)%init()
    end do
    do i=1,org%nPoint
      do j=1,size(org%Point(i)%Prt)
        if(org%Point(i)%Prt(j)/=0)then
          call pushArr(r2o(MAP_POINT,org%Point(i)%Prt(j))%dat,i)
          rst(org%Point(i)%Prt(j))%nPoint=rst(org%Point(i)%Prt(j))%nPoint+1
        end if
      end do
    end do
    do i=1,org%nLine
      do j=1,size(org%Line(i)%Prt)
        if(org%Line(i)%Prt(j)/=0)then
          call pushArr(r2o(MAP_LINE,org%Line(i)%Prt(j))%dat,i)
          rst(org%Line(i)%Prt(j))%nLine=rst(org%Line(i)%Prt(j))%nLine+1
        end if
      end do
    end do
    do i=1,org%nFacet
      do j=1,size(org%Facet(i)%Prt)
        if(org%Facet(i)%Prt(j)/=0)then
          call pushArr(r2o(MAP_FACET,org%Facet(i)%Prt(j))%dat,i)
          rst(org%Facet(i)%Prt(j))%nFacet=rst(org%Facet(i)%Prt(j))%nFacet+1
        end if
      end do
    end do
    do i=1,org%nBlock
      do j=1,size(org%Block(i)%Prt)
        if(overlap)then
          if(org%Block(i)%Prt(j)/=0)then
            call pushArr(r2o(MAP_BLOCK,abs(org%Block(i)%Prt(j)))%dat,i)
            rst(abs(org%Block(i)%Prt(j)))%nBlock=rst(abs(org%Block(i)%Prt(j)))%nBlock+1
          end if
        else
          if(org%Block(i)%Prt(j)>0)then
            call pushArr(r2o(MAP_BLOCK,org%Block(i)%Prt(j))%dat,i)
            rst(org%Block(i)%Prt(j))%nBlock=rst(org%Block(i)%Prt(j))%nBlock+1
          end if
        end if
      end do
    end do
  end subroutine
  
end module
