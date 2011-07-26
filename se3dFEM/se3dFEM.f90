!----------------------------------------------------------------------------- best with 100 columns

!******************************
! declare the global variables
!******************************
module globalvar
  ! mesh related variables
  double precision,save,allocatable::Pos(:,:)
  integer,save,allocatable::lEle(:,:),lPoint(:,:),lLine(:,:),lTri(:,:),&
  &                         lQuad(:,:),lTet(:,:),lHex(:,:)
  integer,save::nNode,nEle,nPoint,nLine,nTri,nQuad,nTet,nHex
  double precision,save::PosRange(3,2)
  ! condition related variables
  double precision,save,allocatable::vNeum(:,:),vDiri(:,:)
  logical,save,allocatable::lNeum(:),lDiri(:)
  ! material related variables
  double precision,save::E,poisson,dens,yieldstr,lambda,mu
  character*3,save::mtlid
  ! frontal solving related variables
  double precision,save::kNode(3,3)
  double precision,save,allocatable::Front(:,:),buffFront(:,:)
  integer,save,allocatable::NodeEliTime(:),NodeInFront(:),N2F(:),F2N(:),seqEli(:),buffF2N(:,:)
  integer,save::FrontWidth,noe
  ! results
  double precision,save,allocatable::disp(:,:),strain(:,:),stress(:,:)
end module

!**************
! main program
!**************
program se3dFEM
  ! declare and initialize variables & functions
  use globalvar
  
  ! read mesh file
  call readmsh()
  if(nTet+nHex==0)then
    write(*,'(a)'),'ERROR: no solid element is found in "grid.msh"'
    stop
  end if
  ! read simulation conditions
  call readcod()
  ! read material properties
  call readmtl()
  
  !---------
  ! prepare
  !---------
  allocate(NodeEliTime(nNode))
  allocate(NodeInFront(nNode))
  allocate(N2F(nNode))
  N2F(:)=0
  
  nEle=nTet+nHex
  allocate(lEle(nEle,2)) ! nodes per element & (e.g.) Tet. index
  lEle(1:nTet,1)=4
  lEle(1:nTet,2)=[(i,i=1,nTet)]
  lEle(nTet+1:nTet+nHex,1)=8
  lEle(nTet+1:nTet+nHex,2)=[(i,i=1,nHex)]
  m=maxloc(PosRange(:,2)-PosRange(:,1),1) ! find the direction to sort along
  call sortEle(m)
  
  ! create node elimination time-table
  NodeEliTime(:)=0
  do i=1,nEle
    if(lEle(i,1)==4)then ! for Tet.s
      NodeEliTime(lTet(lEle(i,2),1:4))=i
    end if
    if(lEle(i,1)==8)then ! for Hex.s
      NodeEliTime(lHex(lEle(i,2),1:8))=i
    end if
  end do
  if(minval(NodeEliTime)==0)then
    write(*,'(a)'),'ERROR: exist node(s) not belong to any Tet. or Hex. element'
    write(*,'(a)'),'hint: please try using "Frontal" 3d meshing algorithm in GMSH'
    stop
  end if
  
  ! find the required size of the frontal matrix
  NodeInFront(:)=0
  FrontWidth=0
  do i=1,nEle
    if(lEle(i,1)==4)then !for Tet.s
      NodeInFront(lTet(lEle(i,2),1:4))=1
      FrontWidth=max(sum(NodeInFront),FrontWidth)
      do j=1,4
        if(NodeEliTime(lTet(lEle(i,2),j))==i)then
          NodeInFront(lTet(lEle(i,2),j))=0
        end if
      end do
    end if
    if(lEle(i,1)==8)then !for Hex.s
      NodeInFront(lHex(lEle(i,2),1:8))=1
      FrontWidth=max(sum(NodeInFront),FrontWidth)
      do j=1,8
        if(NodeEliTime(lHex(lEle(i,2),j))==i)then
          NodeInFront(lHex(lEle(i,2),j))=0
        end if
      end do
    end if
  end do
  allocate(F2N(FrontWidth))
  allocate(Front(FrontWidth*3,FrontWidth*3+1))
  allocate(seqEli(nNode))
  allocate(buffF2N(nNode,FrontWidth))
  allocate(buffFront(nNode*3,FrontWidth*3+1))
  F2N(:)=0
  Front(:,:)=0
  noe=0

  write(*,'(a,i8,a)'),'bandwidth:',FrontWidth,'*3'
  
  !----------
  ! assembly
  !----------
  write(*,'(a,$)'),'progress:       '
  do i=1,nEle
    if(lEle(i,1)==4)then ! for Tet.s
      k=lEle(i,2)
      call asmTet(k,i)
    end if
    if(lEle(i,1)==8)then ! for Hex.s
      k=lEle(i,2)
      call asmHex(k,i)
    end if
    ! progress bar
    call progBar(dble(i),dble(nEle))
  end do
  
  ! free RAM as much as possible
  deallocate(lDiri)
  deallocate(lNeum)
  deallocate(vDiri)
  deallocate(vNeum)
  deallocate(NodeEliTime)
  deallocate(NodeInFront)
  deallocate(N2F)
  deallocate(F2N)
  deallocate(Front)
  
  !-------------------
  ! back substitution
  !-------------------
  allocate(disp(nNode,3))
  
  do i=nNode,1,-1
    n=seqEli(i) ! resolve node n
    do j=1,FrontWidth
      if(buffF2N(i,j)/=0.and.buffF2N(i,j)/=n)then
        buffFront(3*n-2:3*n,FrontWidth*3+1)=buffFront(3*n-2:3*n,FrontWidth*3+1)&
        &                           -matmul(buffFront(3*n-2:3*n,3*j-2:3*j),disp(buffF2N(i,j),:))
      end if
    end do
    disp(n,:)=buffFront(3*n-2:3*n,3*FrontWidth+1)
  end do
  
  !------------------
  ! generate results
  !------------------
  allocate(strain(nEle,6))
  allocate(stress(nEle,6))
  
  do i=1,nEle
    if(lEle(i,1)==4)then ! for Tet.s
      call strainTet(strain(i,:),lEle(i,2))
    end if
    if(lEle(i,1)==8)then ! for Hex.s
      call strainHex(strain(i,:),lEle(i,2))
    end if
    call findstress(i)
  end do
  call writerst
  
  ! deallocate memory
  deallocate(Pos)
  deallocate(lEle)
  deallocate(lPoint)
  deallocate(lLine)
  deallocate(lTri)
  deallocate(lQuad)
  deallocate(lTet)
  deallocate(lHex)
  deallocate(seqEli)
  deallocate(buffF2N)
  deallocate(buffFront)
  deallocate(disp)
  deallocate(strain)
  deallocate(stress)
  
  write(*,'(/,a)'),'all done!'
end program
