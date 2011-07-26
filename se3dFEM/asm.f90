!----------------------------------------------------------------------------- best with 100 columns

!******************************************
! assmble Tet. k which is the i th element
!******************************************
subroutine asmTet(k,i)
  use globalvar
  integer k,i
  integer findSpace
  double precision ATet(12,12),rhsTet(12),P1(3),P2(3),P3(3),P4(3)
  integer TetMask(12)
  ! stiffness matrix for element
  P1(:)=Pos(lTet(k,1),:)
  P2(:)=Pos(lTet(k,2),:)
  P3(:)=Pos(lTet(k,3),:)
  P4(:)=Pos(lTet(k,4),:)
  call stiTet(ATet,P1,P2,P3,P4,lambda,mu)
  ! RHS for element
  rhsTet(:)=0d0
  do j=1,4
    if(lNeum(lTet(k,j)).eqv..true.)then
      rhsTet(j*3-2:j*3)=vNeum(lTet(k,j),:)
      lNeum(lTet(k,j))=.false. ! Neumann BC (and body force) are only applied once
    else
      rhsTet(j*3-2:j*3)=[0d0,0d0,0d0]
    end if
  end do
  
  ! allocate space in frontal matrix
  do j=1,4
    if(N2F(lTet(k,j))==0)then
      N2F(lTet(k,j))=findSpace(F2N,FrontWidth)
      F2N(N2F(lTet(k,j)))=lTet(k,j)
    end if
    TetMask(j*3-2:j*3)=N2F(lTet(k,j))*3+[-2,-1,0]
  end do
  
  ! assembly stiffness matrix
  Front(TetMask,TetMask)=Front(TetMask,TetMask)+ATet
  ! assembly RHS
  Front(TetMask,FrontWidth*3+1)=Front(TetMask,FrontWidth*3+1)+rhsTet
  
  ! eliminate when should
  do j=1,4
    if(NodeEliTime(lTet(k,j))==i)then
      ! eliminate
      kNode(:,:)=Front(3*N2F(lTet(k,j))-2:3*N2F(lTet(k,j)),3*N2F(lTet(k,j))-2:3*N2F(lTet(k,j)))
      call inv(kNode,3)
      call solvePiv(N2F(lTet(k,j)))
      !$omp parallel do schedule(static,1)
      do l=1,FrontWidth
        call eli(l,N2F(lTet(k,j))) ! eliminate line l
      end do
      !$omp end parallel do
      ! save eliminated node
      noe=noe+1
      seqEli(noe)=lTet(k,j)
      buffF2N(noe,:)=F2N(:)
      buffFront(3*lTet(k,j)-2:3*lTet(k,j),:)=Front(3*N2F(lTet(k,j))-2:3*N2F(lTet(k,j)),:)
      ! clean up
      Front(3*N2F(lTet(k,j))-2:3*N2F(lTet(k,j)),:)=0d0
      F2N(N2F(lTet(k,j)))=0
      N2F(lTet(k,j))=0
    end if
  end do
end subroutine

!******************************************
! assmble Hex. k which is the i th element
!******************************************
subroutine asmHex(k,i)
  use globalvar
  integer k,i
  integer findSpace
  double precision AHex(24,24),rhsHex(24),P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3)
  integer HexMask(24)
  ! stiffness matrix for element
  P1(:)=Pos(lHex(k,1),:)
  P2(:)=Pos(lHex(k,2),:)
  P3(:)=Pos(lHex(k,3),:)
  P4(:)=Pos(lHex(k,4),:)
  P5(:)=Pos(lHex(k,5),:)
  P6(:)=Pos(lHex(k,6),:)
  P7(:)=Pos(lHex(k,7),:)
  P8(:)=Pos(lHex(k,8),:)
  call stiHex(AHex,P1,P2,P3,P4,P5,P6,P7,P8,lambda,mu)
  ! RHS for element
  rhsHex(:)=0d0
  do j=1,8
    if(lNeum(lHex(k,j)).eqv..true.)then
      rhsHex(j*3-2:j*3)=vNeum(lHex(k,j),:)
      lNeum(lHex(k,j))=.false. ! Neumann BC (and body force) are only applied once
    else
      rhsHex(j*3-2:j*3)=[0d0,0d0,0d0]
    end if
  end do
  
  ! allocate space in frontal matrix
  do j=1,8
    if(N2F(lHex(k,j))==0)then
      N2F(lHex(k,j))=findSpace(F2N,FrontWidth)
      F2N(N2F(lHex(k,j)))=lHex(k,j)
    end if
    HexMask(j*3-2:j*3)=N2F(lHex(k,j))*3+[-2,-1,0]
  end do
  
  ! assembly stiffness matrix
  Front(HexMask,HexMask)=Front(HexMask,HexMask)+AHex
  ! assembly RHS
  Front(HexMask,FrontWidth*3+1)=Front(HexMask,FrontWidth*3+1)+rhsHex
  
  ! eliminate when should
  do j=1,8
    if(NodeEliTime(lHex(k,j))==i)then
      ! eliminate
      kNode(:,:)=Front(3*N2F(lHex(k,j))-2:3*N2F(lHex(k,j)),3*N2F(lHex(k,j))-2:3*N2F(lHex(k,j)))
      call inv(kNode,3)
      call solvePiv(N2F(lHex(k,j)))
      !$omp parallel do schedule(static,1)
      do l=1,FrontWidth
        call eli(l,N2F(lHex(k,j))) ! eliminate line l
      end do
      !$omp end parallel do
      ! save eliminated node
      noe=noe+1
      seqEli(noe)=lHex(k,j)
      buffF2N(noe,:)=F2N(:)
      buffFront(3*lHex(k,j)-2:3*lHex(k,j),:)=Front(3*N2F(lHex(k,j))-2:3*N2F(lHex(k,j)),:)
      ! clean up
      Front(3*N2F(lHex(k,j))-2:3*N2F(lHex(k,j)),:)=0d0
      F2N(N2F(lHex(k,j)))=0
      N2F(lHex(k,j))=0
    end if
  end do
end subroutine

!****************************
! "exiplicitize" pivot row m
!****************************
subroutine solvePiv(m)
  use globalvar,only:Front,kNode,F2N,FrontWidth,lDiri,vDiri
  integer m
  if(lDiri(F2N(m)).eqv..true.)then
    Front(3*m-2:3*m,:)=0d0
    Front(3*m-2,3*m-2:3*m)=[1d0,0d0,0d0]
    Front(3*m-1,3*m-2:3*m)=[0d0,1d0,0d0]
    Front(3*m-0,3*m-2:3*m)=[0d0,0d0,1d0]
    Front(3*m-2:3*m,3*FrontWidth+1)=vDiri(F2N(m),:)
  else
    do i=1,FrontWidth
      if(F2N(i)/=0)then
        Front(3*m-2:3*m,3*i-2:3*i)=matmul(kNode,Front(3*m-2:3*m,3*i-2:3*i))
      end if
    end do
    Front(3*m-2:3*m,3*FrontWidth+1)=matmul(kNode,Front(3*m-2:3*m,3*FrontWidth+1))
  end if
end subroutine

!********************************
! eliminate n th node from row m
!********************************
subroutine eli(m,n)
  use globalvar,only:Front,F2N,FrontWidth,lDiri,vDiri
  integer m,n
  double precision k(3,3)
  if(F2N(m)/=0.and.m/=n)then
    k=Front(3*m-2:3*m,3*n-2:3*n)
    if(lDiri(F2N(n)).eqv..true.)then ! Dirichlet node
      ! eliminate RHS
      Front(3*m-2:3*m,3*FrontWidth+1)=Front(3*m-2:3*m,3*FrontWidth+1)&
      &                              -matmul(k,vDiri(F2N(n),:))
    else ! not Dirichlet node
      ! eliminate stiffness matrix
      do i=1,FrontWidth
        if(F2N(i)/=0)then
          Front(3*m-2:3*m,3*i-2:3*i)=Front(3*m-2:3*m,3*i-2:3*i)&
          &                         -matmul(k,Front(3*n-2:3*n,3*i-2:3*i))
        end if
      end do
      ! eliminate RHS
      Front(3*m-2:3*m,3*FrontWidth+1)=Front(3*m-2:3*m,3*FrontWidth+1)&
      &                              -matmul(k,Front(3*n-2:3*n,3*FrontWidth+1))
    end if
    Front(3*m-2:3*m,3*n-2:3*n)=0d0
  end if
end subroutine
