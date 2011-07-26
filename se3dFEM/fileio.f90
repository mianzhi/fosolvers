!----------------------------------------------------------------------------- best with 100 columns

!***************
! write results
!***************
subroutine writerst()
  use globalvar,only:nNode,nEle,lEle,lTet,lHex,Pos,disp,strain,stress
  integer ifile
  character*200 tempstring

  ifile=10
  open(ifile,file='rst.msh',status='replace')
  ! write header
  write(ifile,'(a)'),'$MeshFormat'
  write(ifile,'(a)'),'2.2 0 8'
  write(ifile,'(a)'),'$EndMeshFormat'
  ! write nodes
  write(ifile,'(a)'),'$Nodes'
  write(tempstring,*),nNode
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nNode
    write(tempstring,'(i8,1x,g14.7,1x,g14.7,1x,g14.7)'),i,Pos(i,1),Pos(i,2),Pos(i,3)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndNodes'
  ! write elements
  write(ifile,'(a)'),'$Elements'
  write(tempstring,*),nEle
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nEle
    if(lEle(i,1)==4)then ! write Tet.s
      write(tempstring,'(i8,1x,i1,1x,i1,1x,i8,1x,i8,1x,i8,1x,i8)'),i,4,0,&
      &    lTet(lEle(i,2),1),lTet(lEle(i,2),2),lTet(lEle(i,2),3),lTet(lEle(i,2),4)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end if
    if(lEle(i,1)==8)then ! write Hex.s
      write(tempstring,'(i8,1x,i1,1x,i1,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)'),i,5,0,&
      &    lHex(lEle(i,2),1),lHex(lEle(i,2),2),lHex(lEle(i,2),3),lHex(lEle(i,2),4),&
      &    lHex(lEle(i,2),5),lHex(lEle(i,2),6),lHex(lEle(i,2),7),lHex(lEle(i,2),8)
      tempstring=adjustl(tempstring)
      write(ifile,'(a)'),trim(tempstring)
    end if
  end do
  write(ifile,'(a)'),'$EndElements'
  ! write displacement
  write(ifile,'(a)'),'$NodeData'
  write(ifile,'(i1)'),1 ! 1 string tag
  write(ifile,'(a)'),'"disp"' ! name of the data-set
  write(ifile,'(i1)'),1 ! 1 real tag
  write(ifile,'(f6.4)'),0.0 ! time
  write(ifile,'(i1)'),3 ! 3 integer tags
  write(ifile,'(i1)'),0 ! time step index
  write(ifile,'(i1)'),3 ! 3 components vector
  write(tempstring,*),nNode
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nNode
    write(tempstring,'(i8,1x,g14.7,1x,g14.7,1x,g14.7)'),i,disp(i,1),disp(i,2),disp(i,3)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndNodeData'
  ! write strain tensor
  write(ifile,'(a)'),'$ElementData'
  write(ifile,'(i1)'),1 ! 1 string tag
  write(ifile,'(a)'),'"strain"' ! name of the data-set
  write(ifile,'(i1)'),1 ! 1 real tag
  write(ifile,'(f6.4)'),0.0 ! time
  write(ifile,'(i1)'),3 ! 3 integer tags
  write(ifile,'(i1)'),0 ! time step index
  write(ifile,'(i1)'),9 ! 9 components tensor
  write(tempstring,*),nEle
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nEle
    write(tempstring,&
    &'(i8,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7)'),&
    &i,strain(i,1),strain(i,4)/2d0,strain(i,5)/2d0,&
    &  strain(i,4)/2d0,strain(i,2),strain(i,6)/2d0,&
    &  strain(i,5)/2d0,strain(i,6)/2d0,strain(i,3)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElementData'
  ! write stress tensor
  write(ifile,'(a)'),'$ElementData'
  write(ifile,'(i1)'),1 ! 1 string tag
  write(ifile,'(a)'),'"stress"' ! name of the data-set
  write(ifile,'(i1)'),1 ! 1 real tag
  write(ifile,'(f6.4)'),0.0 ! time
  write(ifile,'(i1)'),3 ! 3 integer tags
  write(ifile,'(i1)'),0 ! time step index
  write(ifile,'(i1)'),9 ! 9 components tensor
  write(tempstring,*),nEle
  tempstring=adjustl(tempstring)
  write(ifile,'(a)'),trim(tempstring)
  do i=1,nEle
    write(tempstring,&
    &'(i8,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7,1x,g14.7)'),&
    &i,stress(i,1),stress(i,4),stress(i,5),&
    &  stress(i,4),stress(i,2),stress(i,6),&
    &  stress(i,5),stress(i,6),stress(i,3)
    tempstring=adjustl(tempstring)
    write(ifile,'(a)'),trim(tempstring)
  end do
  write(ifile,'(a)'),'$EndElementData'
  close(ifile)
end subroutine

!******************************************************************************
! read simulation conditions (and translate into Dirichlet and Neumann Tables)
!******************************************************************************
subroutine readcod()
  use globalvar,only:mtlid,vNeum,vDiri,lNeum,lDiri,&
  &                  lPoint,lLine,lTri,lQuad,lTet,lHex,&
  &                  nPoint,nLine,nTri,nQuad,nTet,nHex,nNode,&
  &                  Pos,dens
  
  integer codfile,readerr
  character*100 temp_string
  double precision totdim
  double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),P7(3),P8(3)
  double precision temp_double_vect(3),temp_double
  
  ! called functions
  double precision areaTri,volTet
  
  codfile=13
  
  open(codfile,file='conditions.cod',status='old')
  readerr=0
  
  allocate(lDiri(nNode))
  allocate(lNeum(nNode))
  allocate(vDiri(nNode,3))
  allocate(vNeum(nNode,3))
  
  lDiri(:)=.false.
  lNeum(:)=.false.
  vDiri(:,:)=0d0
  vNeum(:,:)=0d0
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(codfile,'(a100)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    ! read material code
    if(temp_string(1:4)=='$Mat')then
      read(codfile,*,iostat=readerr),mtlid
    end if
    ! read surface displacment
    if(temp_string(1:9)=='$DispSurf')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! set Dirichlet vector to each Dirichlet nodes
        do i=1,nTri
          if(lTri(i,4)==l)then
            lDiri(lTri(i,1))=.true.
            lDiri(lTri(i,2))=.true.
            lDiri(lTri(i,3))=.true.
            vDiri(lTri(i,1),:)=temp_double_vect
            vDiri(lTri(i,2),:)=temp_double_vect
            vDiri(lTri(i,3),:)=temp_double_vect
          end if
        end do
        do i=1,nQuad
          if(lQuad(i,5)==l)then
            lDiri(lQuad(i,1))=.true.
            lDiri(lQuad(i,2))=.true.
            lDiri(lQuad(i,3))=.true.
            lDiri(lQuad(i,4))=.true.
            vDiri(lQuad(i,1),:)=temp_double_vect
            vDiri(lQuad(i,2),:)=temp_double_vect
            vDiri(lQuad(i,3),:)=temp_double_vect
            vDiri(lQuad(i,4),:)=temp_double_vect
          end if
        end do
      end do
    end if
    ! read line displacement
    if(temp_string(1:9)=='$DispLine')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! set Dirichlet vector to each Dirichlet nodes
        do i=1,nLine
          if(lLine(i,3)==l)then
            lDiri(lLine(i,1))=.true.
            lDiri(lLine(i,2))=.true.
            vDiri(lLine(i,1),:)=temp_double_vect
            vDiri(lLine(i,2),:)=temp_double_vect
          end if
        end do
      end do
    end if
    ! read point displacement
    if(temp_string(1:10)=='$DispPoint')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! set Dirichlet vector to each Dirichlet nodes
        do i=1,nPoint
          if(lPoint(i,2)==l)then
            lDiri(lPoint(i,1))=.true.
            vDiri(lPoint(i,1),:)=temp_double_vect
          end if
        end do
      end do
    end if
    ! read body force
    if(temp_string(1:10)=='$ForceBody')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! set body force vector to each body nodes
        ! note that source behaves the same way as Neumann BC
        do i=1,nTet
          if(lTet(i,5)==l)then
            lNeum(lTet(i,1))=.true.
            lNeum(lTet(i,2))=.true.
            lNeum(lTet(i,3))=.true.
            lNeum(lTet(i,4))=.true.
            P1=Pos(lTet(i,1),:)
            P2=Pos(lTet(i,2),:)
            P3=Pos(lTet(i,3),:)
            P4=Pos(lTet(i,4),:)
            temp_double=volTet(P1,P2,P3,P4)/4d0
            vNeum(lTet(i,1),:)=vNeum(lTet(i,1),:)+temp_double_vect*temp_double*dens
            vNeum(lTet(i,2),:)=vNeum(lTet(i,2),:)+temp_double_vect*temp_double*dens
            vNeum(lTet(i,3),:)=vNeum(lTet(i,3),:)+temp_double_vect*temp_double*dens
            vNeum(lTet(i,3),:)=vNeum(lTet(i,3),:)+temp_double_vect*temp_double*dens
          end if
        end do
        do i=1,nHex
          if(lHex(i,9)==l)then
            lNeum(lHex(i,1))=.true.
            lNeum(lHex(i,2))=.true.
            lNeum(lHex(i,3))=.true.
            lNeum(lHex(i,4))=.true.
            lNeum(lHex(i,5))=.true.
            lNeum(lHex(i,6))=.true.
            lNeum(lHex(i,7))=.true.
            lNeum(lHex(i,8))=.true.
            P1=Pos(lHex(i,1),:)
            P2=Pos(lHex(i,2),:)
            P3=Pos(lHex(i,3),:)
            P4=Pos(lHex(i,4),:)
            P5=Pos(lHex(i,5),:)
            P6=Pos(lHex(i,6),:)
            P7=Pos(lHex(i,7),:)
            P8=Pos(lHex(i,8),:)
            temp_double=volHex(P1,P2,P3,P4,P5,P6,P7,P8)/8d0
            vNeum(lHex(i,1),:)=vNeum(lHex(i,1),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,2),:)=vNeum(lHex(i,2),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,3),:)=vNeum(lHex(i,3),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,4),:)=vNeum(lHex(i,4),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,5),:)=vNeum(lHex(i,5),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,6),:)=vNeum(lHex(i,6),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,7),:)=vNeum(lHex(i,7),:)+temp_double_vect*temp_double*dens
            vNeum(lHex(i,8),:)=vNeum(lHex(i,8),:)+temp_double_vect*temp_double*dens
          end if
        end do
      end do
    end if
    ! read surface force
    if(temp_string(1:10)=='$ForceSurf')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! find total area of Neumann surface
        totdim=0d0
        do i=1,nTri
          if(lTri(i,4)==l)then
            P1=Pos(lTri(i,1),:)
            P2=Pos(lTri(i,2),:)
            P3=Pos(lTri(i,3),:)
            totdim=totdim+areaTri(P1,P2,P3)
          end if
        end do
        do i=1,nQuad
          if(lQuad(i,5)==l)then
            P1=Pos(lQuad(i,1),:)
            P2=Pos(lQuad(i,2),:)
            P3=Pos(lQuad(i,3),:)
            P4=Pos(lQuad(i,4),:)
            totdim=totdim+areaTri(P1,P2,P3)+areaTri(P3,P4,P1)
          end if
        end do
        ! set Neumann vector to each Neumann nodes
        do i=1,nTri
          if(lTri(i,4)==l)then
            lNeum(lTri(i,1))=.true.
            lNeum(lTri(i,2))=.true.
            lNeum(lTri(i,3))=.true.
            P1=Pos(lTri(i,1),:)
            P2=Pos(lTri(i,2),:)
            P3=Pos(lTri(i,3),:)
            temp_double=areaTri(P1,P2,P3)/totdim/3d0
            vNeum(lTri(i,1),:)=vNeum(lTri(i,1),:)+temp_double_vect*temp_double
            vNeum(lTri(i,2),:)=vNeum(lTri(i,2),:)+temp_double_vect*temp_double
            vNeum(lTri(i,3),:)=vNeum(lTri(i,3),:)+temp_double_vect*temp_double
          end if
        end do
        do i=1,nQuad
          if(lQuad(i,5)==l)then
            lNeum(lQuad(i,1))=.true.
            lNeum(lQuad(i,2))=.true.
            lNeum(lQuad(i,3))=.true.
            lNeum(lQuad(i,4))=.true.
            P1=Pos(lQuad(i,1),:)
            P2=Pos(lQuad(i,2),:)
            P3=Pos(lQuad(i,3),:)
            P4=Pos(lQuad(i,4),:)
            temp_double=(areaTri(P1,P2,P3)+areaTri(P3,P4,P1))/totdim/4d0
            vNeum(lQuad(i,1),:)=vNeum(lQuad(i,1),:)+temp_double_vect*temp_double
            vNeum(lQuad(i,2),:)=vNeum(lQuad(i,2),:)+temp_double_vect*temp_double
            vNeum(lQuad(i,3),:)=vNeum(lQuad(i,3),:)+temp_double_vect*temp_double
            vNeum(lQuad(i,4),:)=vNeum(lQuad(i,4),:)+temp_double_vect*temp_double
          end if
        end do
      end do
    end if
    ! read line force
    if(temp_string(1:10)=='$ForceLine')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! find total length of Neumann line
        totdim=0d0
        do i=1,nLine
          if(lLine(i,3)==l)then
            P1=Pos(lLine(i,1),:)
            P2=Pos(lLine(i,2),:)
            totdim=totdim+sqrt(dot_product((P1-P2),(P1-P2)))
          end if
        end do
        ! set Neumann vector to each Neumann nodes
        do i=1,nLine
          if(lLine(i,3)==l)then
            lNeum(lLine(i,1))=.true.
            lNeum(lLine(i,2))=.true.
            P1=Pos(lLine(i,1),:)
            P2=Pos(lLine(i,2),:)
            temp_double=sqrt(dot_product((P1-P2),(P1-P2)))/totdim/2d0
            vNeum(lLine(i,1),:)=vNeum(lLine(i,1),:)+temp_double_vect*temp_double
            vNeum(lLine(i,2),:)=vNeum(lLine(i,2),:)+temp_double_vect*temp_double
          end if
        end do
      end do
    end if
    ! read point force
    if(temp_string(1:11)=='$ForcePoint')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! set Neumann vector to each Neumann nodes
        do i=1,nPoint
          if(lPoint(i,2)==l)then
            lNeum(lPoint(i,1))=.true.
            vNeum(lPoint(i,1),:)=vNeum(lPoint(i,1),:)+temp_double_vect
          end if
        end do
      end do
    end if
    ! read surface pressure
    if(temp_string(1:9)=='$PresSurf')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double
        ! set Neumann vector to each Neumann nodes
        do i=1,nTri
          if(lTri(i,4)==l)then
            lNeum(lTri(i,1))=.true.
            lNeum(lTri(i,2))=.true.
            lNeum(lTri(i,3))=.true.
            P1=Pos(lTri(i,1),:)
            P2=Pos(lTri(i,2),:)
            P3=Pos(lTri(i,3),:)
            call normal(temp_double_vect,P1,P2,P3)
            vNeum(lTri(i,1),:)=vNeum(lTri(i,1),:)&
            &                  +temp_double_vect*temp_double*areaTri(P1,P2,P3)/3d0
            vNeum(lTri(i,2),:)=vNeum(lTri(i,2),:)&
            &                  +temp_double_vect*temp_double*areaTri(P1,P2,P3)/3d0
            vNeum(lTri(i,3),:)=vNeum(lTri(i,3),:)&
            &                  +temp_double_vect*temp_double*areaTri(P1,P2,P3)/3d0
          end if
        end do
        do i=1,nQuad
          if(lQuad(i,5)==l)then
            lNeum(lQuad(i,1))=.true.
            lNeum(lQuad(i,2))=.true.
            lNeum(lQuad(i,3))=.true.
            lNeum(lQuad(i,4))=.true.
            P1=Pos(lQuad(i,1),:)
            P2=Pos(lQuad(i,2),:)
            P3=Pos(lQuad(i,3),:)
            P4=Pos(lQuad(i,4),:)
            call normal(temp_double_vect,P1,P2,P3)
            vNeum(lQuad(i,1),:)=vNeum(lQuad(i,1),:)&
            &                   +temp_double_vect*temp_double&
            &                    *(areaTri(P1,P2,P3)+areaTri(P3,P4,P1))/4d0
            vNeum(lQuad(i,2),:)=vNeum(lQuad(i,2),:)&
            &                   +temp_double_vect*temp_double&
            &                    *(areaTri(P1,P2,P3)+areaTri(P3,P4,P1))/4d0
            vNeum(lQuad(i,3),:)=vNeum(lQuad(i,3),:)&
            &                   +temp_double_vect*temp_double&
            &                    *(areaTri(P1,P2,P3)+areaTri(P3,P4,P1))/4d0
            vNeum(lQuad(i,4),:)=vNeum(lQuad(i,4),:)&
            &                   +temp_double_vect*temp_double&
            &                    *(areaTri(P1,P2,P3)+areaTri(P3,P4,P1))/4d0
          end if
        end do
      end do
    end if
    ! read node force
    if(temp_string(1:10)=='$ForceNode')then
      do while(readerr==0)
        read(codfile,'(a100)',iostat=readerr),temp_string
        if(temp_string(1:4)=='$End')then
          exit
        end if
        read(temp_string,*),l,temp_double_vect(:)
        ! set Neumann vector to each Neumann nodes
        lNeum(l)=.true.
        vNeum(l,:)=vNeum(l,:)+temp_double_vect
      end do
    end if
  end do
  
  close(codfile)
end subroutine

!**************************
! read material properties
!**************************
subroutine readmtl()
  use globalvar,only:E,poisson,dens,yieldstr,lambda,mu,mtlid
  
  integer mtlfile,readerr
  character*100 temp_string
  character*20 mtlname
  
  E=0d0
  
  mtlfile=12
  
  open(mtlfile,file='materials.mtl',status='old')
  readerr=0
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(mtlfile,'(a100)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    ! read properties
    if(temp_string(2:4)==mtlid)then
      mtlname=temp_string(6:25)
      read(mtlfile,*,iostat=readerr),E
      read(mtlfile,*,iostat=readerr),poisson
      read(mtlfile,*,iostat=readerr),dens
      read(mtlfile,*,iostat=readerr),yieldstr
    end if
  end do
  
  if(E==0)then
    write(*,'(a,a3,a)'),'ERROR: material ',mtlid,' are not found.'
    stop
  else
    write(*,'(a,a)'),'material:',mtlname
    write(*,'(a,g14.7,a)'),"  Young's modulus:",E/1d9,'[GPa]'
    write(*,'(a,g14.7,a)'),"  Poisson's ratio:",poisson
    write(*,'(a,g14.7,a)'),"  density:        ",dens,'[kg/m^3]'
    write(*,'(a,g14.7,a,/)'),"  yield stress:   ",yieldstr/1d6,'[MPa]'
  end if
  
  lambda=E*poisson/(1+poisson)/(1-2*poisson)
  mu=E/2/(1+poisson)
  
  close(mtlfile)
end subroutine

!****************
! read mesh file
!****************
subroutine readmsh()
  use globalvar,only:Pos,&
  &                  lPoint,lLine,lTri,lQuad,lTet,lHex,&
  &                  nNode,nEle,nPoint,nLine,nTri,nQuad,nTet,nHex,&
  &                  PosRange
  
  integer gridfile,readerr,np,nt
  character*100 temp_string
  integer temp_int_vect1(20)
  integer,allocatable::temp_int_mat1(:,:)
  
  np=0
  
  gridfile=11
  
  open(gridfile,file='grid.msh',status='old')
  readerr=0
  
  do while(readerr==0)
    ! skip the irrelevant lines
    do while(readerr==0)
      read(gridfile,'(a100)',iostat=readerr),temp_string
      if(temp_string(1:1)=='$')then
        exit
      end if
    end do
    ! check if finished
    if(readerr/=0)then
      exit
    end if
    ! read node data
    if(temp_string(1:4)=='$NOD'.or.temp_string(1:4)=='$Nod')then
      read(gridfile,*,iostat=readerr),nNode
      allocate(Pos(nNode,3))
      read(gridfile,*,iostat=readerr),j,Pos(1,:)
      PosRange(:,1)=Pos(1,:)
      PosRange(:,2)=Pos(1,:)
      do i=2,nNode
        read(gridfile,*,iostat=readerr),j,Pos(i,:)
        PosRange(1,1)=min(PosRange(1,1),Pos(i,1))
        PosRange(1,2)=max(PosRange(1,2),Pos(i,1))
        PosRange(2,1)=min(PosRange(2,1),Pos(i,2))
        PosRange(2,2)=max(PosRange(2,2),Pos(i,2))
        PosRange(3,1)=min(PosRange(3,1),Pos(i,3))
        PosRange(3,2)=max(PosRange(3,2),Pos(i,3))
        if(i/=j)then
          write(*,'(a)'),'WARNING: node data may be not in sequence'
        end if
      end do
      cycle
    end if
    ! read element data
    if(temp_string(1:4)=='$ELM'.or.temp_string(1:4)=='$Ele')then
      if(temp_string(1:4)=='$ELM')then
        write(*,'(a)'),'WARNING: not capable with GMSH version 1 format'
        write(*,'(a)'),'  please try using the latest version of GMSH'
      end if
      read(gridfile,*,iostat=readerr),nEle
      allocate(lPoint(nEle,2))
      allocate(lLine(nEle,3))
      allocate(lTri(nEle,4))
      allocate(lQuad(nEle,5))
      allocate(lTet(nEle,5))
      allocate(lHex(nEle,9))
      !note the additional byte to storage geometric entity index
      nPoint=0
      nLine=0
      nTri=0
      nQuad=0
      nTet=0
      nHex=0
      do i=1,nEle
        read(gridfile,'(a)',iostat=readerr),temp_string
        read(temp_string,*),j,k,nt
        if(i/=j)then
          write(*,'(a)'),'WARNING: element data may be not in sequence'
        end if
        select case(k) !type of element
          case(15) !point
            np=1
          case(1) !line
            np=2
          case(2) !tri
            np=3
          case(3) !quad
            np=4
          case(4) !tet
            np=4
          case(5) !hex
            np=8
        end select
        read(temp_string,*),temp_int_vect1(1:3+nt+np)
        select case(temp_int_vect1(2)) !type of element
          case(15) !point
            nPoint=nPoint+1
            lPoint(nPoint,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lPoint(nPoint,np+1)=temp_int_vect1(3+nt)
          case(1) !line
            nLine=nLine+1
            lLine(nLine,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lLine(nLine,np+1)=temp_int_vect1(3+nt)
          case(2) !tri
            nTri=nTri+1
            lTri(nTri,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lTri(nTri,np+1)=temp_int_vect1(3+nt)
          case(3) !quad
            nQuad=nQuad+1
            lQuad(nQuad,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lQuad(nQuad,np+1)=temp_int_vect1(3+nt)
          case(4) !tet
            nTet=nTet+1
            lTet(nTet,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lTet(nTet,np+1)=temp_int_vect1(3+nt)
          case(5) !hex
            nHex=nHex+1
            lHex(nHex,1:np)=temp_int_vect1(3+nt+1:3+nt+np)
            lHex(nHex,np+1)=temp_int_vect1(3+nt)
        end select
      end do
      !trim lPoint
      if(nPoint>0)then
        allocate(temp_int_mat1(nPoint,2))
        temp_int_mat1=lPoint(1:nPoint,:)
        deallocate(lPoint)
        allocate(lPoint(nPoint,2))
        lPoint=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lPoint)
        allocate(lPoint(1,1))
      end if
      !trim lLine
      if(nLine>0)then
        allocate(temp_int_mat1(nLine,3))
        temp_int_mat1=lLine(1:nLine,:)
        deallocate(lLine)
        allocate(lLine(nLine,3))
        lLine=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lLine)
        allocate(lLine(1,1))
      end if
      !trim lTri
      if(nTri>0)then
        allocate(temp_int_mat1(nTri,4))
        temp_int_mat1=lTri(1:nTri,:)
        deallocate(lTri)
        allocate(lTri(nTri,4))
        lTri=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lTri)
        allocate(lTri(1,1))
      end if
      !trim lQuad
      if(nQuad>0)then
        allocate(temp_int_mat1(nQuad,5))
        temp_int_mat1=lQuad(1:nQuad,:)
        deallocate(lQuad)
        allocate(lQuad(nQuad,5))
        lQuad=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lQuad)
        allocate(lQuad(1,1))
      end if
      !trim lTet
      if(nTet>0)then
        allocate(temp_int_mat1(nTet,5))
        temp_int_mat1=lTet(1:nTet,:)
        deallocate(lTet)
        allocate(lTet(nTet,5))
        lTet=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lTet)
        allocate(lTet(1,1))
      end if
      !trim lHex
      if(nHex>0)then
        allocate(temp_int_mat1(nHex,9))
        temp_int_mat1=lHex(1:nHex,:)
        deallocate(lHex)
        allocate(lHex(nHex,9))
        lHex=temp_int_mat1
        deallocate(temp_int_mat1)
      else
        deallocate(lHex)
        allocate(lHex(1,1))
      end if
      cycle
    end if
  end do
  
  write(*,'(a)'),'grid summary:'
  write(*,'(a,i8)'),'  number of nodes: ',nNode
  write(*,'(a,i8)'),'  number of points:',nPoint
  write(*,'(a,i8)'),'  number of lines: ',nLine
  write(*,'(a,i8)'),'  number of tri.:  ',nTri
  write(*,'(a,i8)'),'  number of quad.: ',nQuad
  write(*,'(a,i8)'),'  number of tet.:  ',nTet
  write(*,'(a,i8,/)'),'  number of hex.:  ',nHex
  
  close(gridfile)
end subroutine
