!----------------------------------------------------------------------------- best with 100 columns

!***********
! apply BCs
!***********
subroutine appBCs(flag)
  use globalvar
  use CNaddvar
  character*3 flag
  double precision P1(3),P2(3),P3(3),P4(3),Pc(3),Ps(3),Ts,Dist,Flux,sigma,time
  double precision lookupTab
  
  do i=1,nTri
    do j=1,nTet
      do k=1,4
        if(neighbour(j,k)==-i)then ! the neighbour of ghost cell -i
          P1(:)=Pos(lTet(j,1),:)
          P2(:)=Pos(lTet(j,2),:)
          P3(:)=Pos(lTet(j,3),:)
          P4(:)=Pos(lTet(j,4),:)
          Pc=(P1+P2+P3+P4)/4d0
          Ps(:)=(Pos(lTri(i,1),:)+Pos(lTri(i,2),:)+Pos(lTri(i,3),:))/3d0
          Dist=sqrt(dot_product(Pc-Ps,Pc-Ps)) ! distance between cell center and boundary surface
          
          ! Neumann type BC
          sigma=5.67051d-8
          if(flag=='old')then
            Ts=Temp(j)+dot_product(gradT(j,:),Ps-Pc)
            time=t
            Temp(-i)=Temp(j) ! NOTE: the default BC is adiabatic
          else ! NOTE: the convection & radiation BCs can also be implicit.
            Ts=newTemp(j)+dot_product(gradT(j,:),Ps-Pc)
            time=t+dt
            newTemp(-i)=newTemp(j)
          end if
          Flux=0d0
          ! superposition radiation flux
          if(lRadSurf(i).eqv..true.)then
            Flux=Flux+sigma*vRadSurf(i,2)*(Ts**4-vRadSurf(i,1)**4)
          end if
          ! superposition convection flux
          if(lConvSurf(i).eqv..true.)then
            Flux=Flux+vConvSurf(i,2)*(Ts-vConvSurf(i,1))
          end if
          ! superposition predefined flux
          if(lFluxSurf(i).eqv..true.)then
            Flux=Flux+lookupTab(time,codData([0,vFluxSurf(i)],:),ncodData)
          end if
          ! apply Neumann type BC
          if(flag=='old')then
            Temp(-i)=Temp(j)-4*Dist*Flux/(cond(-i)+cond(j))
          else
            newTemp(-i)=newTemp(j)-4*Dist*Flux/(newcond(-i)+newcond(j))
          end if
          
          ! apply surface temperature
          if(lTempSurf(i).eqv..true.)then
            if(flag=='old')then
              Temp(-i)=2d0*lookupTab(time,codData([0,vTempSurf(i)],:),ncodData)-Temp(j)
            else
              newTemp(-i)=2d0*lookupTab(time,codData([0,vTempSurf(i)],:),ncodData)-newTemp(j)
            end if
          end if
          
          ! also update the property of the ghost cell
          do l=1,size(iVol)
            if(flag=='old')then
              cond(-i)=lookupTab(Temp(-i),kTab(l,:,:),nkTab)
            else
              newcond(-i)=lookupTab(newTemp(-i),kTab(l,:,:),nkTab)
            end if
          end do
          
        end if
      end do
    end do
  end do
end subroutine

!**************************
! initialize the variables
!**************************
subroutine init()
  use globalvar,only:Tinit,Temp,gradT,e,rho,c,cond,&
  &                  iVol,nTri,nTet,lTet,neighbour,&
  &                  rhoTab,nrhoTab,cTab,ncTab,kTab,nkTab,eTab,neTab
  double precision lookupTab
  
  allocate(Temp(-nTri:nTet))
  allocate(rho(1:nTet))
  allocate(c(1:nTet))
  allocate(cond(-nTri:nTet))
  allocate(e(1:nTet)) ! NOTE: we don't care about conservation within ghost cells
  allocate(gradT(nTet,3))
  
  Temp(:)=Tinit
  ! initial value of variable & parameters
  do i=1,nTet
    do j=1,size(iVol)
      if(lTet(i,5)==iVol(j))then
        rho(i)=lookupTab(Tinit,rhoTab(j,:,:),nrhoTab)
        c(i)=lookupTab(Tinit,cTab(j,:,:),ncTab)
        cond(i)=lookupTab(Tinit,kTab(j,:,:),nkTab)
        e(i)=lookupTab(Tinit,eTab(j,:,:),neTab)
      end if
    end do
  end do
  ! initial parameters of ghost cells
  do i=1,nTet
    do j=1,size(iVol)
      if(lTet(i,5)==iVol(j))then
        do k=1,4
          if(neighbour(i,k)<0)then ! if has a ghost neighbour
            cond(neighbour(i,k))=lookupTab(Tinit,kTab(j,:,:),nkTab)
          end if
        end do
      end if
    end do
  end do
end subroutine

!***********************************************
! build the table of specific energy e(T)[J/kg]
!***********************************************
subroutine eTabGen()
  use globalvar,only:cTab,ncTab,eTab,neTab,iVol
  double precision Ta,Tb,ca,cb
  double precision lookupTab
  
  n=size(iVol) ! number of materials involved
  allocate(neTab(n))
  allocate(eTab(n,2,maxval(ncTab)*4+1))
  neTab(:)=1
  eTab(:,:,:)=0d0
  
  do i=1,n ! for each material
    ! integrate with assuming the slope of conductivity is piecewise linear
    Ta=0d0
    ca=lookupTab(Ta,cTab(i,:,:),ncTab)
    do j=1,ncTab(i)
      Tb=cTab(i,1,j)
      cb=cTab(i,2,j)
      if(Tb==Ta)then
        ca=cb
        cycle
      end if
      eTab(i,1,neTab(i)+1)=3d0/4d0*Ta+1d0/4d0*Tb
      eTab(i,1,neTab(i)+2)=1d0/2d0*Ta+1d0/2d0*Tb
      eTab(i,1,neTab(i)+3)=1d0/4d0*Ta+3d0/4d0*Tb
      eTab(i,1,neTab(i)+4)=Tb
      eTab(i,2,neTab(i)+1)=(Tb-Ta)*(0.5d0/16d0*(cb-ca)+1d0/4d0*ca)+eTab(i,2,neTab(i))
      eTab(i,2,neTab(i)+2)=(Tb-Ta)*(0.5d0/4d0*(cb-ca)+1d0/2d0*ca)+eTab(i,2,neTab(i))
      eTab(i,2,neTab(i)+3)=(Tb-Ta)*(0.5d0*9d0/16d0*(cb-ca)+3d0/4d0*ca)+eTab(i,2,neTab(i))
      eTab(i,2,neTab(i)+4)=(Tb-Ta)*(0.5d0*(cb-ca)+ca)+eTab(i,2,neTab(i))
      neTab(i)=neTab(i)+4
      Ta=Tb
      ca=cb
    end do
  end do
end subroutine

!**********************************
! find the auxiliary points offset
!**********************************
! Note: we are finding position offset instead of position
subroutine auxGen()
  use globalvar,only:lTet,nTet,Pos,Paux,faceTet
  double precision P1(3),P2(3),P3(3),Pc(3),normVect(3)
  
  allocate(Paux(1:nTet,4,3))
  
  do i=1,nTet
    Pc(:)=sum(Pos(lTet(i,1:4),:),1)/4d0
    do j=1,4
      P1(:)=Pos(lTet(i,faceTet(j,1)),:)
      P2(:)=Pos(lTet(i,faceTet(j,2)),:)
      P3(:)=Pos(lTet(i,faceTet(j,3)),:)
      call normal(normVect,P1,P2,P3)
      Paux(i,j,:)=(P1+P2+P3)/3d0-normVect*dot_product(normVect,(P1+P2+P3)/3d0-Pc)-Pc
    end do
  end do
end subroutine

!*******************************************
! find neighbour cells & assign ghost cells
!*******************************************
subroutine findNeighbour()
  use globalvar,only:lTet,nTet,lTri,nTri,neighbour,faceTet
  logical mask(4)
  
  allocate(neighbour(nTet,4))
  neighbour(:,:)=0
  
  do i=1,nTet
    if(neighbour(i,1)==0.and.neighbour(i,2)==0.and.neighbour(i,3)==0.and.neighbour(i,4)==0)then
      do j=1,nTet ! find neighbour cells
        if(i/=j)then
          mask(:)=.false.
          do k=1,4
            if(lTet(i,k)==lTet(j,1).or.lTet(i,k)==lTet(j,2).or.&
            &  lTet(i,k)==lTet(j,3).or.lTet(i,k)==lTet(j,4))then
              mask(k)=.true.
            end if
          end do
          do k=1,4
            if(mask(faceTet(k,1)).and.mask(faceTet(k,2)).and.mask(faceTet(k,3)))then
              neighbour(i,k)=j
            end if
          end do
          if(neighbour(i,1)/=0.and.neighbour(i,2)/=0.and.&
          &  neighbour(i,3)/=0.and.neighbour(i,4)/=0)then
            exit
          end if
        end if
      end do
      if(j==nTet+1)then ! find neighbour ghost cells
        do j=1,nTri
          mask(:)=.false.
          do k=1,4
            if(lTet(i,k)==lTri(j,1).or.lTet(i,k)==lTri(j,2).or.lTet(i,k)==lTri(j,3))then
              mask(k)=.true.
            end if
          end do
          do k=1,4
            if(neighbour(i,k)==0.and.&
            &  mask(faceTet(k,1)).and.mask(faceTet(k,2)).and.mask(faceTet(k,3)))then
              neighbour(i,k)=-j
            end if
          end do
          if(neighbour(i,1)/=0.and.neighbour(i,2)/=0.and.&
          &  neighbour(i,3)/=0.and.neighbour(i,4)/=0)then
            exit
          end if
        end do
      end if
      ! NOTE: neighbour via contact surface with another sub-rigion would not be a ghost cell
    end if
  end do
end subroutine

!************************
! find the gradient of T
!************************
subroutine findGradT()
  use globalvar,only:gradT,Temp,lTet,nTet,Pos,neighbour,Vol,faceTet
  double precision Afact,Vfact,Te,nVect(3),P1(3),P2(3),P3(3)
  integer isBoundary
  double precision areaTri
  
  gradT(:,:)=0d0
  do i=1,nTet
    Afact=1d0
    Vfact=1d0
    isBoundary=2-sum(sign(1,neighbour(i,:)))/2 ! number of boundary surfaces out of 4 surfaces
    select case(isBoundary)
      case(0)
        Afact=1d0
        Vfact=1d0
      case(1)
        Afact=9d0/16d0 ! shrink the cell so that only available information would be used
        Vfact=27d0/64d0
      case(2)
        Afact=1d0/4d0
        Vfact=1d0/8d0
      case(3)
        Afact=1d0/16d0
        Vfact=1d0/64d0
      case(4)
        Afact=1d0
        Vfact=1d0
    end select
    do j=1,4 ! "Gauss Theorem"
      if(neighbour(i,j)>0)then
        Te=(Temp(i)*Vol(neighbour(i,j))+Temp(neighbour(i,j))*Vol(i))/(Vol(i)+Vol(neighbour(i,j)))
        ! find T on the surface if possible
      else
        Te=Temp(i) ! if not, shrink the cell, so that T at center would be Te
      end if
      P1(:)=Pos(lTet(i,faceTet(j,1)),:)
      P2(:)=Pos(lTet(i,faceTet(j,2)),:)
      P3(:)=Pos(lTet(i,faceTet(j,3)),:)
      call normal(nVect,P1,P2,P3)
      gradT(i,:)=gradT(i,:)+Afact*Te*areaTri(P1,P2,P3)*nVect/Vfact/Vol(i)
      ! long live The Gauss Theorem
    end do
  end do
end subroutine

!***************
! look up table
!***************
function lookupTab(pos,tab,n)
  double precision lookupTab,pos
  integer n
  double precision tab(2,n)
  lookupTab=-1d0
  do i=2,n-1
    if(pos<tab(1,i))then
      lookupTab=(pos-tab(1,i-1))/(tab(1,i)-tab(1,i-1))*(tab(2,i)-tab(2,i-1))+tab(2,i-1)
      exit
    end if
  end do
  if(i==n)then
    lookupTab=(pos-tab(1,i-1))/(tab(1,i)-tab(1,i-1))*(tab(2,i)-tab(2,i-1))+tab(2,i-1)
  end if
end function

!*************************
! volume of a tetrahedron
!*************************
function volTet(P1,P2,P3,P4)
  double precision P1(3),P2(3),P3(3),P4(3),volTet
  double precision y1z2,y1z3,y1z4,y2z1,y2z3,y2z4,y3z1,y3z2,y3z4,y4z1,y4z2,y4z3
  y1z2=P1(2)*P2(3)
  y1z3=P1(2)*P3(3)
  y1z4=P1(2)*P4(3)
  y2z1=P2(2)*P1(3)
  y2z3=P2(2)*P3(3)
  y2z4=P2(2)*P4(3)
  y3z1=P3(2)*P1(3)
  y3z2=P3(2)*P2(3)
  y3z4=P3(2)*P4(3)
  y4z1=P4(2)*P1(3)
  y4z2=P4(2)*P2(3)
  y4z3=P4(2)*P3(3)
  volTet=((P2(1)*(y3z4-y4z3)-P3(1)*(y2z4-y4z2)+P4(1)*(y2z3-y3z2))&
  &      -(P1(1)*(y3z4-y4z3)-P3(1)*(y1z4-y4z1)+P4(1)*(y1z3-y3z1))&
  &      +(P1(1)*(y2z4-y4z2)-P2(1)*(y1z4-y4z1)+P4(1)*(y1z2-y2z1))&
  &      -(P1(1)*(y2z3-y3z2)-P2(1)*(y1z3-y3z1)+P3(1)*(y1z2-y2z1)))/6d0
  volTet=abs(volTet)
end function

!****************************
! normal vector of a polygon
!****************************
subroutine normal(rst,P1,P2,P3)
  double precision rst(3),P1(3),P2(3),P3(3)
  double precision a(3),b(3)
  a(:)=P2(:)-P1(:)
  b(:)=P3(:)-P2(:)
  rst(1)=a(2)*b(3)-a(3)*b(2)
  rst(2)=a(3)*b(1)-a(1)*b(3)
  rst(3)=a(1)*b(2)-a(2)*b(1)
  rst(:)=rst(:)/sqrt(dot_product(rst(:),rst(:)))
end subroutine

!****************************
! surface area of a triangle
!****************************
function areaTri(P1,P2,P3)
  double precision areaTri,P1(3),P2(3),P3(3)
  double precision a(3),b(3)
  a=P2-P1
  b=P3-P1
  areaTri=((a(2)*b(3)-a(3)*b(2))**2d0&
  &       +(a(3)*b(1)-a(1)*b(3))**2d0&
  &       +(a(1)*b(2)-a(2)*b(1))**2d0)**0.5d0/2d0
end function

!**************
! progress bar
!**************
subroutine progBar(v,vtot,l,dt)
  integer l
  double precision v,vtot,dt
  real prog
  prog=v/vtot*100d0
  write(*,'(58a,a,f5.1,a,i3,a,g10.3,a,$)'),&
  &    (char(8),i=1,58),'progress:',prog,'%;      it. steps:',l,';      dt:',dt,'sec'
end subroutine
