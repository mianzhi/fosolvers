!----------------------------------------------------------------------------- best with 100 columns

!***********
! apply BCs
!***********
subroutine appBCs(flag)
  use moduleGrid
  use globalvar
  use CNaddvar
  use moduleWrite,only:t
  character(3) flag
  double precision Pc(3),Ps(3),Ts,Dist,Flux,sigma,time
  double precision lookupTab
  
  do i=1,nFacet
    j=Facet(i)%NeibEle(1) ! j is the neighbour of ghost element -i
    if(j>0.and.Facet(i)%NeibEle(2)==0)then
      Pc(:)=Ele(j)%PC(:)
      Ps(:)=Facet(i)%PC(:)
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
        if(Ele(j)%GeoEnti==iVol(l))then
          if(flag=='old')then
            cond(-i)=lookupTab(Temp(-i),kTab(l,:,:),nkTab(l))
          else
            newcond(-i)=lookupTab(newTemp(-i),kTab(l,:,:),nkTab(l))
          end if
        end if
      end do
      
    end if
  end do
end subroutine

!***********************************************
! build the table of specific energy e(T)[J/kg]
!***********************************************
subroutine geneTab()
  use globalvar,only:cTab,ncTab,eTab,neTab,iVol
  double precision Ta,Tb,ca,cb
  double precision lookupTab !function
  
  n=size(iVol) ! number of materials involved
  allocate(neTab(n))
  allocate(eTab(n,2,maxval(ncTab)*4+1))
  neTab(:)=1
  eTab(:,:,:)=0d0
  
  do i=1,n ! for each material
    ! integrate with assuming the slope of conductivity is piecewise linear
    Ta=0d0
    ca=lookupTab(Ta,cTab(i,:,:),ncTab(i))
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

!**********************************
! find the auxiliary points offset
!**********************************
! Note: we are finding position offset instead of position
subroutine genPaux()
  use moduleGrid
  use globalvar,only:Paux
  double precision Ps(3),Pc(3),normVect(3)
  
  !Note: 6 is the maximum possible number of surfaces of an element, 3 is the number of components
  !      of the vector
  allocate(Paux(nEle,6,3))
  Paux(:,:,:)=0d0
  
  do i=1,nEle
    Pc(:)=Ele(i)%PC(:)
    do j=1,Ele(i)%SurfNum
      Ps(:)=Ele(i)%SurfPC(j,:)
      normVect(:)=Ele(i)%SurfNorm(j,:)
      Paux(i,j,:)=Ps(:)-normVect(:)*dot_product(normVect,Ps-Pc)-Pc(:)
    end do
  end do
end subroutine

!**************************
! initialize the variables
!**************************
! Note: although facets between geometic regions are saved in Facet, they would not be refered as
!       neigbour ghost elements, thus they would not be applied BC etc..
subroutine initVar()
  use moduleGrid
  use globalvar,only:Tinit,Temp,gradT,e,rho,c,cond,iVol,&
  &                  rhoTab,nrhoTab,cTab,ncTab,kTab,nkTab,eTab,neTab
  double precision lookupTab ! function
  
  allocate(Temp(-nFacet:nEle))
  allocate(rho(1:nEle))
  allocate(c(1:nEle))
  allocate(cond(-nFacet:nEle))
  allocate(e(1:nEle)) ! NOTE: we don't care about conservation within ghost cells
  allocate(gradT(nEle,3))
  
  Temp(:)=Tinit
  ! initial value of variable & parameters
  do i=1,nEle
    do j=1,size(iVol)
      if(Ele(i)%GeoEnti==iVol(j))then
        rho(i)=lookupTab(Tinit,rhoTab(j,:,:),nrhoTab(j))
        c(i)=lookupTab(Tinit,cTab(j,:,:),ncTab(j))
        cond(i)=lookupTab(Tinit,kTab(j,:,:),nkTab(j))
        e(i)=lookupTab(Tinit,eTab(j,:,:),neTab(j))
      end if
    end do
  end do
  ! initial parameters of ghost elements
  do i=1,nEle
    do j=1,size(iVol)
      if(Ele(i)%GeoEnti==iVol(j))then
        do k=1,Ele(i)%SurfNum
          if(Ele(i)%Neib(k)<0)then ! if has a ghost neighbour
            cond(Ele(i)%Neib(k))=lookupTab(Tinit,kTab(j,:,:),nkTab(j))
          end if
        end do
      end if
    end do
  end do
end subroutine

!**********************************
! show progress of Crank-Nicholson
!**********************************
! Note: this is different from the showProg() in libfosolvers.
subroutine showProgCN(v,vtot,l,dt)
  integer,intent(in)::l
  double precision,intent(in)::v,vtot,dt
  double precision prog
  prog=v/vtot*100d0
  write(*,'(58a,a,f5.1,a,i3,a,g10.3,a,$)'),&
  &    (char(8),i=1,58),'progress:',prog,'%;      it. steps:',l,';      dt:',dt,'sec'
end subroutine
