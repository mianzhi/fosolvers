!----------------------------------------------------------------------------- best with 100 columns

!*********
! find dt
!*********
function finddt(l)
  use moduleGrid,only:BoundBox,nEle
  use globalvar,only:rho,c,cond,numIt,dt,tWrite
  integer,intent(in)::l
  double precision finddt
  finddt=dt ! not changed in default
  if(l<0)then
    finddt=minval(rho(1:nEle)*c(1:nEle)/cond(1:nEle))
    finddt=finddt*dot_product(BoundBox(2,:)-BoundBox(1,:),BoundBox(2,:)-BoundBox(1,:))/nEle
    ! note: this is just an initial approximation of time step size, would be adjusted later
  else
    if(l<numIt(1))then
      finddt=dt*1.005d0
    end if
    if(l>numIt(2))then
      finddt=dt*0.99d0
    end if
    if(l>=numIt(3))then
      finddt=dt*0.95d0
    end if
  end if
  finddt=min(finddt,tWrite/2d0)
end function

!*************************************************
! additional shared variables for Crank-Nicholson
!*************************************************
! includes implicit parameters and state variables
module CNaddvar
  double precision,save,allocatable::newTemp(:),newe(:),newcond(:),auxTcorr(:,:)
end module

!********************************
! Crank-Nicholson time advancing
!********************************
subroutine CrankNicolson(l)
! note: return the number of iterations to adjust the time step size
  use moduleGrid
  use globalvar
  use CNaddvar
  use moduleWrite,only:t
  double precision flux,source,error,buff,Dist,conds,newconds,area
  double precision lookupTab ! functions
  allocate(newTemp(-nFacet:nEle))
  allocate(newe(1:nEle))
  allocate(newcond(-nFacet:nEle))
  allocate(auxTcorr(nEle,6))
  ! Note: 6 is the maximum possible number of surfaces (auxiliary points) an element can have.
  
  !----------------------------------------------------------
  ! prepare information for t^n-1 & initialize the iteration
  !----------------------------------------------------------
  ! find grad(T) at centers, save to gradT, and find auxiliary points' temperature corrector
  !$omp parallel do
  do i=1,nEle
    call findEleGradScal(i,Temp(1:nEle),gradT(i,:))
    forall(j=1:Ele(i)%SurfNum)
      auxTcorr(i,j)=dot_product(gradT(i,:),Paux(i,j,:))
    end forall
  end do
  !$omp end parallel do
  ! BC for t^n-1
  call appBCs('old')
  ! copy variables & parameters to t^n
  newTemp(:)=Temp(:)
  newcond(:)=cond(:)
  newe(:)=e(:)
  
  !-----------------------------------------------
  ! iteration with Crank-Nicolson implicit scheme
  !-----------------------------------------------
  do l=1,numIt(3)
    error=0d0
    ! BC for t^n
    call appBCs('new')
    ! solve
    !$omp parallel do &
    !$omp default(shared) &
    !$omp private(flux,source,Dist,area,conds,newconds,buff,i,j,k) &
    !$omp reduction(max:error)
    do i=1,nEle
      flux=0d0
      source=0d0
      do j=1,Ele(i)%SurfNum ! find flux
        Dist=sqrt(dot_product(Ele(i)%PC(:)+Paux(i,j,:)-Ele(i)%SurfPC(j,:),&
        &                     Ele(i)%PC(:)+Paux(i,j,:)-Ele(i)%SurfPC(j,:)))
        area=Ele(i)%SurfArea(j)
        if(Ele(i)%Neib(j)<0)then ! for boundary surface
          Dist=2d0*Dist
          conds=(cond(i)+cond(Ele(i)%Neib(j)))/2d0 ! conductivity at the element surface
          newconds=(newcond(i)+newcond(Ele(i)%Neib(j)))/2d0
          flux=flux+0.5d0*conds*area/Dist*(Temp(i)-Temp(Ele(i)%Neib(j))+auxTcorr(i,j))& ! explicit
          &        +0.5d0*newconds*area/Dist*(newTemp(i)-newTemp(Ele(i)%Neib(j))+auxTcorr(i,j))
                                                                                        ! implicit
        else ! for inner surface
          Dist=(1d0+Ele(Ele(i)%Neib(j))%Vol/Ele(i)%Vol)*Dist
          conds=(Ele(i)%Vol+Ele(Ele(i)%Neib(j))%Vol)&
          &     /(Ele(i)%Vol/cond(i)+Ele(Ele(i)%Neib(j))%Vol/cond(Ele(i)%Neib(j)))
          newconds=(Ele(i)%Vol+Ele(Ele(i)%Neib(j))%Vol)&
          &        /(Ele(i)%Vol/newcond(i)+Ele(Ele(i)%Neib(j))%Vol/newcond(Ele(i)%Neib(j)))
          do k=1,Ele(Ele(i)%Neib(j))%SurfNum !find the correct auxTcorr to apply on the neighbour
            if(Ele(Ele(i)%Neib(j))%Neib(k)==i)then
              exit
            end if
          end do
          flux=flux+0.5d0*conds*area/Dist*(Temp(i)-Temp(Ele(i)%Neib(j))& ! explicit
          &                               +auxTcorr(i,j)-auxTcorr(Ele(i)%Neib(j),k))&
          &        +0.5d0*newconds*area/Dist*(newTemp(i)-newTemp(Ele(i)%Neib(j))& ! implicit
          &                                  +auxTcorr(i,j)-auxTcorr(Ele(i)%Neib(j),k))
        end if
      end do
      
      if(lSource(i).eqv..true.)then ! find source
        source=lookupTab(t+dt/2d0,codData([0,vSource(i)],:),ncodData)
      else
        source=0d0
      end if
      
      ! find new specific energy
      newe(i)=beta*(e(i)*rho(i)*Ele(i)%Vol&
      &            -dt*flux& ! flux; note that area is included
      &            +dt*source*Ele(i)%Vol)& ! source
      &            /rho(i)/Ele(i)%Vol&
      &      +(1-beta)*newe(i)
      
      ! recover new temperature from new specific energy & find new conductivity accordingly
      do j=1,size(iVol)
        if(Ele(i)%GeoEnti==iVol(j))then
          buff=newTemp(i)
          newTemp(i)=lookupTab(newe(i),eTab(j,[2,1],:),neTab(j)) ! using T-e table inversely
          if(abs(newTemp(i)-buff)>error)then
            error=abs(newTemp(i)-buff)
          end if
          newcond(i)=lookupTab(newTemp(i),kTab(j,:,:),nkTab(j))
        end if
      end do
    end do
    !$omp end parallel do
    
    ! check whether to stop
    ! Note: ".not.(abs(newe(:))<huge(0d0))" is an alternative to "isnan(new(e))".
    if(error<=tolerance.or.any(.not.(abs(newe(:))<huge(0d0))).or.any(newe(:)<0d0))then
      exit
    end if
  end do
  
  ! advance only if the iteration converges to a reasonable result
  if(l>1.and.l<=numIt(3).and.(all(abs(newe(:))<huge(0d0))).and.(.not.any(newe(:)<0d0)))then
    e(:)=newe(:)
    Temp(:)=newTemp(:)
    cond(:)=newcond(:)
    t=t+dt
  else
    l=numIt(3)+1 ! decrease the time step size
  end if
  
  !----------
  ! clean up
  !----------
  deallocate(newTemp)
  deallocate(newcond)
  deallocate(auxTcorr)
  deallocate(newe)
end subroutine
