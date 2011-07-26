!----------------------------------------------------------------------------- best with 100 columns

!*********
! find dt
!*********
function dtGen(l)
  use globalvar,only:rho,c,cond,PosRange,nTet,numIt,dt,tWrite
  double precision dtGen
  dtGen=dt ! not changed in default
  if(l<0)then
    dtGen=minval(rho(1:nTet)*c(1:nTet)/cond(1:nTet))
    dtGen=dtGen*dot_product(PosRange(:,2)-PosRange(:,1),PosRange(:,2)-PosRange(:,1))/nTet
    ! note: this is just an initial approximation of time step size, would be adjusted later
  else
    if(l<numIt(1))then
      dtGen=dt*1.005d0
    end if
    if(l>numIt(2))then
      dtGen=dt*0.99d0
    end if
    if(l>=numIt(3))then
      dtGen=dt*0.95d0
    end if
  end if
  dtGen=min(dtGen,tWrite/2d0)
end function

!********************************
! Crank-Nicholson time advancing
!********************************
!=============================
! additional shared variables
!=============================
! includes implicit parameters
module CNaddvar
  double precision,save,allocatable::newTemp(:),newe(:),newcond(:),auxTcorr(:,:)
end module
!=======================================
! main program of Crank-Nicholson Sheme
!=======================================
subroutine CrankNicolson(l)
! note: return the number of iterations to adjust the time step size
  use globalvar
  use CNaddvar
  double precision flux,source,error,buff,Pc(3),Ps(3),Dist,conds,newconds,area
  double precision lookupTab,areaTri
  allocate(newTemp(-nTri:nTet))
  allocate(newe(1:nTet))
  allocate(newcond(-nTri:nTet))
  allocate(auxTcorr(nTet,4))
  
  !----------------------------------------------------------
  ! prepare information for t^n-1 & initialize the iteration
  !----------------------------------------------------------
  ! find grad(T) at centers, save to gradT(nTet,3)
  call findGradT()
  ! BC for t^n-1
  call appBCs('old')
  ! find auxiliary points' temperature corrector
  do i=1,nTet
    do j=1,4 ! 4 auxiliary points per cell
      auxTcorr(i,j)=dot_product(gradT(i,:),Paux(i,j,:))
    end do
  end do
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
    do i=1,nTet
      flux=0d0
      source=0d0
      Pc(:)=sum(Pos(lTet(i,1:4),:),1)/4d0 ! position of the cell center
      
      do j=1,4 ! find flux
        Ps(:)=sum(Pos(lTet(i,faceTet(j,:)),:),1)/3d0 ! position of the surface center
        Dist=sqrt(dot_product(Pc-Ps,Pc-Ps))
        area=areaTri(Pos(lTet(i,faceTet(j,1)),:),Pos(lTet(i,faceTet(j,2)),:),&
        &            Pos(lTet(i,faceTet(j,3)),:))
        if(neighbour(i,j)<0)then ! for boundary surface
          Dist=2d0*Dist
          conds=(cond(i)+cond(neighbour(i,j)))/2d0 ! conductivity at the cell surface
          newconds=(newcond(i)+newcond(neighbour(i,j)))/2d0
          flux=flux+0.5d0*conds*area/Dist*(Temp(i)-Temp(neighbour(i,j))+auxTcorr(i,j))& ! explicit
          &        +0.5d0*newconds*area/Dist*(newTemp(i)-newTemp(neighbour(i,j))+auxTcorr(i,j))
                                                                                        ! implicit
        else ! for inner surface
          Dist=(1d0+Vol(neighbour(i,j))/Vol(i))*Dist
          conds=(Vol(i)+Vol(neighbour(i,j)))&
          &     /(Vol(i)/cond(i)+Vol(neighbour(i,j))/cond(neighbour(i,j)))
          newconds=(Vol(i)+Vol(neighbour(i,j)))&
          &        /(Vol(i)/newcond(i)+Vol(neighbour(i,j))/newcond(neighbour(i,j)))
          do k=1,4 !find the correct auxTcorr to apply on the neighbour cell
            if(neighbour(neighbour(i,j),k)==i)then
              exit
            end if
          end do
          flux=flux+0.5d0*conds*area/Dist*(Temp(i)-Temp(neighbour(i,j))& ! explicit
          &                               +auxTcorr(i,j)-auxTcorr(neighbour(i,j),k))&
          &        +0.5d0*newconds*area/Dist*(newTemp(i)-newTemp(neighbour(i,j))& ! implicit
          &                                  +auxTcorr(i,j)-auxTcorr(neighbour(i,j),k))
        end if
      end do
      
      if(lSource(i).eqv..true.)then ! find source
        source=lookupTab(t+dt/2d0,codData([0,vSource(i)],:),ncodData)
      else
        source=0d0
      end if
      
      ! find new specific energy
      newe(i)=beta*(e(i)*rho(i)*Vol(i)&
      &            -dt*flux& ! flux; note that area is included
      &            +dt*source*Vol(i))& ! source
      &            /rho(i)/Vol(i)&
      &      +(1-beta)*newe(i)
      
      ! recover new temperature from new specific energy & find new conductivity accordingly
      do j=1,size(iVol)
        if(lTet(i,5)==iVol(j))then
          buff=newTemp(i)
          newTemp(i)=lookupTab(newe(i),eTab(j,[2,1],:),neTab(j)) ! using T-e table inversely
          if(abs(newTemp(i)-buff)>error)then
            error=abs(newTemp(i)-buff)
          end if
          newcond(i)=lookupTab(newTemp(i),kTab(j,:,:),nkTab(j))
        end if
      end do
    end do
    
    ! check whether to stop
    if(error<=tolerance.or.any(isnan(newe(:))).or.any(newe(:)<0d0))then
      exit
    end if
  end do
  
  ! advance only if the iteration converges to a reasonable result
  if(l>1.and.l<=numIt(3).and.(.not.any(isnan(newe(:)))).and.(.not.any(newe(:)<0d0)))then
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
