!----------------------------------------------------------------------------- best with 100 columns

!> read CGNS from fname into polyGrid
subroutine readCGNSPolyGrid(fname,grid)
  use modPolyGrid
  use CGNS
  character(*),intent(in)::fname !< file name
  class(polyGrid),intent(inout)::grid !< result
  integer::fid,ier,isize(3)
  character(100)::tmpStr
  integer,allocatable::conn(:),tmpI(:)
  real,allocatable::tmpR(:)
  
  call cg_open_f(trim(adjustl(fname)),CG_MODE_READ,fid,ier)
  if(ier/=CG_OK)then
    write(*,'(a,a)')'[E] readCGNSPolyGrid(): cannot open ',trim(adjustl(fname))
    call cg_error_exit_f()
  end if
  ! assumed 1 base, 1 zone
  call cg_zone_read_f(fid,1,1,tmpStr,isize,ier)
  ! find maximum nNE, save in m
  m=0
  call cg_nsections_f(fid,1,1,n,ier)
  do i=1,n
    call cg_section_read_f(fid,1,1,i,tmpStr,iS,i0,i1,i2,i3,ier)
    if(iS==MIXED)then
      call cg_ElementDataSize_f(fid,1,1,i,lConn,ier)
      allocate(conn(lConn))
      allocate(tmpI(isize(2)))
      call cg_Elements_read_f(fid,1,1,i,conn,tmpI,ier)
      j=1
      do while(j<=n)
        call cg_npe_f(conn(j),m0,ier)
        m=max(m,m0)
        j=j+m+1
      end do
      deallocate(conn)
      deallocate(tmpI)
    else
      call cg_npe_f(iS,m0,ier)
      m=max(m,m0)
    end if
  end do
  ! fill in node positions
  call grid%init(isize(1),isize(2),m)
  do i=1,3
    call cg_coord_info_f(fid,1,1,i,iT,tmpStr,ier)
    if(iT==RealSingle)then
      allocate(tmpR(grid%nN))
      call cg_coord_read_f(fid,1,1,trim(adjustl(tmpStr)),iT,1,grid%nN,tmpR,ier)
      grid%pN(i,:)=dble(tmpR(:))
      deallocate(tmpR)
    else
      call cg_coord_read_f(fid,1,1,trim(adjustl(tmpStr)),iT,1,grid%nN,grid%pN(i,:),ier)
    end if
  end do
  ! fill in element data from all sections
  do i=1,n
    call cg_section_read_f(fid,1,1,i,tmpStr,iS,iStart,iEnd,i2,i3,ier)
    call cg_ElementDataSize_f(fid,1,1,i,lConn,ier)
    allocate(conn(lConn))
    allocate(tmpI(isize(2)))
    call cg_Elements_read_f(fid,1,1,i,conn,tmpI,ier)
    if(iS==MIXED)then
      j=1
      k=iStart
      do while(j<=lConn)
        grid%sE(k)=transShape(conn(j))
        call cg_npe_f(conn(j),m,ier)
        grid%nNE(k)=m
        grid%iNE(1:m,k)=conn(j+1:j+m)
        j=j+m+1
        k=k+1
      end do
    else
      grid%sE(iStart:iEnd)=transShape(iS)
      call cg_npe_f(iS,m,ier)
      grid%nNE(iStart:iEnd)=m
      grid%iNE(1:m,iStart:iEnd)=reshape(conn(:),[m,iEnd-iStart+1])
    end if
    grid%gid(iStart:iEnd)=i ! gid is defined as section number
    deallocate(conn)
    deallocate(tmpI)
  end do
  call cg_close_f(fid,ier)
  
contains
  
  !> translate from the shape type of CGNS
  pure function transShape(s)
    use modPolyMesh,only:TRI,QUAD
    use modPolyGrid,only:TET,HEX
    integer,intent(in)::s !< shape type of CGNS
    integer::transShape !< shape type
    
    select case(s)
    case(TRI_3)
      transShape=TRI
    case(QUAD_4)
      transShape=QUAD
    case(TETRA_4)
      transShape=TET
    case(HEXA_8)
      transShape=HEX
    case default
      transShape=0
    end select
  end function
  
end subroutine
