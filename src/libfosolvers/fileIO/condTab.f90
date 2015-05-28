!----------------------------------------------------------------------------- best with 100 columns

!> read conditions from fid into cond
subroutine readCondTab(fid,cond)
  use modCondition
  integer,intent(in)::fid !< file id
  type(condTab),intent(inout)::cond !< result
  
  read(fid,*),n
  read(fid,*),m
  call cond%init(n,m)
  do i=1,n
    read(fid,*)
    read(fid,*),cond%gid(i)
    read(fid,*),cond%t(i)
    read(fid,*),k
    do j=1,k
      read(fid,*),cond%p(j,i)
    end do
  end do
end subroutine
