!> sorting and reordering
module modSort
  private
  
  !> generic sort
  interface sort
    module procedure::qsortI
  end interface
  
  public::sort

contains
  
  !> quicksort integer array (driver)
  subroutine qsortI(a,perm)
    integer,intent(inout)::a(:) !< the array
    integer,intent(inout),optional::perm(*) !< the permutation array
    
    if(present(perm))then
      perm(1:size(a))=[(i,i=1,size(a))]
      call qsortI_actual(a,perm=perm,first=1,last=size(a))
    else
      call qsortI_actual(a,first=1,last=size(a))
    end if
  end subroutine
  
  !> quicksort integer array (actual)
  recursive subroutine qsortI_actual(a,perm,first,last)
    integer,intent(inout)::a(*) !< the array
    integer,intent(inout),optional::perm(*) !< the permutation array
    integer,intent(in)::first,last !< sort section from first to last
    integer::tmp,pivot
    
    pivot=a((first+last)/2)
    i=first
    j=last
    do while(.true.)
      do while(a(i)<pivot)
        i=i+1
      end do
      do while(pivot<a(j))
        j=j-1
      end do
      if(i>=j) exit
      tmp=a(i)
      a(i)=a(j)
      a(j)=tmp
      if(present(perm))then
        tmp=perm(i)
        perm(i)=perm(j)
        perm(j)=tmp
      end if
      i=i+1
      j=j-1
    end do
    if(present(perm))then
      if(first<i-1) call qsortI_actual(a,perm=perm,first=first,last=i-1)
      if(j+1<last) call qsortI_actual(a,perm=perm,first=j+1,last=last)
    else
      if(first<i-1) call qsortI_actual(a,first=first,last=i-1)
      if(j+1<last) call qsortI_actual(a,first=j+1,last=last)
    end if
  end subroutine
  
end module
