!----------------------------------------------------------------------------- best with 100 columns

!> interface to sparse linear solvers and ILU preconditioners (ITSOL)
module modSparse
  use iso_c_binding
  private
  
  !> generic CSR linear equation object
  type,public::linEq
    integer,private::nEq !< number of equations
    integer,private::maxNNz !< maximum number of non-zeros
    integer,allocatable,private::iA(:) !< CSR index array (size nEq+1)
    integer,allocatable,private::jA(:) !< CSR index array (size nNz)
    double precision,allocatable,private::A(:) !< CSR value array
  contains
    procedure::initLinEq
    procedure,public::clear=>clearLinEq
    procedure,public::setCSR=>setCSRLinEq
    procedure,public::setCOO=>setCOOLinEq
    final::purgeLinEq
  end type
  
  !> CSR linear equation object (incompletely) solved by ILUT of SPARSKIT2/ITSOL
  type,extends(linEq),public::ILU
    integer,allocatable,private::ju(:),jlu(:),jw(:) !< auxiliary storage for ILUT
    double precision,allocatable,private::alu(:),w(:) !< auxiliary storage for ILUT
    integer,private::lFill !< maximum number of non-diagonal L or U elements in a row
    integer,private::iwk !< length of alu, jlu
  contains
    procedure,public::init=>initILU
    procedure,public::fact=>factILU
    procedure,public::solve=>solveILU
  end type
  
  !> interface to SPARSKIT2
  interface
    
    !> convert COO format to CSR format
    subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
      integer,intent(in)::nrow,nnz
      double precision,intent(in)::a(*)
      integer,intent(in)::ir(*),jc(*)
      double precision,intent(inout)::ao(*)
      integer,intent(inout)::jao(*),iao(*)
    end subroutine
    
    !> ILUT factorization
    subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
      integer,intent(in)::n
      double precision,intent(in)::a(*)
      integer,intent(in)::ja(*),ia(n+1)
      integer,intent(in)::lfil
      double precision,intent(in)::droptol
      double precision,intent(inout)::alu(*)
      integer,intent(inout)::jlu(*),ju(n)
      integer,intent(in)::iwk
      double precision,intent(inout)::w(n+1)
      integer,intent(inout)::jw(2*n)
      integer,intent(out)::ierr
    end subroutine
    
    !> LU forward and backward solve
    subroutine lusol(n,y,x,alu,jlu,ju)
       integer,intent(in)::n
       double precision,intent(in)::y(n)
       double precision,intent(inout)::x(n)
       double precision,intent(in)::alu(*)
       integer,intent(in)::jlu(*),ju(*)
    end subroutine
    
  end interface
  
contains
  
  !> generic constructor of linEq
  subroutine initLinEq(this,nEq,maxNNz)
    class(linEq),intent(inout)::this !< this linEq
    integer,intent(in)::nEq !< number of equations
    integer,intent(in)::maxNNZ !< maximum number of non-zeros
    
    this%nEq=nEq
    this%maxNNz=maxNNz
    allocate(this%iA(nEq+1))
    allocate(this%jA(maxNNz))
    allocate(this%A(maxNNz))
  end subroutine
  
  !> clear this generic linEq
  subroutine clearLinEq(this)
    class(linEq),intent(inout)::this !< this linEq
    
    deallocate(this%iA,this%jA,this%A)
  end subroutine
  
  !> set matrix with CSR format for this generic linEq
  subroutine setCSRLinEq(this,iA,jA,A)
    class(linEq),intent(inout)::this !< this linEq
    integer,intent(in)::iA(*) !< CSR index array (size nEq+1)
    integer,intent(in)::jA(*) !< CSR index array (size nNz)
    double precision,intent(in)::A(*) !< CSR value array
    integer::nNz
    
    this%iA(1:this%nEq+1)=iA(1:this%nEq+1)
    nNz=iA(this%nEq+1)-1
    this%jA(1:nNz)=jA(1:nNz)
    this%A(1:nNz)=A(1:nNz)
  end subroutine
  
  !> set matrix with COO format for this generic linEq
  subroutine setCOOLinEq(this,iA,jA,A,nNz)
    class(linEq),intent(inout)::this !< this linEq
    integer,intent(in)::iA(*) !< COO row index array (size nNz)
    integer,intent(in)::jA(*) !< COO column index array (size nNz)
    double precision,intent(in)::A(*) !< COO value array
    integer,intent(in)::nNz !< number of non-zeros
    
    call coocsr(this%nEq,nNz,A,iA,jA,this%A,this%jA,this%iA)
  end subroutine
  
  !> generic destructor of linEq
  subroutine purgeLinEq(this)
    type(linEq),intent(inout)::this !< this linEq
    
    call this%clear()
  end subroutine
  
  !> constructor of ILU
  subroutine initILU(this,nEq,maxNNz,lFill)
    class(ILU),intent(inout)::this !< this linEq
    integer,intent(in)::nEq !< number of equations
    integer,intent(in)::maxNNZ !< maximum number of non-zeros
    integer,intent(in)::lFill !< maximum number of non-diagonal L or U elements in a row
    
    call this%initLinEq(nEq,maxNNz)
    this%lFill=min(lFill,this%nEq-1)
    this%iwk=this%nEq*min(2*lFill+1,this%nEq)
    allocate(this%alu(this%iwk))
    allocate(this%jlu(this%iwk))
    allocate(this%ju(this%nEq))
    allocate(this%jw(2*this%nEq))
    allocate(this%w(this%nEq+1))
  end subroutine
  
  !> ILUT factorize this ILU with drop tolerance eps
  subroutine factILU(this,eps)
    class(ILU),intent(inout)::this !< this linEq
    double precision,intent(in),optional::eps !< drop tolerance
    double precision::droptol
    integer::ierr
    
    if(present(eps))then
      droptol=eps
    else
      droptol=0d0
    end if
    call ilut(this%nEq,this%A,this%jA,this%iA,this%lFill,droptol,this%alu,this%jlu,this%ju,&
    &         this%iwk,this%w,this%jw,ierr)
  end subroutine
  
  !> solve this ILU which is already factorized
  subroutine solveILU(this,rhs,x)
    class(ILU),intent(inout)::this !< this linEq
    double precision,intent(in)::rhs(this%nEq) !< the right-hand-side
    double precision,intent(inout)::x(this%nEq) !< the solution
  
    call lusol(this%nEq,rhs,x,this%alu,this%jlu,this%ju)
  end subroutine
  
end module
