!----------------------------------------------------------------------------- best with 100 columns

!> interface to sparse linear solvers and ILU preconditioners (ITSOL)
module modSparse
  use iso_c_binding
  private
  
  ! constants
  integer,parameter,public::CSR_CLEAN=1 !< remove CSR duplicated/zero entries
  integer,parameter,public::CSR_CLEAN_PSORT=2 !< remove CSR duplicated/zero entries, partial sort
  integer,parameter,public::CSR_CLEAN_SORT=3 !< remove CSR duplicated/zero entries, full sort
  
  integer(C_INT),parameter::UMFPACK_CONTROL=20 !< size of umfpack control array
  integer(C_INT),parameter::UMFPACK_INFO=90 !< size of umfpack info array
  integer(C_INT),parameter::UMFPACK_AT=1 !< umfpack uses CSC format, so need to transpose
  
  !> generic CSR linear equation object
  type,public::linEq
    integer,private::nEq !< number of equations
    integer,private::maxNNz !< maximum number of non-zeros
    integer,allocatable,private::iA(:) !< CSR index array (size nEq+1)
    integer,allocatable,private::jA(:) !< CSR index array (size nNz)
    double precision,allocatable,private::A(:) !< CSR value array
  contains
    procedure,public::initLinEq
    procedure,public::clear=>clearLinEq
    procedure,public::setCSR=>setCSRLinEq
    procedure,public::setCOO=>setCOOLinEq
    final::purgeLinEq
  end type
  
  !> CSR linear equation object solved by multi-front method of umfpack
  type,extends(linEq),public::multiFront
    type(C_PTR),private::symbolic=C_NULL_PTR !< symbolic object
    type(C_PTR),private::numeric=C_NULL_PTR !< numeric object
    integer(C_INT),private,allocatable::ciA(:),cjA(:) !< 0-based index
    real(C_DOUBLE),private::control(UMFPACK_CONTROL) !< options
    real(C_DOUBLE),private::info(UMFPACK_INFO) !< exiting info
  contains
    procedure,public::init=>initMultiFront
    procedure,public::clear=>clearMultiFront
    procedure,public::fact=>factMultiFront
    procedure,public::solve=>solveMultiFront
    final::purgeMultiFront
  end type
  
  !> CSR linear equation object (incompletely) solved by ILUT of SPARSKIT2/ITSOL
  type,extends(linEq),public::ILU
    integer,allocatable,private::ju(:),jlu(:),jw(:) !< auxiliary storage for ILUT
    double precision,allocatable,private::alu(:),w(:) !< auxiliary storage for ILUT
    integer,private::lFill !< maximum number of non-diagonal L or U elements in a row
    integer,private::iwk !< length of alu, jlu
  contains
    procedure,public::init=>initILU
    procedure,public::clear=>clearILU
    procedure,public::fact=>factILU
    procedure,public::solve=>solveILU
    final::purgeILU
  end type
  
  !> interface to umfpack C functions and SPARSKIT2
  interface
    
    !> symbolic factorization of CSC format
    subroutine umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info)&
    &  bind(c,name='umfpack_di_symbolic')
      use iso_c_binding
      integer(C_INT),value::n_row,n_col
      integer(C_INT),intent(in)::Ap(*),Ai(*)
      real(C_DOUBLE),intent(in)::Ax(*)
      type(C_PTR)::Symbolic
      real(C_DOUBLE),intent(in)::Control(*)
      real(C_DOUBLE)::Info(*)
    end subroutine
    
    !> numerical factorization of CSC format
    subroutine umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)&
    &  bind(c,name='umfpack_di_numeric')
      use iso_c_binding
      integer(C_INT),intent(in)::Ap(*),Ai(*)
      real(C_DOUBLE),intent(in)::Ax(*)
      type(C_PTR),value::Symbolic
      type(C_PTR)::Numeric
      real(C_DOUBLE),intent(in)::Control(*)
      real(C_DOUBLE)::Info(*)
    end subroutine
    
    !> solve CSC format
    subroutine umfpack_di_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info)&
    &  bind(c,name='umfpack_di_solve')
      use iso_c_binding
      integer(C_INT),value::sys
      integer(C_INT),intent(in)::Ap(*),Ai(*)
      real(C_DOUBLE),intent(in)::Ax(*)
      real(C_DOUBLE)::X(*)
      real(C_DOUBLE),intent(in)::B(*)
      type(C_PTR),value::Numeric
      real(C_DOUBLE),intent(in)::Control(*)
      real(C_DOUBLE)::Info(*)
    end subroutine
    
    !> get default control array
    subroutine umfpack_di_defaults(Control) bind(c,name='umfpack_di_defaults')
      use iso_c_binding
      real(C_DOUBLE)::Control(*)
    end subroutine
    
    !> free symbolic factorization
    subroutine umfpack_di_free_symbolic(Symbolic) bind(c,name='umfpack_di_free_symbolic')
      use iso_c_binding
      type(C_PTR)::Symbolic
    end subroutine
    
    !> free numeric factorization
    subroutine umfpack_di_free_numeric(Numeric) bind(c,name='umfpack_di_free_numeric')
      use iso_c_binding
      type(C_PTR)::Numeric
    end subroutine
    
    !> convert COO format to CSR format
    subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
      integer,intent(in)::nrow,nnz
      double precision,intent(in)::a(*)
      integer,intent(in)::ir(*),jc(*)
      double precision,intent(inout)::ao(*)
      integer,intent(inout)::jao(*),iao(*)
    end subroutine
    
    !> CSR format cleanup and sort
    subroutine clncsr(job,value2,nrow,a,ja,ia,indu,iwk)
      integer,intent(in)::job,value2,nrow
      integer,intent(inout)::ia(nrow+1),ja(*)
      double precision,intent(inout)::a(*)
      integer,intent(inout)::indu(nrow),iwk(nrow+1)
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
    
    if(allocated(this%iA)) deallocate(this%iA)
    if(allocated(this%jA)) deallocate(this%jA)
    if(allocated(this%A)) deallocate(this%A)
  end subroutine
  
  !> set matrix with CSR format for this generic linEq
  subroutine setCSRLinEq(this,iA,jA,A,job)
    class(linEq),intent(inout)::this !< this linEq
    integer,intent(in)::iA(*) !< CSR index array (size nEq+1)
    integer,intent(in)::jA(*) !< CSR index array (size nNz)
    double precision,intent(in)::A(*) !< CSR value array
    integer,intent(in),optional::job !< additional cleanups
    integer::nNz
    integer,allocatable::indu(:),iwk(:)
    
    this%iA(1:this%nEq+1)=iA(1:this%nEq+1)
    nNz=iA(this%nEq+1)-1
    this%jA(1:nNz)=jA(1:nNz)
    this%A(1:nNz)=A(1:nNz)
    if(present(job))then
      allocate(indu(this%nEq))
      allocate(iwk(this%nEq+1))
      call clncsr(job,1,this%nEq,this%A,this%jA,this%iA,indu,iwk)
      deallocate(indu,iwk)
    end if
  end subroutine
  
  !> set matrix with COO format for this generic linEq
  subroutine setCOOLinEq(this,iA,jA,A,nNz,job)
    class(linEq),intent(inout)::this !< this linEq
    integer,intent(in)::iA(*) !< COO row index array (size nNz)
    integer,intent(in)::jA(*) !< COO column index array (size nNz)
    double precision,intent(in)::A(*) !< COO value array
    integer,intent(in)::nNz !< number of non-zeros
    integer,intent(in),optional::job !< additional cleanups
    integer,allocatable::indu(:),iwk(:)
    
    call coocsr(this%nEq,nNz,A,iA,jA,this%A,this%jA,this%iA)
    if(present(job))then
      allocate(indu(this%nEq))
      allocate(iwk(this%nEq+1))
      call clncsr(job,1,this%nEq,this%A,this%jA,this%iA,indu,iwk)
      deallocate(indu,iwk)
    end if
  end subroutine
  
  !> generic destructor of linEq
  subroutine purgeLinEq(this)
    type(linEq),intent(inout)::this !< this linEq
    
    call this%clear()
  end subroutine
  
  !> constructor of multiFront
  subroutine initMultiFront(this,nEq,maxNNz)
    class(multiFront),intent(inout)::this !< this multiFront
    integer,intent(in)::nEq !< number of equations
    integer,intent(in)::maxNNZ !< maximum number of non-zeros
    
    call this%initLinEq(nEq,maxNNz)
    allocate(this%ciA(nEq+1))
    allocate(this%cjA(maxNNZ))
    call umfpack_di_defaults(this%control)
  end subroutine
  
  !> clear this multiFront
  subroutine clearMultiFront(this)
    use iso_c_binding
    class(multiFront),intent(inout)::this !< this multiFront
    
    call this%linEq%clear()
    call umfpack_di_free_symbolic(this%symbolic)
    call umfpack_di_free_numeric(this%numeric)
    if(allocated(this%ciA)) deallocate(this%ciA)
    if(allocated(this%cjA)) deallocate(this%cjA)
  end subroutine
  
  !> destructor of multiFront
  subroutine purgeMultiFront(this)
    type(multiFront),intent(inout)::this !< this multiFront
    
    call this%clear()
  end subroutine
  
  !> factorize this multiFront with umfpack
  subroutine factMultiFront(this)
    use iso_c_binding
    class(multiFront),intent(inout)::this !< this multiFront
    integer::nnz
    
    nnz=this%iA(this%nEq+1)-1
    this%ciA(:)=this%iA(:)-1
    this%cjA(1:nnz)=this%jA(1:nnz)-1
    call umfpack_di_symbolic(this%nEq,this%nEq,this%ciA,this%cjA,this%A,this%symbolic,&
    &                        this%control,this%info)
    call umfpack_di_numeric(this%ciA,this%cjA,this%A,this%symbolic,this%numeric,&
    &                       this%control,this%info)
  end subroutine
  
  !> solve this multiFront which is already factorized
  subroutine solveMultiFront(this,rhs,x)
    use iso_c_binding
    class(multiFront),intent(inout)::this !< this multiFront
    double precision,intent(in)::rhs(this%nEq) !< the right-hand-side
    double precision,intent(inout)::x(this%nEq) !< the solution
    
    call umfpack_di_solve(UMFPACK_AT,this%ciA,this%cjA,this%A,x,rhs,this%numeric,&
    &                     this%control,this%info)
  end subroutine
  
  !> constructor of ILU
  subroutine initILU(this,nEq,maxNNz,lFill)
    class(ILU),intent(inout)::this !< this ILU
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
  
  !> clear this ILU
  subroutine clearILU(this)
    class(ILU),intent(inout)::this !< this ILU
    
    call this%linEq%clear()
    if(allocated(this%ju)) deallocate(this%ju)
    if(allocated(this%jlu)) deallocate(this%jlu)
    if(allocated(this%jw)) deallocate(this%jw)
    if(allocated(this%alu)) deallocate(this%alu)
    if(allocated(this%w)) deallocate(this%w)
  end subroutine
  
  !> destructor of ILU
  subroutine purgeILU(this)
    type(ILU),intent(inout)::this !< this ILU
    
    call this%clear()
  end subroutine
  
  !> ILUT factorize this ILU with drop tolerance eps
  subroutine factILU(this,eps)
    class(ILU),intent(inout)::this !< this ILU
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
    class(ILU),intent(inout)::this !< this ILU
    double precision,intent(in)::rhs(this%nEq) !< the right-hand-side
    double precision,intent(inout)::x(this%nEq) !< the solution
  
    call lusol(this%nEq,rhs,x,this%alu,this%jlu,this%ju)
  end subroutine
  
end module
