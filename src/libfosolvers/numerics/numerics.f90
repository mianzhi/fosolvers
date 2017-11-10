!----------------------------------------------------------------------------- best with 100 columns

!> algebraic solvers and interface to external algebraic solvers
module modNumerics
  use iso_c_binding
  private
  
  !> generic non-linear equations, i.e. x=f(x), solved by KINSOL
  type::noLinEq
    integer(C_INT),private::nEq=0 !< number of equations
    type(C_FUNPTR),private::f=C_NULL_FUNPTR !< pointer to the RHS function
    type(C_PTR),private::work=C_NULL_PTR !< pointer to KINSOL problem object
    type(C_PTR),private::x=C_NULL_PTR !< pointer to the solution N_Vector
    type(C_PTR),private::xScale=C_NULL_PTR !< pointer to the solution scaling N_Vector
    type(C_PTR),private::rScale=C_NULL_PTR !< pointer to the residual scaling N_Vector
  contains
    procedure::initNoLinEq ! specific class uses different arguments, saved the init() name
    procedure,public::clear=>clearNoLinEq
    procedure,public::setTol=>setTolNoLinEq
    procedure,public::setMaxIt=>setMaxItNoLinEq
    procedure,public::getNIt=>getNItNoLinEq
    final::purgeNoLinEq
  end type
  
  !> non-linear fixed point equations, i.e. x=f(x), solved by KINSOL fixed point iteration
  type,extends(noLinEq),public::fixPt
  contains
    procedure,public::init=>initFixPt
    procedure,public::solve=>solveFixPt
  end type
  
  !> non-linear equations, i.e. f(x)=0, solved by KINSOL Newton-Krylov method
  type,extends(noLinEq),public::NewtonKrylov
    type(C_FUNPTR),private::pSet=C_NULL_FUNPTR !< pointer to the preconditioner setter
    type(C_FUNPTR),private::pSolve=C_NULL_FUNPTR !< pointer to the preconditioner solver
  contains
    procedure,public::init=>initNewtonKrylov
    procedure,public::solve=>solveNewtonKrylov
  end type
  
  !> interface to SUNDIALS C functions
  interface
    
    !> NVECTOR create N_Vector object
    function n_vnew(nEq) bind(c,name='N_VNew_Serial')
      use iso_c_binding
      integer(C_INT),value::nEq !< number of equations
      type(C_PTR)::n_vnew !< value of the N_Vector pointer
    end function
    
    !> NVECTOR destroy N_Vector object, without nullify pointer value
    subroutine n_vdestroy(vector) bind(c,name='N_VDestroy_Serial')
      use iso_c_binding
      type(C_PTR),value::vector !< pointer to the N_Vector to be destroyed
    end subroutine
    
    !> NVECTOR set constant value in entries of the N_Vector object
    subroutine n_vconst(val,vector) bind(c,name='N_VConst_Serial')
      use iso_c_binding
      real(C_DOUBLE),value::val !< the constant value
      type(C_PTR),value::vector !< pointer to the N_Vector
    end subroutine
    
    !> NVECTOR get the vector length of the N_Vector object
    function n_vgetlength(vector) bind(c,name='N_VGetLength_Serial')
      use iso_c_binding
      type(C_PTR),value::vector !< pointer to the N_Vector
      integer(C_LONG)::n_vgetlength !< the length of vector
    end function
    
    !> NVECTOR get the data pointer of the N_Vector object, assuming data is contagious
    function n_vgetarraypointer(vector) bind(c,name='N_VGetArrayPointer_Serial')
      use iso_c_binding
      type(C_PTR),value::vector !< pointer to the N_Vector
      type(C_PTR)::n_vgetarraypointer !< the pointer to the contagious array
    end function
    
    !> KINSOL create memory object
    function kincreate() bind(c,name='KINCreate')
      use iso_c_binding
      type(C_PTR)::kincreate !< value of the memory pointer
    end function
    
    !> KINSOL initialize memory object
    function kininit(mem,func,tmpl) bind(c,name='KINInit')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_FUNPTR),value::func !< function pointer
      type(C_PTR),value::tmpl !< template N_Vector pointer
      integer(C_INT)::kininit !< error code
    end function
    
    !> KINSOL free memory object
    subroutine kinfree(mem) bind(c,name='KINFree')
      use iso_c_binding
      type(C_PTR)::mem !< location of memory pointer
    end subroutine
    
    !> KINSOL solve
    function kinsol(mem,x,method,xScale,rScale) bind(c,name='KINSol')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_PTR),value::x !< solution N_Vector
      integer(C_INT),value::method !< method code
      type(C_PTR),value::xScale !< solution scaling N_Vector
      type(C_PTR),value::rScale !< residual scaling N_Vector
      integer(C_INT)::kinsol !< returns error code
    end function
    
    !> KINSOL set tolerance of the residual
    function kinsetfuncnormtol(mem,tol) bind(c,name='KINSetFuncNormTol')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      real(C_DOUBLE),value::tol !< residual tolerance
      integer(C_INT)::kinsetfuncnormtol !< error code
    end function
    
    !> KINSOL set maximum number of iterations for the non-linear system
    function kinsetnummaxiters(mem,maxit) bind(c,name='KINSetNumMaxIters')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG),value::maxit !< maximum number of iterations
      integer(C_INT)::kinsetnummaxiters !< error code
    end function
    
    !> KINSOL set number of iterations for Anderson acceleration
    function kinsetmaa(mem,maa) bind(c,name='KINSetMAA')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG),value::maa !< number of iterations for Anderson acceleration
      integer(C_INT)::kinsetmaa !< error code
    end function
    
    !> KINSOL initialize direct dense linear solver
    function kindense(mem,nEq) bind(c,name='KINDense')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG),value::nEq !< dimension of the problem
      integer(C_INT)::kindense !< error code
    end function
    
    !> KINSOL initialize GMRES linear solver
    function kinspgmr(mem,maxl) bind(c,name='KINSpgmr')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_INT),value::maxl !< maximum dimension of the Krylov subspace
      integer(C_INT)::kinspgmr !< error code
    end function
    
    !> KINSOL get number of iterations performed for the non-linear system
    function kingetnumnonlinsolviters(mem,nit) bind(c,name='KINGetNumNonlinSolvIters')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG)::nit !< number of iterations
      integer(C_INT)::kingetnumnonlinsolviters !< error code
    end function
    
  end interface
  
  ! public procedures
  public::associateVector
  
contains
  
  !> associate fortran pointer v with N_Vector vector
  subroutine associateVector(vector,v)
    type(C_PTR)::vector !< the N_Vector object
    double precision,pointer,intent(inout)::v(:) !< fortran pointer
    type(C_PTR)::ptr
    integer(C_LONG)::n
    
    n=n_vgetlength(vector)
    ptr=n_vgetarraypointer(vector)
    call c_f_pointer(ptr,v,shape=[n])
  end subroutine
  
  !> generic constructor of noLinEq
  subroutine initNoLinEq(this,nEq,f)
    class(noLinEq),intent(inout)::this !< this noLinEq
    integer,intent(in)::nEq !< number of equations
    external::f !< the RHS subroutine in Fortran
    
    this%nEq=nEq
    this%f=c_funloc(f)
    this%work=kincreate()
    if(.not.c_associated(this%work))then
      write(*,'(a)')"[E] initNoLinEq(): KINSOL memory object not created"
      stop
    end if
    this%x=n_vnew(this%nEq)
    if(.not.c_associated(this%x))then
      write(*,'(a)')"[E] initNoLinEq(): solution N_Vector not allocated"
      stop
    end if
    this%xScale=n_vnew(this%nEq)
    if(.not.c_associated(this%xScale))then
      write(*,'(a)')"[E] initNoLinEq(): solution scaling N_Vector not allocated"
      stop
    end if
    this%rScale=n_vnew(this%nEq)
    if(.not.c_associated(this%rScale))then
      write(*,'(a)')"[E] initNoLinEq(): residual scaling N_Vector not allocated"
      stop
    end if
    ! no scaling by default
    call n_vconst(1d0,this%xScale)
    call n_vconst(1d0,this%rScale)
  end subroutine
  
  !> clear this generic noLinEq
  subroutine clearNoLinEq(this)
    class(noLinEq),intent(inout)::this !< this noLinEq
    
    this%nEq=0
    this%f=C_NULL_FUNPTR
    if(c_associated(this%work)) call kinfree(this%work)
    if(c_associated(this%x))then
      call n_vdestroy(this%x)
      this%x=C_NULL_PTR ! NOTE: NVECTOR does not nullify pointer value
    end if
    if(c_associated(this%xScale))then
      call n_vdestroy(this%xScale)
      this%xScale=C_NULL_PTR ! NOTE: NVECTOR does not nullify pointer value
    end if
    if(c_associated(this%rScale))then
      call n_vdestroy(this%rScale)
      this%rScale=C_NULL_PTR ! NOTE: NVECTOR does not nullify pointer value
    end if
  end subroutine
  
  !> set the residual tolerance for this noLinEq
  subroutine setTolNoLinEq(this,tol)
    class(noLinEq),intent(inout)::this !< this noLinEq
    double precision,intent(in)::tol !< residual tolerance
    real(C_DOUBLE)::c_tol
    integer(C_INT)::info
    
    c_tol=tol
    info=kinsetfuncnormtol(this%work,c_tol)
    if(info/=0)then
      write(*,'(a,i3)')"[E] setTolNoLinEq(): KINSetFuncNormTol error code ",info
      stop
    end if
  end subroutine
  
  !> set the maximum number of iterations for this noLinEq
  subroutine setMaxItNoLinEq(this,maxit)
    class(noLinEq),intent(inout)::this !< this noLinEq
    integer,intent(in)::maxit !< maximum number of iterations
    integer(C_LONG)::c_maxit
    integer(C_INT)::info
    
    c_maxit=maxit
    info=kinsetnummaxiters(this%work,c_maxit)
    if(info/=0)then
      write(*,'(a,i3)')"[E] setMaxItNoLinEq(): KINSetNumMaxIters error code ",info
      stop
    end if
  end subroutine
  
  !> get the number of iterations performed for this noLinEq
  subroutine getNItNoLinEq(this,nit)
    class(noLinEq),intent(inout)::this !< this noLinEq
    integer,intent(out)::nit !< number of iterations
    integer(C_LONG)::c_nit
    integer(C_INT)::info
    
    info=kingetnumnonlinsolviters(this%work,c_nit)
    nit=int(c_nit)
    if(info/=0)then
      write(*,'(a,i3)')"[E] setNItNoLinEq(): KINGetNumNonlinSolvIters error code ",info
      stop
    end if
  end subroutine
  
  !> generic destructor of noLinEq
  subroutine purgeNoLinEq(this)
    type(noLinEq),intent(inout)::this !< this noLinEq
    
    call this%clear()
  end subroutine
  
  !> initialize the problem of fixed point equations
  subroutine initFixPt(this,nEq,f,maa)
    class(fixPt),intent(inout)::this !< this fixPt
    integer,intent(in)::nEq !< number of equations
    external::f !< the RHS subroutine in Fortran
    integer,intent(in),optional::maa !< number of iterations for Anderson acceleration
    integer(C_LONG)::c_maa
    integer(C_INT)::info
    
    call this%noLinEq%initNoLinEq(nEq,f)
    if(present(maa))then
      c_maa=maa
      info=kinsetmaa(this%work,c_maa)
      if(info/=0)then
        write(*,'(a,i3)')"[E] initFixPt(): KINSetMAA error code ",info
        stop
      end if
    end if
    info=kininit(this%work,this%f,this%x) ! use solution vector as template
    if(info/=0)then
      write(*,'(a,i3)')"[E] initFixPt(): KINInit error code ",info
      stop
    end if
  end subroutine
  
  !> initialize the Newton-Krylov problem
  subroutine initNewtonKrylov(this,nEq,f,maxl)
    class(NewtonKrylov),intent(inout)::this !< this NewtonKrylov
    integer,intent(in)::nEq !< number of equations
    external::f !< the RHS subroutine in Fortran
    integer,intent(in),optional::maxl !< maximum dimension of the Krylov subspace
    integer(C_INT)::c_maxl,info
    
    call this%noLinEq%initNoLinEq(nEq,f)
    info=kininit(this%work,this%f,this%x) ! use solution vector as template
    if(info/=0)then
      write(*,'(a,i3)')"[E] initNewtonKrylov(): KINInit error code ",info
      stop
    end if
    if(present(maxl))then
      c_maxl=maxl
    else
      c_maxl=0
    end if
    info=kinspgmr(this%work,c_maxl)
    if(info/=0)then
      write(*,'(a,i3)')"[E] initNewtonKrylov(): KINSpgmr error code ",info
      stop
    end if
  end subroutine
  
  !> solve this system of fixed point equations
  subroutine solveFixPt(this,x,info)
    class(fixPt),intent(inout)::this !< this fixPt
    double precision,intent(inout)::x(:) !< the initial guess and solution
    integer,optional,intent(out)::info !< exit code
    double precision,pointer::xPtr(:)
    integer(C_INT)::c_info
    integer(C_INT),parameter::KIN_FP=3
    
    call associateVector(this%x,xPtr)
    xPtr(1:this%nEq)=x(1:this%nEq)
    c_info=kinsol(this%work,this%x,KIN_FP,this%xScale,this%rScale)
    if(c_info<0)then
      write(*,'(a,i3)')"[W] solveFixPt(): KINSol exit code ",c_info
    end if
    if(present(info))then
      info=c_info
    end if
    x(1:this%nEq)=xPtr(1:this%nEq)
  end subroutine
  
  !> solve this Newton-Krylov problem
  subroutine solveNewtonKrylov(this,x,info)
    class(NewtonKrylov),intent(inout)::this !< this NewtonKrylov
    double precision,intent(inout)::x(:) !< the initial guess and solution
    integer,optional,intent(out)::info !< exit code
    double precision,pointer::xPtr(:)
    integer(C_INT)::c_info
    integer(C_INT),parameter::KIN_LINESEARCH=1
    
    call associateVector(this%x,xPtr)
    xPtr(1:this%nEq)=x(1:this%nEq)
    c_info=kinsol(this%work,this%x,KIN_LINESEARCH,this%xScale,this%rScale)
    if(c_info<0)then
      write(*,'(a,i3)')"[W] solveNewtonKrylov(): KINSol exit code ",c_info
    end if
    if(present(info))then
      info=c_info
    end if
    x(1:this%nEq)=xPtr(1:this%nEq)
  end subroutine
  
end module
