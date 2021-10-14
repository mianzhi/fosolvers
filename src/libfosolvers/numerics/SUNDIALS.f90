!----------------------------------------------------------------------------- best with 100 columns

!> interface to SUNDIALS
module modSUNDIALS
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
    procedure,public::setScale=>setScaleNoLinEq
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
    type(C_PTR),private::ls=C_NULL_PTR !< pointer to the linear solver
    type(C_FUNPTR),private::pSet=C_NULL_FUNPTR !< pointer to the preconditioner setter
    type(C_FUNPTR),private::pSol=C_NULL_FUNPTR !< pointer to the preconditioner solver
  contains
    procedure,public::init=>initNewtonKrylov
    procedure,public::solve=>solveNewtonKrylov
    procedure,public::clear=>clearNewtonKrylov
    procedure,public::setLazyPset=>setLazyPsetNewtonKrylov
    final::purgeNewtonKrylov
  end type
  
  !> generic ODE, i.e. dx/dt=f(x), solved by CVODE
  type::ODE
    integer(C_INT),private::nEq=0 !< number of equations
    type(C_FUNPTR),private::f=C_NULL_FUNPTR !< pointer to the RHS function
    type(C_PTR),private::work=C_NULL_PTR !< pointer to CVODE problem object
    type(C_PTR),private::x=C_NULL_PTR !< pointer to the solution N_Vector
  contains
    procedure::initODE ! specific class uses different arguments, saved the init() name
    procedure,public::clear=>clearODE
    procedure,public::setIV=>setIVODE
    procedure,public::setTol=>setTolODE
    procedure,public::setMaxSteps=>setMaxStepsODE
    procedure,public::getDt=>getDtODE
    final::purgeODE
  end type
  
  !> ODE (likely to be stiff), i.e. dx/dt=f(x) solved by CVODE BDF-Newton-Krylov method
  type,extends(ODE),public::BDFNewtonKrylov
    type(C_PTR),private::ls=C_NULL_PTR !< pointer to the linear solver
    type(C_FUNPTR),private::pSet=C_NULL_FUNPTR !< pointer to the preconditioner setter
    type(C_FUNPTR),private::pSol=C_NULL_FUNPTR !< pointer to the preconditioner solver
  contains
    procedure,public::init=>initBDFNewtonKrylov
    procedure,public::step=>stepBDFNewtonKrylov
    procedure,public::clear=>clearBDFNewtonKrylov
    final::purgeBDFNewtonKrylov
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
    
    !> KINSOL set linear solver
    function kinsetlinearsolver(mem,ls,mat) bind(c,name='KINSetLinearSolver')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_PTR),value::ls !< sparse iterative linear solver pointer
      type(C_PTR),value::mat !< matrix template
      integer(C_INT)::kinsetlinearsolver !< error code
    end function
    
    !> KINSOL set preconditoner setup and solve routines
    function kinsetpreconditioner(mem,pSet,pSol) bind(c,name='KINSetPreconditioner')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_FUNPTR),value::pSet !< preconditioner setup function pointer
      type(C_FUNPTR),value::pSol !< preconditioner solve function pointer
      integer(C_INT)::kinsetpreconditioner !< error code
    end function
    
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
    
    !> KINSOL set whether to skip initial preconditioner setup
    function kinsetnoinitsetup(mem,noInitSetup) bind(c,name='KINSetNoInitSetup')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      logical(C_BOOL),value::noInitSetup !< whether to skip
      integer(C_INT)::kinsetnoinitsetup !< error code
    end function
    
    !> KINSOL set maximum nonlinear iterations between preconditioner setups
    function kinsetmaxsetupcalls(mem,m) bind(c,name='KINSetMaxSetupCalls')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_INT),value::m !< maximum number of nonlinear iterations
      integer(C_INT)::kinsetmaxsetupcalls !< error code
    end function
    
    !> KINSOL get number of iterations performed for the non-linear system
    function kingetnumnonlinsolviters(mem,nit) bind(c,name='KINGetNumNonlinSolvIters')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG)::nit !< number of iterations
      integer(C_INT)::kingetnumnonlinsolviters !< error code
    end function
    
    !> CVODE create memory object
    function cvodecreate(m1,m2) bind(c,name='CVodeCreate')
      use iso_c_binding
      integer(C_INT),value::m1 !< ODE method
      integer(C_INT),value::m2 !< nonlinear equation method
      type(C_PTR)::cvodecreate !< value of the memory pointer
    end function
    
    !> CVODE initialize memory object
    function cvodeinit(mem,func,t0,y0) bind(c,name='CVodeInit')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_FUNPTR),value::func !< function pointer
      real(C_DOUBLE),value::t0 !< initial time value
      type(C_PTR),value::y0 !< initial state value N_Vector pointer
      integer(C_INT)::cvodeinit !< error code
    end function
    
    !> CVODE reinitialize the time and state
    function cvodereinit(mem,t0,y0) bind(c,name='CVodeReInit')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      real(C_DOUBLE),value::t0 !< initial time value
      type(C_PTR),value::y0 !< initial state value N_Vector pointer
      integer(C_INT)::cvodereinit !< error code
    end function
    
    !> CVODE free memory object
    subroutine cvodefree(mem) bind(c,name='CVodeFree')
      use iso_c_binding
      type(C_PTR)::mem !< location of memory pointer
    end subroutine
    
    !> CVODE set linear solver
    function cvodesetlinearsolver(mem,ls,mat) bind(c,name='CVodeSetLinearSolver')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_PTR),value::ls !< sparse iterative linear solver pointer
      type(C_PTR),value::mat !< matrix template
      integer(C_INT)::cvodesetlinearsolver !< error code
    end function
    
    !> CVODE set preconditoner setup and solve routines
    function cvodesetpreconditioner(mem,pSet,pSol) bind(c,name='CVodeSetPreconditioner')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      type(C_FUNPTR),value::pSet !< preconditioner setup function pointer
      type(C_FUNPTR),value::pSol !< preconditioner solve function pointer
      integer(C_INT)::cvodesetpreconditioner !< error code
    end function
    
    !> CVODE solve
    function cvode(mem,tOut,x,tReturn,task) bind(c,name='CVode')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      real(C_DOUBLE),value::tOut !< target time
      type(C_PTR),value::x !< solution N_Vector
      real(C_DOUBLE)::tReturn !< returned time
      integer(C_INT),value::task !< task code
      integer(C_INT)::cvode !< returns error code
    end function
    
    !> CVODE set scalar relative and absolute tolerance
    function cvodesstolerances(mem,rTol,aTol) bind(c,name='CVodeSStolerances')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      real(C_DOUBLE),value::rTol !< relative tolerance
      real(C_DOUBLE),value::aTol !< absolute tolerance
      integer(C_INT)::cvodesstolerances !< error code
    end function
    
    !> CVODE set maximum number of time steps
    function cvodesetmaxnumsteps(mem,maxStep) bind(c,name='CVodeSetMaxNumSteps')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG),value::maxStep !< maximum number of time steps
      integer(C_INT)::cvodesetmaxnumsteps !< error code
    end function
    
    !> CVODE get the size of the last time step
    function cvodegetlaststep(mem,dt) bind(c,name='CVodeGetLastStep')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      real(C_DOUBLE)::dt !< size of the last time step
      integer(C_INT)::cvodegetlaststep !< error code
    end function
    
    !> SUNLinearSolver create GMRES linear solver
    function sunlinsol_spgmr(tmpl,pType,maxl) bind(c,name='SUNLinSol_SPGMR')
      use iso_c_binding
      type(C_PTR),value::tmpl !< template N_Vector pointer
      integer(C_INT),value::pType !< preconditioning type
      integer(C_INT),value::maxl !< maximum dimensions of the Krylov subspace
      type(C_PTR)::sunlinsol_spgmr !< value of the linear solver pointer
    end function
    
    !> SUNLinearSolver free linear solver
    function sunlinsolfree(ls) bind(c,name='SUNLinSolFree')
      use iso_c_binding
      type(C_PTR),value::ls !< linear solver pointer
      integer(C_INT)::sunlinsolfree !< error code
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
    procedure(integer(C_INT))::f !< the RHS subroutine in Fortran
    
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
  
  !> set the solution and residual scaling vectors for this noLinEq
  subroutine setScaleNoLinEq(this,xScale,rScale)
    class(noLinEq),intent(inout)::this !< this noLinEq
    double precision,intent(in)::xScale(:) !< solution scaling vector
    double precision,intent(in)::rScale(:) !< residual scaling vector
    double precision,pointer::xScalePtr(:),rScalePtr(:)
    
    call associateVector(this%xScale,xScalePtr)
    call associateVector(this%rScale,rScalePtr)
    xScalePtr(1:this%nEq)=xScale(1:this%nEq)
    rScalePtr(1:this%nEq)=rScale(1:this%nEq)
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
    procedure(integer(C_INT))::f !< the RHS subroutine in Fortran
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
  subroutine initNewtonKrylov(this,nEq,f,maxl,pSet,pSol)
    class(NewtonKrylov),intent(inout)::this !< this NewtonKrylov
    integer,intent(in)::nEq !< number of equations
    procedure(integer(C_INT))::f !< the RHS subroutine in Fortran
    integer,intent(in),optional::maxl !< maximum dimension of the Krylov subspace
    procedure(integer(C_INT)),optional::pSet !< preconditioner setup function in Fortran
    procedure(integer(C_INT)),optional::pSol !< preconditioner solve function in Fortran
    integer(C_INT)::pOpt,c_maxl,info
    integer(C_INT),parameter::PREC_NONE=0
    integer(C_INT),parameter::PREC_RIGHT=2
    
    call this%noLinEq%initNoLinEq(nEq,f)
    info=kininit(this%work,this%f,this%x) ! use solution vector as template
    if(info/=0)then
      write(*,'(a,i3)')"[E] initNewtonKrylov(): KINInit error code ",info
      stop
    end if
    pOpt=merge(PREC_RIGHT,PREC_NONE,present(pSet).and.present(pSol))
    if(present(maxl))then
      c_maxl=maxl
    else
      c_maxl=0
    end if
    this%ls=sunlinsol_spgmr(this%x,pOpt,c_maxl) ! use solution vector as template
    if(.not.c_associated(this%ls))then
      write(*,'(a)')"[E] initNewtonKrylov(): linear solver object not allocated"
      stop
    end if
    info=kinsetlinearsolver(this%work,this%ls,C_NULL_PTR)
    if(info/=0)then
      write(*,'(a,i3)')"[E] initNewtonKrylov(): KINSetLinearSolver error code ",info
      stop
    end if
    if(present(pSet).and.present(pSol))then
      this%pSet=c_funloc(pSet)
      this%pSol=c_funloc(pSol)
      info=kinsetpreconditioner(this%work,this%pSet,this%pSol)
      if(info/=0)then
        write(*,'(a,i3)')"[E] initNewtonKrylov(): KINSetPreconditioner error code ",info
        stop
      end if
    end if
  end subroutine
  
  !> clear this NewtonKrylov
  subroutine clearNewtonKrylov(this)
    class(NewtonKrylov),intent(inout)::this !< this NewtonKrylov
    integer(C_INT)::info
    
    call this%noLinEq%clear()
    this%pSet=C_NULL_FUNPTR
    this%pSol=C_NULL_FUNPTR
    if(c_associated(this%ls))then
      info=sunlinsolfree(this%ls)
      if(info/=0)then
        write(*,'(a,i3)')"[E] clearNewtonKrylov(): SUNLinSolFree error code ",info
        stop
      end if
      this%ls=C_NULL_PTR
    end if
  end subroutine
  
  !> set this NewtonKrylov to skip initial prec. setup and has low prec. setup frequency
  subroutine setLazyPsetNewtonKrylov(this,skip,nit)
    class(NewtonKrylov),intent(inout)::this !< this NewtonKrylov
    logical,intent(in),optional::skip !< whether to skip initial prec. setup
    integer,intent(in),optional::nit !< maximum number of nonlinear iterations between prec. setups
    integer(C_INT)::info,c_nit
    logical(C_BOOL)::c_skip
    
    if(present(skip))then
      c_skip=skip
    else
      c_skip=.true.
    end if
    if(present(nit))then
      c_nit=nit
    else
      c_nit=100
    end if
    info=kinsetnoinitsetup(this%work,c_skip)
    if(info/=0)then
      write(*,'(a,i3)')"[E] setLazyPsetNewtonKrylov(): KINSetNoInitSetup error code ",info
      stop
    end if
    info=kinsetmaxsetupcalls(this%work,c_nit)
    if(info/=0)then
      write(*,'(a,i3)')"[E] setLazyPsetNewtonKrylov(): KINSetMaxSetupCalls error code ",info
      stop
    end if
  end subroutine
  
  !> destructor of NewtonKrylov
  subroutine purgeNewtonKrylov(this)
    type(NewtonKrylov),intent(inout)::this !< this NewtonKrylov
    
    call this%clear()
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
  
  !> generic constructor of ODE
  subroutine initODE(this,nEq,f)
    class(ODE),intent(inout)::this !< this ODE
    integer,intent(in)::nEq !< number of equations
    procedure(integer(C_INT))::f !< the RHS subroutine in Fortran
    
    this%nEq=nEq
    this%f=c_funloc(f)
    this%x=n_vnew(this%nEq)
    if(.not.c_associated(this%x))then
      write(*,'(a)')"[E] initODE(): solution N_Vector not allocated"
      stop
    end if
  end subroutine
  
  !> clear this generic ODE
  subroutine clearODE(this)
    class(ODE),intent(inout)::this !< this ODE
    
    this%nEq=0
    this%f=C_NULL_FUNPTR
    if(c_associated(this%work)) call cvodefree(this%work)
    if(c_associated(this%x))then
      call n_vdestroy(this%x)
      this%x=C_NULL_PTR ! NOTE: NVECTOR does not nullify pointer value
    end if
  end subroutine
  
  !> set the initial time and state values for this ODE
  subroutine setIVODE(this,t0,x0)
    class(ODE),intent(inout)::this !< this ODE
    double precision,intent(in)::t0 !< initial time
    double precision,intent(in)::x0(:) !< initial state
    real(C_DOUBLE)::c_t0
    double precision,pointer::xPtr(:)
    integer(C_INT)::c_info
    
    c_t0=t0
    call associateVector(this%x,xPtr)
    xPtr(1:this%nEq)=x0(1:this%nEq)
    c_info=cvodereinit(this%work,c_t0,this%x)
    if(c_info/=0)then
      write(*,'(a,i3)')"[E] setIVODE(): CVodeReInit exit code ",c_info
      stop
    end if
  end subroutine
  
  !> set the relative and absolute tolerance for this ODE
  subroutine setTolODE(this,rTol,aTol)
    class(ODE),intent(inout)::this !< this ODE
    double precision,intent(in)::rTol !< relative tolerance
    double precision,intent(in)::aTol !< absolute tolerance
    real(C_DOUBLE)::c_rTol,c_aTol
    integer(C_INT)::c_info
    
    c_rTol=rTol
    c_aTol=aTol
    c_info=cvodesstolerances(this%work,c_rTol,c_aTol)
    if(c_info/=0)then
      write(*,'(a,i3)')"[E] setTolODE(): CVodeSStolerances exit code ",c_info
      stop
    end if
  end subroutine
  
  !> set the maximum number of time steps for this ODE
  subroutine setMaxStepsODE(this,maxSteps)
    class(ODE),intent(inout)::this !< this ODE
    integer,intent(in)::maxSteps !< maximum number of steps
    integer(C_LONG)::c_maxSteps
    integer(C_INT)::c_info
    
    c_maxSteps=maxSteps
    c_info=cvodesetmaxnumsteps(this%work,c_maxSteps)
    if(c_info/=0)then
      write(*,'(a,i3)')"[E] setMaxStepsODE(): CVodeSetMaxNumSteps exit code ",c_info
      stop
    end if
  end subroutine
  
  !> get the size of the last time step for this ODE
  subroutine getDtODE(this,dt)
    class(ODE),intent(inout)::this !< this ODE
    double precision,intent(out)::dt !< the size of the last time step
    real(C_DOUBLE)::c_dt
    integer(C_INT)::c_info
    
    c_info=cvodegetlaststep(this%work,c_dt)
    dt=c_dt
    if(c_info/=0)then
      write(*,'(a,i3)')"[E] getDtODE(): CVodeGetLastStep exit code ",c_info
      stop
    end if
  end subroutine
  
  !> generic destructor of ODE
  subroutine purgeODE(this)
    type(ODE),intent(inout)::this !< this ODE
    
    call this%clear()
  end subroutine
  
  !> initialize this BDFNewtonKrylov
  subroutine initBDFNewtonKrylov(this,nEq,f,maxl,pSet,pSol)
    class(BDFNewtonKrylov),intent(inout)::this !< this NewtonKrylov
    integer,intent(in)::nEq !< number of equations
    procedure(integer(C_INT))::f !< the RHS subroutine in Fortran
    integer,intent(in),optional::maxl !< maximum dimension of the Krylov subspace
    procedure(integer(C_INT)),optional::pSet !< preconditioner setup function in Fortran
    procedure(integer(C_INT)),optional::pSol !< preconditioner solve function in Fortran
    integer(C_INT)::c_maxl,info,pOpt
    real(C_DOUBLE)::t=0d0
    integer(C_INT),parameter::CV_BDF=2
    integer(C_INT),parameter::CV_NEWTON=2
    integer(C_INT),parameter::PREC_NONE=0
    integer(C_INT),parameter::PREC_LEFT=1
    
    call this%ODE%initODE(nEq,f)
    this%work=cvodecreate(CV_BDF,CV_NEWTON)
    if(.not.c_associated(this%work))then
      write(*,'(a)')"[E] initBDFNewtonKrylov(): CVODE memory object not created"
      stop
    end if
    info=cvodeinit(this%work,this%f,t,this%x) ! use solution vector as template
    if(info/=0)then
      write(*,'(a,i3)')"[E] initBDFNewtonKrylov(): CVodeInit error code ",info
      stop
    end if
    pOpt=merge(PREC_LEFT,PREC_NONE,present(pSet).and.present(pSol))
    if(present(maxl))then
      c_maxl=maxl
    else
      c_maxl=0
    end if
    this%ls=sunlinsol_spgmr(this%x,pOpt,c_maxl) ! use solution vector as template
    if(.not.c_associated(this%ls))then
      write(*,'(a)')"[E] initBDFNewtonKrylov(): linear solver object not allocated"
      stop
    end if
    info=cvodesetlinearsolver(this%work,this%ls,C_NULL_PTR)
    if(info/=0)then
      write(*,'(a,i3)')"[E] initBDFNewtonKrylov(): CVodeSetLinearSolver error code ",info
      stop
    end if
    if(present(pSet).and.present(pSol))then
      this%pSet=c_funloc(pSet)
      this%pSol=c_funloc(pSol)
      info=cvodesetpreconditioner(this%work,this%pSet,this%pSol)
      if(info/=0)then
        write(*,'(a,i3)')"[E] initBDFNewtonKrylov(): CVodeSetPreconditioner error code ",info
        stop
      end if
    end if
  end subroutine
  
  !> clear this BDFNewtonKrylov
  subroutine clearBDFNewtonKrylov(this)
    class(BDFNewtonKrylov),intent(inout)::this !< this NewtonKrylov
    integer(C_INT)::info
    
    call this%ODE%clear()
    this%pSet=C_NULL_FUNPTR
    this%pSol=C_NULL_FUNPTR
    if(c_associated(this%ls))then
      info=sunlinsolfree(this%ls)
      if(info/=0)then
        write(*,'(a,i3)')"[E] clearBDFNewtonKrylov(): SUNLinSolFree error code ",info
        stop
      end if
      this%ls=C_NULL_PTR
    end if
  end subroutine
  
  !> destructor of BDFNewtonKrylov
  subroutine purgeBDFNewtonKrylov(this)
    type(BDFNewtonKrylov),intent(inout)::this !< this NewtonKrylov
    
    call this%clear()
  end subroutine
  
  !> step this BDFNewtonKrylov
  subroutine stepBDFNewtonKrylov(this,tOut,t,x,info)
    class(BDFNewtonKrylov),intent(inout)::this !< this BDFNewtonKrylov
    double precision,intent(in)::tOut !< target time
    double precision,intent(out)::t !< returned time
    double precision,intent(inout)::x(:) !< the solution
    integer,optional,intent(out)::info !< exit code
    double precision,pointer::xPtr(:)
    real(C_DOUBLE)::c_tOut,c_t
    integer(C_INT)::c_info
    integer(C_INT),parameter::CV_ONE_STEP=2
    
    c_tOut=tOut
    call associateVector(this%x,xPtr)
    c_info=cvode(this%work,c_tOut,this%x,c_t,CV_ONE_STEP)
    if(c_info<0)then
      write(*,'(a,i3)')"[W] stepBDFNewtonKrylov(): CVode exit code ",c_info
    end if
    if(present(info))then
      info=c_info
    end if
    t=c_t
    x(1:this%nEq)=xPtr(1:this%nEq)
  end subroutine
  
end module
