!----------------------------------------------------------------------------- best with 100 columns

!> algebraic solvers and interface to external algebraic solvers
module modNumerics
  use iso_c_binding
  private
  
  !> non-linear fixed point equations, i.e. x=f(x), solved by KINSOL fixed point iteration
  type,public::fixPtEq
    integer(C_INT),private::nEq=0 !< number of equations
    type(C_FUNPTR),private::f=C_NULL_FUNPTR !< pointer to the RHS function
    type(C_PTR),private::work=C_NULL_PTR !< pointer to KINSOL problem object
    type(C_PTR),private::x=C_NULL_PTR !< pointer to the solution N_Vector
    type(C_PTR),private::xScale=C_NULL_PTR !< pointer to the solution scaling N_Vector
    type(C_PTR),private::rScale=C_NULL_PTR !< pointer to the residual scaling N_Vector
  contains
    procedure,public::init=>initFixPtEq
    procedure,public::clear=>clearFixPtEq
    procedure,public::solve=>solveFixPtEq
    final::purgeFixPtEq
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
    
    !> KINSOL set number of iterations for Anderson acceleration
    function kinsetmaa(mem,maa) bind(c,name='KINSetMAA')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG)::maa !< number of iterations for Anderson acceleration
      integer(C_INT)::kinsetmaa !< error code
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
  
  !> initialize the problem of fixed point equations
  subroutine initFixPtEq(this,nEq,f,maa)
    class(fixPtEq),intent(inout)::this !< this fixPtEq
    integer,intent(in)::nEq !< number of equations
    external::f !< the RHS subroutine in Fortran
    integer,intent(in),optional::maa !< number of iterations for Anderson acceleration
    integer(C_LONG)::c_maa
    integer(C_INT)::info
    
    this%nEq=nEq
    this%f=c_funloc(f)
    this%work=kincreate()
    if(.not.c_associated(this%work))then
      write(*,'(a)')"[E] initFixPtEq(): KINSOL memory object not created"
      stop
    end if
    this%x=n_vnew(this%nEq)
    if(.not.c_associated(this%x))then
      write(*,'(a)')"[E] initFixPtEq(): solution N_Vector not allocated"
      stop
    end if
    this%xScale=n_vnew(this%nEq)
    if(.not.c_associated(this%xScale))then
      write(*,'(a)')"[E] initFixPtEq(): solution scaling N_Vector not allocated"
      stop
    end if
    this%rScale=n_vnew(this%nEq)
    if(.not.c_associated(this%rScale))then
      write(*,'(a)')"[E] initFixPtEq(): residual scaling N_Vector not allocated"
      stop
    end if
    if(present(maa))then
      c_maa=maa
      info=kinsetmaa(this%work,c_maa)
      if(info/=0)then
        write(*,'(a,i3)')"[E] initFixPtEq(): KINSetMAA error code ",info
        stop
      end if
    end if
    info=kininit(this%work,this%f,this%x) ! use solution vector as template
    if(info/=0)then
      write(*,'(a,i3)')"[E] initFixPtEq(): KINInit error code ",info
      stop
    end if
    ! no scaling by default
    call n_vconst(1d0,this%xScale)
    call n_vconst(1d0,this%rScale)
  end subroutine
  
  !> clear this system of fixed point equations
  subroutine clearFixPtEq(this)
    class(fixPtEq),intent(inout)::this !< this fixPtEq
    
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
  
  !> solve this system of fixed point equations
  subroutine solveFixPtEq(this,x)
    class(fixPtEq),intent(inout)::this !< this fixPtEq
    double precision,intent(inout)::x(:) !< the initial guess and solution
    double precision,pointer::xPtr(:)
    integer(C_INT)::info
    integer(C_INT),parameter::KIN_FP=3
    
    call associateVector(this%x,xPtr)
    xPtr(1:this%nEq)=x(1:this%nEq)
    info=kinsol(this%work,this%x,KIN_FP,this%xScale,this%rScale)
    if(info/=0)then
      write(*,'(a,i3)')"[E] solveFixPtEq(): KINSol error code ",info
      stop
    end if
    x(1:this%nEq)=xPtr(1:this%nEq)
  end subroutine
  
  !> destructor of fixPtEq
  subroutine purgeFixPtEq(this)
    type(fixPtEq),intent(inout)::this !< this fixPtEq
    
    call this%clear()
  end subroutine
  
end module
