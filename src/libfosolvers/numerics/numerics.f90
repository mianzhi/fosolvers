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
    type(C_PTR),private::tmpl=C_NULL_PTR !< pointer to a template N_Vector
  contains
    procedure,public::init=>initFixPtEq
    procedure,public::clear=>clearFixPtEq
    final::purgeFixPtEq
  end type
  
  !> interface to SUNDIALS C functions
  interface
    
    !> NVECTOR create N_Vector object
    function n_vnew_serial(nEq) bind(c,name='N_VNew_Serial')
      use iso_c_binding
      integer(C_INT),value::nEq !< number of equations
      type(C_PTR)::n_vnew_serial !< value of the N_Vector pointer
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
    
    !> KINSOL set number of iterations for Anderson acceleration
    function kinsetmaa(mem,maa) bind(c,name='KINSetMAA')
      use iso_c_binding
      type(C_PTR),value::mem !< memory pointer
      integer(C_LONG)::maa !< number of iterations for Anderson acceleration
      integer(C_INT)::kinsetmaa !< error code
    end function
    
  end interface
  
contains
  
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
    this%tmpl=n_vnew_serial(this%nEq)
    if(.not.c_associated(this%tmpl))then
      write(*,'(a)')"[E] initFixPtEq(): template N_Vector not allocated"
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
    info=kininit(this%work,this%f,this%tmpl) ! FIXME should use a dummy function compatible with KINSOL
    if(info/=0)then
      write(*,'(a,i3)')"[E] initFixPtEq(): KINInit error code ",info
      stop
    end if
  end subroutine
  
  !> clear this system of fixed point equations
  subroutine clearFixPtEq(this)
    class(fixPtEq),intent(inout)::this !< this fixPtEq
    
    this%nEq=0
    this%f=C_NULL_FUNPTR
    if(c_associated(this%work)) call kinfree(this%work)
  end subroutine
  
  !> destructor of fixPtEq
  subroutine purgeFixPtEq(this)
    type(fixPtEq),intent(inout)::this !< this fixPtEq
    
    call this%clear()
  end subroutine
  
end module
