!----------------------------------------------------------------------------- best with 100 columns

!> solve nonlinear system
module moduleNonlinearSolve
  private
  
  ! variables
  integer,public::nEq
  
  ! procedures
  public::solveNonlinear
  
  !> model of problem function
  interface
    function modelProblemFunc(r)
      double precision r(:)
      double precision modelProblemFunc(size(r))
    end function
  end interface
  procedure(modelProblemFunc),pointer::ProblemFunc=>null() !< the problem to be solved
  public::ProblemFunc
  
contains
  
  !> solve nonlinear system \f$ F(u)=0 \f$ using sundials-kinsol
  subroutine solveNonlinear(u)
    double precision,intent(inout)::u(:) !< initial guess and final solution
    integer*8 iout(15)
    integer ier
    double precision rout(2),koefScal(size(u))
    
    nEq=size(u)
    koefScal(:)=1d0
    call fnvinits(3,nEq,ier)
    call fkinmalloc(iout,rout,ier)
    call fkinspgmr(50,10,ier)
    call fkinsol(u,1,koefScal,koefScal,ier)
    call fkinfree()
  end subroutine
  
end module

!> the function used by sundials-kinsol \f$ fval=fkfun(u) \f$
subroutine fkfun(u,fval,ier)
  use moduleNonlinearSolve
  double precision::u(*) !< input
  double precision::fval(*) !< output
  integer::ier !< error flag
  
  if(associated(ProblemFunc))then
    fval(1:nEq)=ProblemFunc(u(1:nEq))
    ier=0
  else
    ier=1
  end if
end subroutine
