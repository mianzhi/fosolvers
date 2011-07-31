!----------------------------------------------------------------------------- best with 100 columns

module moduleGrid
  private
  
  type,public::typeNode
    double precision Pos(3)
  end type
  
  type,public::typePoint
    integer NodeInd
    integer GeoEnti
  contains
    procedure,public::findPos=>findPointPos
  end type
  
contains
  subroutine findPointPos(this,rst)
    class(typePoint),intent(in)::this
    double precision,intent(out)::rst(3)
    
    rst=[1d0,1d0,1d0]
  end subroutine
end module
