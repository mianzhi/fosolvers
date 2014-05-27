!----------------------------------------------------------------------------- best with 100 columns

function polyX1()
  use modPolyX
  integer polyX1
  type(polyX)::a
  
  polyX1=0
  call a%init(4,1,4)
  a%pN(:,1)=[0d0,0d0,0d0]
  a%pN(:,2)=[1d0,0d0,0d0]
  a%pN(:,3)=[0d0,2d0,0d0]
  a%pN(:,4)=[0d0,0d0,3d0]
  a%nNE(1)=4
  a%iNE(:,1)=[1,2,3,4]
  if(norm2(a%p(1)-[0.25d0,0.5d0,0.75d0])>tiny(1d0))then
    polyX1=polyX1+1
  end if
  call a%clear()
end function
