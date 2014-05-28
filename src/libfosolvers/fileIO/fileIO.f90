!----------------------------------------------------------------------------- best with 100 columns

!> file I/O module
module modFileIO
  public
  
  interface
    
    !> read GTS from fid into mesh
    subroutine readGTS(fid,mesh)
      use modpolyMesh
      integer,intent(in)::fid !< file id
      type(polyMesh),intent(inout)::mesh !< result
    end subroutine
    
  end interface
  
end module
