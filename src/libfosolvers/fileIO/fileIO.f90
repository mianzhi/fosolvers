!----------------------------------------------------------------------------- best with 100 columns

!> file I/O module
module modFileIO
  public
  
  ! read VTK from fid into a generic object
  interface readVTK
    ! read into polyX
    subroutine readVTKPolyX(fid,poly)
      use modPolyX
      integer,intent(in)::fid !< file id
      class(polyX),intent(inout)::poly !< result
    end subroutine
  end interface
  
  ! not generic interface
  interface
    
    !> read GTS from fid into mesh
    subroutine readGTS(fid,mesh)
      use modPolyMesh
      integer,intent(in)::fid !< file id
      type(polyMesh),intent(inout)::mesh !< result
    end subroutine
    
  end interface
  
end module
