!----------------------------------------------------------------------------- best with 100 columns

!> file I/O module
module modFileIO
  public
  
  !> read CGNS from fname into a generic object
  interface readCGNS
    !> read CGNS from fname into polyGrid
    subroutine readCGNSPolyGrid(fname,grid)
      use modPolyGrid
      character(*),intent(in)::fname !< file name
      class(polyGrid),intent(inout)::grid !< result
    end subroutine
  end interface
  
  !> read VTK from fid into a generic object
  interface readVTK
    !> read VTK from fid into polyX
    subroutine readVTKPolyX(fid,poly)
      use modPolyX
      integer,intent(in)::fid !< file id
      class(polyX),intent(inout)::poly !< result
    end subroutine
  end interface
  
  !> write VTK into fid from a generic object
  interface writeVTK
    !> write VTK into fid from poly
    subroutine writeVTKPolyX(fid,poly)
      use modPolyX
      integer,intent(in)::fid !< file id
      class(polyX),intent(in)::poly !< source
    end subroutine
    !> write VTK into fid from nodal scalar field
    subroutine writeVTKScal(fid,key,a)
      integer,intent(in)::fid !< file id
      character(*),intent(in)::key !< data name
      double precision,intent(in)::a(:) !< data
    end subroutine
    !> write VTK into fid from nodal vector field
    subroutine writeVTKVect(fid,key,a)
      integer,intent(in)::fid !< file id
      character(*),intent(in)::key !< data name
      double precision,intent(in)::a(:,:) !< data
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
