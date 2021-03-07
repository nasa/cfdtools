!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module packages data structures suited to [some of] the legacy serial
!  file formats supported by the Visualization Toolkit (VTK), along with the
!  Fortran 90 routine[s] for reading such files.  It was prompted by program
!  SHADOWGRAPH in order to read outputs from the US3D flow solver.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module VTK_file_header_module  ! Header items common to all legacy VTK files

   type VTK_file_header_type
      character (len=26)  :: VTK_version     ! # vtk DataFile Version x.x
      character (len=128) :: title           ! Embedded blanks permitted
      character (len=6)   :: file_format     ! ASCII | BINARY
      character (len=17)  :: dataset_type    ! E.g., UNSTRUCTURED_GRID
   end type VTK_file_header_type

   end module VTK_file_header_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module VTK_unstructured_grid_module

   integer, parameter :: mxscalar = 8       ! Limit on # scalar functions

   type unstructured_type
      integer           :: npoints          ! # points in the unstructured grid
      integer           :: ncells           ! # cells  in the unstructured grid
      integer           :: nlist            ! # integers in all (n,i1:in) tuples
      integer           :: nscalar          ! # scalar fns. at each (x,y,z) pt.
      logical           :: cell_centered    ! T if function data not vrtx-cntrd.
      character (len=8) :: fname(mxscalar)  ! Avoid pointers here
      real,    pointer  :: x(:), y(:), z(:) ! Coordinates for all grid points
      integer, pointer  :: list(:)          ! Packed list of (n,i1:in) tuples
      integer, pointer  :: cell_type(:)     ! E.g., 12 = VTK_HEXAHEDRON
      real,    pointer  :: f(:,:)           ! Scalar function data for all grid
   end type unstructured_type               ! points (1:npoints,1:nscalar)

   end module VTK_unstructured_grid_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module VTK_IO_module

   public :: read_VTK_unstructured        ! Read one binary unstructured dataset
   public :: VTK_centers_to_vertices      ! Convert function data in-place

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_VTK_unstructured (lundat, lunerr, header, us_grid, ios)

!     Read one file in a VTK binary stream format as written by the US3D flow
!     solver.  The file is expected to contain (x,y,z) coordinates for all
!     grid points (only one unstructured grid in VTK terminology), along with
!     density and its x/y/z partial derivatives, although the routine will read
!     any number of flow variables.  US3D grids may be so large that multiple
!     such files may be required by serial applications like SHADOWGRAPH with
!     limited memory on typical processors.
!
!     The legacy/serial unstructured grid format allows for interleaving cells
!     of more than one type.  These cells are defined by a packed set of tuples,
!     with each tuple containing a cell vertex count n followed by n indices
!     into the 1-D arrays for x, y, z and into the function array (:,1:nscalar).
!
!     While indices in the VTK file being read are expected to start at 0, this
!     routine adds 1 to all point indices returned in us_grid%list(:) in order
!     to match the indexing of %x/y/z(:) and %f(:,:) that always starts at 1.
!
!     See the available VTK File Formats documentation on-line for more details.
!
!     The file being read is expected to be open upon entry, and is not closed
!     here (possibly in a loop over such files).  All storage is allocated here.
!
!     Current limitations:
!
!     o  Only scalar functions are handled at present (no vector functions).
!     o  The number of scalar functions in the file is expected to be input
!        to this routine as us_grid%nscalar.
!     o  Function data may be cell- or vertex-centered, but only one type is
!        presently handled in a given file
!     o  Converting cell-centered function data to cell vertices (only) is an
!        option
!
!     History:
!
!     07/08/10  D.A.Saunders  Initial adaptation of a FieldView reading routine.
!                             Handle a known number of scalar functions only for
!                             now.  If more generality is needed, a prescanning
!                             routine is suggested to count how many functions
!                             are present in the file.  We avoid counting here.
!     07/15/10     "    "     Return us_grid%cell_centered so the application
!                             can call new VTK_centers_to_vertices if necessary.
!
!     Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use VTK_file_header_module
      use VTK_unstructured_grid_module

      implicit none

!     Arguments:

      integer, intent (in)        :: lundat   ! For file being read
      integer, intent (in)        :: lunerr   ! For error messages
      type (VTK_file_header_type) :: header   ! VTK file header derived type
      type (unstructured_type)    :: us_grid  ! Unstructured grid derived type
      integer, intent (out)       :: ios      ! 0 => no error detected

!     Local constants:

      integer, parameter :: mxchars = 128  ! Plenty for likely character data

!     Local variables:

      integer   :: i, icell, il, isize, n, ncells, nchars, nf, nlist, npoints
      character :: buffer*(mxchars), key*17, keyword*17, LF*1, precision*6

!     Execution:
!     ----------

      key    = '#'   ! For error message purposes
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      header%VTK_version = buffer(1:nchars)

      key    = 'Title'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      header%title = buffer(1:nchars)

      key    = 'BINARY'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      keyword(1:6) = buffer(1:6)
      if (keyword(1:6) /= key(1:6)) go to 930  ! Deal with ASCII only if reqd.

      header%file_format = buffer(1:nchars)

      key    = 'DATASET'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      keyword(1:7) = buffer(1:7)
      if (keyword(1:7) /= key(1:7)) go to 930

      key     = 'UNSTRUCTURED_GRID'
      keyword = buffer(9:25)
      if (keyword /= key) go to 930

      header%dataset_type = key

!     Look for this type of record:   POINTS 1234567 double
!     -----------------------------------------------------

      key    = 'POINTS'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      read (buffer(1:nchars), *, iostat=ios) keyword, npoints, precision
      if (ios /= 0) go to 920
      if (keyword(1:6) /= key(1:6)) go to 930

      write (lunerr, '(a, i10)') '# VTK file points to be read:', npoints
      us_grid%npoints = npoints
      isize           = npoints

!     Allocate and read the point (x,y,z) coordinates as triples:

      allocate (us_grid%x(npoints), us_grid%y(npoints), us_grid%z(npoints), &
                stat=ios)
      if (ios /= 0) go to 940

      read (lundat, iostat=ios) &
         (us_grid%x(i), us_grid%y(i), us_grid%z(i), i = 1, npoints)
      if (ios /= 0) go to 950

!     Look for this type of record:   CELLS 1000000 9000000
!     -----------------------------------------------------

      key    = 'CELLS'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      read (buffer(1:nchars), *, iostat=ios) keyword, ncells, nlist
      if (ios /= 0) go to 920
      if (keyword(1:5) /= key(1:5)) go to 930

      write (lunerr, '(a, i10)') '# VTK file cells to be read:', ncells
      us_grid%ncells = ncells
      us_grid%nlist  = nlist
      isize          = nlist

!     Allocate and read the list of tuples defining all cells:

      allocate (us_grid%list(nlist), stat=ios)
      if (ios /= 0) go to 940

      read (lundat, iostat=ios) us_grid%list(:)
      if (ios /= 0) go to 950

!     Make all the grid point indices start at 1 to match %x/y/z(:) and %f(:,:):

      il = 1
      do icell = 1, ncells
         n = us_grid%list(il)
         us_grid%list(il+1:il+n) = us_grid%list(il+1:il+n) + 1
         il = il + n + 1
      end do

!     Look for this type of record:   CELL_TYPES 1000000
!     --------------------------------------------------

      key    = 'CELL_TYPES'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      read (buffer(1:nchars), *, iostat=ios) keyword, isize
      if (ios /= 0) go to 920
      if (keyword(1:5) /= key(1:5)) go to 930
      if (isize /= ncells) go to 960

!     Allocate and read the cell types:

      allocate (us_grid%cell_type(ncells), stat=ios)
      if (ios /= 0) go to 940

      read (lundat, iostat=ios) us_grid%cell_type(:)
      if (ios /= 0) go to 950

!     Allow for function data to be vertex- or cell-centered, but not both.
!     Look for a record like    POINT_DATA 1234567   or   CELL_DATA 1000000
!     ---------------------------------------------------------------------

      key    = 'POINT_DATA'
      nchars = mxchars

      call read_to_line_feed (lundat, buffer, nchars, ios)
      if (ios /= 0) go to 910

      read (buffer(1:nchars), *, iostat=ios) keyword, isize
      if (ios /= 0) go to 920

      select case (keyword(1:10))

         case ('POINT_DATA')

            us_grid%cell_centered = .false.
            if (isize /= npoints) go to 970

         case ('CELL_DATA ')

            us_grid%cell_centered = .true.
            if (isize /= ncells) go to 960

         case default

            write (lunerr, '(/, 2a)') 'Expecting POINT- or CELL_DATA; got ', &
               keyword
            ios = 8
            go to 999

      end select

!     The following could be modularized and called in a loop of length 2 over
!     the above case statement.  For now, only one case is permitted per file,
!     so keep it all together here.

!     Allocate and read the known number of scalar functions.
!     -------------------------------------------------------
!     If the application does not know that number, a preliminary scan routine
!     would be one possibility.  Here, we avoid any such scanning.

      nf  = min (mxscalar, us_grid%nscalar)
      key = keyword  ! For the error message

      allocate (us_grid%f(isize,nf), stat=ios)
      if (ios /= 0) go to 940

      do n = 1, nf

!        Look for this type of record:   SCALARS varname double [1]

         key    = 'SCALARS'
         nchars = mxchars

         call read_to_line_feed (lundat, buffer, nchars, ios)
         if (ios /= 0) go to 910

         read (buffer(1:nchars), *, iostat=ios) keyword, us_grid%fname(n), &
            precision
         if (ios /= 0) go to 920
         if (keyword(1:7) /= key(1:7)) go to 930

!        Look for this type of record:   LOOKUP_TABLE default

         key    = 'LOOKUP_TABLE'
         nchars = mxchars

         call read_to_line_feed (lundat, buffer, nchars, ios)
         if (ios /= 0) go to 910

         read (buffer(1:nchars), *, iostat=ios) keyword  ! Non-default table
         if (ios /= 0) go to 920                         ! not handled
         if (keyword(1:12) /= key(1:12)) go to 930

!        Read this set of scalar function values:

         read (lundat, iostat=ios) us_grid%f(:,n)
         if (ios /= 0) go to 950

      end do  ! Next scalar function

!     Verify EOF?  Don't bother.  Rethink if vector data & cell- and/or vertex-
!     centered data are ever handled, as they are not expected to be.

      go to 999

 910  write (lunerr, '(/, 2a)') 'Error looking for VTK keyword ', trim (key)
      go to 999

 920  write (lunerr, '(/, 4a)') 'Error decoding VTK ', trim (key), &
         ' record:  ', buffer(1:nchars)
      go to 999

 930  write (lunerr, '(/, 4a)') 'Expecting VTK keyword ', trim (key), &
         ' but found ', keyword
      ios = 3
      go to 999

 940  write (lunerr, '(/, 3a, i10)') 'Error allocating for VTK file data ', &
         trim (key), '; size:', isize
      go to 999

 950  write (lunerr, '(/, 2a)') 'Error reading data for keyword ', trim (key)
      go to 999

 960  write (lunerr, '(/, 2a, 2i10)') trim (key), ' count /= ncells:', &
         isize, ncells
      ios = 6
      go to 999

 970  write (lunerr, '(/, 2a, 2i10)') trim (key), ' count /= npoints:', &
         isize, npoints
      ios = 7

 999  return

      end subroutine read_VTK_unstructured

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine VTK_centers_to_vertices (lundat, lunerr, header, us_grid, ios)

!     Convert the function data in a VTK unstructured dataset from cell centers
!     to vertices.  Inverse distance weighted averaging is employed.  Function
!     values are updated in-place, but it is not practical to avoid considerable
!     extra storage temporarily.

!     History:
!
!     07/15/10  D.A.Saunders  US3D function data turned out to be cell-centered.
!                             While copying each function value to all vertices
!                             of each cell produced reasonable results on a
!                             coarse grid test case, we provide this option to
!                             avoid hints of discontinuities in shadowgraph
!                             images (and possibly elsewhere some day).
!
!     Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use VTK_file_header_module
      use VTK_unstructured_grid_module

      implicit none

!     Arguments:

      integer, intent (in)        :: lundat  ! For scratch file I/O, closed here
      integer, intent (in)        :: lunerr  ! For error messages
      type (VTK_file_header_type) :: header  ! VTK file header derived type
      type (unstructured_type)    :: us_grid ! Unstructured grid derived type
      integer, intent (out)       :: ios     ! 0 => no error detected

!     Local constants:

      real, parameter :: one = 1., zero = 0.

!     Local variables:

      integer   :: i, ic, il, iv, ncells, nf, npoints, nv
      integer   :: list(8)
      real      :: rnv, w, xbar, ybar, zbar
      real, allocatable :: sumw(:), sumwf(:,:)

!     Execution:
!     ----------

      npoints = us_grid%npoints
      ncells  = us_grid%ncells
      nf      = us_grid%nscalar

      allocate (sumw(npoints), sumwf(npoints,nf), stat=ios)
      if (ios /= 0) go to 910

      sumw(:)    = zero
      sumwf(:,:) = zero

      il = 1
      do ic = 1, ncells
         nv = us_grid%list(il)
         list(1:nv) = us_grid%list(il+1:il+nv)

         xbar = zero;  ybar = zero;  zbar = zero
         do iv = 1, nv
            xbar = xbar + us_grid%x(list(iv))
            ybar = ybar + us_grid%y(list(iv))
            zbar = zbar + us_grid%z(list(iv))
         end do
         rnv  = one / real (nv)
         xbar = xbar * rnv
         ybar = ybar * rnv
         zbar = zbar * rnv

         do iv = 1, nv
            i = list(iv)
            w = one / sqrt ((us_grid%x(i) - xbar)**2 + &
                            (us_grid%y(i) - ybar)**2 + &
                            (us_grid%z(i) - zbar)**2)
            sumw(i)    = sumw(i)    + w
            sumwf(i,:) = sumwf(i,:) + w*us_grid%f(ic,:)
         end do
         il = il + nv + 1
      end do

      deallocate (us_grid%f)

      do i = 1, npoints
         sumwf(i,:) = sumwf(i,:) / sumw(i)
      end do

      deallocate (sumw)

      allocate (us_grid%f(npoints,nf), stat=ios)
      if (ios /= 0) go to 920

      us_grid%f(:,:) = sumwf(:,:)

      deallocate (sumwf)
      go to 999

 910  write (lunerr, '(/, a, i10, a, i3)') &
         'VTK_centers_to_vertices:  cannot allocate sumw & sumwf workspace:', &
         npoints, '  x', nf + 1
      go to 999

 920  write (lunerr, '(/, a, i10, a, i3)') &
         'VTK_centers_to_vertices:  cannot allocate vertex-centered f:', &
         npoints, '  x', nf

 999  return

      end subroutine VTK_centers_to_vertices

   end module VTK_IO_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
