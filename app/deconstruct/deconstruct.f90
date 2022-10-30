!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program deconstruct
!
!     Description:
!
!        (Originally, for testing unstructured grid forms oF ADT searching:)
!
!        Convert a multiblock structured surface or volume grid to a surface
!     triangulation or tetrahedral volume mesh.  Quadrilaterals become pairs
!     of triangles; hex cells become five tetrahedra each.
!
!        (This version:)
!
!        Preserving the input quad surface cells or hex volume cells is now
!     an option, as needed for use by US3D,  Avoiding duplicate vertices is
!     an added complication, that has initially been handled with a "brute-
!     force" method that will be slow on a large grid.
!
!        Outputs are in Tecplot format, as a single zone.
!
!        This version now allows for an accompanying function file that has
!     to be vertex-centered.  Note that a cell-centered output from DPLR is
!     really the grid formed by the cell centroids along with the cell-cen-
!     tered solution at those centroids,  so the files are actually vertex-
!     centered.   Output many be in either point or block order, the latter
!     being more efficient for Tecplot.
!
!        A function file option has NOT been implemented for the original
!     cell-splitting options.
!
!     Control File Format (standard input):
!
!        surface_grid.g               ! Structured grid
!        none                         ! Matching flow solution
!        n                            ! Subdivide cells?
!        y                            ! Eliminate repeated (x,y,z)s?
!        unstructured-unique-surf.dat ! Unstructured output (Tecplottable)
!        y                            ! Point order? n = block order
!
!     Procedures:
!
!        triangulation_io  I/O utilities for unstructured data
!        xyzq_io package   I/O utilities for PLOT3D grid and function files
!
!     Performance:
!
!        A measure of performance on a single Intel core (2022) is this result:
!
!        Volume data with 1015040 cells and 942590 hex cells and 15 functions:
!        333711623 bytes Oct 22 00:10 6b-hex-unique.dat
!
!        ~20 minutes
!
!     History:
!
!        08/16/04  D.A.Saunders  Initial adaptation of EXTRACT_BLOCKS, to
!                                provide data for testing the tetrahedral mesh
!                                variant of the ADT search package.
!        09/25/2022   "     "    Option to preserve quad/hex cells.
!                                Duplicate vertices are NOT suppressed.
!        10/14/2022   "     "    75 years since Chuck Yeager/Bell X-1/Mach 1+.
!                                Debugging proved unusually difficult.
!                                Performance on a big voolume grid may be
!                                intolerable.  Example: a 6-block volume grid
!                                with 942590 points required about 17 minutes
!                                with no function data.
!        10/18-21/22  "     "    Handled optional function data, starting with
!                                introduction of triangulation_io.  This means
!                                working with xyz(3,npoints) rather than the
!                                original x,y,z vectors.
!        10/23/2022   "     "    An input control file can be used in lieu of
!                                interactive inputs.
!
!     Author:  David Saunders, ELORET, Inc. at NASA Ames Research Center, CA.
!              Now with AMA, Inc. at Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure        ! For one structured grid block
   use xyzq_io_module              ! For structured grid data
   use tri_header_structure        ! Data structures used by triangulation_io
   use tri_zone_structure
   use triangulation_io            ! Tecplot-type unstructured I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd = 5,        &
      luncrt = 6,        &
      lunf   = 1,        &
      lung   = 2,        &
      lunout = 3

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character (11), parameter :: &
      form = 'unformatted'

!  Variables:

   integer :: &
      i, j, k, i1, i2, i3, i4, i5, i6, i7, i8, ib, ie, ios, ipoint, ipt, &
      mi, mj, mk, nblocks, nelements, nf, ni, nij, nj, nk, npoints, npts, offset

   integer, allocatable, dimension (:,:) :: &
      conn

   logical :: &
      cell_centered, centers, formatted_in, formatted_out, no_repeats, &
      point_order, subdivide, vertices

   character (1) :: &
      answer

   character (132) :: &
      filename

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      xyz
   type (tri_header_type) :: &
      header_out
   type (tri_type), pointer, dimension (:) :: &
      zone

!  Execution:

!  Prompt for and open the input grid:

   write (luncrt, '(a)', advance='no') &
      'Enter the input Plot3D surface or volume grid: '
   ios = 1
   do while (ios /= 0)
      read (lunkbd, *) filename
      call determine_grid_form (filename, lung, formatted_in, ios)
      i1 = 1;  if (formatted_in) i1 = 3
      open (lung, file=filename, status='old', form=form(i1:), &
            iostat=ios)
   end do

!! call file_prompt (lung, -luncrt, 'input grid', 'old', true, &
!!                   formatted_in, ios)
!! if (ios /= 0) go to 99

!  Read the input grid header:

   call xyz_header_io (1, lung, formatted_in, nblocks, xyz, ios)
   if (ios /= 0) go to 99

!  A function file is optional:

   nf = 0
   cell_centered = .false.
   write (luncrt, '(a)', advance='no') 'Optional function file, or ''none'': '

   ios = 1
   do while (ios /= 0)
      read (lunkbd, *) filename
      if (filename(1:4) == 'none') exit

      open (lunf, file=filename, status='old', form=form(i1:), &
            iostat=ios)
      nf = 1  ! Unknown, but nonzero
   end do

!  Read the input function file header?

   if (nf /= 0) then
      call q_header_io (1, lunf, formatted_in, nblocks, nf, xyz, ios)
      if (ios /= 0) go to 99

      header_out%numf = nf

!     Check for files that do not have matching dimensions:

      vertices = true;  centers = true

      do ib = 1, nblocks
         mi = xyz(ib)%mi;  mj = xyz(ib)%mj;  mk = xyz(ib)%mk
         ni = xyz(ib)%ni;  nj = xyz(ib)%nj;  nk = xyz(ib)%nk

         vertices = mi == ni   .and. mj == nj   .and. mk == nk   .and. vertices
         centers  = mi == ni-1 .and. mj == nj-1 .and. mk == nk-1 .and. centers

         if (.not. vertices) then
            write (luncrt, '(/, a, i4, /, a, 3i5)') &
               ' Grid/flow field dimension mismatch.  Block:', ib, &
               ' Grid: ', ni, nj, nk, ' Flow: ', mi, mj, mk
            ios = 1
            go to 99
         end if

      end do

      cell_centered = centers  ! Must be vertex-centered
!!    if (centers) header_out%fileform = 2
   end if

   header_out%fileform = 1  ! Vertex-centered
   npoints   = 0
   nelements = 0

   do ib = 1, nblocks
      nk = xyz(ib)%nk
      npoints   = npoints   +  xyz(ib)%ni*xyz(ib)%nj*nk
      nelements = nelements + (xyz(ib)%ni - 1)*(xyz(ib)%nj - 1) * max (nk-1, 1)
   end do

   write (luncrt, '(/, a, 2i9)') &
      'Input # vertices and cells:', npoints, nelements

!  Break up the cells?

   nk = xyz(1)%nk
   if (nk == 1) then
      write (luncrt, '(a)', advance ='no') &
         'Triangulate the cells? [y|n]: '
   else
      write (luncrt, '(a)', advance ='no') &
         'Turn the hex cells into tetrahedra? [y|n]: '
   end if
   read  (lunkbd, *) answer
   subdivide = answer == 'y' .or. answer == 'Y'

   write (luncrt, '(a)', advance = 'no') &
      'Eliminate repeated (x,y,z)s? [y|n]: '
   read  (lunkbd, *) answer
   no_repeats = answer == 'y' .or. answer == 'Y'

!  Display the block dimensions to reassure the user:

   write (luncrt, '(/, a)') 'Block dimensions found:'
   write (luncrt, '(4(i6, 2X, 3i5))') &
      (ib, xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk, ib = 1, nblocks)

!  Prompt for and open the Tecplottable output file:

   formatted_out = true  ! Avoid Tecplot binaries; only Tecplot can read them

!! call file_prompt (lunout, -luncrt, 'output dataset', 'unknown', false, &
!!                   formatted_out, ios)
!! if (ios /= 0) go to 99              !Switch to triangulation_io utilities

   write (luncrt, '(a)', advance = 'no') &
      'Output dataset (formatted, Tecplottable): '
   read  (lunkbd, *) header_out%filename

   write (luncrt, '(a)', advance = 'no') &
      'Point order? n = block order [y|n]: '
   read  (lunkbd, *) answer
   point_order = answer == 'y' .or. answer == 'Y'

   header_out%formatted = formatted_out

   if (nk == 1) then  ! Surface
      header_out%nvertices = 4  ! Quads?
      if (subdivide) header_out%nvertices = 3  ! Triangles
   else  ! Volume data
      header_out%nvertices = 8  ! Hex cells?
      if (subdivide) header_out%nvertices = 4  ! Tetrahedra
   end if

   header_out%numf = nf
   header_out%nzones = 1
   header_out%datapacking = 0  ! Point order?
   if (.not. point_order) header_out%datapacking = 1  ! Block order

   header_out%combine_zones = false
   header_out%centroids_to_vertices = false
   header_out%title = 'Unstructured data in Fluent order'

   allocate (zone(1))
   zone(1)%zone_title    = ' '

   if (subdivide) then
      call break_cells ()
   else
      call preserve_cells ()
   end if

!  Save the unstructured surface or volume grid (Tecplotable):

   allocate (header_out%varname(3+nf))

   header_out%varname(1) = 'x'
   header_out%varname(2) = 'y'
   header_out%varname(3) = 'z'

   do i = 1, nf
      call numbered_name ('f', i, header_out%varname(3+i), j)
   end do

   call tri_write (lunout, header_out, zone, ios)
   if (ios /= 0) go to 99

   close (lunout)

99 continue

! *** stop ! Avoid system dependencies.

!  Internal procedures for program deconstruct:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine break_cells ()
!
!     Original decomposition of surface quad cells into two triangles or hex
!     volume cells into five tetrahedra, for all blocks of a structured grid.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: il

!     Execution:

      if (nk == 1) then ! Triangulate the surface grid

         nelements = nelements*2

         zone(1)%nelements    = nelements
         zone(1)%zone_type    = 'FETRIANGLE'
         zone(1)%element_type = 'TRIANGLE'
         zone(1)%nnodes       = npoints
         zone(1)%nelements    = nelements

         call tri_zone_allocate (header_out, zone(1), ios)
         if (ios /= 0) go to 99

         ipoint = 1  ! For packed x/y/z coordinates
         ipt    = 1  ! For packed function values

         do ib = 1, nblocks

            call xyz_allocate (xyz(ib), ios)
            if (ios /= 0) go to 99

            ni = xyz(ib)%ni
            nj = xyz(ib)%nj
            npts = ni*nj

            call xyz_block_io (1, lung, formatted_in, npts,          &
                               xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)
            if (ios /= 0) go to 99

            if (nf > 0) then

               call q_allocate (xyz(ib), nf, ios)
               if (ios /= 0) go to 99

               call q_block_io (1, lunf, formatted_in, nf, xyz(ib)%mi, &
                                xyz(ib)%mj, 1, xyz(ib)%q, ios)
               if (ios /= 0) go to 99

               do j = 1, nj
                  do i = 1, ni
                     zone(1)%f(:,ipt) = xyz(ib)%q(:,i,j,1)
                     ipt = ipt + 1
                  end do
               end do

               deallocate (xyz(ib)%q)

            end if

            call copy_3d_to_2d (ni, nj, 1, xyz(ib)%x, npoints, 1, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, 1, xyz(ib)%y, npoints, 2, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, 1, xyz(ib)%z, npoints, 3, ipoint, &
                                zone(1)%xyz)
            ipoint = ipoint + npts

            deallocate (xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, 2i5)') &
                  ' Trouble deallocating input grid block #', ib, ios
               go to 99
            end if

         end do

!!       allocate (conn(3,nelements))

         offset = 0 ! For current block
         ie = 0

         do ib = 1, nblocks
            ni = xyz(ib)%ni
            nj = xyz(ib)%nj
            do j = 1, nj - 1
               do i = 1, ni - 1
                  ie = ie + 1
                  i1 = i + (j - 1)*ni + offset ! (i,j)
                  i2 = i1 + 1                  ! (i+1,j)
                  i3 = i1 + ni                 ! (i,j+1)
                  zone(1)%conn(1,ie) = i1
                  zone(1)%conn(2,ie) = i2 
                  zone(1)%conn(3,ie) = i3 
                  ie = ie + 1
                  zone(1)%conn(1,ie) = i2
                  zone(1)%conn(2,ie) = i2 + ni         ! (i+1,j+1)
                  zone(1)%conn(3,ie) = i3
               end do
            end do
            offset = offset + ni*nj
         end do

         if (ie /= nelements) then
            write (luncrt, '(/, a, 2i9)') &          ! Should not happen
               ' ???:  ie, nelements =', ie, nelements
            go to 99
         end if

      else ! Generate tetrahedra from the hex cells

         nelements = nelements * 5

         zone(1)%nelements    = nelements
         zone(1)%zone_type    = 'FETETRAHEDRON'
         zone(1)%element_type = 'TETRAHEDRON'
         zone(1)%nnodes       = npoints
         zone(1)%nelements    = nelements

         call tri_zone_allocate (header_out, zone(1), ios)
         if (ios /= 0) go to 99

         ipoint = 1  ! For packed x/y/z/ coordiantes
         ipt    = 1  ! For packed function values

         do ib = 1, nblocks

            call xyz_allocate (xyz(ib), ios)
            if (ios /= 0) go to 99

            ni = xyz(ib)%ni
            nj = xyz(ib)%nj
            nk = xyz(ib)%nk
            npts = ni*nj*nk

            call xyz_block_io (1, lung, formatted_in, npts, &
                               xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)
            if (ios /= 0) go to 99

            if (nf > 0) then

               call q_allocate (xyz(ib), nf, ios)
               if (ios /= 0) go to 99

               call q_block_io (1, lunf, formatted_in, nf, xyz(ib)%mi, &
                                xyz(ib)%mj, 1, xyz(ib)%q, ios)
               if (ios /= 0) go to 99

               do k = 1, nk
                  do j = 1, nj
                     do i = 1, ni
                        zone(1)%f(:,ipt) = xyz(ib)%q(:,i,j,k)
                        ipt = ipt + 1
                     end do
                  end do
               end do

               deallocate (xyz(ib)%q)

            end if

!           Store packed coordinates for eliminating repeated points:

            call copy_3d_to_2d (ni, nj, nk, xyz(ib)%x, npoints, 1, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, nk, xyz(ib)%y, npoints, 2, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, nk, xyz(ib)%z, npoints, 3, ipoint, &
                                zone(1)%xyz)

            ipoint = ipoint + npts

            deallocate (xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, 2i5)') &
                  ' Trouble deallocating input grid block #', ib, ios
               go to 99
            end if

         end do

         offset = 0  ! For current block
         ie = 0

         do ib = 1, nblocks
            ni  = xyz(ib)%ni
            nj  = xyz(ib)%nj
            nk  = xyz(ib)%nk
            nij = ni*nj
            do k = 1, nk - 1
               do j = 1, nj - 1
                  do i = 1, ni - 1
                     ie = ie + 1
                     i1 = i + (j - 1)*ni + (k - 1)*nij + offset  ! (i,  j,  k)
                     i2 = i1 + 1                                 ! (i+1,j,  k)
                     i3 = i1 + ni                                ! (i,  j+1,k)
                     i4 = i1 + nij                               ! (i,  j,  k+1)
                     zone(1)%conn(1,ie) = i1           ! (i,j,k)
                     zone(1)%conn(2,ie) = i2           ! (i+1,j,k)
                     zone(1)%conn(3,ie) = i3           ! (i,j+1,k)
                     zone(1)%conn(4,ie) = i4           ! (i,j,k+1)
                     ie = ie + 1
                     zone(1)%conn(1,ie) = i3 + 1       ! (i+1,j+1,k)
                     zone(1)%conn(2,ie) = i3           ! (i,j+1,k)
                     zone(1)%conn(3,ie) = i2           ! (i+1,j,k)
                     zone(1)%conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                     ie = ie + 1
                     zone(1)%conn(1,ie) = i2 + nij     ! (i+1,j,k+1)
                     zone(1)%conn(2,ie) = i2           ! (i+1,j,k)
                     zone(1)%conn(3,ie) = i4           ! (i,j,k+1)
                     zone(1)%conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                     ie = ie + 1
                     zone(1)%conn(1,ie) = i3 + nij     ! (i,j+1,k+1)
                     zone(1)%conn(2,ie) = i4           ! (i,j,k+1)
                     zone(1)%conn(3,ie) = i3           ! (i,j+1,k)
                     zone(1)%conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                     ie = ie + 1
                     zone(1)%conn(1,ie) = i2           ! (i+1,j,k)
                     zone(1)%conn(2,ie) = i3           ! (i,j+1,k)
                     zone(1)%conn(3,ie) = i4           ! (i,j,k+1)
                     zone(1)%conn(4,ie) = i4 + 1 + ni  ! (i+1,j+1,k+1)
                  end do
               end do
            end do
            offset = offset + nij*nk
         end do

         if (ie /= nelements) then
            write (luncrt, '(/, a, 3i9)') &          ! Should not happen
               ' ???:  ie, nelements =', ie, nelements
            go to 99
         end if

      end if

      close (lunout)

  99  continue

      end subroutine break_cells

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine preserve_cells ()
!
!     Convert all cells of all structured surface or volume grid blocks to
!     unstructured form.  Duplicate edge point coordinates are suppressed.
!     Hex vertices should follow the Fluent convention: 1234 in circular
!     order on one face and 5678 on the opposite face for volume grids.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: il
      real, allocatable, dimension (:) :: x, y, z   ! Packed coordinates

!     Execution:

      if (nk == 1) then  ! Surface quads

         zone(1)%nelements    = nelements
         zone(1)%zone_type    = 'FEQUADRILATERAL'
         zone(1)%element_type = 'QUADRILATERAL'
         zone(1)%nnodes       = npoints
         zone(1)%nelements    = nelements

         call tri_zone_allocate (header_out, zone(1), ios)
         if (ios /= 0) go to 99

         ipoint = 1  ! For packed x/y/z coords.
         ipt    = 1  ! For packed function values

         do ib = 1, nblocks

            call xyz_allocate (xyz(ib), ios)
            if (ios /= 0) go to 99

            ni = xyz(ib)%ni
            nj = xyz(ib)%nj
            npts = ni*nj

            call xyz_block_io (1, lung, formatted_in, npts, &
                               xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)
            if (ios /= 0) go to 99

            if (ib == 1) then
               write (luncrt, '(/, a, /, (3es16.8))') &
                  'First three input x/y/zs:', &
                  (xyz(1)%x(i,1,1), xyz(1)%y(i,1,1), xyz(1)%z(i,1,1), i = 1, 3)
            else if (ib == nblocks) then
               write (luncrt, '(/, a, /, (3es16.8))') &
                  'Last three input x/y/zs:', &
                  (xyz(ib)%x(ni,j,1), xyz(ib)%y(ni,j,1), xyz(ib)%z(ni,j,1), &
                  j = nj-2, nj)
            end if

            if (nf > 0) then

               call q_allocate (xyz(ib), nf, ios)
               if (ios /= 0) go to 99

               call q_block_io (1, lunf, formatted_in, nf, xyz(ib)%mi, &
                                xyz(ib)%mj, 1, xyz(ib)%q, ios)
               if (ios /= 0) go to 99

               if (ib == 1) then
                  write (luncrt, '(/, a)')'First three input sets of functions:'
                  do i = 1, 3
                     write (luncrt, '(99es16.8)') xyz(1)%q(:,i,1,1)
                  end do
               else if (ib == nblocks) then
                  write (luncrt, '(/, a)') 'Last three input sets of functions:'
                  do j = nj - 2, nj
                     write (luncrt, '(99es16.8)') xyz(ib)%q(:,ni,j,1)
                  end do
               end if

               do j = 1, nj
                  do i = 1, ni
                     zone(1)%f(:,ipt) = xyz(ib)%q(:,i,j,1)
                     ipt = ipt + 1
                  end do
               end do

               deallocate (xyz(ib)%q)

            end if

!           Store packed coordinates for eliminating repeated points:

            call copy_3d_to_2d (ni, nj, 1, xyz(ib)%x, npoints, 1, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, 1, xyz(ib)%y, npoints, 2, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, 1, xyz(ib)%z, npoints, 3, ipoint, &
                                zone(1)%xyz)
            ipoint = ipoint + npts

            deallocate (xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, 2i5)') &
                  ' Trouble deallocating input grid block #', ib, ios
               go to 99
            end if

         end do

!        The following does not suppress common block edge points:

         allocate (zone(1)%conn(4,nelements))

         offset = 0 ! For current block
         ie = 0

         do ib = 1, nblocks
            ni = xyz(ib)%ni
            nj = xyz(ib)%nj
            do j = 1, nj - 1
               do i = 1, ni - 1
                  ie = ie + 1
                  i1 = i + (j - 1)*ni + offset ! (i,  j)
                  i2 = i1 + 1                  ! (i+1,j)
                  i3 = i2 + ni                 ! (i+1,j+1)
                  i4 = i3 - 1                  ! (i,  j+1)
                  zone(1)%conn(1,ie) = i1
                  zone(1)%conn(2,ie) = i2 
                  zone(1)%conn(3,ie) = i3 
                  zone(1)%conn(4,ie) = i4
               end do
            end do
            offset = offset + ni*nj
         end do

         if (ie /= nelements) then  ! Programming error
            write (luncrt, '(/, a, 2i9)') & 
               ' preserve_cells: ???:  ie, nelements =', ie, nelements
            go to 99
         end if

!        Eliminate repeated edge points in place?

         if (no_repeats .and. nblocks > 1) then
            
            call unique_xyz_list (npoints, zone(1)%xyz, nf, zone(1)%f, 4, &
                                  nelements, zone(1)%conn)
            zone(1)%nnodes = npoints
         end if

      else  ! Volume grid: Turn each hex cell into Fluent order:

         zone(1)%nelements    = nelements
         zone(1)%zone_type    = 'FEBRICK'
         zone(1)%element_type = 'BRICK'
         zone(1)%nnodes       = npoints
         zone(1)%nelements    = nelements

         call tri_zone_allocate (header_out, zone(1), ios)
         if (ios /= 0) go to 99

         ipoint = 1  ! For packed x/y/z coords.
         ipt    = 1  ! For packed function values

         do ib = 1, nblocks

            call xyz_allocate (xyz(ib), ios)
            if (ios /= 0) go to 99

            ni  = xyz(ib)%ni
            nj  = xyz(ib)%nj
            nk  = xyz(ib)%nk
            npts = ni*nj*nk

            call xyz_block_io (1, lung, formatted_in, npts, &
                               xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)
            if (ios /= 0) go to 99

            if (ib == 1) then
               write (luncrt, '(/, a, /, (3es16.8))') &
                  'First three input x/y/zs:', &
                  (xyz(1)%x(i,1,1), xyz(1)%y(i,1,1), xyz(1)%z(i,1,1), i = 1, 3)
            else if (ib == nblocks) then
               write (luncrt, '(/, a, /, (3es16.8))') &
                  'Last three input x/y/zs:', &
                  (xyz(ib)%x(ni,nj,k), xyz(ib)%y(ni,nj,k), xyz(ib)%z(ni,nj,k), &
                  k = nk-2, nk)
            end if

            if (nf > 0) then

               call q_allocate (xyz(ib), nf, ios)
               if (ios /= 0) go to 99

               call q_block_io (1, lunf, formatted_in, nf, xyz(ib)%mi, &
                                xyz(ib)%mj, xyz(ib)%mk, xyz(ib)%q, ios)
               if (ios /= 0) go to 99

               if (ib == 1) then
                  write (luncrt, '(/, a)')'First three input sets of functions:'
                  do i = 1, 3
                     write (luncrt, '(99es16.8)') xyz(1)%q(:,i,1,1)
                  end do
               else if (ib == nblocks) then
                  write (luncrt, '(/, a)') 'Last three input sets of functions:'
                  do k = nk - 2, nk
                     write (luncrt, '(99es16.8)') xyz(ib)%q(:,ni,nj,k)
                  end do
               end if

               do k = 1, nk
                  do j = 1, nj
                     do i = 1, ni
                        zone(1)%f(:,ipt) = xyz(ib)%q(:,i,j,k)
                        ipt = ipt + 1
                     end do
                  end do
               end do

               deallocate (xyz(ib)%q)

            end if

!           Store packed coordinates for eliminating repeated points:

            call copy_3d_to_2d (ni, nj, nk, xyz(ib)%x, npoints, 1, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, nk, xyz(ib)%y, npoints, 2, ipoint, &
                                zone(1)%xyz)
            call copy_3d_to_2d (ni, nj, nk, xyz(ib)%z, npoints, 3, ipoint, &
                                zone(1)%xyz)
            ipoint = ipoint + npts

            deallocate (xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, 2i5)') &
                  ' Trouble deallocating input grid block #', ib, ios
               go to 99
            end if

         end do

!        The following does not suppress common block boundary points:

         offset = 0 ! For current block
         ie = 0

         do ib = 1, nblocks
            ni  = xyz(ib)%ni
            nj  = xyz(ib)%nj
            nk  = xyz(ib)%nk
            nij = ni*nj
            do k = 1, nk - 1
               do j = 1, nj - 1
                  do i = 1, ni - 1
                     ie = ie + 1
                     i1 = i + (j - 1)*ni + (k - 1)*nij + offset  ! (i,  j,  k)
                     i2 = i1 + 1                                 ! (i+1,j  ,k)
                     i3 = i2 + ni                                ! (i+1,j+1,k)
                     i4 = i1 + ni                                ! (i  ,j+1,k)
                     i5 = i1 + nij                               ! (i  ,j,  k+1)
                     i6 = i2 + nij                               ! (i+1,j,  k+1)
                     i7 = i3 + nij                               ! (i+1,j+1,k+1)
                     i8 = i4 + nij                               ! (i,  j+1,k+1)
                     zone(1)%conn(1,ie) = i1
                     zone(1)%conn(2,ie) = i2
                     zone(1)%conn(3,ie) = i3
                     zone(1)%conn(4,ie) = i4
                     zone(1)%conn(5,ie) = i5
                     zone(1)%conn(6,ie) = i6
                     zone(1)%conn(7,ie) = i7
                     zone(1)%conn(8,ie) = i8
                  end do
               end do
            end do
            offset = offset + nij*nk
         end do

         if (ie /= nelements) then
            write (luncrt, '(/, a, 2i9)') &  ! This should never happen
               'preserve_cells: ???:  ie, nelements =', ie, nelements
            go to 99
         end if

!        Eliminate repeated edge points in place?

         if (no_repeats .and. nblocks > 1) then

            call unique_xyz_list (npoints, zone(1)%xyz, nf, zone(1)%f, 8, &
                                  nelements, zone(1)%conn)
            zone(1)%nnodes = npoints
         end if

      end if

  99  continue

      end subroutine preserve_cells

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine copy_3d_to_2d (ni, nj, nk, x3d, n, ixyz, i1, x2d)
!
!     Copy one coordinate per call.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: ni, nj, nk, n  ! 3D and 2D array dimensions
      real,    intent (in)  :: x3d(ni,nj,nk)  ! 3D array to copy
      integer, intent (in)  :: ixyz           ! 1 ==> x; 2 ==> y; 3 ==> z
      integer, intent (in)  :: i1             ! First index to copy to
      real,    intent (out) :: x2d(3,n)       ! 2D form of same x y or z coords.

!     Local variables:

      integer :: i, ic, j, k

!     Execution:

      ic = i1

      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               x2d(ixyz,ic) = x3d(i,j,k)
!!             write (77, '(i1, i8, es16.8)') ixyz, ic, x2d(ixyz,ic)
               ic = ic + 1
            end do
         end do
      end do

      end subroutine copy_3d_to_2d

   end program deconstruct
