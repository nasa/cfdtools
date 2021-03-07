!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program shadowgraph
!
!  Description:
!
!     For a real gas CFD solution on a structured or unstructured grid in three-
!  space, construct an image that simulates an experimental shadowgraph. Options
!  for interferograms and schlieren are also [to be] provided.
!
!     [Actually, only the mixture density is treated at present (ns = 1), and it
!  is unlikely that distinguishing the species densities (and their derivatives)
!  would produce perceptible differences in the images.  In the case of FUN3D,
!  where the derivatives must be included with the density, handling multiple
!  species could also involve much larger files for negligible gain.]
!
!     [Interferograms have not been implemented yet.  Simulated schlieren and
!  shadowgraphs can both be obtained via Tecplot gray-scale contouring of the
!  various functions written to quadratures.f[u] - see a sample layout.]
!
!  Control File (Standard Input):
!
!     SHADOWGRAPH Control File
!     1              ! Image type    1|2|3 = shadowgraph|schlieren|interferogram
!     1              ! Grid type 1|2|3|4|5 = P3D|P3D/iblnk|DPLR ovrst|FUN3D|VTK
!     3              ! Dimensions      2|3 = 2D/(x,z)|3D/(x,y,z)
!     mygrid.gu      ! Grid file name|FUN3D/US3D file name|@name -> FUN3D list
!     F              ! Formatted?      T|F = ASCII|unformatted
!     myflow.fu      ! Flow data file name (mixture density or blank, density)
!     F              ! Formatted?      T|F = ASCII|unformatted
!     F              ! ASCII results?  T|F = ASCII|unformatted
!     0.  1.  0.     ! Unit view vector;  0. 1. 0. => no rotation of CFD data
!     2000 2000      ! Image pixel counts
!     999. 999.      ! Image width, height; CFD units; 999. => grid bbox values
!     999. 999. 999. ! Image center (x,y,z) before any rotation; 999. => "   "
!     10.            ! Distance of shadowgraph image from the exiting rays
!     999999.        ! Clip image data further than this # std. devs. from mean
!
!  Notes:
!
!     o  Input flowfield datasets for structured grids need only one function,
!        density, except that DPLR overset cases have blanking information as
!        function 1 and density as function 2.  Partial derivatives are
!        computed by SHADOWGRAPH.
!     o  Unstructured datasets are expected to include density as function 1
!        and the partial x/y/z derivatives of density as functions 2, 3, 4.
!     o  US3D introduced a VTK (Virtual ToolKit) unstructured binary format that
!        allows for cell-centered and/or vertex-centered function data.  This
!        program handles one or the other, but not both at the same time.
!     o  US3D datasets are expected to contain x/y/z coordinates at grid points
!        and function data at the cell centers.  This program detects cell-
!        centered data and interpolates it to the vertices via inverse distance
!        weighted averaging.
!     o  Cart3D datasets as written at NASA Ames are expected to have both the
!        grid coordinates and the function data at cell centers, and this is no
!        different from having both at the vertices as far as SHADOWGRAPH is
!        concerned.
!     o  Both US3D and Cart3D can employ grid type 5 (VTK unstructured binary),
!        which (like the earlier-handled FUN3D format but with different node-
!        numbering conventions) allows for cells that are tetrahedra, hexahedra,
!        prisms, or pyramids.
!
!  History:
!
!     Jan-Mar 2010  D.A.Saunders  Initial implementation of some capabilities
!                                 published for program CISS by Leslie Yates in
!                                 AIAA Journal Vol. 31, No. 10, October 1993:
!                                 "Images Constructed from Computed Flowfields."
!     04/13/10      DAS           Option to read integration data from a prior
!                                 run rather than regenerate it before testing
!                                 changes to the image construction details.
!     04/15/10       "            After too much effort trying to get useful
!                                 images from actual deflections and accumulated
!                                 pixel counts, it appears that contour plots of
!                                 the integrated quantities (1/n)(dn/dx) and
!                                 (1/n)(dn/dz) in gray-scale provide schlieren
!                                 images corresponding to horizontal & vertical
!                                 knife edges respectively, while the magnitude
!                                 of this vector pair simulates a shadowgraph.
!     04/19/10       "            Choosing to do spline integrations along the
!                                 lines of sight comes at a high price: each
!                                 line requires scanning the entire flow field.
!                                 [A 2000 x 2000 case for a Shuttle solution
!                                 with 8,143,360 grid cells takes 4.5 hours.]
!                                 Therefore, an option to use the trapezoidal
!                                 rule (as in CISS) is also provided, allowing
!                                 each cell of the flow field to be visited
!                                 just once for the rays that intersect it.
!                                 The above case now takes 4.5 minutes, and the
!                                 differences are imperceptible everywhere
!                                 except along the top of the OMS pod. Puzzling.
!     04/28/10       "            PLOT3D/iblank option implemented.
!     05/05/10       "            A large OVERFLOW grid shows that the bounding
!                                 box of the grid can be way too big to define
!                                 the imaged region.  Retain that option, but
!                                 provide for confining the image to smaller
!                                 regions.  Also, Bil Kleb suggested using a
!                                 unit vector to allow any view direction.
!     06/30/10       "            The 04/13/10 option for reprocessing the
!                                 integrated data is retained but should not be
!                                 needed following the 04/15/10 realization.
!                                 Entering a non-defaulted image center is best
!                                 done in the input grid coordinate system.
!                                 This means entering three coordinates, even
!                                 though only the rotated x & z are the only
!                                 ones used to define the output image center.
!                                 Note that rotations are about a line through
!                                 the origin of the input grid coords.
!     07/09/10       "            Filled in the FUN3D option (FieldView data
!                                 format, 32-bit, unformatted).  The whole
!                                 program can safely be compiled for 32-bit
!                                 arithmetic.
!     07/28/10       "            PLOT3D input files are now processed one block
!                                 at a time.  One price paid is having to read
!                                 (and possibly rotate) the grid twice in order
!                                 to default the image size, but worse is the
!                                 impact on program modularity.  Now, several
!                                 steps are interleaved rather awkwardly and
!                                 more variables are moved to the highest level
!                                 for inheritance by the internal procedures.
!                                 The spline quadrature option had to be killed,
!                                 and so did the reprocessing option.
!     08/02/10       "            Wondering about the 5.*machine_epsilon by
!                                 Jan-Renee Carlson prompted use of x1 - eps
!                                 for the i2 test, to err on the safe side
!                                 as for i1; likewise for j2 (symmetric now).
!     08/03/10       "            Extended the FUN3D case to allow multiple
!                                 input files, one zone per file, via a list
!                                 of file names indicated by @<list file name>.
!     11/04/10   D. Saunders      Incorporated DPLR overset case, requiring that
!                ERC, Inc./ARC    flow function 1 = blank, function 2 = density,
!                                 and igrid_type = 3.
!     07/07/11       "            Added US3D/VTK unstructured grid option
!                                 (igrid_type = 5).  Cart3D data is expected to
!                                 use the same format as well.
!     07/12/11       "            VTK hexagon and prism cells needed translating
!                                 to FieldView's node-numbering convention.
!     07/15/11       "            If cell-centered VTK function data are found,
!                                 convert to cell vertices in-place.
!     07/19/11       "            Realized the intended (1/n) factor in the
!                                 integrated quantities was unintentionally
!                                 missing.  Results appear identical, though,
!                                 because n is very close to 1, and the partial
!                                 derivatives dominate.
!     07/20/11       "            The documentation above now covers files from
!                                 Cart3D solutions, which can employ the same
!                                 VTK format as US3D but can omit the conversion
!                                 of function data here from cell centers to
!                                 cell vertices because both the coordinates and
!                                 the functions are cell-centered (but can be
!                                 treated as though they are vertex-centered).
!     07/24/11       "            The VTK cases were incrementing the il pointer
!                                 at the end of each cell case, forgetting the
!                                 tests for cycling (skipping the cell) above.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use element_type_module    ! Derived data type for a FieldView cell type
   use FieldView_module       ! Reading routine for FUN3D's FieldView-type data
   use VTK_file_header_module ! Derived data type for a VTK file header
   use VTK_unstructured_grid_module  ! Derived type, VTK unstructured grid data
   use VTK_IO_module          ! I/O utilities for VTK legacy files
   use grid_block_structure   ! Derived data type for one PLOT3D grid block
   use grid_block_utilities   ! Includes one for setting block bounding boxes
   use xyzq_io_module         ! I/O package for PLOT3D files

   implicit none

!  Local constants:

   integer,   parameter :: maxns = 5                 ! Species limit for now
   integer,   parameter :: lungrid = 1, lunflow = 2  ! Input CFD solution
   integer,   parameter :: lunout  = 3, luninfo = 4  ! Image output; other info.
   integer,   parameter :: lunstdin = 5, lunlog = 6  ! Control inputs & run log
   integer,   parameter :: lunlist = 7               ! For list of FUN3D files
   real,      parameter :: default = 999.            ! To use grid bounding box
   logical,   parameter :: false = .false., true = .true.

!  Local variables:

   integer :: &
      igrid_type, image_type, ios, ipixels, jpixels,       &
      n, nblocks, ncells, ndim, nelem, ni, nj, nk, nnodes, npoints, ns, nvars

   real :: &
      angle, camera_distance, clip_factor, dx, dz, grid_bbox(4), height,       &
      Px(3), Qx(3), view_vector(3), width, x1, x2, xcenter, ycenter, zcenter,  &
      z1, z2

   logical :: &
      ASCII_results, formatted_grid, formatted_flow, iblanked, list_of_files,  &
      preprocess

   real, allocatable, dimension (:,:) :: &
      fimage, ximage, zimage               ! Image for gray-scale contouring

   real, allocatable, dimension (:,:) :: &
      deflected_i, deflected_j

   real, allocatable, dimension (:,:,:) :: &
      quadratures                          ! Ray integrals of n, 1/n dn/dx, ...

   real (kind=4), pointer, dimension (:) :: &
      x, y, z                              ! Node coordinates, FUN3D case

   real (kind=4), pointer :: &             ! Flow variables rho and x/y/z
      f(:,:)                               ! derivatives of rho, FUN3D case

   character (len=80) :: &
      nextfilename                         ! For the @listname case

   character (len=80), pointer :: &
      varname(:)                           ! Flow variable names, FUN3D case

   type (grid_type), pointer :: &          ! Grid + species densities + derived
      xyzq(:)                              ! functions, PLOT3D/[no]iblank cases

   type (element_type), pointer :: &       ! For attributes of each cell type
      element(:)                           ! from FUN3D, including node indices

   type (VTK_file_header_type) :: &        ! For any VTK file
      VTK_file_header

   type (unstructured_type) :: &           ! For US3D's VTK unstructured grid
      VTK_us_grid                          ! format

!  Execution:
!  ----------


   call read_controls ()          ! These are all local procedures below

   if (ios /= 0) go to 99         ! Single termination philosophy

   call scan_grid ()              ! Determine any rotation, and if the image
                                  ! size & location are defaulted, scan the grid
   if (ios /= 0) go to 99         ! for its data range, one zone at a time

   call image_grid ()             ! Set up 2-D incident ray/image-related arrays


   call process_cfd_solution ()   ! From the CFD solution, compute refractive-
                                  ! index-based functions associated with the
   if (ios /= 0) go to 99         ! incident rays

   call shadow_image ()           ! Derive shadowgraph-type image data from
                                  ! the refractive-index-based calculations

   call write_results ()          ! Save plottable results


99 continue  ! Avoid possible system "stop" quirks


!  Internal procedures for program shadowgraph in the order they're used:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_controls ()  ! And open input and output files.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local constants:

      character, parameter :: format*11 = 'unformatted'

!     Local variables:

      integer   :: i
      character :: filename*80, trouble*16

!     Execution:

      open (lunlog, file='shadowgraph.log', status='unknown')

!     Read the control file on standard input:

      read (lunstdin, '(a)', iostat=ios)  ! Skip line 1
      if (ios /= 0) then
         trouble = 'first';  go to 90
      end if

      read (lunstdin, *, iostat=ios) image_type
      if (ios /= 0) then
         trouble = 'image type';  go to 90
      end if

      read (lunstdin, *, iostat=ios) igrid_type
      if (ios /= 0) then
         trouble = 'grid type';  go to 90
      end if
      iblanked = igrid_type == 2

      read (lunstdin, *, iostat=ios) ndim
      if (ios /= 0) then
         trouble = '# dimensions';  go to 90
      end if

      if (igrid_type < 4) then  ! Separate grid and function files

         read (lunstdin, *, iostat=ios) filename
         if (ios /= 0) then
            trouble = 'grid file name';  go to 90
         end if

         read (lunstdin, *, iostat=ios) formatted_grid
         if (ios /= 0) then
            trouble = 'grid formatting';  go to 90
         end if
         i = 1;  if (formatted_grid) i = 3

         open (lungrid, file=filename, form=format(i:11), status='old',        &
               iostat=ios)
         if (ios /= 0) then
            trouble = 'grid file';  go to 95
         end if

         read (lunstdin, *, iostat=ios) filename
         if (ios /= 0) then
            trouble = 'flow file name';  go to 90
         end if

         read (lunstdin, *, iostat=ios) formatted_flow
         if (ios /= 0) then
            trouble = 'flow formatting';  go to 90
         end if
         i = 1;  if (formatted_flow) i = 3

         open (lunflow, file=filename, form=format(i:11), status='old',        &
               iostat=ios)
         if (ios /= 0) then
            trouble = 'flow file';  go to 95
         end if

      else  ! igrid_type = 4 (FUN3D) or 5 (VTK)

         ns = 1 ! Multiple species + their derivatives are not handled

         read (lunstdin, *, iostat=ios) filename
         if (ios /= 0) then
            trouble = 'FUN3D or VTK file name';  go to 90
         end if

         ! FUN3D datasets have one zone per file; allow for multiple files

         list_of_files = filename(1:1) == '@'

         if (.not. list_of_files) then

            nblocks = 1

            if (igrid_type == 4) then
               open (lungrid, file=filename, form=format(1:11), status='old',  &
                     iostat=ios)
               if (ios /= 0) then
                  trouble = 'FUN3D file';  go to 95
               end if
            else  ! igrid_type = 5
               open (lungrid, file=filename, form=format(1:11), status='old',  &
                     access='stream', iostat=ios)
               if (ios /= 0) then
                  trouble = 'VTK file';  go to 95
               end if
            end if

         else  ! Count the number of FUN3D zones listed, using nblocks

            open (lunlist, file=filename(2:), form=format(3:11), status='old', &
                  iostat=ios)
            if (ios /= 0) then
               trouble = 'FUN3D file list';  go to 95
            end if

            nblocks = 0
            do while (ios == 0)  ! Till EOF (ios < 0)
               read (lunlist, *, iostat=ios) nextfilename
               if (ios == 0) then  ! Trap any open errors early
                  nblocks = nblocks + 1
                  open (lungrid, file=nextfilename, form=format(1:11),         &
                        status='old', iostat=ios)
                  if (ios /= 0) then
                     write (lunlog, '(/, 2a)') &
                        'Trouble opening ', trim (nextfilename)
                     go to 90
                  end if
                  close (lungrid)
               end if
            end do

            rewind (lunlist)

            write (lunlog, '(/, a, i4)') '# FUN3D datasets specified:', nblocks

         end if

         do i = 1, 3
            read (lunstdin, *)
         end do

      end if

      read (lunstdin, *, iostat=ios) ASCII_results
      if (ios /= 0) then
         trouble = 'results format';  go to 90
      end if

      read (lunstdin, *, iostat=ios) view_vector
      if (ios /= 0) then
         trouble = 'unit view vector';  go to 90
      end if

      read (lunstdin, *, iostat=ios) ipixels, jpixels
      if (ios /= 0) then
         trouble = 'image dimensions';  go to 90
      end if

      read (lunstdin, *, iostat=ios) width, height
      if (ios /= 0) then
         trouble = 'image width and height';  go to 90
      end if

      read (lunstdin, *, iostat=ios) xcenter, ycenter, zcenter
      if (ios /= 0) then
         trouble = 'xcenter, ycenter, zcenter';  go to 90
      end if

      preprocess = xcenter == default  ! No need to read blocks twice if the
                                       ! image size and location are entered
      if (.not. preprocess) then
         if (width == default .or. height == default .or. &
            zcenter == default) then
            write (lunlog, '(/, 2a)') &
               'If xcenter is specified (/= 999.), width, height & zcenter ',  &
               'must be specified as well.'
            ios = 1
            go to 99
         end if
      end if

      read (lunstdin, *, iostat=ios) camera_distance
      if (ios /= 0) then
         trouble = 'camera distance';  go to 90
      end if

      read (lunstdin, *, iostat=ios) clip_factor
      if (ios /= 0) then
         trouble = 'image clipping factor';  go to 90
      end if

      go to 99

  90  continue
      write (lunlog, '(/, 3a)') '*** Error reading ', trim (trouble), ' line.'
      go to 99

  95  write (lunlog, '(/, 4a)') '*** Error opening ', trim (trouble), ': ',    &
                                trim (filename)
  99  return

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scan_grid ()

!     Determine any rotation of the grid coordinates from the look vector,
!     and if the image width/height/center are defaulted, scan the grid one
!     zone at a time, assigning bounding boxes (after any rotation) as we go,
!     else avoid the I/O till the processing phase proper.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use trigd
      implicit none

!     Local constants:

      real, parameter :: big = 1.e+32, one = 1., ninety = 90., zero = 0.

!     Local variables:

      integer :: i, ib, npoints_total
      real    :: length

!     Execution:

      length = sqrt (dot_product (view_vector, view_vector))  ! Normalize it

      view_vector(:) = view_vector(:) / length

      Px(1) = zero;  Px(2) = zero;  Px(3) = zero  ! P is taken to be the origin

      if (abs (view_vector(2)) == one) then       ! (0, 1, 0) => side view, DPLR

         angle = zero

      else if (abs (view_vector(3)) == one) then  ! (0, 0, 1) => top view

         angle = ninety;  Qx(1) = one;  Qx(2) = zero;  Qx(3) = zero

      else if (abs (view_vector(1)) == one) then  ! (1, 0, 0) => front view

         angle = ninety;  Qx(1) = zero;  Qx(2) = zero;  Qx(3) = one

      else  ! More general view: v.u gives angle; v x u gives unit vector => Q

         angle = acosd (view_vector(2))  ! Dot product of v with u = (0, 1, 0)

         Qx(1) = -view_vector(3);  Qx(2) = zero   ! From the cross product
         Qx(3) =  view_vector(1)

      end if

      if (.not. preprocess) then  ! No need to read the grid twice, but ...

         if (angle /= zero) then  ! Rotate the specified image center

            call rotate3d (1, xcenter, ycenter, zcenter, angle, &
                           Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))
         end if

         if (igrid_type /= 5) go to 99  ! For single zone, read it here

      end if


      grid_bbox(1:3:2) = big;  grid_bbox(2:4:2) = -big  ! x and z mins * maxes

!     Read the grid one block/zone at a time and set block bounding boxes
!     so that the image size and location can be defaulted:

      npoints_total = 0

      select case (igrid_type)

         case (1, 2, 3)  ! PLOT3D multiblock, possibly iblanked or DPLR overset

            call xyz_header_io (1, lungrid, formatted_grid, nblocks, xyzq, ios)
            if (ios /= 0) go to 99

            do ib = 1, nblocks

               if (.not. iblanked) then
                  call xyz_allocate (xyzq(ib), ios)
               else
                  call xyzi_allocate (xyzq(ib), ios)
               end if

               if (ios /= 0) go to 99

               npoints = xyzq(ib)%ni * xyzq(ib)%nj * xyzq(ib)%nk
               npoints_total = npoints_total + npoints

               if (.not. iblanked) then

                  call xyz_block_io (1, lungrid, formatted_grid, npoints,      &
                                     xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z, ios)
               else
                  call xyzi_block_io (1, lungrid, formatted_grid, npoints,     &
                                      xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z,      &
                                      xyzq(ib)%iblank, ios)
               end if

               if (ios /= 0) go to 99

               if (angle /= zero) then

                  call rotate3d (npoints, xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z,  &
                                angle, Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))
               end if

               call block_data_range (xyzq(ib))  ! In grid_block_utilities.f90

               grid_bbox(1) = min (grid_bbox(1), xyzq(ib)%xmin)
               grid_bbox(2) = max (grid_bbox(2), xyzq(ib)%xmax)
               grid_bbox(3) = min (grid_bbox(3), xyzq(ib)%zmin)
               grid_bbox(4) = max (grid_bbox(4), xyzq(ib)%zmax)

               deallocate (xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z)
               if (iblanked) deallocate (xyzq(ib)%iblank)

            end do

            rewind (lungrid)

         case (4)  ! FUN3D unstructured (32-bit unformatted FieldView format)

            do ib = 1, nblocks  ! One zone or block per file

               if (list_of_files) then

                  read (lunlist, *, iostat=ios) nextfilename

                  open (lungrid, file=nextfilename, form='unformatted',        &
                        status='old', iostat=ios)

                  if (ios /= 0) then
                     write (lunlog, '(/, 3a, i5)') &
                        'Unable to open ', trim (nextfilename), '.  ios:', ios
                     go to 99
                  end if

               end if

               call read_fieldview (lungrid, lunlog, nvars, nnodes, x, y, z,   &
                                    f, varname, nelem, element, ios)
               if (ios /= 0) go to 99

               if (angle /= zero) then

                  call rotate3d (nnodes, x, y, z, angle, &
                                 Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))

                  call rotate3d (nnodes, f(1,2), f(1,3), f(1,4), angle, & !df/dx
                                 Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3)) !etc.
               end if

               do i = 1, nnodes
                  grid_bbox(1) = min (x(i), grid_bbox(1))
                  grid_bbox(2) = max (x(i), grid_bbox(2))
                  grid_bbox(3) = min (z(i), grid_bbox(3))
                  grid_bbox(4) = max (z(i), grid_bbox(4))
               end do

               npoints_total = npoints_total + nnodes

               close (lungrid)

               if (list_of_files) then
                  if (nblocks > 1) then
                     deallocate (x, y, z, f, varname, element, stat=ios)
                     if (ios /= 0) then
                        write (lunlog, '(/, a, 2i6)') &
                           'Trouble deallocating FUN3D data. Zone, ios:', n, ios
                        go to 99
                     end if
                  end if
               end if

            end do  ! Next block or zone

            if (list_of_files) rewind (lunlist)

         case (5)  ! US3D VTK unstructured grid file (all cells in one file)

            VTK_us_grid%nscalar = 4

            call read_VTK_unstructured (lungrid, lunlog, VTK_file_header, &
                                        VTK_us_grid, ios)
            if (ios /= 0) go to 99

            close (lungrid)

!           Allow for either vertex- or cell-centered function data:

            if (VTK_us_grid%cell_centered) then  ! Convert fn. data in-place

               write (lunlog, '(/, a)') &
                  'Converting cell-centered function values to vertex values.'

               call VTK_centers_to_vertices (lungrid, lunlog, &
                                             VTK_file_header, VTK_us_grid, &
                                             ios)
               if (ios /= 0) go to 99
!!!            write (90, '(1p, 8e16.8)') VTK_us_grid%f
!!!            if (ios <= 0) stop

            end if

            npoints_total = VTK_us_grid%npoints

            if (angle /= zero) then

               call rotate3d (npoints_total, VTK_us_grid%x,           &
                              VTK_us_grid%y, VTK_us_grid%z, angle,    &
                              Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))

               call rotate3d (npoints_total, VTK_us_grid%f(1,2), &  ! Derivs.
                              VTK_us_grid%f(1,3), VTK_us_grid%f(1,4), angle, &
                              Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))
            end if

            do i = 1, npoints_total
               grid_bbox(1) = min (VTK_us_grid%x(i), grid_bbox(1))
               grid_bbox(2) = max (VTK_us_grid%x(i), grid_bbox(2))
               grid_bbox(3) = min (VTK_us_grid%z(i), grid_bbox(3))
               grid_bbox(4) = max (VTK_us_grid%z(i), grid_bbox(4))
            end do

         case default

            ! Unknown dataset type

      end select

      if (igrid_type <= 3) then  ! Display block bounding boxes for possible use

         write (lunlog, '(/, a)') &
            'Grid block dimensions and bounding box limits ([rotated] x & z):'
         write (lunlog, '(4i6, 3x, 1p, 2e17.8, 3x, 2e17.8)')  &
            (ib, xyzq(ib)%ni,   xyzq(ib)%nj,   xyzq(ib)%nk,   &
                 xyzq(ib)%xmin, xyzq(ib)%xmax, xyzq(ib)%zmin, &
                 xyzq(ib)%zmax, ib = 1, nblocks)
      end if

      write (lunlog, '(/, a, i11)') &
         'Total # grid points found:', npoints_total
      write (lunlog, '(/, a, /, 6x, 1p, 2e17.8, 3x, 2e17.8)') &
         'Grid bounding box in view direction, x & z:', grid_bbox(:)

 99   return

      end subroutine scan_grid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine image_grid ()
!
!     Set up a rectangular array of light ray locations.
!     These serve as the image pixels as well since deflections are small.
!     Also allocate storage for the various integrations along light rays,
!     and for the final image, and initialize both arrays.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local constants:

      real, parameter :: half = 0.5, zero = 0.

!     Local variables:

      integer :: i, j, ni, nj

!     Execution:

      allocate (ximage(ipixels,jpixels), zimage(ipixels,jpixels), &
                fimage(ipixels,jpixels), quadratures(3,ipixels,jpixels))
      allocate (deflected_i(ipixels,jpixels), deflected_j(ipixels,jpixels))

      fimage(:,:)      = zero;  quadratures(:,:,:) = zero
      deflected_i(:,:) = zero;  deflected_j(:,:)   = zero

!     Generate the image pixels:

      if (xcenter == default) then
         x1 = grid_bbox(1);  x2 = grid_bbox(2)
         z1 = grid_bbox(3);  z2 = grid_bbox(4)
      else
         x1 = xcenter - half * width;  z1 = zcenter - half * height
         x2 = xcenter + half * width;  z2 = zcenter + half * height
      end if

      dx = (x2 - x1) / real (ipixels - 1)
      dz = (z2 - z1) / real (jpixels - 1)

      do i = 1, ipixels
         ximage(i,:) = x1 + dx * real (i-1)
      end do

      do j = 1, jpixels
         zimage(:,j) = z1 + dz * real (j-1)
      end do

      end subroutine image_grid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine process_cfd_solution ()

!     Read the CFD data one zone at a time (unless it's not multizone data, in
!     which case it's already been read) and process the cells one at a time.
!     Input species densities are converted to n - 1 (n = refractive index) and
!     the functions of n needed for the specified image type.  These grid point
!     functions will be interpolated where each incident ray intersects a cell
!     face, and integrated along the ray.
!
!                     ns
!        n  =  1  +  SUM  kappa  rho         ! kappa  =  Gladstone-Dale constant
!                    i=1       i    i        !      i    for ith species
!
!                1 dn
!        de   =  - --  (partial derivatives) ! changes in angular deflections
!          x     n dx                        ! are integrated along rays to
!                                            ! calculate total deflections,
!                1 dn                        ! used for both schlieren images
!        de   =  - --                        ! and shadowgraphs
!          z     n dz
!
!                 2 pi
!        f(n)  =  ---- (n - n )              ! Integration of f(n) along rays
!                lambda      0               ! gives total phase shifts for
!                                            ! interferograms
!
!     Spline quadrature method (slow, and now gone; can't do 1 block at a time):
!
!        For each incident ray, check the entire flow volume (eliminating as
!        many grid blocks and grid cells as possible via stored bounding boxes);
!        calculate line segments of data, sort them in y, eliminate duplicate
!        end points, perform monotonic spline integrations.
!
!     Trapezoidal rule quadrature (fast):
!
!        For each volume grid cell, find the incident rays that intersect it and
!        increment their numerical quadratures via 2-point averages times the
!        distance between two intersection points.  Results are independent of
!        the order in which the cells are processed.
!
!     Updating pixels in the image plane according to the integrated deflections
!     [as suggested by CISS] appears doomed to a handful of hits at best per
!     pixel.  The integrated quantities themselves are measures of deflections,
!     so contouring them directly makes sense and appears to be all that is
!     needed.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     implicit none

!     Local constants:

      real, parameter :: zero  = 0.
      real, parameter :: GDair = 0.000292/1.293 ! m^3/kg
                                                ! Use table look-up for real gas
!     Local variables:

      integer :: i, ib, j, k, nf
      logical :: lxyz(3)

      real, allocatable :: qtemp(:,:,:)

!     Execution:

      lxyz(:) = true;  lxyz(2) = false  ! Suppress y derivs. (structured cases)

!     Read and process one block/zone at a time:

      select case (igrid_type)

         case (1, 2, 3)  ! PLOT3D multiblock, possibly iblanked, or DPLR overset

            call xyz_header_io (1, lungrid, formatted_grid, nblocks, xyzq, ios)
            if (ios /= 0) go to 99

            call q_header_io (1, lunflow,formatted_flow, nblocks, nf, xyzq, ios)
            if (ios /= 0) go to 99

            do ib = 1, nblocks

               ns = nf

               if (.not. iblanked) then
                  call xyz_allocate (xyzq(ib), ios)
               else  ! For OVERFLOW or DPLR overset
                  call xyzi_allocate (xyzq(ib), ios)
               end if

               if (ios /= 0) go to 99

               ni = xyzq(ib)%ni;  nj = xyzq(ib)%nj;  nk = xyzq(ib)%nk
               npoints = ni * nj * nk

               if (.not. iblanked) then
                  call xyz_block_io (1, lungrid, formatted_grid, npoints, &
                                     xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z, ios)
               else
                  call xyzi_block_io (1, lungrid, formatted_grid, npoints, &
                                     xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z,   &
                                     xyzq(ib)%iblank, ios)
               end if

               if (ios /= 0) go to 99

               if (angle /= zero) then

                  call rotate3d (npoints, xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z,  &
                                angle, Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))
               end if

               call q_allocate (xyzq(ib), ns, ios)
               if (ios /= 0) go to 99

               call q_block_io (1, lunflow, formatted_flow, ns, ni, nj, nk,    &
                                xyzq(ib)%q, ios)
               if (ios /= 0) go to 99

!              DPLR overset files are assumed to have function 1 = iblank.
!              Shuffle such data to make it look like a PLOT3D/iblank case:

               if (igrid_type == 3) then

                  allocate (xyzq(ib)%iblank(ni,nj,nk), qtemp(ni,nj,nk))

                  do k = 1, nk
                     do j = 1, nj
                        do i = 1, ni
                           xyzq(ib)%iblank(i,j,k) = nint (xyzq(ib)%q(1,i,j,k))
                           qtemp(i,j,k) = xyzq(ib)%q(2,i,j,k)
                        end do
                     end do
                  end do

                  deallocate (xyzq(ib)%q)

                  allocate (xyzq(ib)%q(1,ni,nj,nk))

                  xyzq(ib)%q(1,:,:,:) = qtemp(:,:,:)

                  deallocate (qtemp)

                  ns = ns - 1

               end if

!              Derive refractive index n from species densities.  Actually,
!              n - 1 is assigned since its derivatives match those of n.

               allocate (xyzq(ib)%r(ni,nj,nk))

               if (ns == 1) then
                  xyzq(ib)%r(:,:,:) = GDair * xyzq(ib)%q(1,:,:,:)  ! n - 1
               else  ! Multiple species case will probably never be completed
               end if

               deallocate (xyzq(ib)%q)

               allocate (xyzq(ib)%fx(ni,nj,nk), xyzq(ib)%fz(ni,nj,nk))

               allocate (xyzq(ib)%fy(1,1,1))  ! Avoid unassociated pointer

               call flow_gradients (ni, nj, nk, lxyz, xyzq(ib)%x,              &
                                    xyzq(ib)%y,  xyzq(ib)%z,  xyzq(ib)%r,      &
                                    xyzq(ib)%fx, xyzq(ib)%fy, xyzq(ib)%fz)

               call one_cell_all_rays (ib)  ! Trace rays through each cell of
                                            ! this block & accumulate integrals

               deallocate (xyzq(ib)%x,  xyzq(ib)%y,  xyzq(ib)%z,  xyzq(ib)%r,  &
                           xyzq(ib)%fx, xyzq(ib)%fy, xyzq(ib)%fz)

               if (iblanked .or. igrid_type == 3) deallocate (xyzq(ib)%iblank)

            end do  ! Next block

         case (4)  ! FUN3D unstructured

            do n = 1, nblocks  ! One zone or block per file

               if (list_of_files .and. nblocks > 1) then  ! Errors already chkd.

                  read (lunlist, *, iostat=ios) nextfilename

                  open (lungrid, file=nextfilename, form='unformatted',        &
                        status='old', iostat=ios)
               end if

               if (nblocks > 1 .or. .not. preprocess) then

                  call read_fieldview (lungrid, lunlog, nvars, nnodes,         &
                                       x, y, z, f, varname, nelem, element, ios)
                  if (ios /= 0) go to 99

                  if (angle /= zero) then

                     call rotate3d (nnodes, x, y, z, angle, &
                                    Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))

                     call rotate3d (nnodes, f(1,2), f(1,3), f(1,4), angle,     &
                                    Px(1), Px(2), Px(3), Qx(1), Qx(2), Qx(3))
                  end if

                  close (lungrid)

               end if

               f(:,:) =  GDair * f(:,:)  ! n - 1 overwrites density
                                         ! Density derivatives are also scaled

               call one_cell_all_rays (1)

               deallocate (x, y, z, f, element, varname)

            end do  ! Next zone if any

         case (5)  ! US3D/VTK unstructured grid format (single zone)

            VTK_us_grid%f(:,:) =  GDair * VTK_us_grid%f(:,:) ! n - 1 <-- density
                                           ! Density derivatives are also scaled

            call one_cell_all_rays (1)

            deallocate (VTK_us_grid%x, VTK_us_grid%y,    VTK_us_grid%z, &
                        VTK_us_grid%f, VTK_us_grid%list, VTK_us_grid%cell_type)

         case default

      end select

 99   return

      end subroutine process_cfd_solution

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine one_cell_all_rays (iblock)

!     For each cell of the indicated block, determine which incident rays may
!     pass through it, and look for at most 2 line/face intersections.
!     If the numerical integrations are performed with the trapezoidal rule,
!     the outer loop can be over grid cells, not incident rays.  Simply sum
!     the contributions to any ray that encounters the current cell, which
!     needs to be visited only once.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) :: iblock  ! Grid block/zone number

!     Local constants:

      real, parameter :: half = 0.5, one = 1.

!     Local variables:

      integer :: i, i1, i2, ib, ielem, il, in1, in2, in3, in4, in5, in6, it
      integer :: j, j1, j2, jl, jt, k, ncell, npts
      integer :: list(8)
      real    :: dy, fpts(3,2), ypts(2)
      real, dimension (8) :: xcell, ycell, zcell, fcell, fxcell, fzcell

      real, save :: dxr, dzr, eps, x1meps, x1peps, xcmax, xcmin,               &
                    z1meps, z1peps, zcmax, zcmin

!     Execution:

      ib = iblock  ! Avoid repeated references to the argument

      if (ib == 1) then
         eps    = 5.*(x2 - x1)*epsilon (eps)
         x1peps = x1 + eps                 ! Guard against i1 = 2 when x = x1
         x1meps = x1 - eps                 ! Likewise, i2 /= ipixels + 1 at xmax
         z1peps = z1 + eps                 ! We err on the side of including a
         z1meps = z1 - eps                 ! ray that barely misses a cell, not
                                           ! omitting a ray that barely hits it
         dxr = one / dx
         dzr = one / dz

         write (lunlog, '(a)')
      end if

      select case (igrid_type)

         case (1)  ! Plain PLOT3D multiblock

            ni = xyzq(ib)%ni
            nj = xyzq(ib)%nj
            nk = xyzq(ib)%nk

            write (lunlog, '(a, i5, 3x, 3i5)') &
               'Grid block being processed:', ib, ni, nj, nk

            do k = 1, nk - 1
               do j = 1, nj - 1
                  do i = 1, ni - 1

                     xcmin = min &
                            (xyzq(ib)%x(i,j,k),     xyzq(ib)%x(i+1,j,k),       &
                             xyzq(ib)%x(i,j+1,k),   xyzq(ib)%x(i+1,j+1,k),     &
                             xyzq(ib)%x(i,j,k+1),   xyzq(ib)%x(i+1,j,k+1),     &
                             xyzq(ib)%x(i,j+1,k+1), xyzq(ib)%x(i+1,j+1,k+1))
                     if (xcmin > x2) cycle ! Cell is right of image subregion
                     xcmax = max &
                            (xyzq(ib)%x(i,j,k),     xyzq(ib)%x(i+1,j,k),       &
                             xyzq(ib)%x(i,j+1,k),   xyzq(ib)%x(i+1,j+1,k),     &
                             xyzq(ib)%x(i,j,k+1),   xyzq(ib)%x(i+1,j,k+1),     &
                             xyzq(ib)%x(i,j+1,k+1), xyzq(ib)%x(i+1,j+1,k+1))
                     if (xcmax < x1) cycle ! Cell is left of image subregion
                     zcmin = min &
                            (xyzq(ib)%z(i,j,k),     xyzq(ib)%z(i+1,j,k),       &
                             xyzq(ib)%z(i,j+1,k),   xyzq(ib)%z(i+1,j+1,k),     &
                             xyzq(ib)%z(i,j,k+1),   xyzq(ib)%z(i+1,j,k+1),     &
                             xyzq(ib)%z(i,j+1,k+1), xyzq(ib)%z(i+1,j+1,k+1))
                     if (zcmin > z2) cycle ! Cell is above image subregion
                     zcmax = max &
                            (xyzq(ib)%z(i,j,k),     xyzq(ib)%z(i+1,j,k),       &
                             xyzq(ib)%z(i,j+1,k),   xyzq(ib)%z(i+1,j+1,k),     &
                             xyzq(ib)%z(i,j,k+1),   xyzq(ib)%z(i+1,j,k+1),     &
                             xyzq(ib)%z(i,j+1,k+1), xyzq(ib)%z(i+1,j+1,k+1))
                     if (zcmax < z1) cycle ! Cell is below image subregion

!                    Pixel range that may encounter this cell (which may
!                    be only partially inside the image subregion):

                     i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                     i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                     j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                     j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                  if (ib == 1) then
!!!                     write (50, '(3i4, 1p, 4e15.7, 4i4)') &
!!!                        i,j,k,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2
!!!                  else
!!!                     stop
!!!                  end if

                     do jt = j1, j2
                        do it = i1, i2

                           call process_hex_cell (ni, nj, nk, i, j, k,         &
                                         xyzq(ib)%x, xyzq(ib)%y,  xyzq(ib)%z,  &
                                         xyzq(ib)%r, xyzq(ib)%fx, xyzq(ib)%fz, &
                                         ximage(it,jt), zimage(it,jt),         &
                                         npts, ypts, fpts)

                           if (npts == 2) then  ! Trapezoid rule increment

!                             The first quadrature here is of n - 1, but it
!                             should be that of n - n0 for interferograms.

                              dy = half * abs ((ypts(2) - ypts(1)))
                              quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                     dy*(fpts(1,1) + fpts(1,2))
                              fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                              fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                              quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                     dy*(fpts(2,1) + fpts(2,2))
                              quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                     dy*(fpts(3,1) + fpts(3,2))
                           end if  ! Else npts = 0

                        end do
                     end do

                  end do  ! Next block i
               end do  ! Next block j
            end do  ! Next block k

         case (2, 3)  ! PLOT3D/iblanked or DPLR overset;
                      ! for now, just skip if a vertex is blank

            ni = xyzq(ib)%ni
            nj = xyzq(ib)%nj
            nk = xyzq(ib)%nk

            write (lunlog, '(a, i5, 3x, 3i5)') &
               'Grid block being processed:', ib, ni, nj, nk

            do k = 1, nk - 1
               do j = 1, nj - 1
                  do i = 1, ni - 1

                     if (xyzq(ib)%iblank(i,j,k)       <= 0) cycle
                     if (xyzq(ib)%iblank(i+1,j,k)     <= 0) cycle
                     if (xyzq(ib)%iblank(i,j+1,k)     <= 0) cycle
                     if (xyzq(ib)%iblank(i+1,j+1,k)   <= 0) cycle
                     if (xyzq(ib)%iblank(i,j,k+1)     <= 0) cycle
                     if (xyzq(ib)%iblank(i+1,j,k+1)   <= 0) cycle
                     if (xyzq(ib)%iblank(i,j+1,k+1)   <= 0) cycle
                     if (xyzq(ib)%iblank(i+1,j+1,k+1) <= 0) cycle

                     xcmin = min &
                            (xyzq(ib)%x(i,j,k),     xyzq(ib)%x(i+1,j,k),       &
                             xyzq(ib)%x(i,j+1,k),   xyzq(ib)%x(i+1,j+1,k),     &
                             xyzq(ib)%x(i,j,k+1),   xyzq(ib)%x(i+1,j,k+1),     &
                             xyzq(ib)%x(i,j+1,k+1), xyzq(ib)%x(i+1,j+1,k+1))
                     if (xcmin > x2) cycle ! Cell is right of image subregion
                     xcmax = max &
                            (xyzq(ib)%x(i,j,k),     xyzq(ib)%x(i+1,j,k),       &
                             xyzq(ib)%x(i,j+1,k),   xyzq(ib)%x(i+1,j+1,k),     &
                             xyzq(ib)%x(i,j,k+1),   xyzq(ib)%x(i+1,j,k+1),     &
                             xyzq(ib)%x(i,j+1,k+1), xyzq(ib)%x(i+1,j+1,k+1))
                     if (xcmax < x1) cycle ! Cell is left of image subregion
                     zcmin = min &
                            (xyzq(ib)%z(i,j,k),     xyzq(ib)%z(i+1,j,k),       &
                             xyzq(ib)%z(i,j+1,k),   xyzq(ib)%z(i+1,j+1,k),     &
                             xyzq(ib)%z(i,j,k+1),   xyzq(ib)%z(i+1,j,k+1),     &
                             xyzq(ib)%z(i,j+1,k+1), xyzq(ib)%z(i+1,j+1,k+1))
                     if (zcmin > z2) cycle ! Cell is above image subregion
                     zcmax = max &
                            (xyzq(ib)%z(i,j,k),     xyzq(ib)%z(i+1,j,k),       &
                             xyzq(ib)%z(i,j+1,k),   xyzq(ib)%z(i+1,j+1,k),     &
                             xyzq(ib)%z(i,j,k+1),   xyzq(ib)%z(i+1,j,k+1),     &
                             xyzq(ib)%z(i,j+1,k+1), xyzq(ib)%z(i+1,j+1,k+1))
                     if (zcmax < z1) cycle ! Cell is below image subregion

!                    Pixel range that may encounter this cell (which may
!                    be only partially inside the image subregion):

                     i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                     i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                     j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                     j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                  if (ib == 47) then
!!!                     write (50, '(3i4, 1p, 4e15.7, 4i4)') &
!!!                        i,j,k,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2
!!!                  end if

                     do jt = j1, j2
                        do it = i1, i2

                           call process_hex_cell (ni, nj, nk, i, j, k,         &
                                         xyzq(ib)%x, xyzq(ib)%y,  xyzq(ib)%z,  &
                                         xyzq(ib)%r, xyzq(ib)%fx, xyzq(ib)%fz, &
                                         ximage(it,jt), zimage(it,jt),         &
                                         npts, ypts, fpts)

                           if (npts == 2) then  ! Trapezoid rule increment

!                             The first quadrature here is of n - 1, but it
!                             should be that of n - n0 for interferograms.

                              dy = half * abs ((ypts(2) - ypts(1)))
                              quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                     dy*(fpts(1,1) + fpts(1,2))
                              fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                              fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                              quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                     dy*(fpts(2,1) + fpts(2,2))
                              quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                     dy*(fpts(3,1) + fpts(3,2))
                           end if  ! Else npts = 0

                        end do
                     end do

                  end do  ! Next block i
               end do  ! Next block j
            end do  ! Next block k

         case (4)  ! FUN3D unstructured

            do ielem = 1, nelem

               ncell = element(ielem)%ncell

               select case (element(ielem)%cell_type)

                  case ('tet')

                     do i = 1, ncell

                        in1 = element(ielem)%c2n(1,i)
                        in2 = element(ielem)%c2n(2,i)
                        in3 = element(ielem)%c2n(3,i)
                        in4 = element(ielem)%c2n(4,i)

                        xcmin = min (x(in1), x(in2), x(in3), x(in4))
                        if (xcmin > x2) cycle ! Cell is right of image subregion

                        xcmax = max (x(in1), x(in2), x(in3), x(in4))
                        if (xcmax < x1) cycle ! Cell is left of image subregion

                        zcmin = min (z(in1), z(in2), z(in3), z(in4))
                        if (zcmin > z2) cycle ! Cell is above image subregion

                        zcmax = max (z(in1), z(in2), z(in3), z(in4))
                        if (zcmax < z1) cycle ! Cell is below image subregion

!                       Pixel range that may encounter this cell (which may
!                       be only partially inside the image subregion):

                        i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                        i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                        j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                        j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                     write (50, '(i4, 1p, 4e15.7, 4i4)') &
!!!                        i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                        do jt = j1, j2
                           do it = i1, i2

                              call process_tet_cell (nnodes, &
                                                     element(ielem)%c2n(:,i),  &
                                                     x, y, z,                  &
                                                     f(:,1), f(:,2), f(:,4),   &
                                                     ximage(it,jt),            &
                                                     zimage(it,jt),            &
                                                     npts, ypts, fpts)

                              if (npts == 2) then  ! Trapezoid rule increment

!                                The first quadrature here is of n - 1, but it
!                                should be that of n - n0 for interferograms.

                                 dy = half * abs ((ypts(2) - ypts(1)))
                                 quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                      dy*(fpts(1,1) + fpts(1,2))
                                 fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                                 fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                                 quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                      dy*(fpts(2,1) + fpts(2,2))
                                 quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                      dy*(fpts(3,1) + fpts(3,2))
                              end if  ! Else npts = 0

                           end do
                        end do

                        if (mod (i, 100000) == 0 .or. i == ncell) then
                           write (lunlog, '(a, i10, a, i10)') &
                              '   Finished tet cell', i, ' of', ncell
                        end if

                     end do  ! Next tet cell

                  case ('hex')

                     do i = 1, ncell

!                       Repacking allows the structured-grid routine to be used:

                        xcell(:) = x(element(ielem)%c2n(:,i))

                        xcmin = minval (xcell)
                        if (xcmin > x2) cycle ! Cell is right of image subregion

                        xcmax = maxval (xcell)
                        if (xcmax < x1) cycle ! Cell is left of image subregion

                        zcell(:) = z(element(ielem)%c2n(:,i))

                        zcmin = minval (zcell)
                        if (zcmin > z2) cycle ! Cell is above image subregion

                        zcmax = maxval (zcell)
                        if (zcmax < z1) cycle ! Cell is below image subregion

!                       Pixel range that may encounter this cell (which may
!                       be only partially inside the image subregion):

                        i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                        i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                        j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                        j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                     write (50, '(i7, 1p, 4e15.7, 4i6)') &
!!!                        i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                        ycell(:)  = y(element(ielem)%c2n(:,i))
                        fcell(:)  = f(element(ielem)%c2n(:,i),1)
!!!                     write (51, '(i7, 1p, 8e15.7)') i,fcell

                        fxcell(:) = f(element(ielem)%c2n(:,i),2)
!!!                     write (52, '(i7, 1p, 8e15.7)') i,fxcell

                        fzcell(:) = f(element(ielem)%c2n(:,i),4)
!!!                     write (54, '(i7, 1p, 8e15.7)') i,fzcell

                        do jt = j1, j2
                           do it = i1, i2

                              call process_hex_cell (2, 2, 2, 1, 1, 1,         &
                                                     xcell,  ycell,  zcell,    &
                                                     fcell, fxcell, fzcell,    &
                                                     ximage(it,jt),            &
                                                     zimage(it,jt),            &
                                                     npts, ypts, fpts)

                              if (npts == 2) then  ! Trapezoid rule increment

!                                The first quadrature here is of n - 1, but it
!                                should be that of n - n0 for interferograms.

                                 dy = half * abs ((ypts(2) - ypts(1)))
                                 quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                      dy*(fpts(1,1) + fpts(1,2))
                                 fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                                 fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                                 quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                      dy*(fpts(2,1) + fpts(2,2))
                                 quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                      dy*(fpts(3,1) + fpts(3,2))
                              end if  ! Else npts = 0

                           end do
                        end do

                        if (mod (i, 100000) == 0 .or. i == ncell) then
                           write (lunlog, '(a, i10, a, i10)') &
                              '   Finished hex cell', i, ' of', ncell
                        end if

                     end do  ! Next hexagonal cell

                  case ('prz')

                     do i = 1, ncell

                        in1 = element(ielem)%c2n(1,i)
                        in2 = element(ielem)%c2n(2,i)
                        in3 = element(ielem)%c2n(3,i)
                        in4 = element(ielem)%c2n(4,i)
                        in5 = element(ielem)%c2n(5,i)
                        in6 = element(ielem)%c2n(6,i)

                        xcmin = min (x(in1),x(in2),x(in3),x(in4),x(in5),x(in6))
                        if (xcmin > x2) cycle ! Cell is right of image subregion

                        xcmax = max (x(in1),x(in2),x(in3),x(in4),x(in5),x(in6))
                        if (xcmax < x1) cycle ! Cell is left of image subregion

                        zcmin = min (z(in1),z(in2),z(in3),z(in4),z(in5),z(in6))
                        if (zcmin > z2) cycle ! Cell is above image subregion

                        zcmax = max (z(in1),z(in2),z(in3),z(in4),z(in5),z(in6))
                        if (zcmax < z1) cycle ! Cell is below image subregion

!                       Pixel range that may encounter this cell (which may
!                       be only partially inside the image subregion):

                        i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                        i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                        j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                        j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                     write (50, '(i4, 1p, 4e15.7, 4i4)') &
!!!                        i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                        do jt = j1, j2
                           do it = i1, i2

                              call process_prism (nnodes, &
                                                  element(ielem)%c2n(:,i),     &
                                                  x, y, z,                     &
                                                  f(:,1), f(:,2), f(:,4),      &
                                                  ximage(it,jt),               &
                                                  zimage(it,jt),               &
                                                  npts, ypts, fpts)

                              if (npts == 2) then  ! Trapezoid rule increment

!                                The first quadrature here is of n - 1, but it
!                                should be that of n - n0 for interferograms.

                                 dy = half * abs ((ypts(2) - ypts(1)))
                                 quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                      dy*(fpts(1,1) + fpts(1,2))
                                 fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                                 fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                                 quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                      dy*(fpts(2,1) + fpts(2,2))
                                 quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                      dy*(fpts(3,1) + fpts(3,2))
                              end if  ! Else npts = 0

                           end do
                        end do

                        if (mod (i, 100000) == 0 .or. i == ncell) then
                           write (lunlog, '(a, i10, a, i10)') &
                              '   Finished prism cell', i, ' of', ncell
                        end if

                     end do  ! Next prism cell

                  case ('pyr')

                     do i = 1, ncell

                        in1 = element(ielem)%c2n(1,i)
                        in2 = element(ielem)%c2n(2,i)
                        in3 = element(ielem)%c2n(3,i)
                        in4 = element(ielem)%c2n(4,i)
                        in5 = element(ielem)%c2n(5,i)

                        xcmin = min (x(in1), x(in2), x(in3), x(in4), x(in5))
                        if (xcmin > x2) cycle ! Cell is right of image subregion

                        xcmax = max (x(in1), x(in2), x(in3), x(in4), x(in5))
                        if (xcmax < x1) cycle ! Cell is left of image subregion

                        zcmin = min (z(in1), z(in2), z(in3), z(in4), z(in5))
                        if (zcmin > z2) cycle ! Cell is above image subregion

                        zcmax = max (z(in1), z(in2), z(in3), z(in4), z(in5))
                        if (zcmax < z1) cycle ! Cell is below image subregion

!                       Pixel range that may encounter this cell (which may
!                       be only partially inside the image subregion):

                        i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                        i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                        j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                        j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                     write (50, '(i4, 1p, 4e15.7, 4i4)') &
!!!                        i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                        do jt = j1, j2
                           do it = i1, i2

                              call process_pyramid (nnodes, &
                                                    element(ielem)%c2n(:,i),   &
                                                    x, y, z,                   &
                                                    f(:,1), f(:,2), f(:,4),    &
                                                    ximage(it,jt),             &
                                                    zimage(it,jt),             &
                                                    npts, ypts, fpts)

                              if (npts == 2) then  ! Trapezoid rule increment

!                                The first quadrature here is of n - 1, but it
!                                should be that of n - n0 for interferograms.

                                 dy = half * abs ((ypts(2) - ypts(1)))
                                 quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                      dy*(fpts(1,1) + fpts(1,2))
                                 fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                                 fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                                 quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                      dy*(fpts(2,1) + fpts(2,2))
                                 quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                      dy*(fpts(3,1) + fpts(3,2))
                              end if  ! Else npts = 0

                           end do
                        end do

                        if (mod (i, 100000) == 0 .or. i == ncell) then
                           write (lunlog, '(a, i10, a, i10)') &
                              '   Finished pyramid cell', i, ' of', ncell
                        end if

                     end do  ! Next pyramid cell

                  case default

               end select  ! Cell type

               deallocate (element(ielem)%c2n)

            end do  ! Next FUN3D cell

         case (5)  ! US3D/VTK unstructured grid format

            npoints = VTK_us_grid%npoints
            ncell   = VTK_us_grid%ncells
            il      = 1                ! Index into packed list of point indices

            do i = 1, ncell

               select case (VTK_us_grid%cell_type(i))

                  case (10)  ! VTK_TETRA

                     in1 = VTK_us_grid%list(il+1)
                     in2 = VTK_us_grid%list(il+2)
                     in3 = VTK_us_grid%list(il+3)
                     in4 = VTK_us_grid%list(il+4)
                     jl  = il ! In case we cycle
                     il  = il + 5

                     xcmin = min (VTK_us_grid%x(in1), VTK_us_grid%x(in2), &
                                  VTK_us_grid%x(in3), VTK_us_grid%x(in4))
                     if (xcmin > x2) cycle ! Cell is right of image subregion

                     xcmax = max (VTK_us_grid%x(in1), VTK_us_grid%x(in2), &
                                  VTK_us_grid%x(in3), VTK_us_grid%x(in4))
                     if (xcmax < x1) cycle ! Cell is left of image subregion

                     zcmin = min (VTK_us_grid%z(in1), VTK_us_grid%z(in2), &
                                  VTK_us_grid%z(in3), VTK_us_grid%z(in4))
                     if (zcmin > z2) cycle ! Cell is above image subregion

                     zcmax = max (VTK_us_grid%z(in1), VTK_us_grid%z(in2), &
                                  VTK_us_grid%z(in3), VTK_us_grid%z(in4))
                     if (zcmax < z1) cycle ! Cell is below image subregion

!                    Pixel range that may encounter this cell (which may
!                    be only partially inside the image subregion):

                     i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                     i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                     j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                     j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                  write (50, '(i4, 1p, 4e15.7, 4i4)') &
!!!                     i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                     do jt = j1, j2
                        do it = i1, i2

                           call process_tet_cell (npoints,                     &
                                                  VTK_us_grid%list(jl+1:jl+4), &
                                                  VTK_us_grid%x,               &
                                                  VTK_us_grid%y,               &
                                                  VTK_us_grid%z,               &
                                                  VTK_us_grid%f(:,1),          &
                                                  VTK_us_grid%f(:,2),          &
                                                  VTK_us_grid%f(:,4),          &
                                                  ximage(it,jt),               &
                                                  zimage(it,jt),               &
                                                  npts, ypts, fpts)

                           if (npts == 2) then  ! Trapezoid rule increment

!                             The first quadrature here is of n - 1, but it
!                             should be that of n - n0 for interferograms.

                              dy = half * abs ((ypts(2) - ypts(1)))
                              quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                     dy*(fpts(1,1) + fpts(1,2))
                              fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                              fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                              quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                     dy*(fpts(2,1) + fpts(2,2))
                              quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                     dy*(fpts(3,1) + fpts(3,2))
                           end if  ! Else npts = 0

                        end do
                     end do

                  case (12)  ! VTK_HEXAHEDRON

!                    Repacking allows the structured-grid routine to be used,
!                    but the node numbering convention isn't right yet:

                     list(1) = VTK_us_grid%list(il+1)
                     list(2) = VTK_us_grid%list(il+2)
                     list(4) = VTK_us_grid%list(il+3)
                     list(3) = VTK_us_grid%list(il+4)
                     list(5) = VTK_us_grid%list(il+5)
                     list(6) = VTK_us_grid%list(il+6)
                     list(8) = VTK_us_grid%list(il+7)
                     list(7) = VTK_us_grid%list(il+8)
                     il      = il + 9

                     xcell(:) = VTK_us_grid%x(list(:))

                     xcmin = minval (xcell)
                     if (xcmin > x2) cycle ! Cell is right of image subregion

                     xcmax = maxval (xcell)
                     if (xcmax < x1) cycle ! Cell is left of image subregion

                     zcell(:) = VTK_us_grid%z(list(:))

                     zcmin = minval (zcell)
                     if (zcmin > z2) cycle ! Cell is above image subregion

                     zcmax = maxval (zcell)
                     if (zcmax < z1) cycle ! Cell is below image subregion

!                    Pixel range that may encounter this cell (which may
!                    be only partially inside the image subregion):

                     i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                     i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                     j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                     j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                  write (50, '(i7, 0p, 4f6.2, 1p, 6e15.7, 4i6)') &
!!!                    i,x1,x2,z1,z2,dxr,dzr,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                     ycell(:)  = VTK_us_grid%y(list(:))
                     fcell(:)  = VTK_us_grid%f(list(:),1)
!!!                  write (51, '(i7, 1p, 8e15.7)') i,fcell

                     fxcell(:) = VTK_us_grid%f(list(:),2)
!!!                  write (52, '(i7, 1p, 8e15.7)') i,fxcell

                     fzcell(:) = VTK_us_grid%f(list(:),4)
!!!                  write (54, '(i7, 1p, 8e15.7)') i,fzcell

                     do jt = j1, j2
                        do it = i1, i2

                           call process_hex_cell (2, 2, 2, 1, 1, 1,         &
                                                  xcell,  ycell,  zcell,    &
                                                  fcell, fxcell, fzcell,    &
                                                  ximage(it,jt),            &
                                                  zimage(it,jt),            &
                                                  npts, ypts, fpts)

                           if (npts == 2) then  ! Trapezoid rule increment

!                             The first quadrature here is of n - 1, but it
!                             should be that of n - n0 for interferograms.

                              dy = half * abs ((ypts(2) - ypts(1)))
                              quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                     dy*(fpts(1,1) + fpts(1,2))
                              fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                              fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                              quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                     dy*(fpts(2,1) + fpts(2,2))
                              quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                     dy*(fpts(3,1) + fpts(3,2))
                           end if  ! Else npts = 0

                        end do
                     end do

                  case (13)  ! VTK_WEDGE (prism)

!                    We translate from VTK convention to FieldView ordering:

                     in2 = VTK_us_grid%list(il+1)
                     in3 = VTK_us_grid%list(il+2)
                     in5 = VTK_us_grid%list(il+3)
                     in1 = VTK_us_grid%list(il+4)
                     in4 = VTK_us_grid%list(il+5)
                     in6 = VTK_us_grid%list(il+6)
                     jl  = il ! In case we cycle
                     il  = il + 7

                     xcmin = min (VTK_us_grid%x(in1), VTK_us_grid%x(in2), &
                                  VTK_us_grid%x(in3), VTK_us_grid%x(in4), &
                                  VTK_us_grid%x(in5), VTK_us_grid%x(in6))
                     if (xcmin > x2) cycle ! Cell is right of image subregion

                     xcmax = max (VTK_us_grid%x(in1), VTK_us_grid%x(in2), &
                                  VTK_us_grid%x(in3), VTK_us_grid%x(in4), &
                                  VTK_us_grid%x(in5), VTK_us_grid%x(in6))
                     if (xcmax < x1) cycle ! Cell is left of image subregion

                     zcmin = min (VTK_us_grid%z(in1), VTK_us_grid%z(in2), &
                                  VTK_us_grid%z(in3), VTK_us_grid%z(in4), &
                                  VTK_us_grid%z(in5), VTK_us_grid%z(in6))
                     if (zcmin > z2) cycle ! Cell is above image subregion

                     zcmax = max (VTK_us_grid%z(in1), VTK_us_grid%z(in2), &
                                  VTK_us_grid%z(in3), VTK_us_grid%z(in4), &
                                  VTK_us_grid%z(in5), VTK_us_grid%z(in6))
                     if (zcmax < z1) cycle ! Cell is below image subregion

!                    Pixel range that may encounter this cell (which may
!                    be only partially inside the image subregion):

                     i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                     i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                     j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                     j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                  write (50, '(i4, 1p, 4e15.7, 4i4)') &
!!!                     i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                     do jt = j1, j2
                        do it = i1, i2

                           call process_prism (npoints,                        &
                                               VTK_us_grid%list(jl+1:jl+6),    &
                                               VTK_us_grid%x,                  &
                                               VTK_us_grid%y,                  &
                                               VTK_us_grid%z,                  &
                                               VTK_us_grid%f(:,1),             &
                                               VTK_us_grid%f(:,2),             &
                                               VTK_us_grid%f(:,4),             &
                                               ximage(it,jt),                  &
                                               zimage(it,jt),                  &
                                               npts, ypts, fpts)

                           if (npts == 2) then  ! Trapezoid rule increment

!                             The first quadrature here is of n - 1, but it
!                             should be that of n - n0 for interferograms.

                              dy = half * abs ((ypts(2) - ypts(1)))
                              quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                     dy*(fpts(1,1) + fpts(1,2))
                              fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                              fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                              quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                     dy*(fpts(2,1) + fpts(2,2))
                              quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                     dy*(fpts(3,1) + fpts(3,2))
                           end if  ! Else npts = 0

                        end do
                     end do

                  case (14)  ! VTK_PYRAMID

                     in1 = VTK_us_grid%list(il+1)
                     in2 = VTK_us_grid%list(il+2)
                     in3 = VTK_us_grid%list(il+3)
                     in4 = VTK_us_grid%list(il+4)
                     in5 = VTK_us_grid%list(il+5)
                     jl  = il ! In case we cycle
                     il  = il + 6

                     xcmin = min (VTK_us_grid%x(in1), VTK_us_grid%x(in2), &
                                  VTK_us_grid%x(in3), VTK_us_grid%x(in4), &
                                  VTK_us_grid%x(in5))
                     if (xcmin > x2) cycle ! Cell is right of image subregion

                     xcmax = max (VTK_us_grid%x(in1), VTK_us_grid%x(in2), &
                                  VTK_us_grid%x(in3), VTK_us_grid%x(in4), &
                                  VTK_us_grid%x(in5))
                     if (xcmax < x1) cycle ! Cell is left of image subregion

                     zcmin = min (VTK_us_grid%z(in1), VTK_us_grid%z(in2), &
                                  VTK_us_grid%z(in3), VTK_us_grid%z(in4), &
                                  VTK_us_grid%z(in5))
                     if (zcmin > z2) cycle ! Cell is above image subregion

                     zcmax = max (VTK_us_grid%z(in1), VTK_us_grid%z(in2), &
                                  VTK_us_grid%z(in3), VTK_us_grid%z(in4), &
                                  VTK_us_grid%z(in5))
                     if (zcmax < z1) cycle ! Cell is below image subregion

!                    Pixel range that may encounter this cell (which may
!                    be only partially inside the image subregion):

                     i1 = max (int ((xcmin - x1peps) * dxr) + 1, 1)
                     i2 = min (int ((xcmax - x1meps) * dxr) + 1, ipixels)
                     j1 = max (int ((zcmin - z1peps) * dzr) + 1, 1)
                     j2 = min (int ((zcmax - z1meps) * dzr) + 1, jpixels)

!!!                  write (50, '(i4, 1p, 4e15.7, 4i4)') &
!!!                     i,xcmin,xcmax,zcmin,zcmax,i1,i2,j1,j2

                     do jt = j1, j2
                        do it = i1, i2

                           call process_pyramid (npoints,                      &
                                                 VTK_us_grid%list(jl+1:jl+5),  &
                                                 VTK_us_grid%x,                &
                                                 VTK_us_grid%y,                &
                                                 VTK_us_grid%z,                &
                                                 VTK_us_grid%f(:,1),           &
                                                 VTK_us_grid%f(:,2),           &
                                                 VTK_us_grid%f(:,4),           &
                                                 ximage(it,jt),                &
                                                 zimage(it,jt),                &
                                                 npts, ypts, fpts)

                           if (npts == 2) then  ! Trapezoid rule increment

!                             The first quadrature here is of n - 1, but it
!                             should be that of n - n0 for interferograms.

                              dy = half * abs ((ypts(2) - ypts(1)))
                              quadratures(1,it,jt) = quadratures(1,it,jt) + &
                                                     dy*(fpts(1,1) + fpts(1,2))
                              fpts(2:3,1) = fpts(2:3,1)/(fpts(1,1) + one)
                              fpts(2:3,2) = fpts(2:3,2)/(fpts(1,2) + one)
                              quadratures(2,it,jt) = quadratures(2,it,jt) + &
                                                     dy*(fpts(2,1) + fpts(2,2))
                              quadratures(3,it,jt) = quadratures(3,it,jt) + &
                                                     dy*(fpts(3,1) + fpts(3,2))
                           end if  ! Else npts = 0

                        end do
                     end do

                  case default

               end select  ! US3D/VTK cell type

               if (mod (i, 100000) == 0 .or. i == ncell) then
                  write (lunlog, '(a, i10, a, i10)') &
                     '   Finished cell', i, ' of', ncell
               end if

            end do  ! Next US3D/VTK cell

         case default

      end select  ! Grid type

      end subroutine one_cell_all_rays

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine shadow_image ()

!     Derive the indicated image type from the integrated functions of
!     refractive index.  Use the same pixel coordinates as the array of
!     incident rays, since deflections are tiny.  [Is this naive?]
!     Note that tan (eps) ~ eps for these small deflections.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local constants:

      real, parameter :: angle_limit = 0.07  ! Radians, ~ 4 degrees
      real, parameter :: one = 1., zero = 0.

!     Local variables:

      integer :: i, ic, j, jc, npixels
      real    :: epsx_max, epsx_min, epsz_max, epsz_min, eps_scale
      real    :: x1, x_deflected, z1, z_deflected
      real    :: dlimit, dmax, raw, xmean, xsd, zmean, zsd
      real    :: xlimit1, xlimit2, zlimit1, zlimit2

!     Execution:

      npixels = ipixels * jpixels
      x1 = ximage(1,1);  z1 = zimage(1,1)

      select case (image_type)

         case (1, 2)  ! Shadowgraph or schlieren

!           Caution:  Heuristics ahead!
!           Shift/scale the integrated deflection angles so that they don't
!           exceed the range of validity of tan (eps) ~ eps, which is about
!           +/4 degrees ~ +/- 0.07 radians.  Then the camera distance can
!           vary the absolute distance deflections.

            epsx_max = quadratures(2,1,1);  epsx_min = epsx_max
            epsz_max = quadratures(3,1,1);  epsz_min = epsz_max

            do j = 1, jpixels
               do i = 1, ipixels
                  epsx_max = max (quadratures(2,i,j), epsx_max)
                  epsx_min = min (quadratures(2,i,j), epsx_min)
                  epsz_max = max (quadratures(3,i,j), epsz_max)
                  epsz_min = min (quadratures(3,i,j), epsz_min)
               end do
            end do

            dmax = max (abs (epsx_max), abs (epsx_min), &
                        abs (epsz_max), abs (epsz_min))
            eps_scale = angle_limit / dmax

!           However ...
!           Propulsion plumes can give extreme data ranges where the large
!           values swamp the smaller values of interest.  Therefore, scale
!           the deflections based on their standard deviation, not their max.
!           In this call to a one-pass method, 3 is the stride:

            call stdev1p (npixels, quadratures(2,1,1), 3, xmean, xsd)
            call stdev1p (npixels, quadratures(3,1,1), 3, zmean, zsd)

            write (lunlog, '(/, a, 1p, 2e14.6)') &
               'Mean x & z deflections (raw/unscaled):      ', xmean, zmean
            write (lunlog, '(a, 1p, 2e14.6)') &
               'Corresponding standard deviations:          ', xsd, zsd
            write (lunlog, '(/, a, 1p, e14.6)') &
               'Max/min-based scale for deflection angles:  ', eps_scale

            dlimit = clip_factor * max (xsd, zsd)  ! Large clip has no effect

            if (dlimit < dmax) then

               eps_scale = angle_limit / dlimit

               write (lunlog, '(/, a, 1p, e14.6)') &
                  'Standard-deviation-based scale (initial):   ', eps_scale

               xlimit1 = xmean - dlimit;  xlimit2 = xmean + dlimit
               zlimit1 = zmean - dlimit;  zlimit2 = zmean + dlimit

               do j = 1, jpixels
                  do i = 1, ipixels
                     if (quadratures(1,i,j) == zero) cycle ! Ray missed CFD grid
                     quadratures(2,i,j) = max (xlimit1, min (xlimit2, &
                                               quadratures(2,i,j)))
                     quadratures(3,i,j) = max (zlimit1, min (zlimit2, &
                                               quadratures(3,i,j)))
                  end do
               end do

               call stdev1p (npixels, quadratures(2,1,1), 3, xmean, xsd)
               call stdev1p (npixels, quadratures(3,1,1), 3, zmean, zsd)

               dlimit    = clip_factor * max (xsd, zsd)
               eps_scale = angle_limit / dlimit

               write (lunlog, '(/, a, 1p, 2e14.6)') &
                  'Mean x & z deflections (adjusted/unscaled): ', xmean, zmean
               write (lunlog, '(a, 1p, 2e14.6)') &
                  'Corresponding standard deviations:          ', xsd, zsd
               write (lunlog, '(/, a, 1p, e14.6)') &
                  'Standard-deviation-based scale (final):     ', eps_scale

            end if

            eps_scale = camera_distance * eps_scale

!           Contouring the deflected offsets from the "lower left" pixel shows
!           an interesting view of the flow features when the two sets of
!           contours are overlaid.

            do j = 1, jpixels
               do i = 1, ipixels
                  if (quadratures(1,i,j) == zero) cycle ! No intersection w/ CFD

                  x_deflected = eps_scale*quadratures(2,i,j) + ximage(i,j)
                  deflected_i(i,j) = (x_deflected - x1)/dx + one
                  ic = max (1, min (ipixels, nint ((x_deflected - x1)/dx + 1)))

                  z_deflected = eps_scale*quadratures(3,i,j) + zimage(i,j)
                  deflected_j(i,j) = (z_deflected - z1)/dz + one
                  jc = max (1, min (jpixels, nint ((z_deflected - z1)/dz + 1)))

                  fimage(ic,jc) = fimage(ic,jc) + one
               end do
            end do

         case (3)  ! Interferogram (uncertain at this stage)

         case default

      end select

      end subroutine shadow_image

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_results ()

!     Various other quantities are saved for plotting as well as the shadow-type
!     image data.  Tecplot can apply equations and do the gray-scale contouring.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local variables:

      integer :: ib

!     Execution:

!     Incident ray grid = image pixel grid:

      if (ASCII_results) then
         open  (luninfo, file='pixels.g', form='formatted', status='unknown')
         write (luninfo, '(i1)') 1
         write (luninfo, '(2i5)') ipixels, jpixels
         write (luninfo, '(1p, 8e16.7)') ximage, zimage
      else
         open  (luninfo, file='pixels.gu', form='unformatted', status='unknown')
         write (luninfo) 1
         write (luninfo) ipixels, jpixels
         write (luninfo) ximage, zimage
      end if
      close (luninfo)

!     Partial dn/dx and dn/dz on flow volume grid:

!!!   open  (luninfo, file='fx.fu', form='unformatted', status='unknown')
!!!   write (luninfo) nblocks
!!!   write (luninfo) &
!!!      (xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, 1, ib = 1, nblocks)
!!!   do ib = 1, nblocks
!!!      write (luninfo) xyzq(ib)%fx
!!!   end do
!!!   close (luninfo)

!!!   open  (luninfo, file='fz.fu', form='unformatted', status='unknown')
!!!   write (luninfo) nblocks
!!!   write (luninfo) &
!!!      (xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, 1, ib = 1, nblocks)
!!!   do ib = 1, nblocks
!!!      write (luninfo) xyzq(ib)%fz
!!!   end do
!!!   close (luninfo)

!     Integrated quantities (before possible reprocessing):

      if (ASCII_results) then
         open  (luninfo, file='quadratures.f', form='formatted', &
                status='unknown')
         write (luninfo, '(i1)') 1
         write (luninfo, '(3i5)') ipixels, jpixels, 3
         write (luninfo, '(1p, 8e16.7)') &
            quadratures(1,:,:), quadratures(2,:,:), quadratures(3,:,:)
      else
         open  (luninfo, file='quadratures.fu', form='unformatted', &
                status='unknown')
         write (luninfo) 1
         write (luninfo) ipixels, jpixels, 3
         write (luninfo) &
            quadratures(1,:,:), quadratures(2,:,:), quadratures(3,:,:)
      end if
      close (luninfo)

!     Image intensities at each pixel:

      if (ASCII_results) then
         open  (luninfo, file='image.f', form='formatted', status='unknown')
         write (luninfo, '(i1)') 1
         write (luninfo, '(3i5)') ipixels, jpixels, 1
         write (luninfo, '(1p, 8e16.7)') fimage
      else
         open  (luninfo, file='image.fu', form='unformatted', status='unknown')
         write (luninfo) 1
         write (luninfo) ipixels, jpixels, 1
         write (luninfo) fimage
      end if
      close (luninfo)

!     Fractional i and j pixel indices of deflected rays:

      if (ASCII_results) then
         open  (luninfo, file='ic.f', form='formatted', status='unknown')
         write (luninfo, '(i1)') 1
         write (luninfo, '(3i5)') ipixels, jpixels, 1
         write (luninfo, '(10f10.3)') deflected_i
      else
         open  (luninfo, file='ic.fu', form='unformatted', status='unknown')
         write (luninfo) 1
         write (luninfo) ipixels, jpixels, 1
         write (luninfo) deflected_i
      end if
      close (luninfo)

      if (ASCII_results) then
         open  (luninfo, file='jc.f', form='formatted', status='unknown')
         write (luninfo, '(i1)') 1
         write (luninfo, '(3i5)') ipixels, jpixels, 1
         write (luninfo, '(10f10.3)') deflected_j
      else
         open  (luninfo, file='jc.fu', form='unformatted', status='unknown')
         write (luninfo) 1
         write (luninfo) ipixels, jpixels, 1
         write (luninfo) deflected_j
      end if
      close (luninfo)

      end subroutine write_results

   end program shadowgraph
