! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program cross_sections
!
!  Description:
!
!     For a CFD surface dataset, structured or unstructured in three-space,
!  perform slicing at the indicated stations and process the slice data in
!  some way.  Common usage is expected to involve just one function, so no
!  provision was made for more than one originally.  EXCEPTION:  Slicing at
!  multiple azimuthal angles in a front view has been retrofitted, and there
!  it makes sense to allow any number of surface functions, all output along
!  slices in Tecplot form, as opposed to tabulating statistics for one function
!  at many X stations.  (Later:  This version can now read Tecplot datasets as
!  well as PLOT3D datasets.  If the number of functions is greater than 1, the
!  original tabulation of max/min/mean data is simply not an option - the slices
!  can only be output in Tecplot form, whether at X or angular stations.)
!
!     The initial application is to surface heating data, where maxima and
!  minima and associated quantities are required at each of many stations.
!  Means and standard deviations are calculated via arc-length-based monotonic
!  spline quadrature.
!
!     The slicing is implemented for X stations only, with the option to rotate
!  the surface coordinates as needed for such a simplifying approach.
!  EXCEPTION:  Slicing at angular stations assumes the slices are parallel to
!  the X axis with no rotation option to ensure this, at least initially.
!  Program ADJUST_GRID is available if necessary to satisfy this assumption.
!  The common origin of the angular slices need not be on the X axis.
!
!     Two-point line segments are determined from triangular elements met by a
!  given slice plane, and all such segments from any one slice are organized
!  into one or more contiguous curves with no duplicate points (except possibly
!  the end points of a closed curve), as would be needed for line plotting
!  purposes - especially if more than one geometry component is represented by
!  the surface grid.
!
!     The option for slicing at angular stations can handle both half and whole
!  bodies in a way appropriate to an axisymmetric blunt body.  E.g., "spokes"
!  of a half body all start with arc-length s = 0 at the origin of the slices
!  (if only one contiguous curve is present in each slice).  Likewise, s starts
!  at 0 in a consistent way for the angular slices of a whole body (either at
!  or the other, depending on the ordering of the first slice, since each slice
!  after the first is arranged to start at the end nearest to the start of the
!  previous slice).
!
!  Control File (Standard Input):
!
!     CROSS_SECTIONS Control File
!     1          ! Analysis type:      1|2 = max/min/mean table|Tecplotable cuts
!     1          ! Grid type:      1|2|3|4 = PLOT3D|PLOT3D/iblank|Tecplot|FUN3D
!     3          ! Dimensions:         2|3 = 2D/(x,z)(inactive)|3D/(x,y,z)
!     mygrid.gu  ! Grid file name | Tecplot or FUN3D dataset name
!     F          ! Formatted?          T|F = ASCII|unformatted
!     myflow.fu  ! Flow data file name if PLOT3D data
!     F          ! Formatted?          T|F = ASCII|unformatted
!     slices.dat ! Output file:        Tabulation|Tecplot slices
!     X          ! slice coordinate:   x|y|z|t; t[heta] => angular stations
!     1000       ! Number of slices          (uniform for now)
!     0. 0. 0.   ! Origin for clocking angle calculations or angular slices
!     T          ! Merge subslices?    T|F = 1|1+ tabulation lines/slice
!
!  Output Table Format (some digits suppressed, with merging turned off):
!
!  # slice station       fmax    clock       fmin    clock       mean    st.dev.
!
!      2  -11.1657  3.580E+04  -19.450  3.142E+04  -90.000  3.345E+04  1.372E+03
!      2  -11.1657  4.072E+04   58.846  3.881E+04   17.308  4.035E+04  4.290E+02
!      3  -11.1257  2.603E+04 -167.549  1.402E+04  -90.000  1.962E+04  3.768E+03
!      3  -11.1257  3.228E+04   58.846  2.958E+04  164.087  3.182E+04  6.079E+02
!      4  -11.0858  3.088E+04  124.365  5.343E+03 -103.596  2.282E+04  9.238E+03
!      5  -11.0458  3.106E+04  136.212  3.523E+03  -90.000  2.571E+04  7.127E+03
!      :     :       :            :      :            :      :          :
!
!      [An extra arc length column is now included, with the option to merge
!       subslices or not.  Merging is performed by storing subtotals that can
!       be combined into single totals for the means and standard deviations.
!       Mins. and maxs. are determined across all subslices.]
!
!  Output Tecplot Slices Format (analysis type 2, one zone per [sub]slice):
!
!     TITLE = ""
!     VARIABLES = "x"
!     "y"
!     "z"
!     "s"
!     "v1"
!      :
!     "vn"
!     ZONE T="Theta=  0.0000"
!     I=123, J=1, K=1, ZONETYPE=Ordered
!     DATAPACKING=POINT
!     1.234 0.000 0.000 0.000 5.678 9....
!      :     :     :     :     :     :
!     ZONE T="Theta=  1.0000"
!     I=...
!      :     :     :     :     :     :
!
!     Actually, if there is more than one contiguous curve in a slice, that
!     number of zones is output for that slice, all with the same zone title.
!
!  Origins:
!
!        Todd White requested the initial capability to analyze a database of
!     Launch Abort Vehicle flow solutions with TPS thickness in mind.  David
!     Saunders had already written a couple of surface slicing applications
!     employing quadrilateral/triangle utilities originated by Scott Thomas.
!     The framework of the more recent SHADOWGRAPH program was considered
!     more flexible than the earlier MULTICUT and MULTIPLOT.
!
!  History:
!
!     Feb., 2010  D. A. Saunders  Initial adaptation of the SHADOWGRAPH
!                                 framework, which already allowed for
!                                 choices of grid types and outputs.
!     03/03/10       "     "      Added option to handle iblanked grids;
!                                 # slices < 0 turns on printing of slices
!                                 (for analysis-type 1 cases).
!     04/12/10       "     "      Todd asked for a single line of tabulated
!                                 output where two or more have been produced
!                                 if more than one subslice is identified.
!                                 The option to preserve the subslices is
!                                 retained with the addition of an extra column
!                                 showing total arc length for each subslice.
!     05/18/10       "     "      Storing slices before analyzing them was using
!                                 lots of memory for Todd's extreme cases.
!                                 This version fully processes each slice before
!                                 proceeding to the next, at a cost in clarity.
!     09/14/10       "     "      Todd's use of SORT_SURFACE_SLICE suggested
!                                 extending CROSS_SECTIONS to do angular slices.
!     09/15/10       "     "      Half bodies require special treatment of the
!                                 symmetry-plane cut, which is really 2 "spokes"
!                                 where the other cuts are each just 1 spoke.
!                                 Arc lengths should all start at the origin.
!     09/17/10       "     "      Whole bodies also require ensuring that arc
!                                 lengths start at consistent ends of slices as
!                                 angular stations are processed.
!     10/08/10       "     "      Incorporated the option to read Tecplot
!                                 surface datasets (ASCII, nf >= 1 with names).
!                                 The original tabulation is available for any
!                                 type of dataset if (and only if) nf = 1.
!     10/13/10       "     "      Allowed for iblanking if nf >= 1: cutsurfi_nf.
!     11/05/14       "     "      Flow file name = 'none' wasn't handled for
!                                 PLOT3D-file input.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_header_structure  ! See Tecplot_io.f90 for these modules
   use grid_block_structure
   use tecplot_io_module
   use xyzq_io_module         ! I/O package for PLOT3D files

   implicit none

!  Local constants:

   integer,   parameter :: maxseg = 5000          ! Max. # line segments/slice
   integer,   parameter :: maxend = maxseg * 2    ! Max. # line segment end pts.
   integer,   parameter :: naddim = 500           ! Max. # separate curves/cut
   integer,   parameter :: lungrid = 1, lunflow = 2, lunout = 3
   integer,   parameter :: lunstdin = 5, lunlog = 6
   integer,   parameter :: name_limit = 32         ! Must match Tecplot_io.f90
   real,      parameter :: eps = 1.e-12           ! Applied to data range -> tol
   real,      parameter :: one = 1., zero = 0.
   logical,   parameter :: false = .false., true = .true.

!  Local variables:

   integer   :: analysis_type, igrid_type, ndim, nslices
   integer   :: ier, ios, islice, lunf, n, nblocks, ncells, nf, ni, nj

   real      :: grid_bbox(6)              ! Data ranges of entire surface grid
   real      :: xcenter, ycenter, zcenter ! Origin for clock angles/angular sl.
   real      :: xcut                      ! Moved here after not storing slices

   logical   :: formatted_grid, formatted_flow, iblanked, merge_subslices,     &
                print_header, testing, xcuts, ysymmetry, zsymmetry

   character :: filename*128, slice_coordinate*1

   character, allocatable, dimension (:) :: &
      temp_name*(name_limit)              ! For inserting arc length when a
                                          ! Tecplot dataset is read and Tecplot-
                                          ! able output is specified

   integer   :: nad(naddim)               ! nad(i) points to curve i of a slice
                                          ! within x/y/z/fslice(1:nf,:)
   integer, allocatable, dimension (:) :: &
      ic_block_root                       ! Allows skipping entire grid blocks

   real, allocatable, dimension (:) :: &
      tslice, xslice, yslice, zslice,  &  ! Packed coordinates of one slice
      xstation                            ! and x stations for all slices

   real, allocatable, dimension (:,:) :: &
      block_bbox, cell_bbox, fslice       ! Bounding boxes and packed functions

   type (grid_header) :: &
      header                              ! Input structured Tecplot data info

   type (grid_type), pointer :: &         ! Input structured surface data proper
      block(:)


!  Execution:
!  ----------


   call read_controls ()          ! These are all local procedures below

   if (ios /= 0) go to 99         ! Single termination philosophy


   call read_cfd_solution ()      ! Read the CFD grid and the flow field data

   if (ios /= 0) go to 99

   call rotate_grid ()            ! So that slicing can be done at X stations
                                  ! (or parallel to Ox for angular stations)
   if (ios /= 0) go to 99


   call compute_bounding_boxes () ! Set up grid/block/cell bounding boxes

   if (ios /= 0) go to 99


   call slice_stations ()         ! Set up slice stations [after any rotation]


   call process_cfd_solution ()   ! From the surface dataset, calculate
                                  ! slice data as 1+ curves per slice and
                                  ! analyze each slice as it's found to avoid
                                  ! storing potentially voluminous data

!! call process_slice_data ()     ! Analyze the slices in whatever way is
!!                                ! specified by analysis_type
!!                                ! No - it's at the lower level now

   call save_results ()           ! If not done already


99 continue  ! Avoid possible system "stop" quirks


!  Local procedures for program CROSS_SECTIONS in the order they're used:

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
      character :: trouble*16

!     Execution:

!     Open a log file for diagnostics:

      open (lunlog, file='cross_sections.log', status='unknown')

!     Read the control file on standard input:

      read (lunstdin, '(a)', iostat=ios)  ! Skip line 1
      if (ios /= 0) then
         trouble = 'first';  go to 90
      end if

      read (lunstdin, *, iostat=ios) analysis_type
      if (ios /= 0) then
         trouble = 'analysis type';  go to 90
      end if

      read (lunstdin, *, iostat=ios) igrid_type
      if (ios /= 0) then
         trouble = 'grid type';  go to 90
      end if

      read (lunstdin, *, iostat=ios) ndim
      if (ios /= 0) then
         trouble = '# dimensions';  go to 90
      else if (ndim == 2) then
         write (lunlog, '(/, a)') 'Slices of 2-D data are not implemented yet.'
         ios = 1;  go to 99
      end if

      read (lunstdin, *, iostat=ios) filename
      if (ios /= 0) then
         trouble = 'grid file name';  go to 90
      end if

      read (lunstdin, *, iostat=ios) formatted_grid
      if (ios /= 0) then
         trouble = 'grid formatting';  go to 90
      end if
      i = 1;  if (formatted_grid) i = 3

      select case (igrid_type)

         case (1, 2) ! PLOT3D

            open (lungrid, file=filename, form=format(i:11), status='old',     &
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

            if (filename(1:4) == 'none') then
               lunf = -lunflow
               nf = 0
            else
               lunf = lunflow
               open (lunflow, file=filename, form=format(i:11), status='old',  &
                     iostat=ios)
               if (ios /= 0) then
                  trouble = 'flow file';  go to 95
               end if
            end if

         case (3) ! Tecplot structured surface

            header%filename  = filename
            header%formatted = true
            header%ndim      = ndim

            read (lunstdin, *, iostat=ios) ! Skip flow filename
            read (lunstdin, *, iostat=ios) ! Skip flow formatting

         case (4) ! FUN3D unstructured surface

            write (lunlog, '(/, a)') 'FUN3D datasets are not supported yet.'
            ios = 1;  go to 99

         case default

            write (lunlog, '(/, a, i5)') 'Unknown grid type:', igrid_type
            ios = 1;  go to 99

      end select

      read (lunstdin, *, iostat=ios) filename
      if (ios /= 0) then
         trouble = 'output file name';  go to 90
      end if

      open (lunout, file=filename, form=format(3:11), status='unknown',        &
            iostat=ios)
      if (ios /= 0) then
         trouble = 'output file';  go to 95
      end if

      read (lunstdin, *, iostat=ios) slice_coordinate
      if (ios /= 0) then
         trouble = 'slice coordinate';  go to 90
      end if

      read (lunstdin, *, iostat=ios) nslices
      if (ios /= 0) then
         trouble = 'number of slices';  go to 90
      end if
      testing = nslices < 0;  nslices = abs (nslices)

      read (lunstdin, *, iostat=ios) xcenter, ycenter, zcenter
      if (ios /= 0) then
         trouble = 'clock angle origin';  go to 90
      end if

      read (lunstdin, *, iostat=ios) merge_subslices
      if (ios /= 0) then
         trouble = 'merging of subslices';  go to 90
      end if

      go to 99

  90  write (lunlog, '(/, 3a)') '*** Error reading ', trim (trouble), ' line.'
      go to 99

  95  write (lunlog, '(/, 4a)') '*** Error opening ', trim (trouble), ': ',    &
                                trim (filename)
  99  return

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_cfd_solution ()

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     implicit none

!     Local variables:

      integer :: iv
      logical :: cell_centered  ! Not used, but needed by xyzq_read, xyziq_read

!     Execution:

      iblanked = false

      select case (igrid_type)

         case (1)  ! Plain PLOT3D multiblock (DPLR)

            call xyzq_read (lungrid, lunf, formatted_grid, nblocks, nf,        &
                            cell_centered, block, ios)

         case (2)  ! PLOT3D/iblanked (DPLR-overset/OVERFLOW)

            iblanked = true

            call xyziq_read (lungrid, lunf, formatted_grid, nblocks, nf,       &
                             cell_centered, block, ios)

         case (3)  ! Tecplot structured surface dataset

            call Tecplot_read (lungrid, header, block, ios)

            if (ios /= 0) go to 99

            nblocks = header%nblocks
            nf      = header%numq

         case (4)  ! FUN3D unstructured surface dataset

            ! TBD (see SHADOWGRAPH)

      end select

      if (analysis_type == 1) then
         if (nf > 1) then
            write (lunlog, '(/, a, i5, a)') '*** Warning: # functions found:', &
               nf, '; Tecplottable slices will be output, not a tabulation.'
            analysis_type = 2
         end if
      end if

      if (analysis_type == 2) then  ! Tecplottable output

         close (lunout)  ! Tec_header_write does the open

         header%filename    = filename  ! Output name now
         header%title       = filename
         header%formatted   = true
         header%ndim        = ndim
         header%numq        = 1 + nf    ! Includes arc length
         header%nblocks     = nslices   ! ?? (could be more zones)
         header%ndatasetaux = 0

!        Set up the Tecplot output header:

         select case (igrid_type)

            case (1, 2, 4)  ! PLOT3D or FUN3D dataset was read

               allocate (header%varname(ndim+1+nf))  ! Insert arc length

               header%varname(1) = 'x'
               header%varname(2) = 'y'
               if (ndim == 3) header%varname(3) = 'z'
               header%varname(ndim+1) = 's'

               do iv = ndim + 1, ndim + nf
                  header%varname(iv+1) = 'f'
               end do

               do iv = ndim + 1, min (ndim + nf, 9)
                  write (header%varname(iv+1)(2:2), '(i1)') iv - ndim
               end do

               if (nf >= 10) then
                  do iv = 10, ndim + nf
                     write (header%varname(iv+1)(2:3), '(i2)') iv - ndim
                  end do
               end if

            case (3)  ! Tecplot dataset was read; need to insert arc length

               allocate (temp_name(ndim+1+nf))

               temp_name(1:ndim) = header%varname(1:ndim)
               temp_name(ndim+1) = 's'
               temp_name(ndim+2:ndim+1+nf) = header%varname(ndim+1:ndim+nf)

               deallocate (header%varname)
               allocate   (header%varname(ndim+1+nf))

               header%varname(:) = temp_name(:)
               deallocate (temp_name)

            case default

         end select

      end if

 99   return

      end subroutine read_cfd_solution

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rotate_grid ()

!     Rotate the CFD data if necessary so that the viewing direction is normal
!     to the y = 0 plane.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local constants:

      real, parameter :: one = 1., ninety = 90., zero = 0.

!     Local variables:

      integer :: ib
      real    :: angle, Px, Py, Pz, Qx, Qy, Qz

!     Execution:

!     If necessary, rotate the CFD grid about the line joining P (Px, Py, Pz)
!     and Q (Qx, Qy, Qz) in-place to give the specified view:

      select case (slice_coordinate)

         case ('x', 'X')  ! Most likely case - stream stations;  no rotation

            angle = zero

         case ('y', 'Y')  ! Span stations

            angle = ninety
            Px = zero;  Py = zero;  Pz = zero
            Qx = zero;  Qy = zero;  Qz = one

         case ('z', 'Z')  ! Vertical stations/horizontal slices

            angle = ninety
            Px = zero;  Py = zero;  Pz = zero
            Qx = zero;  Qy = one;   Qz = zero

         case ('t', 'T')  ! Angular slices parallel to Ox over 180 deg

            angle = zero  ! For now
            slice_coordinate = 'T' ! Whole body should still just slice over 180

         case default

            write (lunlog, '(/, 2a, /, a)') &
               'Bad input slice definition: ', slice_coordinate, &
               'Enter x, y, z, or t[heta] for angular slices.'
            ios = 1

      end select

      if (angle /= zero) then

         select case (igrid_type)

            case (1, 2, 3)  ! Structured multiblock data (iblanked or not)

               do ib = 1, nblocks

                  ni = block(ib)%ni;  nj = block(ib)%nj

                  call rotate3d (ni*nj, block(ib)%x, block(ib)%y, &
                                 block(ib)%z, angle, Px, Py, Pz, Qx, Qy, Qz)
               end do

               call rotate3d (1, xcenter, ycenter, zcenter, angle, &
                              Px, Py, Pz, Qx, Qy, Qz)

            case (4)  ! Unstructured FUN3D dataset

         end select

      end if

      end subroutine rotate_grid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine compute_bounding_boxes ()
!
!     Determine the grid bounding box in the x direction only (after any
!     rotation), since we're slicing at x stations only.
!     If appropriate, set up grid block bounding boxes as well.
!     (Later: the angular slice case, possibly on half bodies, needs the grid
!     bounding box in all directions.)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local variables:

      integer :: i, ib, ic, j
      real    :: xmax, xmin, ymax, ymin, zmax, zmin

!     Execution:

      select case (igrid_type)

         case (1, 2, 3)  ! Structured multiblock surface, blanked or not

            allocate (block_bbox(2,nblocks))

!!!         if (use_cell_bboxes) then

!!!            allocate (ic_block_root(nblocks+1))  ! To increment ic properly

!!!            ncells = 0;  ic_block_root(1) = ncells

!!!            do ib = 1, nblocks
!!!               ni = block(ib)%ni;  nj = block(ib)%nj
!!!               ncells = (ni - 1) * (nj - 1) + ncells
!!!               ic_block_root(ib+1) = ncells
!!!            end do

!!!            allocate (cell_bbox(2,ncells), stat=ios)

!!!            if (ios /= 0) then
!!!               write (lunlog, '(/, a, i12)') &
!!!                  '*** Trouble allocating CFD cell bounding boxes:', ncells
!!!               go to 99
!!!            end if

!!!            write (lunlog, '(/, a, i12)') 'Total # CFD grid cells:', ncells
!!!            ic = 0  ! Cell counter
!!!         end if

            do ib = 1, nblocks

               xmin = block(ib)%x(1,1,1);  xmax = xmin
               ymin = block(ib)%y(1,1,1);  ymax = ymin
               zmin = block(ib)%z(1,1,1);  zmax = zmin

               if (ib == 1) then
                  grid_bbox(1) = xmin;  grid_bbox(2) = xmax
                  grid_bbox(3) = ymin;  grid_bbox(4) = ymax
                  grid_bbox(5) = zmin;  grid_bbox(6) = zmax
               end if

               ni = block(ib)%ni;  nj = block(ib)%nj

               do j = 1, nj
                  do i = 1, ni
                     xmin = min (block(ib)%x(i,j,1), xmin)
                     xmax = max (block(ib)%x(i,j,1), xmax)
                     ymin = min (block(ib)%y(i,j,1), ymin)
                     ymax = max (block(ib)%y(i,j,1), ymax)
                     zmin = min (block(ib)%z(i,j,1), zmin)
                     zmax = max (block(ib)%z(i,j,1), zmax)
                  end do
               end do

               block_bbox(1,ib) = xmin;  block_bbox(2,ib) = xmax

               grid_bbox(1) = min (grid_bbox(1), xmin)
               grid_bbox(2) = max (grid_bbox(2), xmax)
               grid_bbox(3) = min (grid_bbox(3), ymin)
               grid_bbox(4) = max (grid_bbox(4), ymax)
               grid_bbox(5) = min (grid_bbox(5), zmin)
               grid_bbox(6) = max (grid_bbox(6), zmax)

!              Set up cell bounding boxes for the whole grid:

!!!            do j = 1, nj - 1
!!!               do i = 1, ni - 1
!!!                  ic = ic + 1
!!!                  cell_bbox(1,ic) = min &
!!!                     (block(ib)%x(i,j,1),     block(ib)%x(i+1,j,1),         &
!!!                      block(ib)%x(i,j+1,1),   block(ib)%x(i+1,j+1,1))
!!!                  cell_bbox(2,ic) = max &
!!!                     (block(ib)%x(i,j,1),     block(ib)%x(i+1,j,1),         &
!!!                      block(ib)%x(i,j+1,1),   block(ib)%x(i+1,j+1,1))
!!!               end do
!!!            end do

            end do  ! Next grid block

            write (lunlog, '(/, a)') 'Block bounding box x limits:'
            write (lunlog, '(i6, 1p, 2e24.15)') &
               (ib, block_bbox(:,ib), ib = 1, nblocks)

!!!         write (31, '(i9, 1p, 2e15.6)') (ic, cell_bbox(:,ic), ic = 1, ncells)

         case (4)  ! FUN3D, unstructured, TBD

      end select

      write (lunlog, '(/, a, /)') 'Grid bounding box x/y/z limits:'
      write (lunlog, '(1p, 2e20.11, e23.11, e20.11, e23.11, e20.11)')          &
         grid_bbox(:)

 99   return

      end subroutine compute_bounding_boxes

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine slice_stations ()
!
!     Set up a sequence of slice stations along the geometry in the X direction
!     (or over a half or whole rotation in the case of angular stations).
!     Choose slice_coordinate = 'Y' or 'Z' for slicing in other directions, and
!     the grid will be rotated before slicing.  X is the working direction.
!     (For angular cases, it is assumed that no rotation is needed, and "x"
!     will refer to angle theta between 0 and 180 or 360.)
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Local constants:

      real, parameter :: half = 0.5, one80 = 180.

!     Local variables:

      integer :: i
      real    :: d1, d2, dtheta, dx

!     Execution:

      allocate (xstation(nslices))

      select case (slice_coordinate)

         case ('T')  ! T = theta: slices over 180 deg for half or whole body

            xcuts = false

!           ysymmetry is assumed true if only a left or right half is present

            d1 = ycenter - grid_bbox(3)
            d2 = grid_bbox(4) - ycenter
            ysymmetry = abs (d2 - d1) > (grid_bbox(4) - grid_bbox(3)) * half
            d1 = zcenter - grid_bbox(5)
            d2 = grid_bbox(6) - zcenter
            zsymmetry = abs (d2 - d1) > (grid_bbox(6) - grid_bbox(5)) * half

            write (lunlog, '(/, a, 2l3)') &
               'Y and Z symmetry determined:', ysymmetry, zsymmetry

            dtheta = one80 / real (nslices-1)
            do i = 1, nslices
               xstation(i) = dtheta * real (i-1)
            end do

            xstation(nslices) = one80  ! Exactly
            xcenter = grid_bbox(1)     ! Presumed nose location at xmin

         case default  ! Original scheme at X cuts

            xcuts     = true
            ysymmetry = false;  zsymmetry = false
            dx        = (grid_bbox(2) - grid_bbox(1)) / real (nslices - 1)

            xstation(1) = grid_bbox(1) + eps  ! Fudge at the (likely) nose

            do i = 2, nslices - 1
               xstation(i) = grid_bbox(1) + dx * real (i-1)
            end do

            xstation(nslices) = grid_bbox(2) - eps  ! Fudge at (likely) aft end

      end select

      write (lunlog, '(/, a, /)') 'Slice stations:'

      if (nslices <= 100) then
         write (lunlog, '(i5, 1p, e21.12)') (i, xstation(i), i = 1, nslices)
      else
         write (lunlog, '(i5, 1p, e21.12)') (i, xstation(i), i = 1, 10)
         write (lunlog, '(a)') '    :       :'
         write (lunlog, '(i5, 1p, e21.12)') &
            (i, xstation(i), i = nslices-9, nslices)
      end if

      end subroutine slice_stations

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine process_cfd_solution ()

!     Perform the slicing of the given surface dataset, processing each slice
!     as it is made.  (Originally, slices were stored, but not now.)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     implicit none

!     Local variables:

      integer :: i, ib, ic, ier, nnode, nseg
      real    :: epsf, epsxyz, p(3), v(3)

      integer :: ndgcut(2,maxseg)       ! Indices in xyzfcut(*,:) for the end
                                        ! points of each 2-pt. line segment

      real    :: xyzfcut(ndim+nf,maxend)! Packed (x,y,z,f) coordinates of end
                                        ! points of 2-pt. line segments for
                                        ! one slice (no duplicates)
!     Execution:

      epsxyz = (grid_bbox(2) - grid_bbox(1)) * eps ! Tolerance for duplicate pts
      epsf   = epsxyz * 1.e+30          ! Questionable tests in cutsurf/cutorder

      if (xcuts) then                   ! Set a vector normal to each cut plane
         p(:) = zero                    ! Point on a cut; p(1) is set to xcut
         v(:) = zero;  v(1) = one       ! Vector from the pt. defines cut plane
      else
         p(1) = xcenter
         p(2) = ycenter
         p(3) = zcenter                  ! v(2) and v(3) depend on the angular
         v(1) = zero                     ! station, with rotation axis || Ox
      end if

      allocate (xslice(maxseg), yslice(maxseg), zslice(maxseg), &
                tslice(maxseg), fslice(nf,maxseg))

      print_header = true

      select case (igrid_type)

         case (1, 2, 3)  ! Structured multiblock surface data, iblanked or not

            do islice = 1, nslices

               xcut  = xstation(islice)
               nseg  = 0
               nnode = 0

               do ib = 1, nblocks

                  if (xcuts) then
                     if (xcut < block_bbox(1,ib)) cycle
                     if (xcut > block_bbox(2,ib)) cycle

                     p(1) = xcut
                  else
                     v(2) = -cosd (xcut)  ! If ysymmetry;
                     v(3) = -sind (xcut)  ! zsymmetry unlikely and not handled
                  end if

                  ni = block(ib)%ni
                  nj = block(ib)%nj

!!!               if (use_cell_bboxes) ic = ic_block_root(ib)

!                 Cut the current grid block as a set of 2-point line segments:

                  if (.not. iblanked) then

                     if (nf == 1) then  ! Slightly more efficient

                        call cutsurf (block(ib)%x, block(ib)%y, block(ib)%z,   &
                                      block(ib)%q, 1, ni, 1, nj, ni, nj,       &
                                      xyzfcut, nnode, maxend, ndgcut, nseg,    &
                                      maxseg, p, v, epsxyz, epsf, ier)
                     else  ! nf > 1

                        call cutsurf_nf (ni, nj, 1, ni, 1, nj, nf,             &
                                         maxend, maxseg,                       &
                                         block(ib)%x, block(ib)%y,             &
                                         block(ib)%z, block(ib)%q,             &
                                         p, v, epsxyz, nnode, nseg, ndgcut,    &
                                         xyzfcut, ier)
                     end if

                  else

                     if (nf == 1) then  ! Slightly more efficient

                        call cutsurfi (block(ib)%x, block(ib)%y, block(ib)%z,  &
                                       block(ib)%q, block(ib)%iblank,          &
                                       1, ni, 1, nj, ni, nj, xyzfcut,          &
                                       nnode, maxend, ndgcut, nseg,            &
                                       maxseg, p, v, epsxyz, epsf, ier)
                     else

                        call cutsurfi_nf (ni, nj, 1, ni, 1, nj, nf,            &
                                          maxend, maxseg,                      &
                                          block(ib)%x, block(ib)%y,            &
                                          block(ib)%z, block(ib)%q,            &
                                          p, v, epsxyz, nnode, nseg, ndgcut,   &
                                          xyzfcut, ier)
                     end if

                  end if

                  write (lunlog, '(a, 4i7, f12.6)') &
                     '   islice,ib,nseg,nnode,xcut: ', &
                         islice,ib,nseg,nnode,xcut

                  if (ier /= 0) then
                     write (lunlog, '(a, i3, a, i5, a, f12.6, a, i5, /, a)')   &
                        'Cutsurf error:', ier, '  block #:', ib,               &
                        '  Xcut:', xcut, '  # 2-pt. segments so far:', nseg,   &
                        '  Proceeding.'
                  end if

               end do  ! Next grid block

!              Organize the 2-pt. line segments as contiguous curve(s):

               nad(1) = 0  ! # lines for this slice (may be all blanked)

               if (nseg > 0) then

!!!               if (islice <=5) then
!!!                  write (50, '(a, 3i7, f12.6)') &
!!!                     'islice,nseg,nnode,xcut: ', islice,nseg,nnode,xcut
!!!                  write (50, '(i4, i6, i4)') (i, ndgcut(:,i), i = 1, nseg)

!!!                  write (51, *) 'epsx/f: ', epsxyz, epsf
!!!                  if (nf == 1) then
!!!                     write (51, '(i4, 1p, 4e24.15)') &
!!!                              (i, xyzfcut(:,i), i = 1, nnode)
!!!                  else if (nf == 2) then
!!!                     write (51, '(i4, 1p, 5e24.15)') &
!!!                              (i, xyzfcut(:,i), i = 1, nnode)
!!!                  end if
!!!               end if

!!!               call cutorder (nad, naddim, xslice, yslice, zslice, fslice,  &
!!!                              maxseg, xyzfcut, nnode, maxend, ndgcut, nseg, &
!!!                              maxseg, epsxyz, epsf, ier)

                  call cutorder_nf (nf, maxend, maxseg, maxseg, naddim,        &
                                    nnode, nseg, ndgcut, xyzfcut, epsxyz,      &
                                    xslice, yslice, zslice, fslice, nad, ier)

                  if (ier /= 0) then
                     write (lunlog, '(a, i3, a, i5, a, f12.6, a, i5, /, a)')   &
                        'Cutorder error:', ier, '  slice #:', islice,          &
                        '  Xcut:', xcut, '  # 2-pt. segments:', nseg,          &
                        '  Aborting this slice.'
                     exit
                  end if

!                 Analyze this slice now so we don't have to store it:

                  call process_slice_data (islice)

                  if (ios /= 0) then
                     write (lunlog, '(/, a)') &
                        'Bad return from process_slice_data. Aborting.'
                     go to 99
                  end if

               end if

            end do  ! Next slice

         case (4)  ! FUN3D, unstructured, TBD

      end select

 99   return

      end subroutine process_cfd_solution

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine process_slice_data (islice)

!     Analyze or process the indicated slice as specified by analysis_type.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) :: islice  ! For the full-body angular stations case

!     Local constants:

      character, parameter :: spline_method*1 = 'M'  ! Monotonic spline fits

!     Local variables:

      integer :: i, isqmin, l, l1, l2, line, lmax, lmin, n, nlines, npts
      real    :: dsq, dsq1, dsq2, dsqmin, fl, fmax, fmin, fmean, fstdev,       &
                 stotal, thetamax, thetamin, ymax, ymin, zmax, zmin
      real, save :: x1prev, y1prev, z1prev
      real, allocatable :: arcs(:), fsq(:), fsqsum(:), fsum(:), ssum(:)

!     Execution:

      select case (analysis_type)

         case (1)  ! Surface heating statistics

!!!         write (lunlog, '(a, 3i5)') 'islice, nad12: ', islice, nad(1:2)

            if (print_header) then
                print_header = false
                write (lunout, '(/, 2a, /)') &
                   '# slice   station         fmax    clock         fmin',     &
                   '    clock         mean      st.dev.    total arc'
            end if

            nlines = nad(1)

            if (nlines > 0) then

               nad(1) = 1  ! So the loop works for line 1

               if (nlines == 1 .or. .not. merge_subslices) then  ! This is more
                                                                 ! efficient
                  do line = 1, nlines

                     l1   = nad(line)
                     l2   = nad(line+1) - 1
                     npts = l2 - l1 + 1

!!!                  write (lunlog, '(a, 5i7)') &
!!!                     'islice,nlines,l1,l2,npts: ', &
!!!                      islice,nlines,l1,l2,npts

!                    Max. and min.:

                     fmin = fslice(1,l1);  fmax = fmin
                     lmin = l1;  lmax = l1

                     do l = l1, l2
                        fl = fslice(1,l)
                        if (fl < fmin) then
                           fmin = fl;  lmin = l
                        else if (fl > fmax) then
                           fmax = fl;  lmax = l
                        end if
                     end do

                     if (testing) then
                        write (lunlog, '(1p, 4e24.15)') &
                           (xslice(l), yslice(l), &
                            zslice(l), fslice(1,l), l = l1, l2)
                     end if

!                    Corresponding clocking angles:

                     thetamin = atan2d (zslice(lmin) - zcenter, &
                                        yslice(lmin) - ycenter)
                     thetamax = atan2d (zslice(lmax) - zcenter, &
                                        yslice(lmax) - ycenter)

!                    Mean function value along slice, and standard deviation:

                     allocate (arcs(npts), fsq(npts))

                     call chords3d (npts, xslice(l1), yslice(l1), zslice(l1),  &
                                    false, stotal, arcs)

                     if (stotal > zero) then
                        call lcsquad (npts, arcs, fslice(1,l1), zero, stotal,  &
                                      spline_method, fmean)
                        fmean = fmean / stotal
                        i = 1
                        do l = l1, l2
                           fsq(i) = (fslice(1,l) - fmean)**2
                           i = i + 1
                        end do
                        call lcsquad (npts, arcs, fsq, zero, stotal,           &
                                      spline_method, fstdev)
                        fstdev = sqrt (fstdev / stotal)
                     else
                        fmean  = fslice(1,l1)
                        fstdev = zero
                     end if

                     deallocate (arcs, fsq)

                     write (lunout, &
                        '(i5, f12.6, 1p, e13.5, 0p, f9.3, 1p, e13.5, 0p, f9.3, &
                          1p, 3e13.5)') &
                        islice, xcut, fmax, thetamax, fmin, thetamin, fmean,   &
                        fstdev, stotal

                  end do  ! Next subslice (if any)

               else  ! More than one subslice, and we're merging them

!                 We need the (arc-length-based) mean before we can do the
!                 mean-distance-from-the-mean calculation.  It's not worth
!                 trying to avoid calculating arc lengths more than once.

                  lmin = nad(1);          lmax = lmin
                  fmin = fslice(1,lmin);  fmax = fmin

                  allocate (ssum(nlines), fsum(nlines))

                  do line = 1, nlines

                     l1   = nad(line)
                     l2   = nad(line+1) - 1
                     npts = l2 - l1 + 1

!!!                  write (lunlog, '(a, 5i7)') &
!!!                     'islice,nlines,l1,l2,npts: ', &
!!!                      islice,nlines,l1,l2,npts

                     do l = l1, l2  ! Max. & min over all subslices
                        fl = fslice(1,l)
                        if (fl < fmin) then
                           fmin = fl;  lmin = l
                        else if (fl > fmax) then
                           fmax = fl;  lmax = l
                        end if
                     end do

                     if (testing) then
                        write (lunlog, '(1p, 4e24.15)') &
                           (xslice(l), yslice(l), &
                            zslice(l), fslice(1,l), l = l1, l2)
                     end if

                     allocate (arcs(npts))

                     call chords3d (npts, xslice(l1), yslice(l1), zslice(l1),  &
                                    false, stotal, arcs)
                     ssum(line) = stotal

!                    Integral of f ds along subslice:

                     if (stotal > zero) then
                        call lcsquad (npts, arcs, fslice(1,l1), zero, stotal,  &
                                      spline_method, fmean)
                        fsum(line) = fmean
                     else
                        fsum(line) = fslice(1,l1)
                     end if

                     deallocate (arcs)

                  end do  ! Next subslice for the overall mean calculation

!                 Mean across all subslices:

                  stotal = zero
                  fmean  = zero
                  do line = 1, nlines
                     stotal = stotal + ssum(line)
                     fmean  = fmean  + fsum(line)
                  end do
                  fmean = fmean / stotal

                  deallocate (fsum)

!                 Integrate the squared deviations from the overall mean:

                  allocate (fsqsum(nlines))

                  do line = 1, nlines

                     l1   = nad(line)
                     l2   = nad(line+1) - 1
                     npts = l2 - l1 + 1

                     allocate (arcs(npts), fsq(npts))

                     call chords3d (npts, xslice(l1), yslice(l1), zslice(l1),  &
                                    false, stotal, arcs)

                     if (stotal > zero) then
                        i = 1
                        do l = l1, l2
                           fsq(i) = (fslice(1,l) - fmean)**2
                           i = i + 1
                        end do
                        call lcsquad (npts, arcs, fsq, zero, stotal,           &
                                      spline_method, fstdev)
                        fsqsum(line) = fstdev
                     else
                        fsqsum(line) = zero
                     end if

                     deallocate (arcs, fsq)

                  end do  ! Next subslice

!                 Overall mean deviation from overall mean:

                  stotal = zero
                  fstdev = zero
                  do line = 1, nlines
                     stotal = stotal + ssum(line)
                     fstdev = fstdev + fsqsum(line)
                  end do
                  fstdev = sqrt (fstdev / stotal)

                  deallocate (ssum, fsqsum)

!                 Clocking angles for overal min. and max.:

                  thetamin = atan2d (zslice(lmin) - zcenter, &
                                     yslice(lmin) - ycenter)
                  thetamax = atan2d (zslice(lmax) - zcenter, &
                                     yslice(lmax) - ycenter)

                  write (lunout, &
                     '(i5, f12.6, 1p, e13.5, 0p, f9.3, 1p, e13.5, 0p, f9.3,    &
                       1p, 3e13.5)') &
                     islice, xcut, fmax, thetamax, fmin, thetamin, fmean,      &
                     fstdev, stotal

               end if  ! Merge/don't-merge branch

               nad(1) = nlines

            end if  ! No-line/at-least-one-line-at-this-station branch

         case (2)  ! Write one Tecplot zone per [sub]slice

            if (print_header) then

                print_header = false
                header%datapacking = 0  ! POINT order

!               Open file and use the Tecplot_io utility:

                call Tec_header_write (lunout, header, ios)

                if (ios /= 0) then
                   write (lunlog, '(/, a)') &
                      'Bad return from Tec_header_write.  Aborting.'
                   go to 99
                end if

            end if

            nlines = nad(1)  ! >= 0
            nad(1) = 1       ! So the loop works for line 1

            do line = 1, nlines

               l1   = nad(line)
               l2   = nad(line+1) - 1
               npts = l2 - l1 + 1

               write (*, '(a, i4, 4i5)') &
                  'slice, # lines, l1, l2, npts:', islice, nlines, l1, l2, npts

               allocate (arcs(npts))

!              Specialized tweak of a symmetry-plane slice.  We do the slice
!              twice, and output the upper or lower portion only to get one
!              spoke worth as for other slices.

               if (nlines == 1) then
                  if (ysymmetry) then  ! Assume z is up (sorry)
                     if (xcut == zero .or. xcut == 180.) then

!                       Find the point nearest the slice origin at the nose:

                        dsqmin = 1.e+6;  isqmin = 1

                        do i = 1, npts
                           dsq = (xslice(i) - xcenter)**2 + &
                                 (yslice(i) - ycenter)**2 + &
                                 (zslice(i) - zcenter)**2
                           if (dsq < dsqmin) then
                               dsqmin = dsq;  isqmin = i
                           end if
                        end do

                        if (zslice((1 + isqmin) / 2) > zcenter) then
                           if (xcut == zero) then
                              l2 = isqmin
                           else  ! 180.
                              l1 = isqmin
                           end if
                        else
                           if (xcut == zero) then
                              l1 = isqmin
                           else
                              l2 = isqmin
                           end if
                        end if

                        npts = l2 - l1 + 1

                     end if
                  end if
               end if

!              Start half-body arc lengths from the origin of angular slices:

               if (nlines == 1) then
                  if (ysymmetry .or. zsymmetry) then  ! Make s = 0 at center

                     if (xslice(l1) > xslice(l2)) then
                        call rverse (npts, xslice(l1), xslice(l1))
                        call rverse (npts, yslice(l1), yslice(l1))
                        call rverse (npts, zslice(l1), zslice(l1))
                        do n = 1, nf
                           tslice(1:npts)         = fslice(n,l1:l1+npts-1)
                           call rverse (npts, tslice, tslice)
                           fslice(n,l1:l1+npts-1) = tslice(1:npts)
                        end do
                     end if

                  else

                     if (xcuts) then  ! What?

                     else  ! Start s = 0 at end nearest previous start

                        if (islice > 1) then
                           dsq1 = (x1prev - xslice(1))**2 + &
                                  (y1prev - yslice(1))**2 + &
                                  (z1prev - zslice(1))**2
                           dsq2 = (x1prev - xslice(npts))**2 + &
                                  (y1prev - yslice(npts))**2 + &
                                  (z1prev - zslice(npts))**2
                           if (dsq1 > dsq2) then
                              call rverse (npts, xslice, xslice)
                              call rverse (npts, yslice, yslice)
                              call rverse (npts, zslice, zslice)
                              do n = 1, nf
                                 tslice(1:npts)   = fslice(n,1:npts)
                                 call rverse (npts, tslice, tslice)
                                 fslice(n,1:npts) = tslice(1:npts)
                              end do
                           end if
                        end if

                        x1prev = xslice(1)
                        y1prev = yslice(1)
                        z1prev = zslice(1)

                     end if

                  end if
               end if

               call chords3d (npts, xslice(l1), yslice(l1), zslice(l1), false, &
                              stotal, arcs)

               write (lunout, '(a, f8.4, a)') 'ZONE T="Theta=', xcut, '"'
               write (lunout, '(a, i5, a, /, a)') &
                  'I=', npts, ', J=1, K=1, ZONETYPE=ORDERED', &
                  'DATAPACKING=POINT'
               do l = l1, l2
                  write (lunout, '(1p, 20e15.7)') &
                     xslice(l), yslice(l), zslice(l), arcs(l-l1+1), fslice(:,l)
               end do

               deallocate (arcs)

            end do  ! Next [sub]slice

      end select

 99   return

      end subroutine process_slice_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine save_results ()  ! If not done already

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Execution:

      end subroutine save_results

   end program cross_sections

