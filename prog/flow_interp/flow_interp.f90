!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module face_info_module ! Module for program flow_interp (but unused following option to suppress out-of-range solution blocks).

   public :: face_type

   type face_type
      integer, dimension (6) :: i1, i2, j1, j2, k1, k2             ! Indices defining faces 1-6 in the usual imin, imax, ... order
      real,    dimension (6) :: xmin, xmax, ymin, ymax, zmin, zmax ! Data ranges for the faces
      logical, dimension (6) :: exterior                           ! T|F; F means the face appears to be interior because it
   end type face_type                                              ! matches another face

   end module face_info_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program flow_interp

!  This program interpolates a flow field solution from one multiblock structured grid to another.  The solution may be either
!  cell-centered or vertex-centered.  The same underlying geometry is assumed to be represented by both grids, but the outer
!  boundaries and/or block counts and point distributions may differ.
!
!  Intended application:
!
!     >  Ideally, all target grid cells will be contained by the grid associated with the given flow solution.
!     >  For target cells outside the solution grid, the nearest solution cell is used.  Thus, extrapolation is explicitly avoided.
!     >  The standard approach has long been to use the ADT (Alternating Digital Tree) bounding-box-based technique to ensure
!        efficiency no matter how many of the target points are outside the solution grid.
!     >  Alternatively, KDTREE searches may be specified: these determine just the nearest grid point (probably a cell centroid)
!        much more efficiently than determining the best point within or on a cell and interpolation coefficients to go with it.
!        The search tree requires less memory as well, which may make this the only option on very large grids to be searched.
!     >  Most recently, a hybrid option may be specified: KDTREE searches (not on the input grid vertices, but on the associated
!        cell centroids generated here) followed by refinement with one call per search to find the nearest point to the target
!        point that is not outside the indicated closest-centroid cell.  This is an attempt to combine the best features of the
!        KDTREE and ADT methods.  Unfortunately, on typical hypersonic grids, the centroid closest to a target point may not be
!        that of the best cell that actually contains the point.  (This appears to be possible for high-aspect-ratio cells where
!        the grid spacing in some direction is varying.)  This hybrid method may still be the best compromise for interpolating
!        within very large grids.  Nearest-point (method 2) or refined-within-nearest-centroid-cell (method 3) flow values may be
!        adequate for starting new flow solutions, but the best-possible interpolated values of the standard ADT method are
!        recommended for, say, line-of-sight interpolations intended for radiation calculations.
!     >  The refinements of the hybrid method 3 cost 2.5 - 3 times the cost of the plain KDTREE method 2 results, while the ADT
!        method 1 may be as much as 100 times slower than plain KDTREE on large grids, and it requires more memory.
!     >  Note that application to different geometries (even slightly perturbed ones) is NOT appropriate:  boundary layer point
!        distributions are all-important, and the real-space interpolations performed here will not preserve the boundary layer
!        flow variables (or gradients, especially) anywhere the two surface grids differ more than infinitesimally.
!     >  Further, if the new grid points resolve surface curvature better than the grid being searched, not even the ADT method
!        can guarantee properly varying boundary layer interpolations.  Only for local problems (where the in-flow at interpolated
!        boundaries has to become a fixed boundary condition) should this be an issue, though; otherwise, the flow solver should
!        clean up the interpolated starting guess.
!     >  For interpolating just surface solutions, see RADIAL_INTERP or SURFACE_INTERP.  (RADIAL_INTERP does correct for added
!        resolution at the wall, but it assumes only one layer of grid blocks as is common for atmospheric entry vehicles, while
!        FLOW_INTERP makes no such assumption.)
!     >  This version tabulates results for the special case of single-radial-line target block(s), as needed for comparisons with
!        wind tunnel boundary layer measurements or lines of sight for radiation calculations.  Target ni = nj = 1 invokes this
!        option.  However, hemisphere lines of sight are constructed as 1 x 1 x nk blocks and produce excessive output, so if the
!        number of target blocks is excessive, we suppress this option.
!
!  Assumptions:
!
!     >  All files are PLOT3D-type, multiblock, 3-D (not 2-D), formatted or unformatted.
!     >  If the flow solution block dimensions match the grid block dimensions, the solution is assumed to be vertex-centered.
!        If they differ by 1, the solution must be cell-centered.  Any other finding means the files are inconsistent.
!
!  Clarification for DPLR Users:
!
!        DPLR-related usage normally involves cell-centered solutions with "halo" cells included.  The associated grids should
!        also be cell-centered with halo cells, in which case there is no apparent difference from vertex-centered data as far
!        as FLOW_INTERP is concerned.
!
!  Algorithm:
!
!     >  If the option to attempt suppression of input flow blocks that won't ever be needed is invoked, do the following pre-
!        processing.  (Otherwise, simply read the entire input flow solution.)
!
!        >  Read one block of the input flow grid at at time and determine the data ranges.  Rewind the file.
!
!        >  Do likewise for the target grid block(s), and determine their overall data range.  Rewind the file.
!
!        >  Determine which flow blocks to include in the search tree using a small safety margin.  (The main worry is with target
!           points that might be outside the input flow grid.  The searching then produces the point with the shortest orthogonal
!           projected distance to an outer boundary cell face.  We cannot be certain that a block suppressed via data range
!           information doesn't actually contain the best choice of cell for an "orphan" point.)
!
!        >  Set up a "grid_type" array suited to the blocks to be included, then read only those blocks.  The remaining steps are
!           the same as for the no-suppression case.
!
!     >  For each block in the target grid:
!
!        >  Read the target grid block.
!
!        >  If this is the first target block, perform some initialization:
!
!           >  If the solution is cell-centered but the grid is not, interpolate it to the vertices of the solution grid (in-place).
!              (See the above clarification about halo cells.)
!
!           >  Build the solution grid search tree from all volume cells of all (unsuppressed) blocks.  (Method 1 uses bounding
!              box techniques for the searching; method 2 builds its tree from the input solution grid points and works with just
!              distances, while method 3 builds its tree from the centroids of the input solution grid cells.)
!
!        >  For each point in the target block:
!
!           >  Locate the nearest solution cell point by searching the ADT, or just the nearest grid point if KDTREE method 2 is
!              specified, or the nearest cell centroid if KDTREE method 3 is specified.  For method 3, refine the search by
!              calculating the nearest point of the closest-centroid cell to the target point (but this cell may not be the best
!              cell - working with cell centroids is not bullet-proof the way the ADT method is).
!
!           >  Interpolate the flow solution to the current target point using the solution cell found and the interpolation
!              coefficients produced by the search (unless method = 2, in which case the flow associated with the nearest data
!              point is simply copied).
!
!              Note that the interpolations are "trilinear" within a single cell, but this familiar formulation is not really
!              linear because it contains nonlinear terms, but the effect is (tri)linear if the cell is perfectly rectangular.
!
!        >  If the input flow was cell-centered, meaning dimension differences of 1, convert the interpolated flow to be likewise.
!
!        >  Output the interpolated solution for the current block.
!
!  Control file format ('flow_interp.inp')
!
!     FLOW_INTERP controls for case xxx
!     ------------- INPUT SOLUTION GRID -------------
!     baseline.xyz            Solution grid file name
!     T                       Formatted? [T|F]
!     ------------- INPUT FLOW SOLUTION -------------
!     baseline.f              Flow solution file name
!     T                       Formatted? [T|F]
!     ----------------- TARGET GRID -----------------
!     densified.xyz           Input target grid file name
!     T                       Formatted? [T|F]
!     ---------- INTERPOLATED FLOW SOLUTION ---------
!     densified.f             Output flow file name
!     T                       Formatted? [T|F]
!     ------------ MISCELLANEOUS CONTROLS -----------
!     0.0001                  Distance tolerance for determining whether a target point is inside the solution grid or not
!     --------------- OPTIONAL CONTROLS -------------
!     T                       Suppress solution blocks if possible? [T|F]
!     1                       Method:  1 = plain ADT; 2 = plain KDTREE; 3 = KDTREE + nearest-centroid-cell refinement
!
!  Sponsor:
!
!     Aerothermodynamics Branch, NASA Ames Research Center
!
!  History:
!
!     03/24/04  DAS  Initial implementation, styled after GRID_MORPH_1, with argument-driven INTERP_Q portion reusable.
!     03/29/04   "   Introduced face_type structure to help identify exterior faces of blocks.
!     04/04/04   "   Made use of xyzq_io package.
!     06/15/04   "   Replaced RIPPLE3D with a volume grid version of the Alternating Digital Tree package from Stanford.
!     02/01/05   "   Added special treatment of boundary layer profiles defined as single radial line(s) in the target block(s).
!     05/12/06   "   The profile tabulation format now handles any number of functions, not 6 or fewer.
!     08/07/08   "   Todd White suggested an option to suppress solution blocks in the hope of speeding things up in some cases.
!                    He also has parallelization in mind.
!     08/08/08   "   A case involving 1.5 million target points and 2 x 8.7 million points from 2 x 47 blocks (full Shuttle) takes
!                    9.5 minutes the standard way and 8 minutes if only the relevant 21 blocks are searched.  This represents only
!                    a modest 16% reduction on a demanding case, so the standard approach has been quite efficient all along.
!                    The printout at the end of each k plane has been commented out now.
!     08/22/08   "   The format flag for the target file was being passed as that for the solution file.
!     06/06/10   "   The 08/08/08 history comment above had typos.  No other change was made.
!     06/24/13   "   Option to use KDTREE searches for nearest points (centroids) only.  This may be adequate for some purposes
!                    such as starting some flow solutions, and uses less memory.  It's about 2 orders of magnitude faster.
!     07/26/13   "   Hybrid method refines each KDTREE result by computing the best point (+ associated interpolation coefficients)
!                    within the indicated cell.  This is 2.5 - 3 x slower than plain KDTREE, but still much faster than plain ADT
!                    on large grids.  Note that, for DPLR applications, the input grid and flow are normally at cell centers,
!                    treated as vertices here.  But the hybrid KDTREE + refinement method requires the search to identify a best
!                    cell, not vertex, so we actually generate and search the centroids of the input grid cells.  This does take
!                    more memory, but not as much as the ADT search tree requires (which is 18 reals + 8 integers per cell).
!     07/30/13   "   Replaced nearest_hex_point (intended for unstructured grids) with the analogous nearest_brick_point (suited to
!                    structured grids).  Clarified the potential imprecision of the hybrid method 3.  As explained above, the
!                    hybrid method is not guaranteed to find the best cell for each search - it merely finds the cell whose
!                    centroid is closest to the target, but some other cell(s) may actually contain the target.
!     08/05/13   "   All ADT variants have been combined into a module now for distribution reasons (generic build & search calls).
!     02/22/17   "   Suppress profile tabulation (to standard output) if the target grid seems to be a set of hemisphere lines of
!                    sight (more than 100 blocks, all 1 x 1 x nk).
!     10/15/19   "   Mysterious interpolations along lines of sight for radiation calculations revealed that the hybrid method 3 is
!                    potentially seriously flawed.  At a shock envelope boundary, volume grid cells are typically close together
!                    near the shock in the off-wall direction, yet relatively large in the surface index directions. This means that
!                    cell centroids can be a poor measure of the correct cells within which to refine interpolation coefficients.
!                    Method 1 (ADT searching) has been adjusted to cope with occasional matrix singularity observed in boundary
!                    layer regions, and is strongly recommended as the preferred method where practical. The advice presented above
!                    under "Intended application" remains valid.
!
!  Author:           David Saunders, ELORET Corporation at NASA Ames Research Center, Moffett Field, CA
!                    (later with ERC, Inc. then AMA, Inc. at NASA ARC)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use grid_block_structure  ! For on structured grid block
   use xyzq_io_module        ! For PLOT3D files

   implicit none

!  Local constants:
!  !!!!!!!!!!!!!!!!

   integer, parameter ::   &
      lunctl = 1,          & ! Control file
      lunxyz = 2,          & ! Input solution grid
      luntar = 3,          & ! Input target grid
      lunq   = 4,          & ! Input flow solution; output flow solution
      luncrt = 6             ! Screen diagnostics

   character, parameter :: &
      format * 11 = 'unformatted'

!  Derived data types:
!  !!!!!!!!!!!!!!!!!!!

   type (grid_type), pointer, dimension(:) :: &
      solution, target

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: &
      i1, ib, ios, method, mi, mj, mk, n, nblocks_solution, nblocks_target, ni, nj, nk, npts, num_q

   real :: &
      dtol, time_end, time_start

   logical :: &
      cell_centered, cleanup, formatted_out, formatted_sol, formatted_tar, initialize, try_to_suppress

   character :: &
      filename * 80

!  Execution:
!  !!!!!!!!!!

!  Open the control file:

   open (lunctl, file='flow_interp.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Unable to open flow_interp.inp.'
      go to 999
   end if

!  Echo it to the output log, using filename as a buffer:

   do
      read (lunctl, '(a)', iostat=ios) filename
      if (ios < 0) exit

      write (luncrt, '(1x, a)') filename(1:len_trim(filename))
   end do

   rewind (lunctl)

!  Open the input grid and flow solution files:
!  --------------------------------------------

   read (lunctl, *)                  ! Case description
   read (lunctl, *)                  ! Header
   read (lunctl, *) filename         ! Input grid name
   read (lunctl, *) formatted_sol

   i1 = 1; if (formatted_sol) i1 = 3

   open (lunxyz, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open input solution grid: ', filename(1:len_trim(filename))
      go to 999
   end if

   read (lunctl, *)
   read (lunctl, *) filename         ! Input flow solution name
   read (lunctl, *) formatted_sol    ! Assume the same as the grid in order to use the I/O package

   open (lunq, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open input flow solution file: ', filename(1:len_trim(filename))
      go to 999
   end if

!  Open the target grid file:
!  --------------------------

   read (lunctl, *)
   read (lunctl, *) filename         ! Target grid name
   read (lunctl, *) formatted_tar

   i1 = 1; if (formatted_tar) i1 = 3

   open (luntar, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open target grid: ', filename(1:len_trim(filename))
      go to 999
   end if

!  There is no need to store the entire target grid.  Read its header:

   call xyz_header_io (1, luntar, formatted_tar, nblocks_target, target, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble allocating the target grid blocks or reading the header records.'
      go to 999
   end if

!  Read remaining controls, including any optional ones:
!  -----------------------------------------------------

   read (lunctl, *)
   read (lunctl, *) filename         ! Output flow solution name
   read (lunctl, *) formatted_out
   read (lunctl, *)
   read (lunctl, *) dtol             ! Distance tolerance for diagnostic purposes only

   try_to_suppress = .false.
   method = 1  ! Plain ADT

   read (lunctl, *, iostat=ios)      ! Optional inputs title line
   if (ios == 0) then
      read (lunctl, *, iostat=ios) try_to_suppress
      if (ios == 0) then
         read (lunctl, *, iostat=ios) method
      end if
   end if

   close (lunctl)

!  The input solution files may need preprocessing:
!  ------------------------------------------------

   if (.not. try_to_suppress) then   ! Allocate work-space and read all blocks of the solution grid and flow field

      call xyzq_read (lunxyz, lunq, formatted_sol, nblocks_solution, num_q, cell_centered, solution, ios)

      if (ios /= 0) go to 999

   else 

      call read_only_relevant_solution_blocks (ios)  ! Local procedure below

      cell_centered = solution(1)%mi == solution(1)%ni - 1  ! For DPLR, dimensions are equal - both grid & flow are cell-centered

      if (ios /= 0) go to 999

   end if

!  Set up the interpolated flow solution file:
!  -------------------------------------------

   i1 = 1; if (formatted_out) i1 = 3

   open (lunq, file=filename, status='unknown', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open output flow file: ', filename(1:len_trim(filename))
      go to 999
   end if

   do ib = 1, nblocks_target  ! Set up the interpolated flow block dimensions

      if (cell_centered) then
         target(ib)%mi = target(ib)%ni - 1;  target(ib)%mj = target(ib)%nj - 1;  target(ib)%mk = target(ib)%nk - 1
      else
         target(ib)%mi = target(ib)%ni;      target(ib)%mj = target(ib)%nj;      target(ib)%mk = target(ib)%nk
      end if

   end do

!  As with the target grid, there is no need to store the entire interpolated flow.  Just manage one block at a time.

   call q_header_io (2, lunq, formatted_out, nblocks_target, num_q, target, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the interpolated flow header records.'
      go to 999
   end if

!  Perform the flow field interpolation one target block at a time:
!  ----------------------------------------------------------------

   call cpu_time (time_start)

   do ib = 1, nblocks_target

      initialize = ib == 1
      cleanup    = ib == nblocks_target

!     Set up for and read the current target grid block:

      call xyz_allocate (target(ib), ios)

      ni = target(ib)%ni
      nj = target(ib)%nj
      nk = target(ib)%nk

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating target grid.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 999
      end if

      npts = ni * nj * nk

      call xyz_block_io (1, luntar, formatted_tar, npts, target(ib)%x, target(ib)%y, target(ib)%z, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading target grid.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 999
      end if

!     Allocate space for the corresponding interpolated flow block:

      call q_allocate (target(ib), num_q, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating interpolated flow space.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 4i5)') ' Dimensions: ', ni, nj, nk, num_q
         go to 999
      end if

!     Interpolate the flow at the points (or cell centers) of the current target block.
!     Some preprocessing of the input solution blocks is performed when ib = 1.

      call interp_q (nblocks_solution, nblocks_target, num_q, solution, cell_centered, dtol, method, ib, target(ib), &
                     initialize, cleanup, luncrt)
!     -------------

!     Save results for this target block:

      mi = target(ib)%mi
      mj = target(ib)%mj
      mk = target(ib)%mk

      call q_block_io (2, lunq, formatted_out, num_q, mi, mj, mk, target(ib)%q, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble writing interpolated flow solution.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 4i5)') ' Dimensions and # functions: ', mi, mj, mk, num_q
         go to 999
      end if

      deallocate (target(ib)%x, target(ib)%y, target(ib)%z, target(ib)%q)

   end do ! Next target block

   call cpu_time (time_end)

   write (luncrt, '(/, a, f9.2)') ' Total run time, seconds:', time_end - time_start

   close (luntar)
   close (lunq)

!  Done.

   999 continue  ! Avoid system-dependent STOP behavior

!  Local procedures for program flow_interp:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine find_block_data_range (gridname, ib, block)

!     Determine the data range of the given grid block.  Tabulate results.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      character, intent (in)           :: gridname * (*)  ! Grid name, for tabulation purposes
      integer, intent (in)             :: ib              ! Block number, for tabulation purposes
      type (grid_type), intent (inout) :: block           ! Input with x,y,z fields defined; output with %xmax, %xmin, ... defined.

!     Local constants:

      real, parameter :: big = 1.e+32

!     Local variables:

      integer :: i, j, k

!     Execution:

      block%xmax = -big;  block%xmin = big
      block%ymax = -big;  block%ymin = big
      block%zmax = -big;  block%zmin = big

      do k = 1, block%nk
         do j = 1, block%nj
            do i = 1, block%ni
               block%xmax = max (block%x(i,j,k), block%xmax);  block%xmin = min (block%x(i,j,k), block%xmin)
               block%ymax = max (block%y(i,j,k), block%ymax);  block%ymin = min (block%y(i,j,k), block%ymin)
               block%zmax = max (block%z(i,j,k), block%zmax);  block%zmin = min (block%z(i,j,k), block%zmin)
            end do
         end do
      end do

      if (ib == 1) then
         i = len_trim (gridname)
         write (luncrt, '(/, 3a)') ' Data ranges of the ', gridname(1:i), ' grid blocks:', &
            ' Block            Xmin            Xmax              Ymin            Ymax              Zmin            Zmax'
      end if

      write (luncrt, '(i4, es18.7, es16.7, es18.7, es16.7, es18.7, es16.7)') &
         ib, block%xmin, block%xmax, block%ymin, block%ymax, block%zmin, block%zmax

      end subroutine find_block_data_range

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_only_relevant_solution_blocks (ios)

!     Scan the target and solution grids in the hope of suppressing some of the solution blocks; return only the relevant blocks.
!     The files are assumed to be open already, and the target grid header has been read.  The target file is rewound and left in
!     the same state upon return as it was upon entry.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Argument:

      integer, intent (out) :: ios  ! 0 means no problem

!     Local constants:

      real, parameter :: big = 1.e+32

!     Local variables:

      integer :: i, ib, ibkeep, ibomit, j, k;  integer, allocatable :: ib_suppress(:)
      real    :: dmargin, xmax_tar, xmin_tar, ymax_tar, ymin_tar, zmax_tar, zmin_tar
      logical :: suppress

!     Execution:

      dmargin  = dtol * 100.   ! Err on the safe side for testing overlapping of blocks

      xmax_tar = -big;  xmin_tar = big
      ymax_tar = -big;  ymin_tar = big
      zmax_tar = -big;  zmin_tar = big

!     Determine the overall data range of the entire target grid:

      do ib = 1, nblocks_target

         call xyz_allocate (target(ib), ios)

         ni = target(ib)%ni;  nj = target(ib)%nj;  nk = target(ib)%nk

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating target grid block #', ib, ';  status:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 999
         end if

         npts = ni * nj * nk

         call xyz_block_io (1, luntar, formatted_tar, npts, target(ib)%x, target(ib)%y, target(ib)%z, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading target grid block #', ib, ';  I/O status:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 999
         end if

         call find_block_data_range ('target', ib, target(ib))

         deallocate (target(ib)%x, target(ib)%y, target(ib)%z)

         xmax_tar = max (target(ib)%xmax, xmax_tar);  xmin_tar = min (target(ib)%xmin, xmin_tar)
         ymax_tar = max (target(ib)%ymax, ymax_tar);  ymin_tar = min (target(ib)%ymin, ymin_tar)
         zmax_tar = max (target(ib)%zmax, zmax_tar);  zmin_tar = min (target(ib)%zmin, zmin_tar)

      end do

      rewind (luntar)  ! Ready for processing one block at a time (no need to store all target blocks)

      deallocate (target)

      call xyz_header_io (1, luntar, formatted_tar, nblocks_target, target, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' *** read_only_relevant_blocks:  Trouble reallocating the target grid blocks or reading the header records.'
         go to 999
      end if

!     -----------------------------------------------------------------------------------------------------------

!     Read the solution grid blocks one at a time to determine which ones are outside the target grid data range:

      call xyz_header_io (1, lunxyz, formatted_sol, nblocks_solution, solution, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' *** read_only_relevant_blocks:  Trouble allocating the solution grid blocks or reading the header records.'
         go to 999
      end if

      call q_header_io (1, lunq, formatted_sol, nblocks_solution, num_q, solution, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' *** read_only_relevant_blocks:  Trouble allocating the solution flow blocks or reading the header records.'
         go to 999
      end if

      ibkeep = 0;  ibomit = 0;  allocate (ib_suppress(nblocks_solution))

      do ib = 1, nblocks_solution

         call xyz_allocate (solution(ib), ios)

         ni = solution(ib)%ni;  nj = solution(ib)%nj;  nk = solution(ib)%nk

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating solution grid block #', ib, ';  status:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 999
         end if

         npts = ni * nj * nk

         call xyz_block_io (1, lunxyz, formatted_sol, npts, solution(ib)%x, solution(ib)%y, solution(ib)%z, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading solution grid block #', ib, ';  I/O status:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 999
         end if

         call find_block_data_range ('solution', ib, solution(ib))

         suppress = (xmin_tar - solution(ib)%xmax > dmargin) .or. (solution(ib)%xmin - xmax_tar > dmargin) .and. &
                    (ymin_tar - solution(ib)%ymax > dmargin) .or. (solution(ib)%ymin - ymax_tar > dmargin) .and. &
                    (zmin_tar - solution(ib)%zmax > dmargin) .or. (solution(ib)%zmin - zmax_tar > dmargin)

         if (suppress) then
            deallocate (solution(ib)%x, solution(ib)%y, solution(ib)%z)
            ibomit = ibomit + 1;  ib_suppress(ibomit) = ib
         else
            ibkeep = ibkeep + 1
            if (ibkeep < ib) then  ! Shift the grid block up in the solution array
               solution(ibkeep)%ni = ni;  solution(ibkeep)%nj = nj;  solution(ibkeep)%nk = nk
               call xyz_allocate (solution(ibkeep), ios)
               if (ios /= 0) then
                  write (luncrt, '(/, a, i4, a, i4)') ' Trouble reallocating solution grid block #', ib, ';  I/O status:', ios
                  write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
                  go to 999
               end if
               solution(ibkeep)%x = solution(ib)%x;  solution(ibkeep)%y = solution(ib)%y;  solution(ibkeep)%z = solution(ib)%z
               deallocate (solution(ib)%x, solution(ib)%y, solution(ib)%z)
            end if
         end if

!        The solution block flow data have to be read whether the block is being retained or not.

         call q_allocate (solution(ib), num_q, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating solution flow space for block #', ib, ';  I/O status:', ios
            write (luncrt, '(a, 4i5)') ' Dimensions: ', ni, nj, nk, num_q
            go to 999
         end if

         mi = solution(ib)%mi;  mj = solution(ib)%mj;  mk = solution(ib)%mk

         call q_block_io (1, lunq, formatted_sol, num_q, mi, mj, mk, solution(ib)%q, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading flow solution block #', ib, ';  I/O status:', ios
            write (luncrt, '(a, 4i5)') ' Dimensions and # functions: ', mi, mj, mk, num_q
            go to 999
         end if

         if (suppress) then
            deallocate (solution(ib)%q)
         else
            if (ibkeep < ib) then
               solution(ibkeep)%mi = mi; solution(ibkeep)%mj = mj; solution(ibkeep)%mk = mk
               call q_allocate (solution(ibkeep), num_q, ios)
               if (ios /= 0) then
                  write (luncrt, '(/, a, i4, a, i4)') &
                     ' Trouble reallocating solution flow space for block #', ibkeep, ';  I/O status:', ios
                  write (luncrt, '(a, 4i5)') ' Dimensions: ', ni, nj, nk, num_q
                  go to 999
               end if
               solution(ibkeep)%q  = solution(ib)%q
               deallocate (solution(ib)%q)
            end if
         end if

      end do  ! Next solution block

      if (ibkeep == nblocks_solution) then
         write (luncrt, '(/, a, i4, a)') 'All solution blocks are being retained for the searching:', ibkeep, ' of them.'
      else
         nblocks_solution = ibkeep
         write (luncrt, '(/, a, i4, a)') 'Number of blocks to be retained for the searching:', ibkeep, '   Block numbers omitted:'
         write (luncrt, '(20i4)') ib_suppress(1:ibomit)
      end if

      deallocate (ib_suppress)
      close (lunxyz)
      close (lunq)

 999  return

      end subroutine read_only_relevant_solution_blocks

   end program flow_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine interp_q (nblocks_solution, nblocks_target, num_q, solution, cell_centered, dtol, method, iblock_target, &
                        target_block, initialize, cleanup, luncrt)
!
!  This routine interpolates the given flow solution onto the indicated block of a target grid.
!  It is the argument-driven, reusable portion of program flow_interp, q.v.  Briefly, it handles cell-centered or
!  vertex-centered flow solutions, assumes the surface grids match very closely, and treats one target block at a time
!  so the calling program can avoid storing more than that if there is no real need to do so.
!
!  03/23/04  DAS  Initial implementation.
!  03/25/04   "   Let TRILINT (called by RIPPLE3D) iterate more by raising its TOOBIG measure on p,q,r greatly.
!                 Also, introduced close_point in place of nearest_face_point (n, n), using 4 subdivisions.
!  03/26/04   "   Use the result from the previous i if possible.
!  03/29/04   "   Replaced nearest_face_point with nearest_boundary_point (exterior faces only).
!  04/01/04   "   If the search result is barely outside the solution grid, RIPPLE3D's result isn't easily usable because it
!                 points to the lower left cell vertex and cannot point to an outer boundary point.
!  04/07/04   "   The previous (ic,jc,kc) was being lost in the loop over solution blocks when the new target is outside a block.
!  06/15/04   "   Replaced RIPPLE3D with a volume grid version of the Alternating Digital Tree package from Stanford.
!  02/01/05   "   Added special treatment of boundary layer profiles or lines of sight defined as single radial line(s) in the
!                 target block(s).
!  08/07/08   "   The data range calculations have been suppressed following the option at the higher level to look for solution
!                 blocks that are out of the target block range and can therefore be omitted from the search tree.
!  06/24/13   "   Method 2 option to use KDTREE searches (nearest grid point only).
!  07/26/13   "   Method 3 option to refine each KDTREE result with a call to nearest_hex_point on the relevant cell.
!                 This requires searching cell centroids (generated here), not cell vertices.
!  07/30/13   "   The nearest centroid isn't always that of the best cell.  Switch from nearest_hex_point to a nearest_brick_point
!                 variant, with possible application to more than one nearest neighbor in mind.
!  08/05/13   "   All ADT variants have been combined into a module now for distribution reasons (generic build & search calls).
!  02/22/17   "   Added nblocks_target argument in order to suppress profile tabulations for excessive numbers of profiles that
!                 are probably the many lines in a hemisphere lines of sight dataset.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA (later with ERC and AMA, Inc./NASA ARC).
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use grid_block_structure  ! Derived data type for one structured grid block
   use adt_utilities         ! ADT search package (all variants)
   use xyzq_io_module        ! For PLOT3D file I/O
   use face_info_module      ! Local module
   use kdtree2_module        ! Locates nearest 1 or more neighbors in k-dimensional space; kdtree2 indicates (1:k,1:npts) indexing
                             ! of the list of data points, where the original kdtree package assumed C-appropriate order (1:n,1:k).
   implicit none

!  Arguments:
!  !!!!!!!!!!

   integer, intent (in)              :: nblocks_solution            ! Number of blocks in solution grid
   integer, intent (in)              :: nblocks_target              ! Needed to decide whether to suppress profile tabulations
   integer, intent (in)              :: num_q                       ! Number of flow variables
   type (grid_type), intent (inout)  :: solution(nblocks_solution)  ! Solution grid and flow blocks input;
                                                                    ! block data ranges are calculated here if initialize = true;
                                                                    ! flow values may be changed (to vertex values if they are
                                                                    ! input as cell-centered); they are NOT transferred back if
                                                                    ! cleanup is true.
   logical, intent (in)              :: cell_centered               ! True means solution flow is initially transferred to the
                                                                    ! solution vertices, and the output flow is transferred from
                                                                    ! the target vertices to the target cell centers
   real,    intent (in)              :: dtol                        ! Distance tolerance for considering a solution cell close
   integer, intent (in)              :: method                      ! 1 => plain ADT (bounding-block-type searching + refinement);
                                                                    ! 2 => KDTREE (nearest data point only/no refinement);
                                                                    ! 3 => KDTREE + refinement of each search for interpln. coefs.
   integer, intent (in)              :: iblock_target               ! Current target block number, for the running log
   type (grid_type), intent (inout)  :: target_block                ! Current target block, input with its grid defined (only),
                                                                    ! and output with its flow field assigned via interpolation
                                                                    ! or from nearest boundary point
   logical, intent (in)              :: initialize                  ! True means allocate solution-specific work-space (normally
                                                                    ! when iblock_target = 1, but we allow for more than one set
                                                                    ! of solution blocks in the same run)
   logical, intent (in)              :: cleanup                     ! True means deallocate solution-specific work-space before
                                                                    ! returning (but the flow is always left at the vertices
                                                                    ! even if it was originally at the cell centers)
   integer, intent (in)              :: luncrt                      ! Logical unit for running log and error messages

!  Local constants:
!  ----------------

   integer, parameter :: nsuppress = 100   ! For detecting too many target profiles to bother tabulating (presumably hemisphere LOS)

   real,    parameter :: big  = 1.e+32, &
                         dpqr = 1.0,    &  ! Tolerance for sum of (p,q,r) deviations from 0. or 1.; play safe (sum isn't precise)
                         half = 0.5,    &
                         one  = 1.0,    &
                         ten  = 10.0,   &
                         two  = 2.0,    &
                         xtra = 1.e-15, &  ! Scaled by block data range; used as a round-off margin for inside/outside tests
                         zero = 0.0

   logical, parameter :: true = .true.

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: i, ib, ic, icell, ihex, iv, j, jc, jv, k, kc, kv, mi, mj, mk, n, ni, ninside, nj, nk, nout
   logical :: profile
   logical, save :: suppress
   real    :: dsqmin, dsqmax, dtolsq, extra, p, pm1, q, qm1, r, rm1, time1, time2, unused, wall_distance
   real    :: coefs(8), xvert(3,8)
   real, dimension (3) :: delta_coords, interp_coords, target_coords, wall_coords

   type (kdtree2_result), dimension (1) :: kdresult  ! Each result is 1 distance**2/index pair; package allows for 1+ nearest pts.

!  Variables saved for further calls:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   character (42),   save :: pformat
   integer,          save :: nhex, nxyz
   integer,          dimension (:,:), allocatable, save :: conn
   type (face_type), dimension (:),   allocatable, save :: faces
   real,             dimension (:),   allocatable, save :: xmax, xmin, ymax, ymin, zmax, zmin  ! Include slight round-off margins

   real,             dimension (:,:), allocatable, save :: xyz_list  ! Repacked coordinates to be searched
   type (kdtree2), pointer, save :: kdtree

!  Execution:
!  !!!!!!!!!!

!  Pre-process the solution blocks?

   if (initialize) then

      n = nblocks_solution

!     Handle any number of functions in profile tabulations:

      pformat = '(i4, i9, 4x, 3i4, 3f12.8, es15.7, nfe12.4)'
      write (pformat(35:36), '(i2)') num_q

!     But suppress tabulations if we seem to have a hemisphere lines of sight target dataset:

      suppress = nblocks_target > nsuppress

      if (cell_centered) then ! Transfer all solution blocks to the vertices

         do ib = 1, n
            call centers_vertices (1, num_q, solution(ib))  ! Local procedure below
         end do

      end if

!!!   Determine the solution block data ranges, etc.:

!!!   allocate (faces(n))

!!!   call find_data_ranges (n, solution, faces)  ! Local procedure below

!!!   Pad the solution block data ranges slightly to allow for round-off?  (This may not gain anything.  Only on symmetry plane?)

!!!   allocate (xmax(n), xmin(n), ymax(n), ymin(n), zmax(n), zmin(n))

!!!   do ib = 1, n
!!!      extra    = (solution(ib)%xmax - solution(ib)%xmin) * xtra
!!!      xmax(ib) =  solution(ib)%xmax + extra
!!!      xmin(ib) =  solution(ib)%xmin - extra
!!!      extra    = (solution(ib)%ymax - solution(ib)%ymin) * xtra
!!!      ymax(ib) =  solution(ib)%ymax + extra
!!!      ymin(ib) =  solution(ib)%ymin - extra
!!!      extra    = (solution(ib)%zmax - solution(ib)%zmin) * xtra
!!!      zmax(ib) =  solution(ib)%zmax + extra
!!!      zmin(ib) =  solution(ib)%zmin - extra
!!!   end do
 
!!!   Mark solution block faces as exterior or interior (now informative only):

!!!   call find_exterior_faces (n, solution, faces)  ! Local procedure below

!     Build an Alternating Digital Tree for all cells of all solution blocks.
!     First, count the cells so the "connectivity" array (now just block number and (i,j,k) for each cell) can be allocated:

      nhex = 0;  nxyz = 0
      do ib = 1, n
         ni = solution(ib)%ni;  nj = solution(ib)%nj;  nk = solution(ib)%nk
         nhex = nhex + (ni - 1)*(nj - 1)*(nk - 1)
         nxyz = nxyz + ni*nj*nk
      end do

      call cpu_time (time1)

      if (method == 1) then  ! ADT searches & trilinear interpolations

         allocate (conn(4,nhex))

         call build_adt (n, solution, nhex, conn, unused)

      else  ! KDTREE searches/no interpolations

!        Set up all search grid points (or search grid cell centroids if method = 3) as a list of 3-vectors:

         call deconstruct_search_grid ()  ! Internal procedure below

         kdtree => kdtree2_create (xyz_list, sort=.false., rearrange=.false.)

      end if

      call cpu_time (time2)

      write (luncrt, '(/, a, 2i10, f8.2)') ' Search tree # points, # cells, and tree build time, s:', nxyz, nhex, time2 - time1

   end if        ! End of initialization


   dtolsq  = dtol * dtol
   ninside = 0   ! # target block points found inside a solution cell
   nout    = 0   ! # such points found more than the distance tolerance from the nearest solution cell

   ni = target_block%ni
   nj = target_block%nj
   nk = target_block%nk

   profile = ni == 1 .and. nj == 1 ! A single radial line gets special treatment for wind tunnel data comparison or lines of sight
   profile = profile .and. .not. suppress  ! ... unless there are so many lines of sight that the output is excessive

   if (profile) then
      write (luncrt, '(/, a, i4, /, a)') '! Profile #:', iblock_target, &
         '!  k  Soln. block  ic  jc  kc           p           q           r  Wall distance  Interpolated flow'
   else
      write (luncrt, '(/, a, i4, a, 3i5, a, i8)') &
         ' Beginning target grid block', iblock_target, '.  Dimensions:', ni, nj, nk, '  Total:', ni * nj * nk
   end if

   dsqmax = zero

   call cpu_time (time1)

!  For each point of the current target block:

   if (method == 1) then  ! ADT search and interpolations

      do k = 1, nk
         do j = 1, nj
            do i = 1, ni

               target_coords(1) = target_block%x(i,j,k)
               target_coords(2) = target_block%y(i,j,k)
               target_coords(3) = target_block%z(i,j,k)

               call search_adt (target_coords, ihex, p, q, r, dsqmin, true, nblocks_solution, solution, nhex, conn, interp_coords)

!!!            delta_coords(:) = target_coords(:) - interp_coords(:)

!!!            WRITE (6, '(a, 4i4, a, i6, 3f12.8, a, 3es15.7)') &
!!!               ' ib,i,j,k:', iblock_target, i, j, k, ' hex,pqr:', ihex, p, q, r, ' dxyz:', delta_coords

               dsqmax = max (dsqmin, dsqmax)

               if (dsqmin < dtolsq) then ! The best solution cell found was within the distance tolerance
                  ninside = ninside + 1
               else
                  nout = nout + 1
!!!               write (luncrt, '(a, 4i4, a, es13.5)') &
!!!                  ' Tolerance exceeded for target (i,j,k,ib) =', i, j, k, iblock_target, '  min. distance:', sqrt (dsqmin)
               end if

!              Use the fractional coordinates to interpolate the flow solution to the target point:

               n  = conn(1,ihex) ! Block #
               ic = conn(2,ihex) ! Lower left cell indices
               jc = conn(3,ihex)
               kc = conn(4,ihex)

               pm1 = one - p;  qm1 = one - q;  rm1 = one - r

               target_block%q(:,i,j,k) = &
                  rm1 * (qm1 * (pm1 * solution(n)%q(:,ic,jc,  kc  ) + p * solution(n)%q(:,ic+1,jc,  kc  ))  + &
                           q * (pm1 * solution(n)%q(:,ic,jc+1,kc  ) + p * solution(n)%q(:,ic+1,jc+1,kc  ))) + &
                    r * (qm1 * (pm1 * solution(n)%q(:,ic,jc,  kc+1) + p * solution(n)%q(:,ic+1,jc,  kc+1))  + &
                           q * (pm1 * solution(n)%q(:,ic,jc+1,kc+1) + p * solution(n)%q(:,ic+1,jc+1,kc+1)))

               if (profile) then
                  if (k == 1) then
                     wall_coords(1) = qm1 * (pm1 * solution(n)%x(ic,jc,  1) + p * solution(n)%x(ic+1,jc,  1)) + &
                                        q * (pm1 * solution(n)%x(ic,jc+1,1) + p * solution(n)%x(ic+1,jc+1,1))
                     wall_coords(2) = qm1 * (pm1 * solution(n)%y(ic,jc,  1) + p * solution(n)%y(ic+1,jc,  1)) + &
                                        q * (pm1 * solution(n)%y(ic,jc+1,1) + p * solution(n)%y(ic+1,jc+1,1))
                     wall_coords(3) = qm1 * (pm1 * solution(n)%z(ic,jc,  1) + p * solution(n)%z(ic+1,jc,  1)) + &
                                        q * (pm1 * solution(n)%z(ic,jc+1,1) + p * solution(n)%z(ic+1,jc+1,1))
                  end if
                  wall_distance = sqrt ((wall_coords(1) - target_coords(1))**2 + &
                                        (wall_coords(2) - target_coords(2))**2 + &
                                        (wall_coords(3) - target_coords(3))**2)
                  write (luncrt, pformat) k, n, ic, jc, kc, p, q, r, wall_distance, target_block%q(:,1,1,k) 
               end if

            end do ! Next i for this target block
         end do ! Next j for this target block

!!!      if (.not. profile) then
!!!         write (luncrt, '(a, i4, a, 3i5, a, 2i9)') &
!!!            ' Target block:', iblock_target, '   point:', i-1,j-1,k, '  in, out:', ninside, nout
!!!      end if

      end do ! Next k for this target block

!     Transfer interpolated flow back to the target cell centers?

      if (cell_centered) call centers_vertices (2, num_q, target_block)

   else  ! KDTREE search for nearest point (vertex of a cell-centered grid for method 2, or centroid of same for method 3)

      do k = 1, nk  ! Loop over target points for this target block
         do j = 1, nj
            do i = 1, ni

               target_coords(1) = target_block%x(i,j,k)
               target_coords(2) = target_block%y(i,j,k)
               target_coords(3) = target_block%z(i,j,k)

               call kdtree2_n_nearest (tp=kdtree, qv=target_coords, nn=1, results=kdresult)

               n = kdresult(1)%idx  ! Best data point (method 2) or centroid/cell (method 3)

               ib = conn(1,n) ! Block #
               ic = conn(2,n) ! Lower-left cell indices of best point found in the grid searched.
               jc = conn(3,n) ! (Vertex-centered grid for method 2; cell-centered grid for method 3.)
               kc = conn(4,n) ! These are indices in the vertex-centered grid either way.

               if (method == 2) then  ! We accept just the flow solution at the point nearest the target:

                  target_block%q(:,i,j,k) = solution(ib)%q(:,ic,jc,kc)
!!!               delta_coords(:) = target_coords(:) - xyz_list(:,n)
                  dsqmin = kdresult(1)%dis

               else  ! Method 3 refines the search to match the ADT method's result;  n is now the centroid or cell number

!                 WARNING:  This hybrid method is no longer recommended for such applications as setting up up line-of-sight data
!                 for radiation calculations.  High aspect ratio cells near the shock envelope can mean that cell centroids are a
!                 poor measure of the correct cells withini which to refine the interpolation coefficients.
!                 Conceivably, picking the best cell of the eight (or fewer) cells common to a best vertex found via a KDTREE
!                 search would still prove competitive with ADT searching for large grids, but such a possibility has not been
!                 implemented here.

!!!               mi = solution(ib)%ni - 1  ! Solution centroid grid dimensions
!!!               mj = solution(ib)%nj - 1
!!!               mk = solution(ib)%nk - 1
!!!               icell = n                                    ! Cell centroid number in list searched by kdtree.
!!!               kv = (icell - 1)/(mj*mi) + 1                 ! These are the lower left corner indices of that cell,
!!!               jv = (icell - (kv - 1)*(mj*mi))/mi + 1       ! but we would still need to store the block # for each cell.
!!!               iv =  icell - (kv - 1)*(mj*mi) - (jv - 1)*mi ! Store all the indices instead, as originally done.

                  call nearest_brick_point (solution(ib)%ni, solution(ib)%nj, solution(ib)%nk, &
                                            solution(ib)%x,  solution(ib)%y,  solution(ib)%z, ic, jc, kc, &
                                            target_coords, interp_coords, coefs, dsqmin)
!!!               call nearest_hex_point (xvert, target_coords, interp_coords, coefs, dsqmin)

                  target_block%q(:,i,j,k) = coefs(1)*solution(ib)%q(:,ic,jc,kc)     + coefs(2)*solution(ib)%q(:,ic+1,jc,kc)   + &
                                            coefs(3)*solution(ib)%q(:,ic,jc+1,kc)   + coefs(4)*solution(ib)%q(:,ic+1,jc+1,kc) + &
                                            coefs(5)*solution(ib)%q(:,ic,jc,kc+1)   + coefs(6)*solution(ib)%q(:,ic+1,jc,kc+1) + &
                                            coefs(7)*solution(ib)%q(:,ic,jc+1,kc+1) + coefs(8)*solution(ib)%q(:,ic+1,jc+1,kc+1)
                  delta_coords(:) = target_coords(:) - interp_coords(:)
               end if

!!!            write (*, '(a, 4i4, a, 3es15.7, a, i10, a, 4i4, a, es13.6, a, 9es15.7)') 'Target ib,i,j,k:', iblock_target,i,j,k, &
!!!               ' xyztarg:', target_block%x(i,j,k), target_block%y(i,j,k), target_block%z(i,j,k), ' cell #:', n, &
!!!               ' ib,iv,jv,kv:', ib,iv,jv,kv, ' dsq:', dsqmin, ' dxyz, txyz, ixyz:', delta_coords, target_coords, interp_coords

               dsqmax = max (dsqmin, dsqmax)

               if (dsqmin < dtolsq) then ! The best solution point (or cell for method 3) found was within the distance tolerance
                  ninside = ninside + 1
               else
                  nout = nout + 1
                  write (luncrt, '(a, 4i4, a, es13.5)') &
                     ' Tolerance exceeded for target (i,j,k,ib) =', i, j, k, iblock_target, '  min. distance:', sqrt (dsqmin)
               end if

            end do ! Next i for this target block
         end do ! Next j for this target block

!!!      if (.not. profile) then
!!!         write (luncrt, '(a, i4, a, 3i5, a, 2i9)') &
!!!            ' Target block:', iblock_target, '   point:', i-1,j-1,k, '  in, out:', ninside, nout
!!!      end if

      end do ! Next k for this target block

      if (cell_centered .and. method == 3) call centers_vertices (2, num_q, target_block)  ! ? Cell-centered is unlikely

   end if  ! Branch over search method

   if (.not. profile) then
      call cpu_time (time2)
      write (luncrt, '(a, i4, a, i9, i8, a, es13.5, a, f8.2)') &
         ' Finished  target grid block', iblock_target,  ';  in, out:', ninside, nout, &
         ';  largest distance from target:', sqrt (dsqmax), '  CPU seconds:', time2 - time1
   end if

!  Clean up?  (This allows for another set of solution blocks in the same run.)

   if (cleanup) then

!!!   deallocate (faces, xmax, xmin, ymax, ymin, zmax, zmin)

   end if

   return

!  Local procedures for subroutine interp_q:
!  -----------------------------------------

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine centers_vertices (mode, num_q, block)

!     Transfer flow variables between cell centers and cell vertices for one grid block.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: mode  ! 1 = centers -> vertices;
                                    ! 2 = vertices -> centers

      integer, intent (in) :: num_q ! # flow variables

      type (grid_type)     :: block ! The %q field is updated in-place

!     Local constants:

      real, parameter :: eighth = 0.125

!     Local variables:

      integer :: i, ia, ib, j, ja, jb, k, ka, kb, m, mi, mj, mk, ni, nj, nk

      real, allocatable, dimension (:,:,:,:) :: t

!     Execution:

      ni = block%ni;  nj = block%nj;  nk = block%nk
      mi = ni - 1;    mj = nj - 1;    mk = nk - 1

      select case (mode)

         case (1) ! Centers -> vertices

            allocate (t(num_q,ni,nj,nk))
 
            do k = 1, nk
               ka = max (1, k - 1);  kb = min (k, mk)
               do j = 1, nj
                  ja = max (1, j - 1);  jb = min (j, mj)
                  do i = 1, ni
                     ia = max (1, i - 1);  ib = min (i, mi)
                     t(:,i,j,k) = (block%q(:,ia,ja,ka) + block%q(:,ia,ja,kb)  + &
                                   block%q(:,ia,jb,ka) + block%q(:,ia,jb,kb)  + &
                                   block%q(:,ib,ja,ka) + block%q(:,ib,ja,kb)  + &
                                   block%q(:,ib,jb,ka) + block%q(:,ib,jb,kb)) * eighth
                  end do
               end do
            end do

            deallocate (block%q)

            allocate (block%q(num_q,ni,nj,nk))

            block%q = t

            deallocate (t)

!           Updating block%mi/mj/mk is redundant - the interpolation will use ni/nj/nk

         case (2) ! Vertices -> centers

            allocate (t(num_q,mi,mj,mk))

            deallocate (block%q)

            allocate (block%q(num_q,mi,mj,mk))

            do k = 1, mk
               ka = k;  kb = k + 1
               do j = 1, mj
                  ja = j;  jb = j + 1
                  do i = 1, mi
                     ia = i;  ib = i + 1
                     t(:,i,j,k) = (block%q(:,ia,ja,ka) + block%q(:,ia,ja,kb)  + &
                                   block%q(:,ia,jb,ka) + block%q(:,ia,jb,kb)  + &
                                   block%q(:,ib,ja,ka) + block%q(:,ib,ja,kb)  + &
                                   block%q(:,ib,jb,ka) + block%q(:,ib,jb,kb)) * eighth
                  end do
               end do
            end do

            block%q = t

            deallocate (t)

!           Updating block%mi/mj/mk is redundant - they should already be correct

         case default

            write (luncrt, '(a, i6)') ' Centers_vertices mode argument is bad:', mode

      end select

      end subroutine centers_vertices

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine deconstruct_search_grid ()

!     Reorganize the structured multiblock grid being searched as a list of coordinates for the KDTREE package.
!     Initially, at least, we make no attempt to eliminate duplicates on block boundaries.
!     If method = 3, we actually need to search for best cells, not best data points, so we work with cell centroids.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: eighth = 0.125

!     Local variables:

      integer :: i, ia, ib, j, ja, jb, k, ka, kb, m, n

!     Execution:

      n = 0

      if (method == 2) then  ! Input vertices are searched

         allocate (xyz_list(3,nxyz), conn(4,nxyz))  ! We still need conn to grab the nearest flow data from its structured block
 
         do m = 1, nblocks_solution
            do k = 1, solution(m)%nk
               do j = 1, solution(m)%nj
                  do i = 1, solution(m)%ni
                     n = n + 1
                     xyz_list(1,n) = solution(m)%x(i,j,k)
                     xyz_list(2,n) = solution(m)%y(i,j,k)
                     xyz_list(3,n) = solution(m)%z(i,j,k)
                     conn(1,n) = m
                     conn(2,n) = i
                     conn(3,n) = j
                     conn(4,n) = k
                  end do
               end do
            end do
         end do

      else  ! Method = 3 requires searching for the nearest cell, not vertex, so we work with cell centroids

         allocate (xyz_list(3,nhex), conn(4,nhex))

         do m = 1, nblocks_solution
            do k = 1, solution(m)%nk - 1
               ka = k;  kb = k + 1
               do j = 1, solution(m)%nj - 1
                  ja = j;  jb = j + 1
                  do i = 1, solution(m)%ni - 1
                     n = n + 1
                     ia = i;  ib = i + 1
                     xyz_list(1,n) = (solution(m)%x(ia,ja,ka) + solution(m)%x(ia,ja,kb)  + &
                                      solution(m)%x(ia,jb,ka) + solution(m)%x(ia,jb,kb)  + &
                                      solution(m)%x(ib,ja,ka) + solution(m)%x(ib,ja,kb)  + &
                                      solution(m)%x(ib,jb,ka) + solution(m)%x(ib,jb,kb)) * eighth
                     xyz_list(2,n) = (solution(m)%y(ia,ja,ka) + solution(m)%y(ia,ja,kb)  + &
                                      solution(m)%y(ia,jb,ka) + solution(m)%y(ia,jb,kb)  + & 
                                      solution(m)%y(ib,ja,ka) + solution(m)%y(ib,ja,kb)  + & 
                                      solution(m)%y(ib,jb,ka) + solution(m)%y(ib,jb,kb)) * eighth
                     xyz_list(3,n) = (solution(m)%z(ia,ja,ka) + solution(m)%z(ia,ja,kb)  + &
                                      solution(m)%z(ia,jb,ka) + solution(m)%z(ia,jb,kb)  + & 
                                      solution(m)%z(ib,ja,ka) + solution(m)%z(ib,ja,kb)  + & 
                                      solution(m)%z(ib,jb,ka) + solution(m)%z(ib,jb,kb)) * eighth
                     conn(1,n) = m
                     conn(2,n) = i
                     conn(3,n) = j
                     conn(4,n) = k
                  end do
               end do
            end do
         end do

      end if

      end subroutine deconstruct_search_grid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine find_data_ranges (nblocks, grid, faces)
 
!     Determine the (x,y,z) data ranges for the blocks of the given multiblock structured grid, and for their faces.
!     Also, set up the (i,j,k) index ranges defining all faces of the blocks.  Doing this first helps.
!
!     03/29/04  David Saunders, ELORET Corporation.  Initial implementation.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:
!     !!!!!!!!!!

      integer, intent (in)              :: nblocks        ! Number of blocks in the grid

      type (grid_type), intent (inout)  :: grid(nblocks)  ! Grid blocks, input with (x,y,z) fields defined, output with overall
                                                          ! data ranges defined
      type (face_type), intent (out)    :: faces(nblocks) ! Face info. for the grid blocks

!     Local constants:
!     !!!!!!!!!!!!!!!!

      real, parameter :: big = 1.e+32

!     Local variables:
!     !!!!!!!!!!!!!!!!

      integer :: i, j, k, ib, if, ni, nj, nk

!     Execution:

!     Establish the index ranges defining faces 1-6:

      do ib = 1, nblocks
         ni = grid(ib)%ni
         nj = grid(ib)%nj
         nk = grid(ib)%nk

         faces(ib)%i1(:) = 1;   faces(ib)%i1(2) = ni
         faces(ib)%i2(:) = ni;  faces(ib)%i2(1) = 1
         faces(ib)%j1(:) = 1;   faces(ib)%j1(4) = nj
         faces(ib)%j2(:) = nj;  faces(ib)%j2(3) = 1
         faces(ib)%k1(:) = 1;   faces(ib)%k1(6) = nk
         faces(ib)%k2(:) = nk;  faces(ib)%k2(5) = 1

         grid(ib)%xmin     =  big;  grid(ib)%ymin     =  big;  grid(ib)%zmin     =  big
         grid(ib)%xmax     = -big;  grid(ib)%ymax     = -big;  grid(ib)%zmax     = -big
         faces(ib)%xmin(:) =  big;  faces(ib)%ymin(:) =  big;  faces(ib)%zmin(:) =  big
         faces(ib)%xmax(:) = -big;  faces(ib)%ymax(:) = -big;  faces(ib)%zmax(:) = -big
      end do 

!     Calculate the data ranges of the faces:

      do ib = 1, nblocks
         do if = 1, 6
            do k = faces(ib)%k1(if), faces(ib)%k2(if)
               do j = faces(ib)%j1(if), faces(ib)%j2(if)
                  do i = faces(ib)%i1(if), faces(ib)%i2(if)
                     faces(ib)%xmin(if) = min (grid(ib)%x(i,j,k), faces(ib)%xmin(if))
                     faces(ib)%xmax(if) = max (grid(ib)%x(i,j,k), faces(ib)%xmax(if))
                     faces(ib)%ymin(if) = min (grid(ib)%y(i,j,k), faces(ib)%ymin(if))
                     faces(ib)%ymax(if) = max (grid(ib)%y(i,j,k), faces(ib)%ymax(if))
                     faces(ib)%zmin(if) = min (grid(ib)%z(i,j,k), faces(ib)%zmin(if))
                     faces(ib)%zmax(if) = max (grid(ib)%z(i,j,k), faces(ib)%zmax(if))
                  end do
               end do
            end do
         end do
      end do

      write (*, '(/, a)') &
         ' Data ranges of the solution grid blocks:', &
         ' Block  Face            Xmin            Xmax              Ymin            Ymax              Zmin            Zmax'
      do ib = 1, nblocks
         write (*, '(i4, i6, es18.7, es16.7, es18.7, es16.7, es18.7, es16.7)') &
            (ib, if, faces(ib)%xmin(if), faces(ib)%xmax(if), faces(ib)%ymin(if), faces(ib)%ymax(if), &
                     faces(ib)%zmin(if), faces(ib)%zmax(if), if = 1, 6)
      end do

!     Deduce the data ranges of the blocks:

      do ib = 1, nblocks
         do if = 1, 6
            grid(ib)%xmin = min (faces(ib)%xmin(if), grid(ib)%xmin)
            grid(ib)%xmax = max (faces(ib)%xmax(if), grid(ib)%xmax)
            grid(ib)%ymin = min (faces(ib)%ymin(if), grid(ib)%ymin)
            grid(ib)%ymax = max (faces(ib)%ymax(if), grid(ib)%ymax)
            grid(ib)%zmin = min (faces(ib)%zmin(if), grid(ib)%zmin)
            grid(ib)%zmax = max (faces(ib)%zmax(if), grid(ib)%zmax)
         end do
      end do

      write (*, '(/, a)') &
         ' Block            Xmin            Xmax              Ymin            Ymax              Zmin            Zmax'
      do ib = 1, nblocks
         write (*, '(i4, 3(es18.7, es16.7))') &
            ib, grid(ib)%xmin, grid(ib)%xmax, grid(ib)%ymin, grid(ib)%ymax, grid(ib)%zmin, grid(ib)%zmax
      end do

      end subroutine find_data_ranges

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine find_exterior_faces (nblocks, grid, faces)

!     Determine whether the faces of the given multiblock structured grid are interior or exterior.
!     The (x,y,z) extremes of each face of each block are compared with the face extremes of the other blocks.
!
!     03/29/04  David Saunders, ELORET Corporation.  Initial implementation.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:
!     !!!!!!!!!!

      integer, intent (in)              :: nblocks        ! Number of blocks in the grid

      type (grid_type), intent (in)     :: grid(nblocks)  ! Grid blocks, input with x/y/z/min/max fields defined

      type (face_type), intent (inout)  :: faces(nblocks) ! Face info. for the grid blocks; all fields except %exterior are input

!     Local constants:
!     !!!!!!!!!!!!!!!!

      real, parameter :: eps = 1.e-6                      ! (Overall data scale) * eps is used as an (x,y,z) comparison tolerance

!     Local variables:
!     !!!!!!!!!!!!!!!!

      integer :: ib, ic, if, ig
      real    :: tol, xmin, xmax, ymin, ymax, zmin, zmax

!     Execution:
!     !!!!!!!!!!

!     Determine the overall data range:

      xmin = grid(1)%xmin;  xmax = grid(1)%xmax
      ymin = grid(1)%ymin;  ymax = grid(1)%ymax
      zmin = grid(1)%zmin;  zmax = grid(1)%zmax

      do ib = 1, nblocks
         xmin = min (grid(ib)%xmin, xmin);  xmax = max (grid(ib)%xmax, xmax)
         ymin = min (grid(ib)%ymin, ymin);  ymax = max (grid(ib)%ymax, ymax)
         zmin = min (grid(ib)%zmin, zmin);  zmax = max (grid(ib)%zmax, zmax)
         faces(ib)%exterior(:) = .true.
      end do

      write (*, '(/, (2x, a, es16.7, es18.7))') &
         'Overall Xmin & Xmax:', xmin, xmax, '        Ymin & Ymax:', ymin, ymax, '        Zmin & Zmax:', zmin, zmax

      tol = eps * max (xmax - xmin, ymax - ymin, zmax - zmin)

      write (*, '(/, a, es16.7)') ' Common face data-range tolerance:', tol

!     Compare block faces to find interior faces:

      do ib = 1, nblocks
         do if = 1, 6
            do ic = 1, nblocks
               if (ic == ib) cycle

               do ig = 1, 6
                  if (abs (faces(ic)%xmin(ig) - faces(ib)%xmin(if)) < tol .and. &
                      abs (faces(ic)%xmax(ig) - faces(ib)%xmax(if)) < tol .and. &
                      abs (faces(ic)%ymin(ig) - faces(ib)%ymin(if)) < tol .and. &
                      abs (faces(ic)%ymax(ig) - faces(ib)%ymax(if)) < tol .and. &
                      abs (faces(ic)%zmin(ig) - faces(ib)%zmin(if)) < tol .and. &
                      abs (faces(ic)%zmax(ig) - faces(ib)%zmax(if)) < tol) then
                     faces(ib)%exterior(if) = .false.
                  end if
               end do

            end do
         end do
      end do

      write (*, '(/, a)') ' Block face summary:  T = exterior;  F = interior', ' Block  Face 1  2  3  4  5  6'
      write (*, '(i4, 7x, 6l3)') (ib, faces(ib)%exterior(:), ib = 1, nblocks)
      write (*, '(a)')

      end subroutine find_exterior_faces

   end subroutine interp_q
