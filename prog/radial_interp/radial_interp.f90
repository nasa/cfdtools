!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module edge_info_module     ! Module for program Radial_Interp

   type edge_connectivity_type ! Data structure for the patch number and its edge number that abuts the current patch edge
      integer :: ip            ! Neighboring patch #, or 0 for a symmetry plane edge, or < 0 for some other outer edge
      integer :: ie            ! Edge # of the neighboring patch if ip > 0; ie < 0 means neighboring edge runs in opposite direction
   end type edge_connectivity_type

   type edge_type                                                  ! Data structure for one structured surface patch
      integer, dimension (4) :: i1, i2, j1, j2                     ! Indices defining edges 1-4 in the usual imin, imax, ... order
      real,    dimension (4) :: xmin, xmax, ymin, ymax, zmin, zmax ! Data ranges for the edges
      real,    dimension (4) :: xcorner, ycorner, zcorner          ! Numbered counterclockwise from (i1,j1) vertex
      logical, dimension (4) :: exterior                           ! T|F; F means the edge appears to be interior
      type (edge_connectivity_type), dimension (4) :: edge_con     ! Connectivity info. for each of the 4 edges
   end type edge_type

   type pq_map_type                                                ! Date structure to save old patch (p,q) & (i,j) for new patches.
      real,    dimension (:,:,:), pointer :: pq                    ! The map is used here for edge points only, to ensure that
      integer, dimension (:,:,:), pointer :: ijn                   ! common edge results match exactly.  Trying to store just the
   end type pq_map_type                                            ! edge data invites trouble at the corners, so work with (i,j).

   end module edge_info_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program Radial_Interp

!     This program performs specialized interpolations for hypersonic grids and/or flow fields.  The radial grid lines of volume
!     grids are assumed to extend from the body to the outer boundary for all blocks, and there should be only one layer of blocks.
!     This assumption allows rapid generation of volume grids from a given volume grid and a new surface grid.  The interpolation
!     coefficients readily determined for the new surface grid points are reused to grow the corresponding volume quantities out
!     from the surface parallel to the original radial lines.
!
!     All files are PLOT3D grid and function files, formatted or unformatted.
!
!  Options:
!
!     (1) Interpolate a new volume grid from a given volume grid and a new surface grid.
!         A flow field may optionally accompany the given grid and be interpolated to the new volume grid.
!         The given grid and optional flow solution may each be either vertex-centered or cell-centered.  Three of the implied
!         permutations are relevant:
!
!            (a) Vertex-centered grid + vertex-centered solution:  flow interpolations are vertex-centered
!            (b) Vertex-centered grid + cell  -centered solution:  input flow is moved to vertices; cell-centered flow is output
!            (c) Cell  -centered grid + cell  -centered solution:  treated as for (a) but output flow is really cell-centered
!
!         The new surface grid must be input as vertex-centered.  For case (c), its cell-centered form is derived here.
!         The output blocks contain the same number of radial points as the initial grid.
!         The degenerate case (nk = 1) is meaningful if a flow solution is present.
!
!     (2) Interpolate a flow solution from one grid to another ([largely] consistent) grid.
!         The grids and flow solutions may contain either volume data or surface data (only).
!         The case of surfaces is the same as for option (1) with a flow solution present.
!         In fact, option (2) is distinguished from option (1) by nk > 1 in the new target grid.
!         In the case of volume grids, the old and new numbers of radial points may differ.
!         Two sub-options are provided:
!
!            (i)  is suited to fully consistent radial lines;
!            (ii) is more general and is quite efficient if the radial lines are reasonably consistent.
!
!         Option (2) is intended for situations such as where the new mesh has had its outer boundary moved closer to the body.
!         Again, the new target grid must be input as vertex-centered.  If necessary, its cell-centered form is derived here.
!
!  Assumptions:
!
!     >  At solid wall boundaries (k = 1), all values are surface-cell-centered (no below-the-wall values).
!     >  The same underlying geometry is represented by both surface grids.  Any local discrepancies are small and corrected for.
!     >  The new surface may be a reduced region (as for rapid interpolation of damage/repair volume grids and/or solutions).
!
!  Usage for Surface Damage/Repair Applications:
!
!     Two runs of Radial_Interp are required if a new surface grid in the damage region is the starting point:
!
!     (1) Generate a new volume grid from the old volume grid and new surface (vertex-centered; no flow interpolation), followed by
!     (2) Interpolate a new solution from the cell-centered initial flow and grid, using the same new vertex-centered input surface.
!
!     Treating the cell-centered flow separately (with the cell-centered form of the grid) is preferable to combining runs (1) and
!     (2) because the runs are almost identical on different data.  A second reason is that the ADT search package cannot handle
!     more than one surface to be searched simultaneously.
!
!     Apart from avoiding unnecessary smearing of a finite volume flow solution during interpolation, the intended inclusion
!     of "halo" cells in the cell-centered grids and solutions enables proper handling of boundary conditions, both at the
!     solid wall and at the other outer boundaries of what could well be a reduced set of new grid blocks modeling localized
!     surface damage or repair.  These halo cells can define the (frozen) boundary conditions needed for updating the flow in the
!     damage region more efficiently than treating the entire configuration, without significant loss of accuracy.
!
!  Algorithm:
!
!     The assumptions allow use of surface searches (and possibly radial line searches) in place of volume searches, although
!     the more general option (ii) for flow field interpolation still uses 3-D searches along the new radial lines.
! 
!     >  Read the full original volume (or surface) grid and associated flow solution if present.
!
!     >  If an output grid is implied, read the full new surface grid (vertex-centered) and establish common edge info.
!
!     >  If the solution is cell-centered but the grid is not, interpolate it to the grid vertices via averaging of neighbors.
!
!     >  If the grid and solution are both cell-centered, convert the target grid to cell-centered form.
!
!     >  Generate the ADT search tree from all patches of the new structured surface.
!
!     >  For each block of the new grid:
!
!        >  Read the new grid block (or copy it to the k = 1 plane if it is a surface already fully-read).
!
!        >  For each point (i,j) in the new k = 1 surface face:
!
!           >  Find the nearest original surface grid cell (ADT search).
!
!           >  Any target point further than the distance tolerance from the original surface grid may constitute erroneous usage,
!              but may also be the result of minor damage/repair; proceed anyway with the best surface cell found.
!
!           >  Interpolate the original grid and/or flow solution to the current new surface point (bilinear within a single cell).
!
!           >  If a volume grid is being generated,
!
!                 >  Apply the same interpolation coefficients to the successive k planes to produce the new radial grid line and
!                    corresponding flow if present.
!
!              Otherwise, if a flow field is being interpolated (not just a surface flow),
!
!                 >  If the radial lines are consistent,
!
!                       >  Build the new radial line using the old grid as a template.
!
!                    Otherwise (a little more generally),
!
!                       >  For each k > 1, use the old grid cell indices found for the previous k as good starting guesses for a
!                          RIPPLE3D search for the current k.
!
!                       >  Simple use of TRILINT-type volume interpolation will not suffice in the boundary layer.  Consider the
!                          case of a new grid simply twice as dense as the old grid:  a new in-between point at the surface could
!                          be significantly far away from the foot of the projection to the nearest surface cell in the old grid.
!                          This offset must be taken into account as we move off the new surface.  Thus, if x0 is the foot of the
!                          projected normal for the k = 1 point x1, then for k > 1 we really search for the adjusted point xk where
!                          xk = x(k) - (x1 - x0), not for x(k).
!
!        >  If the optional input flow was cell-centered and the grid was not, convert the interpolated flow to be cell-centered for
!           the new volume block.
!
!        >  Output the interpolated grid block and/or flow solution block.
!
!  Control file format ('radial_interp.inp')
!
!     RADIAL_INTERP controls for case xxx
!     ------------ INPUT VOLUME (or SURFACE) GRID --------------
!     baseline.g              Initial grid
!     T                       Formatted? [T|F]
!     F                       Cell-centered?  [T|F]
!     ----------- INPUT FLOW SOLUTION (IF PRESENT) -------------
!     baseline.f              Associated flow solution, or none.
!     T                       Formatted? [T|F]
!     F                       Cell-centered?  [T|F]
!     ------------ TARGET SURFACE or VOLUME GRID ---------------
!     denser.g                Input new volume grid
!     T                       Formatted? [T|F]
!     -------------- INTERPOLATED SURFACE GRID -----------------
!     interpolated.surf.g     Output surface grid (to check surface searches)
!     T                       Formatted? [T|F]
!     -------------- INTERPOLATED VOLUME GRID ------------------
!     interpolated.vol.g      Output volume grid, or none.
!     T                       Formatted? [T|F]
!     ------------- INTERPOLATED FLOW SOLUTION -----------------
!     denser.f                Output flow field, or none.
!     T                       Formatted? [T|F]
!     --------------- MISCELLANEOUS CONTROLS -------------------
!     1                       Flow interpolation method: 1 assumes consistent radial lines; 2 relaxes this
!     T                       T = allow surface cell extrapolation; F = force all surface (p,q)s into the unit square
!     0.0001                  dtol = projected distance tolerance for surface cell searches (in the same units as the grid)
!     0.01  0.15              Optional inputs ds1 and ds2fraction:  ds1 > 0. means redistribute each radial line a la OUTBOUND.
!
!     NOTE:  dtol < 0 means decay surface mismatches to zero at the outer boundary.
!
!  Sponsor:
!
!     TSA Reacting Flow Environments Branch, NASA Ames Research Center (now Aerothermodynamics Branch)
!
!  History:
!
!     04/13/04  DAS  Initial RADIAL_INTERP design:  James suggested adapting GRID_INTERP (also James's idea) as a specialized
!                    alternative to the fully general FLOW_INTERP.  It takes advantage of single-layer-of-blocks topologies.
!     04/14/04   "   Added output of interpolated surface grid as a check on the surface searching behavior.
!     04/16/04 DS/JR RIPPLE_SURF has been refined to balance the (p,q) and distance tolerances.
!     04/19/04   "   Count (p,q) deviations from the unit square and projected distances above tolerance separately;
!                    make the distance tolerance an input instead of some hard-coded fraction of the data range.
!     04/20/04  DAS  Merged the GRID_INTERP function into RADIAL_INTERP, which now produces grids and/or flow fields.
!     05/03/04   "   Added blanking capability: points near the open shuttle trailing edge were being mapped to the wrong surface.
!     05/20/04   "   Made extensions to ensure exact matches in common block faces of interpolated grids.  Very grubby!
!     05/25/04   "   Surface patch x/y/z extrema aren't always at an edge!  (Consider a wing tip cap.)
!     05/28/04   "   Replaced RIPPLE_SURF approach with ADT search scheme from Stanford.  (But RIPPLE3D is still used elsewhere.)
!     06/07/04   "   The ADT scheme now works directly with the structured grid (no need to repack the coordinates).
!     10/20/04   "   Any left-handed output blocks are now made right-handed by reversing the j indices.
!     11/24/04   "   Two new "cell-centered" inputs allow application to cell-centered flow interpolation (in a separate run
!                    from the grid interpolation) without too many complications.
!     01/19/05   "   Interpolating a flow solution with halo cells failed at the start of block 2.  The common edge point data
!                    apply to the new surface BEFORE recentering.  Not matching common points exactly for flow interpolation is OK.
!     03/24/05   "   Preserve outer boundaries precisely by decaying the adjustment of each new radial line at the surface.
!     05/02/05   "   Decaying can hurt at the outer boundary when higher resolution is specified at the surface.  Carrying any
!                    nonlinear effect of curvature at the surface all the way out can actually be preferable for smooth surfaces.i
!                    Quick fix for now; do it more carefully at the outer boundary some day.
!     06/21/06   "   dtol = -0.0001 (say - i.e., negative) can now be used to decay surface adjustments to zero at the outer
!                    boundary as is prefereable for shallow cavities treated as part of the main volume blocks.
!     01/28/11   "   Chun Tang wondered if an option could be provided to specify the spacing off the wall in the output grid.
!                    He had thin inner blocks whose outer spacing needed to be matched, but more conventional use on a single
!                    layer of blocks may benefit from such an option.  Therefore, look for ds1 and ds2fraction on an optional
!                    last line of the control file.  Optional blanking inputs expected there at one time have been removed.
!                    Initially, redistribute a vertex-centered grid only.  If the need arises for treating a flow solution
!                    similarly, it appears that (for each radial line) the cell-centered grid and flow would need to be converted
!                    to vertices, redistributed, then converted back.  Do that if it's ever called for.
!     08/07/13   "   All ADT variants are now in a single module with generic build_adt and search_adt interfaces.
!
!  Author:
!
!     David Saunders, ELORET Corporation/NASA Ames Research Center, Moffett Field, CA  (now ERC, Inc./ARC)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use grid_block_structure  ! Derived data type for one grid block
   use edge_info_module      ! Local module for ensuring matches at block edges
   use xyzq_io_module        ! PLOT3D file I/O package

   implicit none

!  Local constants:
!  !!!!!!!!!!!!!!!!

   integer, parameter ::   &
      lunctl  = 1,         & ! Control file
      lunxyz  = 2,         & ! Initial grid; input new grid
      lunq    = 3,         & ! Initial flow solution; output flow solution if present
      lunsrf  = 4,         & ! Output interpolated surface grid for checking searches
      luncrt  = 6,         & ! Screen diagnostics
      lunvol  = 7            ! Output volume grid if present

   real,    parameter ::   &
      zero    = 0.

   logical, parameter ::   &
      false   = .false.,   &
      true    = .true.

   character, parameter :: &
      format * 11 = 'unformatted', &
      none   *  4 = 'none'

!  Local data structures:
!  !!!!!!!!!!!!!!!!!!!!!!

   type (grid_type), pointer, dimension(:) :: &
      old_grid, new_grid, new_surf, interp_surf

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: &
      i, i1, ib, ios, j, lunflow, method, mi, mj, mk, n, nblocks_old, nblocks_new, ni, nj, nk, npts, num_q

   real :: &
      ds1, ds2fraction, dtol, overall_data_range

   logical :: &
      cell_centered, centered_flow, centered_grid, cleanup, extrapolate, formatted, formatted_flow, formatted_grid, formatted_srf, &
      initialize, input_f, output_f, output_g, recenter_f, recenter_g

   character :: &
      filename * 80, filename_target * 80

   type (edge_type), allocatable, dimension(:) :: &
      new_edges

!  Execution:
!  !!!!!!!!!!

!  Open the control file:

   open (lunctl, file='radial_interp.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Unable to open radial_interp.inp.'
      go to 999
   end if

!  Echo it to the output log, using filename as a buffer:

   do
      read (lunctl, '(a)', iostat=ios) filename
      if (ios < 0) exit

      write (luncrt, '(1x, a)') trim (filename)
   end do

   rewind (lunctl)

!  Deal with the old grid and (optional) flow solution files:
!  ----------------------------------------------------------

   read (lunctl, *)                  ! Case description
   read (lunctl, *)                  ! Header
   read (lunctl, *) filename         ! Input grid name
   read (lunctl, *) formatted
   read (lunctl, *) centered_grid

   i1 = 1; if (formatted) i1 = 3

   open (lunxyz, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open the initial grid file: ', trim (filename)
      go to 999
   end if

   read (lunctl, *)                  ! Header
   read (lunctl, *) filename         ! Input flow solution name, or none
   read (lunctl, *) ! formatted      ! Assume the same as the grid in order to use the I/O package
   read (lunctl, *) centered_flow

   input_f = filename(1:4) /= none

   if (input_f) then
      lunflow = lunq
      open (lunq, file=filename, status='old', form=format(i1:11), iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open the initial flow solution file: ', trim (filename)
         go to 999
      end if
   else
      lunflow    = -lunq
      num_q      = 0
      recenter_f = false ! Should be set by xyzq_read, but isn't
   end if

   if (centered_grid .and. .not. centered_flow) then
      write (luncrt, '(/, a)') ' A cell-centered flow must accompany a cell-centered grid.'
      go to 999
   end if

   recenter_g = centered_grid .and. centered_flow ! The new target grid will be moved from vertices to cell centers

!  Allocate and read the initial grid and optional flow field:

   call xyzq_read (lunxyz, lunflow, formatted, nblocks_old, num_q, recenter_f, old_grid, ios)

   if (ios /= 0) go to 999

   if (recenter_f) then ! An accompanying flow solution was found to be cell-centered (dimensions one less than the grid)
      if (centered_grid .or. .not. centered_flow) then
         write (luncrt, '(/, a)') ' Flawed controls.  Vertex-centered grid/cell-centered flow assumed from input block dimensions.'
         centered_grid = false;  centered_flow = true
      end if
   end if

   nk = old_grid(1)%nk ! nk = 1 is OK
   do ib = 2, nblocks_old
      if (old_grid(ib)%nk /= nk) then
         write (luncrt, '(/, a, 3i5)') ' Input grid blocks do not have a common number of radial lines:', ib, old_grid(ib)%nk, nk
         ios = 1
      end if
   end do
   if (ios /= 0) go to 999

!  Deal with opening the target grid file:
!  ---------------------------------------

   read (lunctl, *)
   read (lunctl, *) filename_target
   read (lunctl, *) formatted

   i1 = 1; if (formatted) i1 = 3

   open (lunxyz, file=filename_target, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open the target grid: ', trim (filename_target)
      go to 999
   end if

!  Pitfall:  Before the surface connectivity complication, we read the target header here.
!  Blanking requires its number of blocks.  Blanking is no longer an option, but leave this read for now.

   if (formatted) then
      read (lunxyz, *) nblocks_new
   else
      read (lunxyz) nblocks_new
   end if

   rewind (lunxyz)

!  Deal with opening the output interpolated surface grid file (which is for diagnostic use):
!  ------------------------------------------------------------------------------------------

   read (lunctl, *)
   read (lunctl, *) filename
   read (lunctl, *) formatted_srf

   i1 = 1; if (formatted_srf) i1 = 3

   open (lunsrf, file=filename, status='unknown', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open the interpolated surface grid file: ', trim (filename)
      go to 999
   end if

!  Deal with opening the output volume grid (which may have to be suppressed):
!  ---------------------------------------------------------------------------

   read (lunctl, *)
   read (lunctl, *) filename
   read (lunctl, *) formatted_grid

   i1 = 1; if (formatted_grid) i1 = 3

   output_g = filename(1:4) /= none
   if (nk == 1) output_g = .false.  ! I.e., if the old grid is just a surface

   if (output_g) then
      open (lunvol, file=filename, status='unknown', form=format(i1:11), iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open the output volume grid file: ', trim (filename)
         go to 999
      end if
   end if

!  Deal with opening the output interpolated flow file (possibly suppressed):
!  --------------------------------------------------------------------------

   read (lunctl, *)
   read (lunctl, *) filename
   read (lunctl, *) formatted_flow

   i1 = 1; if (formatted_flow) i1 = 3

   output_f = input_f .and. filename(1:4) /= none

   if (output_f) then
      open (lunq, file=filename, status='unknown', form=format(i1:11), iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open the output flow file: ', trim (filename)
         go to 999
      end if
   end if

!  Read miscellaneous control inputs:
!  ----------------------------------

   read (lunctl, *)
   read (lunctl, *) method
   read (lunctl, *) extrapolate ! No longer an option with the ADT search scheme
   read (lunctl, *) dtol        ! For surface search metrics; negate the value if surface adjustments are to be decayed to 0 at nk

!  Read optional redistribution controls:
!  --------------------------------------

   ds1 = zero
   read (lunctl, *, iostat=ios) ds1, ds2fraction

   close (lunctl)

!  Deal with reading the target grid (header or full surface or volume, which may have to be moved to cell centers):
!  -----------------------------------------------------------------------------------------------------------------

   if (output_g) then ! A volume grid is being output, and we want its interior common block faces to match exactly.

!     In order to calculate the new surface connectivity information, we have to read the entire new surface grid.
!     We can no longer read it one block at a time into the k = 1 plane of the new volume grid, as was done originally.
!     Therefore, read the entire surface into new_surf(*), then copy it one block at a time during the interpolation.

      call xyzq_read (lunxyz, -1, formatted, nblocks_new, num_q, cell_centered, new_surf, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading the target surface grid.'
         go to 999
      end if

!     We still need to transcribe the header info. to the new volume grid.
!     It is simplest to rewind and reread the new surface header with the utility that also allocates the block array.

      rewind (lunxyz)

      i1 = 1; if (formatted) i1 = 3

      open (lunxyz, file=filename_target, status='old', form=format(i1:11), iostat=ios)

      call xyz_header_io (1, lunxyz, formatted, nblocks_new, new_grid, ios)

      close (lunxyz)

!     Determine the new surface grid connectivity.  This was added to existing edge-info routines that have been merged:

      allocate (new_edges(nblocks_new))

      call get_surface_grid_info ('new', nblocks_new, new_surf, new_edges, overall_data_range)

!     However, the edge data aren't right if we're going to move the target grid to cell centers.
!     We will pass recenter_g to interp_block as a flag to skip the common-edge-point feature.

      if (recenter_g) then ! We'll be interpolating at the new surface cell centers
         do ib = 1, nblocks_new
            new_grid(ib)%ni = new_grid(ib)%ni + 1 ! Cell-centered with halos
            new_grid(ib)%nj = new_grid(ib)%nj + 1
         end do
      end if 

      new_grid(:)%nk = nk  ! Same as "old" grid when just a new surface is provided

   else ! No output volume grid.  We're reading one as a target.  Forget about exactness in the interpolated flow at common faces:

!     Just read the target grid header here, then one target block at a time later as needed for interpolation:

      call xyz_header_io (1, lunxyz, formatted, nblocks_new, new_grid, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble allocating or reading the target grid header arrays.'
         go to 999
      end if

      if (recenter_g) then ! We'll be interpolating at the new volume cell centers
         do ib = 1, nblocks_new
            new_grid(ib)%ni = new_grid(ib)%ni + 1 ! Cell-centered with halos
            new_grid(ib)%nj = new_grid(ib)%nj + 1
            new_grid(ib)%nk = new_grid(ib)%nk + 1
         end do
      end if

      allocate (new_edges(1)) ! To avoid passing an unallocated argument to interp_block, though it won't be used

   end if

!  Set up the interpolated surface grid block dimensions:
!  ------------------------------------------------------

   allocate (interp_surf(nblocks_new))

   do ib = 1, nblocks_new
      interp_surf(ib)%ni = new_grid(ib)%ni
      interp_surf(ib)%nj = new_grid(ib)%nj
      interp_surf(ib)%nk = 1
   end do

!  Write its header records:

   call xyz_header_io (2, lunsrf, formatted_srf, nblocks_new, interp_surf, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the interpolated surface grid header records.'
      go to 999
   end if

!  Write the interpolated volume grid header records?
!  --------------------------------------------------

   if (output_g) then

      call xyz_header_io (2, lunvol, formatted_grid, nblocks_new, new_grid, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the interpolated volume grid header records.'
         go to 999
      end if

   end if

!  Set up the interpolated flow field block dimensions?
!  ----------------------------------------------------

   if (output_f) then

      do ib = 1, nblocks_new
         if (recenter_f) then
            new_grid(ib)%mi = new_grid(ib)%ni - 1
            new_grid(ib)%mj = new_grid(ib)%nj - 1
            new_grid(ib)%mk = max (1, new_grid(ib)%nk - 1)
         else ! recenter_g or not, but if recenter_g, we'll have to read the new blocks into a buffer first
            new_grid(ib)%mi = new_grid(ib)%ni
            new_grid(ib)%mj = new_grid(ib)%nj
            new_grid(ib)%mk = new_grid(ib)%nk
         end if
      end do

!     Write its header records:

      call q_header_io (2, lunq, formatted_flow, nblocks_new, num_q, new_grid, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the interpolated flow header records.'
         go to 999
      end if

   end if


!  --------------------------------------------------------------------------
!  Perform the flow field and/or grid interpolations one new block at a time:
!  --------------------------------------------------------------------------

   do ib = 1, nblocks_new

      initialize = ib == 1
      cleanup    = ib == nblocks_new

      ni = new_grid(ib)%ni ! Final values, not necessarily those in the input file
      nj = new_grid(ib)%nj
      nk = new_grid(ib)%nk

      call xyz_allocate (new_grid(ib), ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating new grid block #', ib, '.  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 999
      end if

      if (output_g) then ! The new surface grid has already been read, but it may need moving to cell centers.

         if (recenter_g) then

            call recenter_surface (new_surf(ib), new_grid(ib))

         else
            do j = 1, nj
               do i = 1, ni
                  new_grid(ib)%x(i,j,1) = new_surf(ib)%x(i,j,1)
                  new_grid(ib)%y(i,j,1) = new_surf(ib)%y(i,j,1)
                  new_grid(ib)%z(i,j,1) = new_surf(ib)%z(i,j,1)
               end do
            end do
         end if

         deallocate (new_surf(ib)%x, new_surf(ib)%y, new_surf(ib)%z)

      else ! Read the current target grid block (surface or volume), but it may need moving to cell centers:

         if (recenter_g)  then

            call recenter_volume (lunxyz, formatted, new_grid(ib)) ! Reads into a local buffer

         else
            npts = ni * nj * nk

            call xyz_block_io (1, lunxyz, formatted, npts, new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z, ios)

            if (ios /= 0) then
               write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading new grid block #', ib, '.  I/O status:', ios
               write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
               go to 999
            end if
         end if

      end if

!     Allocate space for the corresponding interpolated surface grid block:

      call xyz_allocate (interp_surf(ib), ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating interpolated surface block #', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', interp_surf(ib)%ni, interp_surf(ib)%nj, interp_surf(ib)%nk
         go to 999
      end if

!     Allocate space for the optional interpolated flow block:

      if (output_f) then

         call q_allocate (new_grid(ib), num_q, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating interpolated flow space.  Block #:', ib, '  I/O status:', ios
            write (luncrt, '(a, 4i5)') ' Dimensions: ', new_grid(ib)%mi, new_grid(ib)%mj, new_grid(ib)%mk, num_q
            go to 999
         end if

      end if

!     Perform the indicated interpolations for this new grid block.
!     -------------------------------------------------------------

!     Some preprocessing of the input grid surface faces at the body is performed when ib = 1.

      call interp_block (output_g, output_f, method, dtol, ds1, ds2fraction, nblocks_old, num_q, old_grid, recenter_f, recenter_g, &
                         nblocks_new, ib, new_grid(ib), new_edges, interp_surf(ib), initialize, cleanup, luncrt)
!     -----------------

!     Save results for this new block.
!     First, the interpolated surface grid for checking:

      call xyz_block_io (2, lunsrf, formatted_srf, ni * nj, interp_surf(ib)%x, interp_surf(ib)%y, interp_surf(ib)%z, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble writing interpolated surface grid block #', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5, i9)') ' Dimensions and # pts.: ', interp_surf(ib)%ni, interp_surf(ib)%nj, interp_surf(ib)%nk, npts
         go to 999
      end if

      deallocate (interp_surf(ib)%x, interp_surf(ib)%y, interp_surf(ib)%z)

      if (output_g) then

         npts = ni * nj * nk

         call xyz_block_io (2, lunvol, formatted_grid, npts, new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble writing interpolated volume grid block #', ib, '  I/O status:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 999
         end if

         deallocate (new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z)

      end if

      if (output_f) then

         mi = new_grid(ib)%mi
         mj = new_grid(ib)%mj
         mk = new_grid(ib)%mk

         call q_block_io (2, lunq, formatted_flow, num_q, mi, mj, mk, new_grid(ib)%q, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') ' Trouble writing interpolated flow solution block #', ib, '  I/O status:', ios
            write (luncrt, '(a, 4i5)') ' Dimensions and # functions: ', mi, mj, mk, num_q
            go to 999
         end if

         deallocate (new_grid(ib)%q)

      end if

   end do ! Next new block

   if (output_g) deallocate (new_edges)


!  Done.

   999 continue

   end program Radial_Interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine get_surface_grid_info (case, nblocks, grid, edges, overall_data_range)
!
!  Determine the (x,y,z) data ranges for the k = 1 surface faces of the given multiblock structured grid, and for their edges.
!  Also, set up the (i,j) index ranges defining all edges of these faces.  Doing this first helps.
!  Determine whether the surface patch edges are interior or exterior.  This was carried over from flow_interp and may not be used.
!  In addition, establish edge connectivity information as needed to ensure identical results for matching edge point calculations.
!
!  04/09/04  David Saunders, ELORET.  Edge adaptation of block face routines in flow_interp.
!  05/19/05    "      "               Extensions to allow enforcing of exact matches on common block edges and hence common faces.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Data structures:

   use grid_block_structure  ! External module
   use edge_info_module      ! Local module

   implicit none

!  Arguments:

   character, intent (in)            :: case * 3           ! 'new' is really just to get the connectivity data;
                                                           ! 'old' is the earlier requirement to determine data range & edge status
   integer, intent (in)              :: nblocks            ! Number of blocks in the grid

   type (grid_type), intent (inout)  :: grid(nblocks)      ! Grid blocks, input with (x,y,z) fields defined, output with overall
                                                           ! data ranges defined for the k = 1 faces
   type (edge_type), intent (out)    :: edges(nblocks)     ! Edge info. for the k = 1 faces of the grid blocks

   real, intent (out)                :: overall_data_range ! For scaling tolerances

!  Local constants:

   real, parameter :: big = 1.e+32,  &
                      eps = 1.e-6                          ! (Overall data scale) * eps is used as an edge comparison tolerance
!  Local variables:

   integer :: i, i1a, i1b, i2a, j, j1a, j1b, j2a, ib, ic, ie, ig, ni, nj

   real    :: dsq1, dsq2, tol, xmin, xmax, ymin, ymax, zmin, zmax

!  Execution:

   if (case == 'new') then
      write (6, '(/, a)') ' Ascertain the new surface grid block connectivity:'
   else
      write (6, '(/, a)') ' Ascertain the old surface grid block data ranges and edge status:'
   end if

!  Establish the index ranges defining edges 1-4:

   do ib = 1, nblocks
      ni = grid(ib)%ni
      nj = grid(ib)%nj

      edges(ib)%i1(:) = 1;   edges(ib)%i1(2) = ni
      edges(ib)%i2(:) = ni;  edges(ib)%i2(1) = 1
      edges(ib)%j1(:) = 1;   edges(ib)%j1(4) = nj
      edges(ib)%j2(:) = nj;  edges(ib)%j2(3) = 1

      grid(ib)%xmin     =  big;  grid(ib)%ymin     =  big;  grid(ib)%zmin     =  big
      grid(ib)%xmax     = -big;  grid(ib)%ymax     = -big;  grid(ib)%zmax     = -big
      edges(ib)%xmin(:) =  big;  edges(ib)%ymin(:) =  big;  edges(ib)%zmin(:) =  big
      edges(ib)%xmax(:) = -big;  edges(ib)%ymax(:) = -big;  edges(ib)%zmax(:) = -big
   end do

!  Calculate the data ranges of the edges:

   do ib = 1, nblocks
      do ie = 1, 4
         do j = edges(ib)%j1(ie), edges(ib)%j2(ie)
            do i = edges(ib)%i1(ie), edges(ib)%i2(ie)
               edges(ib)%xmin(ie) = min (grid(ib)%x(i,j,1), edges(ib)%xmin(ie))
               edges(ib)%xmax(ie) = max (grid(ib)%x(i,j,1), edges(ib)%xmax(ie))
               edges(ib)%ymin(ie) = min (grid(ib)%y(i,j,1), edges(ib)%ymin(ie))
               edges(ib)%ymax(ie) = max (grid(ib)%y(i,j,1), edges(ib)%ymax(ie))
               edges(ib)%zmin(ie) = min (grid(ib)%z(i,j,1), edges(ib)%zmin(ie))
               edges(ib)%zmax(ie) = max (grid(ib)%z(i,j,1), edges(ib)%zmax(ie))
            end do
         end do
      end do
   end do

!  Calculate the data ranges of the k = 1 surface faces.  CAUTION:  Extremes aren't always at an edge.  Consider a wing tip cap.

   do ib = 1, nblocks
      grid(ib)%xmin = big;  grid(ib)%xmax = -big
      grid(ib)%ymin = big;  grid(ib)%ymax = -big
      grid(ib)%zmin = big;  grid(ib)%zmax = -big
      do j = 1, grid(ib)%nj
         do i = 1, grid(ib)%ni
            grid(ib)%xmin = min (grid(ib)%x(i,j,1), grid(ib)%xmin)
            grid(ib)%xmax = max (grid(ib)%x(i,j,1), grid(ib)%xmax)
            grid(ib)%ymin = min (grid(ib)%y(i,j,1), grid(ib)%ymin)
            grid(ib)%ymax = max (grid(ib)%y(i,j,1), grid(ib)%ymax)
            grid(ib)%zmin = min (grid(ib)%z(i,j,1), grid(ib)%zmin)
            grid(ib)%zmax = max (grid(ib)%z(i,j,1), grid(ib)%zmax)
         end do
      end do
   end do

   write (6, '(/, a)') &
      ' Data ranges of the surface grid blocks:', &
      ' k=1 Face         Xmin            Xmax              Ymin            Ymax              Zmin            Zmax'

   do ib = 1, nblocks
      write (6, '(i4, 1p, 3(e18.7, e16.7))') &
         ib, grid(ib)%xmin, grid(ib)%xmax, grid(ib)%ymin, grid(ib)%ymax, grid(ib)%zmin, grid(ib)%zmax
   end do

!  Deduce further the overall data range:

   xmin = grid(1)%xmin;  xmax = grid(1)%xmax
   ymin = grid(1)%ymin;  ymax = grid(1)%ymax
   zmin = grid(1)%zmin;  zmax = grid(1)%zmax

   do ib = 1, nblocks
      xmin = min (grid(ib)%xmin, xmin);  xmax = max (grid(ib)%xmax, xmax)
      ymin = min (grid(ib)%ymin, ymin);  ymax = max (grid(ib)%ymax, ymax)
      zmin = min (grid(ib)%zmin, zmin);  zmax = max (grid(ib)%zmax, zmax)
      edges(ib)%exterior(:) = .true.  ! Initially; reset below for matching edges; no longer used, but retained just in case
   end do

   overall_data_range = max (xmax - xmin, ymax - ymin, zmax - zmin)
   tol = eps * overall_data_range

   write (6, '(/, (2x, a, 1p, e16.7, e18.7))') &
      'Overall Xmin & Xmax:', xmin, xmax, '        Ymin & Ymax:', ymin, ymax, '        Zmin & Zmax:', zmin, zmax
   write (6, '(/, a, 1p, e16.7)') ' Overall data range:', overall_data_range

!  Compare surface block edges to find interior edges:

   do ib = 1, nblocks

      do ie = 1, 4
         edges(ib)%edge_con(ie)%ip = -1 ! Mark this edge as having no neighbor initially
         edges(ib)%edge_con(ie)%ie = -1

         do ic = 1, nblocks
            if (ic == ib) cycle

            do ig = 1, 4
               if (abs (edges(ic)%xmin(ig) - edges(ib)%xmin(ie)) < tol .and. &
                   abs (edges(ic)%xmax(ig) - edges(ib)%xmax(ie)) < tol .and. &
                   abs (edges(ic)%ymin(ig) - edges(ib)%ymin(ie)) < tol .and. &
                   abs (edges(ic)%ymax(ig) - edges(ib)%ymax(ie)) < tol .and. &
                   abs (edges(ic)%zmin(ig) - edges(ib)%zmin(ie)) < tol .and. &
                   abs (edges(ic)%zmax(ig) - edges(ib)%zmax(ie)) < tol) then
                  edges(ib)%exterior(ie)    = .false.
                  edges(ib)%edge_con(ie)%ip = ic
                  edges(ib)%edge_con(ie)%ie = ig
                  i1a = edges(ib)%i1(ie); j1a = edges(ib)%j1(ie)       ! (i,j) indices of the ends of edge ie of surface patch ib
                  i2a = edges(ib)%i2(ie); j2a = edges(ib)%j2(ie)
                  i1b = edges(ic)%i1(ig); j1b = edges(ic)%j1(ig)       ! (i,j) indices of one end  of edge ig of surface patch ic
                  dsq1 = (grid(ib)%x(i1a,j1a,1) - grid(ic)%x(i1b,j1b,1)) ** 2 + &
                         (grid(ib)%y(i1a,j1a,1) - grid(ic)%y(i1b,j1b,1)) ** 2 + &
                         (grid(ib)%z(i1a,j1a,1) - grid(ic)%z(i1b,j1b,1)) ** 2
                  dsq2 = (grid(ib)%x(i2a,j2a,1) - grid(ic)%x(i1b,j1b,1)) ** 2 + &
                         (grid(ib)%y(i2a,j2a,1) - grid(ic)%y(i1b,j1b,1)) ** 2 + &
                         (grid(ib)%z(i2a,j2a,1) - grid(ic)%z(i1b,j1b,1)) ** 2
                  if (dsq2 < dsq1) edges(ib)%edge_con(ie)%ie = -ig     ! Matching edge indices run in opposite directions
               end if
            end do
         end do
      end do
   end do

!  Flag any surface edge on a symmetry plane:

   do ib = 1, nblocks
      do ie = 1, 4
         if (edges(ib)%edge_con(ie)%ip == -1) then
            if      (abs (edges(ib)%xmin(ie)) < tol .and. abs (edges(ib)%xmax(ie)) < tol) then
               edges(ib)%edge_con(ie)%ip = 0
            else if (abs (edges(ib)%ymin(ie)) < tol .and. abs (edges(ib)%ymax(ie)) < tol) then
               edges(ib)%edge_con(ie)%ip = 0
            else if (abs (edges(ib)%zmin(ie)) < tol .and. abs (edges(ib)%zmax(ie)) < tol) then
               edges(ib)%edge_con(ie)%ip = 0
            end if
         end if
      end do
   end do

   write (6, '(/, a)') ' Surface face edge summary:  T = exterior; F = interior', ' Block  Edge 1  2  3  4'
   write (6, '(i4, 7x, 4l3)') (ib, edges(ib)%exterior(:), ib = 1, nblocks)
   write (6, '(/, a)') ' Surface face connectivity:', ' Block  Edge      1      2      3      4'
   write (6, '(i4, 8x, i4, i3, i4, i3, i4, i3, i4, i3)') &
     (ib, (edges(ib)%edge_con(ie)%ip, edges(ib)%edge_con(ie)%ie, ie = 1, 4), ib = 1, nblocks)

   end subroutine get_surface_grid_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine interp_block (output_g, output_f, method, dtol, ds1, ds2fraction, nblocks_old, num_q, old_grid,  &
                            recenter_f, recenter_g, nblocks_new, iblock_new, new_grid, new_edges, interp_surf, &
                            initialize, cleanup, luncrt)
!
!  This routine interpolates an existing grid and/or flow solution to [one block of] a new surface or volume grid, under certain
!  assumptions listed in the Radial_Interp main program, q.v.  Briefly, a hypersonic flow grid is assumed, with all radial lines
!  extending from the k = 1 surface to the outer boundary (one layer of blocks only).
!
!  As with the earlier flow_interp and interp_q, this scheme handles cell-centered or vertex-centered flow solutions, assumes the
!  surface grids match very closely, and treats one new block at a time so the calling program can avoid storing more than that
!  if there is no real need to do so.
!
!  04/13/04  DAS  Initial adaptation of ideas from grid_interp and flow_interp.
!  04/14/04   "   Added interpolated surface grid output.
!  04/15/04   "   RIPPLE_SURF has DWEIGHT and ERROR arguments now.
!  04/19/04   "   Distinguish the counts of surface searches out of tolerance w.r.t. (p,q) from those exceeding dtol;
!                 make dtol an input rather than some hard-coded fraction of the data range.
!  04/20/04   "   Merged the GRID_INTERP functionality in with the flow interpolation options, at some cost in clarity.
!  05/03/04   "   Added blanking capability.
!  05/24/04   "   Made use of connectivity data to ensure that common edge results match exactly for a new volume grid.
!  05/28/04   "   Replaced RIPPLE_SURF with the ADT surface search scheme from Stanford.
!  06/02/04   "   Introduced x/y/zcorner(*) to overcome the fact that any number of patches may have a common corner point.
!  06/07/04   "   The ADT scheme now works directly with the structured grid (no need to repack the coordinates).
!  01/19/05   "   Added recenter_g argument to enable skipping use of common edge data when flows with halo cells are treated.
!  03/24/05   "   Preserve outer boundaries precisely by decaying the adjustment at the surface.
!  05/02/05   "   Decaying can hurt at the outer boundary when higher resolution is specified at the surface.  Carrying any
!                 nonlinear effect of curvature at the surface all the way out can actually be preferable.  Quick fix for now;
!                 do it more carefully at the outer boundary some day.
!  06/21/06   "   dtol = -0.0001 (say - i.e., negative) can now be used to decay surface adjustments to zero at the outer
!                 boundary as is prefereable for shallow cavities treated as part of the main volume blocks.
!  01/28/11   "   Remnants of inactive blanking have been removed.  An option to redistribute the points along each new radial line
!                 has been added (extra step for each radial line).
!  08/07/13   "   All ADT variants are now in a single module with generic build_adt and search_adt interfaces.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, CA  (now ERC, Inc./ARC)
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use grid_block_structure  ! Derived data type for one grid block
   use edge_info_module      ! Local module for ensuring matches at block boundaries
   use adt_utilities         ! All ADT variants (generic build_adt and search_adt interfaces)

   implicit none

!  Arguments:
!  !!!!!!!!!!

   logical, intent (in)              :: output_g                   ! T means a volume grid is being generated

   logical, intent (in)              :: output_f                   ! T means a flow solution is being interpolated

   integer, intent (in)              :: method                     ! 1 means we assume the old and new radial lines are consistent;
                                                                   ! 2 means we relax this assumption (only)

   real,    intent (in)              :: dtol                       ! Distance tolerance for considering a surface quad close enough
                                                                   ! -dtol input means decay any discrepancies toward outer bndry.

   real,    intent (in)              :: ds1, ds2fraction           ! If ds1 > 0., all lines of the new grid block are redistributed
                                                                   ! according to these controls as in DPLR/OUTBOUND shock aligning

   integer, intent (in)              :: nblocks_old                ! Number of blocks in the old grid

   integer, intent (in)              :: num_q                      ! Number of variables in the associated flow solution, if any

   type (grid_type), intent (inout)  :: old_grid(nblocks_old)      ! Old grid and flow field;  k = 1 surface face data ranges
                                                                   ! are calculated here if initialize = T;
                                                                   ! flow values may be changed (to vertex values if they are
                                                                   ! input as cell-centered; they are NOT transferred back if
                                                                   ! cleanup = true (only the interpolated flow is)

   logical, intent (in)              :: recenter_f                 ! True means the flow field is initially transferred to the
                                                                   ! old grid vertices, and the output flow is transferred from
                                                                   ! the new grid vertices to the new cell centers (no halos)

   logical, intent (in)              :: recenter_g                 ! True means the starting grid and flow are both cell-centered
                                                                   ! with halos, and the target vertices have to be moved to cell
                                                                   ! centers for flow interpolation.  Skip the mechanism for
                                                                   ! identical results at common edge points in this case.

   integer, intent (in)              :: nblocks_new                ! Needed only for the surface edge connectivity info.

   integer, intent (in)              :: iblock_new                 ! Current new grid block number, for the running log

   type (grid_type), intent (inout)  :: new_grid                   ! Current new grid block, input with its header & x,y,z defined
                                                                   ! for k = 1 (at least) and output with its volume and possibly
                                                                   ! its flow field assigned via interpolation

   type (edge_type), intent (inout)  :: new_edges(nblocks_new)     ! New surface grid connectivity info. to ensure matching results
                                                                   ! at common block faces; mostly input but corner info. saved here

   type (grid_type), intent (inout)  :: interp_surf                ! Current interpolated surface block, input with its header
                                                                   ! defined; output with its (x,y,z)s assigned via interpolation

   logical, intent (in)              :: initialize                 ! True means allocate old-grid-specific work-space (normally
                                                                   ! when iblock_new = 1, but we allow for more than one set
                                                                   ! of old blocks in the same run)

   logical, intent (in)              :: cleanup                    ! True means deallocate old-grid-specific work-space before
                                                                   ! returning (but the flow is always left at the vertices
                                                                   ! even if it was originally at the cell centers)

   integer, intent (in)              :: luncrt                     ! Logical unit for running log and error messages

!  Local constants:
!  !!!!!!!!!!!!!!!!

   real,      parameter :: big  = 1.e+32,   &
                           dwt  = 0.5,      &  ! Weight for RIPPLE_SURF's distance part of the measure of closeness
                           dxyz = 1.e-6,    &  ! Scaled by overall data range; used for padding surface patch data ranges
                           half = 0.5,      &
                           one  = 1.0,      &
                           two  = 2.0,      &
                           zero = 0.0

   logical,   parameter :: false = .false., &
                           true  = .true.

   character, parameter :: lcs_method * 1 = 'M' ! Monotonic spline fits safeguard the symmetry plane during redistribution

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: i, i_use, ib, ic, icorner, idim, ie, ier, ie_use, iedge, ip, ip_use, iquad, &
              j, j_use, jc, jdim, k, kdim, keval, len, n, ni, ninside, nj, nk, nout_r, nout_s

   real    :: cell_volume, decay, ds2, dsqmin, dtolsq, dx, dy, dz, extra, overall_data_range, p, pm1, q, qm1, r, rm1, total,       &
              xfoot, xtarg, yfoot, ytarg, zfoot, ztarg

   real, dimension (3) :: foot_coords, target_coords

   logical :: new, same, save_edge_results, washout

!  Variables saved for further calls:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer, save :: nquad, ngeometric

   integer, allocatable, dimension (:,:), save :: conn

   real,    save :: derivs(3), growth_rate

   logical, save :: redistribute

   type (edge_type),   allocatable, dimension (:), save :: old_edges

   type (pq_map_type), allocatable, dimension (:), save :: pq_map

   real, dimension (:),   allocatable, save :: tline, xline, yline, zline, xmax, xmin, ymax, ymin, zmax, zmin
   real, dimension (:),   allocatable, save :: tnew,  xnew,  ynew,  znew
   real, dimension (:,:), allocatable, save :: arcs


!  Execution:
!  !!!!!!!!!!

   save_edge_results = output_g .and. .not. recenter_g

!  Pre-process the initial grid blocks?

   if (initialize) then

      n = nblocks_old

      if (output_f) then
         if (recenter_f) then ! Transfer all input flow blocks to the vertices

            do ib = 1, n
               call centers_vertices (1, num_q, old_grid(ib))  ! Local procedure below
            end do

         end if
      end if

!     Determine the data ranges of the surface faces of the initial grid.
!     Also, mark the initial grid block edges as exterior or interior.  (This info. is no longer used, though.)
!     The old-grid connectivity info (added most recently for the new grid) is not needed, but is tabulated anyway.

      allocate (old_edges(n))

      call get_surface_grid_info ('old', n, old_grid, old_edges, overall_data_range)

!     Pad the surface face data ranges slightly to allow for round-off.

      allocate (xmax(n), xmin(n), ymax(n), ymin(n), zmax(n), zmin(n))

      extra = overall_data_range * dxyz  ! Margin for being inside or outside a surface face of the grid being searched

      do ib = 1, n                       ! This may not gain anything except at the symmetry plane
         xmax(ib) = old_grid(ib)%xmax + extra
         xmin(ib) = old_grid(ib)%xmin - extra
         ymax(ib) = old_grid(ib)%ymax + extra
         ymin(ib) = old_grid(ib)%ymin - extra
         zmax(ib) = old_grid(ib)%zmax + extra
         zmin(ib) = old_grid(ib)%zmin - extra
      end do

      kdim = old_grid(1)%nk
      do ib = 2, n
         if (old_grid(ib)%nk /= kdim) then
            write (luncrt, '(/, a, i4, a, 2i4)') &
               ' Inconsistent radial dimensions in initial grid block', ib, ':  ', kdim, old_grid(ib)%nk
            stop
         end if
      end do

      if (save_edge_results) allocate (pq_map(nblocks_new)) ! Some new blocks will use a map for an edge of a previous new block

      allocate (xline(kdim), yline(kdim), zline(kdim), tline(kdim)) ! For interim new radial grid lines

      tline(1) = zero

      redistribute = ds1 > zero .and. .not. output_f  ! Don't mess with the cell-centered grid accompanying a flow soln. for now

      if (redistribute) then  ! Option to control the wall spacing
         ngeometric  = 2      ! I.e., no geometric portion
         growth_rate = one
         derivs(1)   = -999.  ! Suppresses unwanted derivative calculations
         allocate (xnew(kdim), ynew(kdim), znew(kdim), tnew(kdim))
      end if

      if (output_f) allocate (arcs(kdim,5)) ! 1-4 for corners, 5 for current new spoke during flow field interpolation

!     Build the ADT:
!     --------------

!     conn(1:3,iquad) is used to store the block # and grid point (i,j) for the cell iquad of the surface to be searched,
!     as part of building a tree from all cells of all blocks.

      nquad = 0
      do ib = 1, n
         nquad = (old_grid(ib)%ni - 1) * (old_grid(ib)%nj - 1) + nquad
      end do

      allocate (conn(3,nquad))

      call build_adt (n, old_grid, nquad, conn)

   end if ! End of initialization


   ninside = 0  ! # points of this new surface grid block meeting the distance tolerance for some old surface grid cell
   nout_s  = 0  ! # such points violating the distance tolerance at the surface
   nout_r  = 0  ! # new volume points outside the old outer boundary

   ni = new_grid%ni
   nj = new_grid%nj
   nk = new_grid%nk

   dtolsq  = dtol * dtol
   washout = dtol < zero

!  Check whether edge (p,q) results for this current new surface patch need to be saved for a subsequent new surface patch:

!! if (output_g) then
!!    save_edge_results = false
!!    do ie = 1, 4
!!       if (new_edges(iblock_new)%edge_con(ie)%ip > iblock_new) then
!!          allocate (pq_map(iblock_new)%pq (2,ni,nj))
!!          allocate (pq_map(iblock_new)%ijn(3,ni,nj))
!!          save_edge_results = .true.
!!          exit
!!       end if
!!    end do
!! end if

!  No:  just checking for two matching edges is insufficient.  Any number of patches may have a common corner.

   if (save_edge_results) then
      allocate (pq_map(iblock_new)%pq (2,ni,nj)) ! Avoiding interior points here would raise other problems
      allocate (pq_map(iblock_new)%ijn(3,ni,nj))
   end if

   write (luncrt, '(/, a, i4, a, 3i5, a, i8, /)') &
      ' Beginning new grid block', iblock_new, '.  Dimensions:', ni, nj, nk, '  Total:', ni * nj * nk

!  For each surface point of the current new grid block:

   do j = 1, nj

      do i = 1, ni

         target_coords(1) = new_grid%x(i,j,1);  xtarg = target_coords(1)
         target_coords(2) = new_grid%y(i,j,1);  ytarg = target_coords(2)
         target_coords(3) = new_grid%z(i,j,1);  ztarg = target_coords(3)

!        Check for reusing the edge results of an earlier new surface patch:

         if (save_edge_results) then

            if (i == 1) then ! Are we at an edge?
               ie = 1
            else if (i == ni) then
               ie = 2
            else if (j == 1)  then
               ie = 3
            else if (j == nj) then
               ie = 4
            else
               ie = 0
            end if

            if (ie > 0) then ! Does this edge match some earlier patch edge?
               ip_use = new_edges(iblock_new)%edge_con(ie)%ip  ! Patch number of neighbor to current edge
               if (ip_use > 0) then
                  if (ip_use < iblock_new) then ! The rest of this is really ugly!
                     ie_use = new_edges(iblock_new)%edge_con(ie)%ie  ! Relevant edge number of neighboring patch
                     iedge  = abs (ie_use)
                     same = (ie <= 2 .and. iedge <= 2) .or. (ie > 2 .and. iedge > 2)

                     if (iedge <= 2) then ! The neighboring edge is an i = 1 or i = imax edge
                        i_use = new_edges(ip_use)%i1(iedge) ! See get_surface_grid_info
                        if (same) then
                           j_use = j;  len = nj
                        else
                           j_use = i;  len = ni
                        end if
                        if (ie_use < 0) j_use = len + 1 - j_use  ! Common edge indices run in opposite directions
                     else ! The neighboring edge is a j = 1 or j = jmax edge
                        j_use = new_edges(ip_use)%j1(iedge)
                        if (same) then
                           i_use = i;  len = ni
                        else
                           i_use = j;  len = nj
                        end if
                        if (ie_use < 0) i_use = len + 1 - i_use
                     end if

                     p  = pq_map(ip_use)%pq (1,i_use,j_use);  q  = pq_map(ip_use)%pq (2,i_use,j_use)
                     ic = pq_map(ip_use)%ijn(1,i_use,j_use);  jc = pq_map(ip_use)%ijn(2,i_use,j_use)
                     n  = pq_map(ip_use)%ijn(3,i_use,j_use)
                     ninside = ninside + 1

!!!                  WRITE (6, '(a,2i4,a,2i4,a,3i4,a,1p,2e19.11)') &
!!!                     ' Reuse at i,j:',i,j, '  ip_use,ie_use:',ip_use,ie_use, '  ic,jc,n:',ic,jc,n, '  p,q:', p, q

                     go to 500  ! Very hard to avoid without excessive indenting - a rare case where it makes the most sense
                  end if
               end if

!              Alternatively, we may be at a corner, where any number of patches may all meet.
!              Look for a corner processed already that matches the current corner.
!              Note that we require exactness, because any discrepancy at the surface gets amplified along the radial line.
!              Observed example:  4.16634E-09 at k = 2 became 1.05211E-04 at the outer boundary.

               do ip_use = 1, iblock_new - 1
                  do icorner = 1, 4
                     if (xtarg == new_edges(ip_use)%xcorner(icorner)) then
                        if (ytarg == new_edges(ip_use)%ycorner(icorner)) then
                           if (ztarg == new_edges(ip_use)%zcorner(icorner)) then

                              select case (icorner) ! Counterclockwise order
                              case (1)
                                 i_use = 1
                                 j_use = 1
                              case (2)
                                 i_use = new_edges(ip_use)%i2(2) ! or other choices from edges 2, 3 or 4
                                 j_use = 1
                              case (3)
                                 i_use = new_edges(ip_use)%i2(2)
                                 j_use = new_edges(ip_use)%j2(2)
                              case (4)
                                 i_use = 1
                                 j_use = new_edges(ip_use)%j2(2)
                              end select

                              p  = pq_map(ip_use)%pq (1,i_use,j_use);  q  = pq_map(ip_use)%pq (2,i_use,j_use)
                              ic = pq_map(ip_use)%ijn(1,i_use,j_use);  jc = pq_map(ip_use)%ijn(2,i_use,j_use)
                              n  = pq_map(ip_use)%ijn(3,i_use,j_use)
                              ninside = ninside + 1

                              WRITE (6, '(a,2i4,a,2i4,a,3i4,a,1p,2e19.11)') &
                                 ' Reuse at i,j:',i,j, '  ip_use,icorner:',ip_use,icorner, '  ic,jc,n:',ic,jc,n, '  p,q:', p, q

                              go to 500 ! Very hard to avoid without excessive indenting
                           end if
                        end if
                     end if
                  end do
               end do

            end if ! End of test for being on an edge
         end if ! End of test for whether to match edges exactly

!        Search for the surface quad closest to the target point:
!        --------------------------------------------------------

         call search_adt (target_coords, iquad, p, q, dsqmin, true, nblocks_old, old_grid, nquad, conn, foot_coords)

!!!      dx = xtarg - foot_coords(1);  dy = ytarg - foot_coords(2);  dz = ztarg - foot_coords(3)
!!!      WRITE (6, '(a, 3i4, a, i7, 2f12.8, a, 1p, 3e16.8)') &
!!!         ' ib,i,j:', iblock_new,i,j, '  iquad,p,q:', iquad,p,q, '  dx/y/z:',dx,dy,dz

         if (dsqmin < dtolsq) then ! The best surface quad found was within the distance tolerance
            ninside = ninside + 1
         else
            nout_s = nout_s + 1
            write (luncrt, '(a, 3i4, a, 1p, e13.5)') &
               ' Tolerance exceeded for new (i,j,n) =', i, j, iblock_new, '  min. distance:', sqrt (dsqmin)
!!!         dx = xtarg - foot_coords(1);  dy = ytarg - foot_coords(2);  dz = ztarg - foot_coords(3)
!!!         write (6, '(a, 3i4, a, i7, 2f12.8, a, 1p, 3e16.8)') &
!!!            ' ib,i,j:', iblock_new,i,j, '  iquad,p,q:', iquad,p,q, '  dx/y/z:',dx,dy,dz
!!!         write (6, '(a, 1p, 3e16.8)') ' xtarg, ytarg, ztarg: ', xtarg, ytarg, ztarg
         end if

         n  = conn(1,iquad) ! Block #
         ic = conn(2,iquad) ! Lower left quad. indices
         jc = conn(3,iquad)

!!!      !!!!!!!!!!
500      continue ! Come here directly if we're reusing edge results from an earlier new surface patch to ensure identical results
!!!      !!!!!!!!!!

!        Interpolate within the best cell found.  The ADT search does not permit extrapolation.
!        --------------------------------------------------------------------------------------

         pm1 = one - p
         qm1 = one - q

!        Save the interpolated surface (x,y,z) for checking problem points.
!        We could probably use foot_coords(*) from SEARCH_ADT, but reuse of results from earlier edge points precludes this
!        unless we store them.  Reevaluating keeps the k = 1 results consistent with those for k = 2:nk (which use dx, dy, dz).

         xfoot = qm1 * (pm1 * old_grid(n)%x(ic,jc,  1) + p * old_grid(n)%x(ic+1,jc,  1)) + &
                   q * (pm1 * old_grid(n)%x(ic,jc+1,1) + p * old_grid(n)%x(ic+1,jc+1,1))

         interp_surf%x(i,j,1) = xfoot;  xline(1) = xfoot;  dx = xtarg - xfoot

         yfoot = qm1 * (pm1 * old_grid(n)%y(ic,jc,  1) + p * old_grid(n)%y(ic+1,jc,  1)) + &
                   q * (pm1 * old_grid(n)%y(ic,jc+1,1) + p * old_grid(n)%y(ic+1,jc+1,1))

         interp_surf%y(i,j,1) = yfoot;  yline(1) = yfoot;  dy = ytarg - yfoot

         zfoot = qm1 * (pm1 * old_grid(n)%z(ic,jc,  1) + p * old_grid(n)%z(ic+1,jc,  1)) + &
                   q * (pm1 * old_grid(n)%z(ic,jc+1,1) + p * old_grid(n)%z(ic+1,jc+1,1))

         interp_surf%z(i,j,1) = zfoot;  zline(1) = zfoot;  dz = ztarg - zfoot

!!!      WRITE (6, '(a, 1p, 6e15.5)') ' x/y/z/foot & dx/y/z:', xfoot, yfoot, zfoot, dx, dy, dz

         if (output_g) then ! Interpolate the rest of the new radial grid line off the surface.

!           The interim new line is needed to calculate arc-length-based decay of the correction at the surface.

            do k = 2, nk
               xline(k) = qm1 * (pm1 * old_grid(n)%x(ic,jc,  k) + p * old_grid(n)%x(ic+1,jc,  k)) + &
                            q * (pm1 * old_grid(n)%x(ic,jc+1,k) + p * old_grid(n)%x(ic+1,jc+1,k))
               yline(k) = qm1 * (pm1 * old_grid(n)%y(ic,jc,  k) + p * old_grid(n)%y(ic+1,jc,  k)) + &
                            q * (pm1 * old_grid(n)%y(ic,jc+1,k) + p * old_grid(n)%y(ic+1,jc+1,k))
               zline(k) = qm1 * (pm1 * old_grid(n)%z(ic,jc,  k) + p * old_grid(n)%z(ic+1,jc,  k)) + &
                            q * (pm1 * old_grid(n)%z(ic,jc+1,k) + p * old_grid(n)%z(ic+1,jc+1,k))
               tline(k) = tline(k-1) + sqrt ((xline(k) - xline(k-1))**2 + (yline(k) - yline(k-1))**2 + (zline(k) - zline(k-1))**2)
            end do

!           The k = 1 point is already in place.

            if (washout) then

!              Adjust the outer k points by the DECAYED difference in surface xyzs, so the original outer boundary is preserved.

               total = one / tline(nk)
               do k = 2, nk
                  decay = (tline(nk) - tline(k)) * total
                  new_grid%x(i,j,k) = xline(k) + decay * dx
                  new_grid%y(i,j,k) = yline(k) + decay * dy
                  new_grid%z(i,j,k) = zline(k) + decay * dz
               end do

            else ! For smooth surfaces with curvature and added resolution, not decaying can actually help the outer boundary.

               do k = 2, nk
                  new_grid%x(i,j,k) = xline(k) + dx
                  new_grid%y(i,j,k) = yline(k) + dy
                  new_grid%z(i,j,k) = zline(k) + dz
               end do

            end if

!           Control the grid spacing at the wall?  (It may actually be input to match spacing from an inner layer of blocks.)

            if (redistribute) then  ! As in OUTBOUND and DPLR's alignment of the outer grid boundary with the shock

               do k = 1, nk
                  xline(k) = new_grid%x(i,j,k)
                  yline(k) = new_grid%y(i,j,k)
                  zline(k) = new_grid%z(i,j,k)
               end do

               call chords3d (nk, xline, yline, zline, false, total, tline)

               call expdis5 (1, zero, total, ds1, nk, tnew, -luncrt)  ! One-sided Vinokur-type stretching, safeguarded

               ds2 = (total - tnew(nk-1)) * ds2fraction

               call blgrid (nk, ds1, ds2, ngeometric, growth_rate, tnew, luncrt, ier)  ! 2-sided Vinokur stretching

               keval = 1;  new = true  ! Arc-length-based local cubic spline interpolation

               do k = 1, nk
                  call plscrv3d (nk, xline, yline, zline, tline, lcs_method, new, false, tnew(k), keval, &
                                 xnew(k), ynew(k), znew(k), derivs)
                  new = false
                  new_grid%x(i,j,k) = xnew(k)
                  new_grid%y(i,j,k) = ynew(k)
                  new_grid%z(i,j,k) = znew(k)
               end do

            end if

            if (output_f) then ! Interpolate the flow solution similarly for this radial line.
               do k = 1, nk
                  new_grid%q(:,i,j,k) = qm1 * (pm1 * old_grid(n)%q(:,ic,jc,  k) + p * old_grid(n)%q(:,ic+1,jc,  k)) + &
                                          q * (pm1 * old_grid(n)%q(:,ic,jc+1,k) + p * old_grid(n)%q(:,ic+1,jc+1,k))
               end do
            end if

            if (save_edge_results) then
               if (i == 1 .or. i == ni .or. j == 1 .or. j == nj) then
                  pq_map(iblock_new)%pq (1,i,j) = p;   pq_map(iblock_new)%pq (2,i,j) = q
                  pq_map(iblock_new)%ijn(1,i,j) = ic;  pq_map(iblock_new)%ijn(2,i,j) = jc;  pq_map(iblock_new)%ijn(3,i,j) = n

!!!               WRITE (6, '(a,2i4,a,3i4,a,1p,2e19.11)') ' Save at i,j:',i,j, '  ic,jc,n:',ic,jc,n, '  p,q:', p, q

                  icorner = 0
                  if (i == 1) then
                     if (j == 1) then
                        icorner = 1
                     else if (j == nj) then
                        icorner = 4
                     end if
                  else if (i == ni) then
                     if (j == 1) then
                        icorner = 2
                     else if (j == nj) then
                        icorner = 3
                     end if
                  end if

                  if (icorner > 0) then
                     new_edges(iblock_new)%xcorner(icorner) = xtarg
                     new_edges(iblock_new)%ycorner(icorner) = ytarg
                     new_edges(iblock_new)%zcorner(icorner) = ztarg
                  end if

               end if
            end if

         else ! We must be doing flow interpolation only

            new_grid%q(:,i,j,1)  = qm1 * (pm1 * old_grid(n)%q(:,ic,jc,  1) + p * old_grid(n)%q(:,ic+1,jc,  1)) + &
                                     q * (pm1 * old_grid(n)%q(:,ic,jc+1,1) + p * old_grid(n)%q(:,ic+1,jc+1,1))

            if (nk > 1) then

               kdim = old_grid(n)%nk

!              Carry the surface interpolation into the volume for this new radial line:

               select case (method)

                  case (1) ! Assume the old and new radial lines are consistent

                     call consistent_radial_lines    ! Local procedure below

                  case (2) ! Assume the new radial lines can diverge from the old ones

                     idim = old_grid(n)%ni
                     jdim = old_grid(n)%nj

                     call more_general_radial_lines  ! Local procedure

                  case default

                     write (luncrt, '(/, a, i6)') ' Invalid method:', method
                     stop

               end select

            end if

         end if

!        The current surface point has been dealt with.

      end do ! Next i for this new surface grid block

      write (luncrt, '(a, i4, a, 2i5, a, i9, 2i8)') &
         ' Block:', iblock_new, '   pt.', i-1, j, ':  in, above_dtol, out_r:', ninside, nout_s, nout_r

   end do ! Next j for this new surface block


!  Interpolations are complete for this new block.
!  -----------------------------------------------

!  Is the block left-handed?

   call check_handedness (new_grid, cell_volume)                  ! Checks the first cell

   if (cell_volume < zero) call fix_handedness (new_grid, num_q)  ! Reverses the j indices in-place without affecting the header


!  Transfer interpolated flow to the cell centers?

   if (output_f) then
      if (recenter_f) call centers_vertices (2, num_q, new_grid)
   end if


   write (luncrt, '(/, a, i4, a, i9, 2i8)') &
      ' Finished block', iblock_new, ':  in, above_dtol, out_r:', ninside, nout_s, nout_r


!  Clean up?  (This allows for another set of input grid blocks in the same run -- not a likely requirement, though.)
!  ---------

   if (cleanup) then

      deallocate (old_edges, xmax, xmin, ymax, ymin, zmax, zmin)

      if (save_edge_results) then

         do ib = 1, nblocks_new
            do ie = 1, 4
               if (new_edges(ib)%edge_con(ie)%ip > ib) then
                  deallocate (pq_map(ib)%pq)
                  deallocate (pq_map(ib)%ijn)
                  exit
               end if
            end do
         end do

         deallocate (pq_map)

      end if

      if (output_f) deallocate (arcs)

      deallocate (conn, xline, yline, zline, tline)

      if (redistribute) deallocate (xnew, ynew, znew, tnew)

   end if

   continue

!  Local procedures for subroutine interp_block:
!  ---------------------------------------------

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

      subroutine consistent_radial_lines

!     Interpolate a flow field along a new radial line that is assumed to be consistent with the old radial lines, given the
!     details of the initial surface interpolation at new surface point (i,j,1).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: k, kc
      real    :: sk

!     Execution:

!     Calculate the radial arc lengths from the old grid for the corners of the current surface cell (ic,kc,1):

      arcs(1,:) = zero

      do k = 2, kdim ! Use the counterclockwise corner numbering convention of PROJECT4
         arcs(k,1) = arcs(k-1,1) + sqrt ((old_grid(n)%x(ic,jc,k)     - old_grid(n)%x(ic,jc,k-1))**2     + &
                                         (old_grid(n)%y(ic,jc,k)     - old_grid(n)%y(ic,jc,k-1))**2     + &
                                         (old_grid(n)%z(ic,jc,k)     - old_grid(n)%z(ic,jc,k-1))**2)
         arcs(k,2) = arcs(k-1,2) + sqrt ((old_grid(n)%x(ic+1,jc,k)   - old_grid(n)%x(ic+1,jc,k-1))**2   + &
                                         (old_grid(n)%y(ic+1,jc,k)   - old_grid(n)%y(ic+1,jc,k-1))**2   + &
                                         (old_grid(n)%z(ic+1,jc,k)   - old_grid(n)%z(ic+1,jc,k-1))**2)
         arcs(k,3) = arcs(k-1,3) + sqrt ((old_grid(n)%x(ic+1,jc+1,k) - old_grid(n)%x(ic+1,jc+1,k-1))**2 + &
                                         (old_grid(n)%y(ic+1,jc+1,k) - old_grid(n)%y(ic+1,jc+1,k-1))**2 + &
                                         (old_grid(n)%z(ic+1,jc+1,k) - old_grid(n)%z(ic+1,jc+1,k-1))**2)
         arcs(k,4) = arcs(k-1,4) + sqrt ((old_grid(n)%x(ic,jc+1,k)   - old_grid(n)%x(ic,jc+1,k-1))**2   + &
                                         (old_grid(n)%y(ic,jc+1,k)   - old_grid(n)%y(ic,jc+1,k-1))**2   + &
                                         (old_grid(n)%z(ic,jc+1,k)   - old_grid(n)%z(ic,jc+1,k-1))**2)
      end do

!     Bidirectional interpolation of the corner arc lengths at the current (p,q) location of each k layer:

      do k = 2, kdim
         arcs(k,5) = qm1 * (pm1 * arcs(k,1) +   p * arcs(k,2)) + &
                       q * (  p * arcs(k,3) + pm1 * arcs(k,4))
      end do

!!!   WRITE (6, '(a, 4i4, a)') ' i,j,kdim,nk:', i,j,kdim,nk, '   interpolated arcs(:,5):'
!!!   WRITE (6, '(15i9)') (k, k = 1, kdim)
!!!   WRITE (6, '(15f9.3)') arcs(1:kdim,5)

      kc = 1
      sk = zero

!!!   k = 1
!!!   write (6,'(a, 3i5, a, 1p, 3e16.8)') ' i,j,k: ', i,j,k, '  x,y,z: ', new_grid%x(i,j,k), new_grid%y(i,j,k), new_grid%z(i,j,k)

      do k = 2, nk

!!!      write (6,'(a, 3i5, a, 1p, 3e16.8)') ' i,j,k: ', i,j,k, '  x,y,z: ', new_grid%x(i,j,k), new_grid%y(i,j,k), new_grid%z(i,j,k)

!        Current new radial line arc-length at point k:

         sk = sk + sqrt ((new_grid%x(i,j,k) - new_grid%x(i,j,k-1))**2 + &
                         (new_grid%y(i,j,k) - new_grid%y(i,j,k-1))**2 + &
                         (new_grid%z(i,j,k) - new_grid%z(i,j,k-1))**2)

         call interval (kdim, arcs(1,5), sk, one, kc) ! Locate the k interval of sk in the underlying radial line

         r = (sk - arcs(kc,5)) / (arcs(kc+1,5) - arcs(kc,5))

         if (r <= one) then
            ninside = ninside + 1
         else
            nout_r = nout_r + 1
         end if

         r = min (r, one)
         rm1 = one - r

!!!      WRITE (6, '(a, i3, a, f9.3, a, 3i4, a, 3f12.8)') ' k:', k, '  sk:', sk, '  ic,jc,kc: ', ic,jc,kc, '  p,q,r: ', p,q,r

         new_grid%q(:,i,j,k) = rm1 * (qm1 * (pm1 * old_grid(n)%q(:,ic,jc,  kc  ) + p * old_grid(n)%q(:,ic+1,jc,  kc  ))  + &
                                        q * (pm1 * old_grid(n)%q(:,ic,jc+1,kc  ) + p * old_grid(n)%q(:,ic+1,jc+1,kc  ))) + &
                                 r * (qm1 * (pm1 * old_grid(n)%q(:,ic,jc,  kc+1) + p * old_grid(n)%q(:,ic+1,jc,  kc+1))  + &
                                        q * (pm1 * old_grid(n)%q(:,ic,jc+1,kc+1) + p * old_grid(n)%q(:,ic+1,jc+1,kc+1)))
      end do

      end subroutine consistent_radial_lines

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine more_general_radial_lines

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: eps = 1.e-5  ! Tolerance for RIPPLE3D/TRILINT tests on p,q,r in [0,1] and for convergence tests

!     Local variables:

      integer :: icell, jcell, kcell, istatus, k
      real    :: p, pm1, q, qm1, r, rm1, xt, yt, zt

!     Execution:

!     We have a reasonable starting guess for RIPPLE3D at each new k:

      icell = ic;  jcell = jc;  kcell = 1

!!!   WRITE (luncrt, '(a, 1p, 3e15.5)') ' dx, dy, dz:', dx, dy, dz

      do k = 2, nk

         xt = new_grid%x(i,j,k)
         yt = new_grid%y(i,j,k)
         zt = new_grid%z(i,j,k)

         if (k < nk) then ! This should handle difficulties in the boundary layer
            xt = xt - dx
            yt = yt - dy
            zt = zt - dz
         end if

!!!      WRITE (luncrt, '(a, 3i4, a, 1p, 3e15.5, a, 3i4)') &
!!!         ' ripple3d at ijk =', i,j,k, '  x/y/znew:', xt, yt, zt, '  i/j/kcell:', icell, jcell, kcell

         call ripple3d (idim, jdim, kdim, 1, idim, 1, jdim,1, kdim,  &
                        old_grid(n)%x, old_grid(n)%y, old_grid(n)%z, &
                        xt, yt, zt, icell, jcell, kcell, eps,        &
                        p, q, r, istatus)

!!!      WRITE (luncrt, '(a, 4i4, a, 4i4, a, 3f15.8)') &
!!!         ' Block #, ijk:', iblock_new, i, j, k, ' using old #, ijk:', n, icell, jcell, kcell, '; p,q,r: ', p, q, r

         if (istatus > 0) then ! No enclosing cell found
            write (luncrt, '(a, 4i4, a, 4i4, a, 3f15.8)') &
               ' Block #, ijk:', iblock_new, i, j, k, ' using old #, ijk:', n, icell, jcell, kcell, '; p,q,r: ', p, q, r
            nout_r = nout_r + 1
         else
            ninside = ninside + 1
         end if

         pm1 = one - p
         qm1 = one - q
         rm1 = one - r

         new_grid%q(:,i,j,k) = &
            rm1 * (qm1 * (pm1 * old_grid(n)%q(:,icell,jcell,  kcell  ) + p * old_grid(n)%q(:,icell+1,jcell,  kcell  ))  + &
                     q * (pm1 * old_grid(n)%q(:,icell,jcell+1,kcell  ) + p * old_grid(n)%q(:,icell+1,jcell+1,kcell  ))) + &
              r * (qm1 * (pm1 * old_grid(n)%q(:,icell,jcell,  kcell+1) + p * old_grid(n)%q(:,icell+1,jcell,  kcell+1))  + &
                     q * (pm1 * old_grid(n)%q(:,icell,jcell+1,kcell+1) + p * old_grid(n)%q(:,icell+1,jcell+1,kcell+1)))
      end do

      end subroutine more_general_radial_lines

   end subroutine interp_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine check_handedness (block, volume)

!  Calculate the volume of the "lower left" cell in a grid block.  If it is
!  negative, the block may be left-handed.
!
!  10/15/04   David Saunders     Copy of "test_volume" from Peter Gnoffo's LaRC
!             ELORET/NASA Ames   program for modeling Shuttle surface damage.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! External module

   implicit none

!  Arguments:

   type (grid_type), intent (in)  :: block  ! One grid block with nk > 1
   real,             intent (out) :: volume

!  Local constants:

   real, parameter :: sixty_fourth = 0.015625

!  Local variables:

   real :: xxi, yxi, zxi, xeta, yeta, zeta, xzeta, yzeta, zzeta

!  Execution:

!  The original multiplies by 1/4 have been combined as 1/64 below.

   xxi =   ( block%x(2,1,1) - block%x(1,1,1)            &
           + block%x(2,2,1) - block%x(1,2,1)            &
           + block%x(2,1,2) - block%x(1,1,2)            &
           + block%x(2,2,2) - block%x(1,2,2) )
   yxi =   ( block%y(2,1,1) - block%y(1,1,1)            &
           + block%y(2,2,1) - block%y(1,2,1)            &
           + block%y(2,1,2) - block%y(1,1,2)            &
           + block%y(2,2,2) - block%y(1,2,2) )
   zxi =   ( block%z(2,1,1) - block%z(1,1,1)            &
           + block%z(2,2,1) - block%z(1,2,1)            &
           + block%z(2,1,2) - block%z(1,1,2)            &
           + block%z(2,2,2) - block%z(1,2,2) )

   xeta  = ( block%x(1,2,1) - block%x(1,1,1)            &
           + block%x(2,2,1) - block%x(2,1,1)            &
           + block%x(1,2,2) - block%x(1,1,2)            &
           + block%x(2,2,2) - block%x(2,1,2) )
   yeta  = ( block%y(1,2,1) - block%y(1,1,1)            &
           + block%y(2,2,1) - block%y(2,1,1)            &
           + block%y(1,2,2) - block%y(1,1,2)            &
           + block%y(2,2,2) - block%y(2,1,2) )
   zeta  = ( block%z(1,2,1) - block%z(1,1,1)            &
           + block%z(2,2,1) - block%z(2,1,1)            &
           + block%z(1,2,2) - block%z(1,1,2)            &
           + block%z(2,2,2) - block%z(2,1,2) )

   xzeta = ( block%x(1,1,2) - block%x(1,1,1)            &
           + block%x(2,1,2) - block%x(2,1,1)            &
           + block%x(1,2,2) - block%x(1,2,1)            &
           + block%x(2,2,2) - block%x(2,2,1) )
   yzeta = ( block%y(1,1,2) - block%y(1,1,1)            &
           + block%y(2,1,2) - block%y(2,1,1)            &
           + block%y(1,2,2) - block%y(1,2,1)            &
           + block%y(2,2,2) - block%y(2,2,1) )
   zzeta = ( block%z(1,1,2) - block%z(1,1,1)            &
           + block%z(2,1,2) - block%z(2,1,1)            &
           + block%z(1,2,2) - block%z(1,2,1)            &
           + block%z(2,2,2) - block%z(2,2,1) )

   volume = ( xxi * ( yeta * zzeta  -  yzeta * zeta )   &
            - yxi * ( xeta * zzeta  -  xzeta * zeta )   &
            + zxi * ( xeta * yzeta  -  xzeta * yeta ) ) * sixty_fourth

   end subroutine check_handedness

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine fix_handedness (block, nvars)

!  Change the handedness of a grid block by reversing the j indices.
!  The block is updated in place so that subsequent I/O can remain efficient.
!
!  10/15/04   David Saunders, ELORET/NASA Ames   Initial implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! External module

   implicit none

!  Arguments:

   type (grid_type), intent (inout) :: block
   integer,          intent (in)    :: nvars ! nvars > 0 means flip "q" too

!  Local variables:

   integer                          :: i, j, k, l, njp1, njby2
   real                             :: t
   real, allocatable, dimension (:) :: s

!  Execution:

   njp1  = block%nj + 1
   njby2 = njp1 / 2

   do k = 1, block%nk
      do i = 1, block%ni
         do j = 1, njby2
            l = njp1 - j
            t = block%x(i,j,k)
            block%x(i,j,k) = block%x(i,l,k)
            block%x(i,l,k) = t
            t = block%y(i,j,k)
            block%y(i,j,k) = block%y(i,l,k)
            block%y(i,l,k) = t
            t = block%z(i,j,k)
            block%z(i,j,k) = block%z(i,l,k)
            block%z(i,l,k) = t
         end do
      end do
   end do

   if (nvars > 0) then

      allocate (s(nvars))

      do k = 1, block%nk
         do i = 1, block%ni
            do j = 1, njby2
               l = njp1 - j
               s = block%q(:,i,j,k)
               block%q(:,i,j,k) = block%q(:,i,l,k)
               block%q(:,i,l,k) = s
            end do
         end do
      end do

      deallocate (s)

   end if

   end subroutine fix_handedness

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine recenter_surface (block_in, block_out)

!  For one surface grid block, convert vertex-centered coordinates to cell-centered coordinates including one-layer halo cells.
!  The input and output blocks are assumed to be in different work-space.
!
!  The original version extrapolated to the centers of halo cells, but DPLR uses boundary values for halo cell coordinates, as here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! External module

   implicit none

!  Arguments:

   type (grid_type), intent (in)    :: block_in   ! Vertex-centered surface patch
   type (grid_type), intent (inout) :: block_out  ! Cell-centered equivalent in k = 1 layer, with halos; dimensions are already set

!  Local constants:

   real, parameter :: fourth = 0.25

!  Local variables:

   integer :: i, j, il, ir, jl, jr, ni, nj

!  Execution:

   ni = block_in%ni;  nj = block_in%nj

   do j = 1, nj + 1
      jl = max (1, j - 1);  jr = min (nj, j)
      do i = 1, ni + 1
         il = max (1, i - 1);  ir = min (ni, i)
         block_out%x(i,j,1) = (block_in%x(il,jl,1) + block_in%x(ir,jl,1) + block_in%x(il,jr,1) + block_in%x(ir,jr,1)) * fourth
         block_out%y(i,j,1) = (block_in%y(il,jl,1) + block_in%y(ir,jl,1) + block_in%y(il,jr,1) + block_in%y(ir,jr,1)) * fourth
         block_out%z(i,j,1) = (block_in%z(il,jl,1) + block_in%z(ir,jl,1) + block_in%z(il,jr,1) + block_in%z(ir,jr,1)) * fourth
      end do
   end do

   end subroutine recenter_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine recenter_volume (lun, formatted, block_out)

!  For one volume grid block, convert vertex-centered coordinates to cell-centered coordinates including one-layer halo cells.
!  The vertex-centered block is read into a local buffer of size determined from the output block dimensions.
!
!  The original version extrapolated to the centers of halo cells, but DPLR uses boundary values for halo cell coordinates, as here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! External module
   use xyzq_io_module        ! External module

   implicit none

!  Arguments:

   integer,          intent (in)    :: lun        ! Logical unit to read the next block from
   logical,          intent (in)    :: formatted  ! F = unformatted
   type (grid_type), intent (inout) :: block_out  ! Cell-centered volume block with halo cells to be output; dimensions are input

!  Local constants:

   real,      parameter :: eighth    = 0.125
   character, parameter :: name * 18 = ' recenter_volume: '

!  Local variables:

   integer              :: i, j, k, il, ir, jl, jr, kl, kr, ni, nj, nk, ios
   type (grid_type)     :: block_in ! For reading one vertex-centered grid block

!  Execution:

   ni = block_out%ni - 1;  block_in%ni = ni  ! Vertex-centered point counts (no halos)
   nj = block_out%nj - 1;  block_in%nj = nj
   nk = block_out%nk - 1;  block_in%nk = nk

   call xyz_allocate (block_in, ios)

   if (ios /= 0) then
     write (*, '(/, 2a)') name, 'allocation error'
     go to 99
   end if

   call xyz_block_io (1, lun, formatted, ni*nj*nk, block_in%x, block_in%y, block_in%z, ios)

   if (ios /= 0) then
     write (*, '(/, 2a)') name, 'block read error'
     go to 99
   end if

   do k = 1, nk + 1
      kl = max (1, k - 1);  kr = min (nk, k)
      do j = 1, nj + 1
         jl = max (1, j - 1);  jr = min (nj, j)
         do i = 1, ni + 1
            il = max (1, i - 1);  ir = min (ni, i)
            block_out%x(i,j,k) = eighth * (                                                                &
               block_in%x(il,jl,kl) + block_in%x(ir,jl,kl) + block_in%x(il,jr,kl) + block_in%x(ir,jr,kl) + &
               block_in%x(il,jl,kr) + block_in%x(ir,jl,kr) + block_in%x(il,jr,kr) + block_in%x(ir,jr,kr) )
            block_out%y(i,j,k) = eighth * (                                                                &
               block_in%y(il,jl,kl) + block_in%y(ir,jl,kl) + block_in%y(il,jr,kl) + block_in%y(ir,jr,kl) + &
               block_in%y(il,jl,kr) + block_in%y(ir,jl,kr) + block_in%y(il,jr,kr) + block_in%y(ir,jr,kr) )
            block_out%z(i,j,k) = eighth * (                                                                &
               block_in%z(il,jl,kl) + block_in%z(ir,jl,kl) + block_in%z(il,jr,kl) + block_in%z(ir,jr,kl) + &
               block_in%z(il,jl,kr) + block_in%z(ir,jl,kr) + block_in%z(il,jr,kr) + block_in%z(ir,jr,kr) )
         end do
      end do
   end do

   deallocate (block_in%x, block_in%y, block_in%z)

99 if (ios /= 0) then
      write (*, '(a, 3i5)') ' Vertex-centered block dimensions: ', ni, nj, nk
      stop
   end if

   end subroutine recenter_volume
