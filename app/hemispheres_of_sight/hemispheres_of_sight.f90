!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program hemispheres_of_sight

!     HEMISPHERES_OF_SIGHT is a variant of CONES_OF_SIGHT, itself an adaptation
!     of LINES_OF_SIGHT, q.v.  Originally, it worked with structured discretiz-
!     ations of lines of latitude and longitude, and retained the functionality
!     of LINES_OF_SIGHT if the lat/long steps were specified as zero.  However,
!     such discretizations produce clustering of the lines of site towards the
!     north pole of each hemisphere.  Therefore, this version triangulates the
!     underlying unit hemisphere (or rather one quadrant of it) to produce more
!     uniform (but not perfectly uniform) spacing of the lines associated with
!     one target body point.  The option to treat more than one body point at a
!     time has been retained, but generating only the primary line(s) of sight
!     is no longer an option - use LINES_OF_SIGHT for those.
!
!     Note that a full 3D volume grid, not half the volume, is all that makes
!     sense here, so be sure to reflect the half grid likely used for the
!     associated flow solution, unless the target points are all on the (y = 0)
!     symmetry plane, in which case the y < 0 half of the output hemisphere
!     data is suppressed.*  This volume grid is assumed to contain a single
!     layer of grid blocks with k in the body normal direction, as is typical
!     for hypersonic flow over atmospheric entry vehicles.  Ox points down-
!     stream, y is positive on the starboard side, and Oz is "up."
!
!     * In this case (all centerline body points), the input volume grid is
!       assumed to be the starboard (right) half.  Trying to reflect an input
!       left half proves more trouble than it is worth, as later steps in the
!       radiation calculations that this is presumably a part of are affected.
!       Therefore, the program endeavors to detect an input left half, and
!       stops early if necessary.
!
!     (Later:)  This version allows for suppressing lines to be processed by a
!     radiation solver as might be needed for a radiometer with a viewing cone
!     angle of less than 90 degrees.  The strategy adopted is this:
!
!     > Look for an optional cone angle with the body pt. x/y/z coordinates.
!     > Calculate all hemisphere lines and save all files as for the 90 degree
!       case.  (Suppressing parts of hemispherical surface files is not
!       practical.)
!     > Count the number of lines of sight (hemisphere vertices) within any
!       cone angle < 90, and write that many to an additional PLOT3D-type file
!       named cone.<angle>.<BP #>.lines.g file.
!     > Companion utility NEQAIR_DATA will process the lines in the latter file,
!       numbered 1:m, say. These LOS-n lines are then processed by NEQAIR, each
!       in its own /LINE-n subdirectory (n = 1:m).
!     > NEQAIR_INTEGRATION uses the same cone-angle calculation to read NEQAIR
!       results for these lines only and enter zero values for the rest before
!       integrating over all cells of the underlying hemisphere quadrants.
!
!  Further details (cone angle = 90 deg case):
!
!     For radiative heat flux calculations on an atmospheric entry vehicle, at
!     a specified body surface point, the "right" answer should account for all
!     possible viewing directions.  For a typical convex body, this means a
!     full hemisphere of lines of sight at each body point, not just the line
!     of sight normal to the body at each point as used with the tangent-slab
!     approximation in a typical radiation solver.  Here, the solver is expected
!     to provide a value of radiance (W.sr^-1.m^-2) for each hemisphere line
!     (vastly more expensive for even moderate discretizations), and a companion
!     utility will integrate those values with respect to solid angle subtended
!     at the body point to produce a truer estimate of radiative heat flux
!     (W.m^-2).
!
!     The primary line of sight is taken to be normal to the k = 1 surface at
!     the indicated target point.  It extends to the outer grid boundary with a
!     point distribution derived from (but not necessarily the same as) that of
!     the nearest local radial grid line.  The secondary lines have the same
!     relative distribution as the primary line.  All discretized lines are
!     saved in PLOT3D multiblock form (one line of sight per block) compatible
!     with the earlier FLOW_INTERP, which can perform the flow interpolations
!     and tabulations that are normally what are really desired for application
!     to hypersonic flows, in this case for radiation calculations involving
!     full angular integration.
!
!     The hemisphere discretization is determined by the prompted-for number of
!     points (ne) along a great circle from pole to equator (and along a 4th
!     of the equator) defining a certain triangulation.  A quarter of the
!     underlying unit hemisphere discretization is produced by slices parallel
!     to the equator through the uniform arc points, starting from near the
!     north pole.  The rest of the hemisphere is obtained by geometric trans-
!     formations.  Common edge points are presently NOT suppressed:  they will
!     be needed for eventual quadratures with respect to solid angle.
!
!     For the ne edge points specified, the number of node points on a quarter
!     hemisphere is ne*(ne + 1)/2, and the number of triangular elements is
!     (ne - 1)**2.  The primary line (along the surface normal through the
!     body point) is line number 1 in the <ID>.<bp #>.lines.g output file(s).
!
!     The triangulated surfaces, before and after the intersections-with-outer-
!     grid-boundary calculations, are output as multizone triangulations for
!     visualization with Tecplot.
!
!  Strategy:
!
!     o  Read the entire volume grid and extract the inner and outer boundaries
!        as multiblock surface grids.
!
!     o  For all body points (expected to be one, but more are allowed for)
!        defining (primary) lines of sight, search the inner boundary and save
!        the relevant patch number, cell indices, and fractional interpolation
!        coefficients (which are already in hand if surface grid indices are
!        input instead of coordinates).
!
!     o  Build a [new?] search tree from the outer boundary, for the inter-
!        section calculations.
!
!     o  Construct a triangulation of a quadrant of a unit hemisphere, as
!        indicated by the input number of edge points, ne.
!
!     o  For each inner surface target point:
!
!        > Construct a two-point line normal to the wall with length that of
!          a local radial grid line.  This should be at least as long as the
!          straight line distance to the outer boundary.
!
!        > This line determines the transformations to be applied to the
!          discretized unit hemisphere quadrant.
!
!        > For each line defined by the two or four transformed quadrants,
!          intersect the line with the outer boundary and impose the indicated
!          point distribution.
!
!        > Save all lines of sight in multiblock form, with each block sized
!          1 x 1 x npts, where npts is the prompted-for number.
!
!  Input body point format (read to EOF):
!
!     Either                                      or
!
!     n   i   j                                   x   y   z
!    [n   i   j                                  [x   y   z
!     n   i   j                                   x   y   z
!     :   :   :]                                  :   :   :]
!
!     where n = block number and k = 1 is implied.
!
!     The presence or absence of a decimal point in line 1 determines whether
!     indices or coordinates are being entered.
!
!  Outputs:
!
!     For each target body point:
!
!         (a) Discretized hemisphere (or half-hemisphere if all target points
!             are on the symmetry plane), tangent to the body at the body pts.
!             for visualization purposes (<ID>.<BP #>.hemi.dat):
!
!             Tecplot triangulation (1 hemisphere quadrant/zone; 2|4 zones)
!
!         (b) The topologically similar triangulated surface resulting from
!             intersecting this hemisphere with the roughly-uniformly-spaced
!             lines of sight, also for visualization purposes
!             (<ID>.<BP #>.boundary.dat):
!
!             Similar Tecplot triangulation
!
!         (c) Lines of sight in PLOT3D multiblock form, 1 line per block,
!             all blocks 1 x 1 x npts, where npts is prompted for, ASCII
!             (<ID>.<BP #>.lines.g):
!
!             nnodes x nquadrants
!             1 1 npts
!             1 1 npts
!               : : :
!
!  History:
!
!     10/07/05  D.A.Saunders  Initial design of LINES_OF_SIGHT.
!     11/22/05    "     "     Expanded the range of t for INTSEC6 to 2, not 1.1;
!                             Added the option to read a structured surface grid
!                             rather than a list of indices or coordinates.
!     03/06/06    "     "     Added the cone angle option.
!     03/07/06    "     "     Renamed it as CONES_OF_SIGHT.
!     03/08/06    "     "     HEMISPHERES_OF_SIGHT variant, using discretized
!                             latitude and longitude.
!     03/20/06    "     "     Users requested separate delta lat/long angles.
!     03/27/06    "     "     Move origin of longitude from xy plane to xz;
!                             suppress half a hemisphere if ytarget = 0.
!     03/29/06    "     "     Specify multiples of the first and last grid
!                             spacings, and allow the # points to vary from nk.
!     03/30/06    "     "     Reversed the storage of latitudes to go from 90 to
!                             zero, and output one line-per-block form with no
!                             replication at the "north pole."
!     02/01/14    "     "     All variants of the ADT routines have been merged
!                             into a module (generic build & search calls).
!     02/12/14    "     "     The search range for secondary-line intersections
!                             has been raised from 0.5 L : 10 L to 0.5 L : 30 L
!                             where L is the length of the initial primary line
!                             that is rotated to all the secondary positions.
!     02/21/14    "     "     A body point on the aft-body of a sphere-cone
!                             showed that 0.5 L is too large.  Use 0.1 L now.
!     03/17/14 -  "     "     Replaced latitude/longitude discretization with
!     03/24/14                triangulation of a unit hemisphere quadrant and
!                             transformations of that, to obtain roughly uniform
!                             solid-angle spacing of the lines of sight defined
!                             by the primary line and its hemisphere for each
!                             target body point.
!     03/26/14    "     "     Realized that the underlying unit hemisphere
!                             quadrant needs to be rotated about its original
!                             Oz axis to be properly axisymmetric, because the
!                             triangulation is done with constant-z slices, not
!                             constant-x (which in retrospect would have been
!                             preferable).
!     03/31/14    "     "     The last fix happily means the primary line is
!                             now the first line, not nnodes - (ne - 1).
!     01/26/15    "     "     An aft body point near the shoulder for MSL
!                             with a longer wake needed a larger search range
!                             for some of the intersections:  double the
!                             default range for aft body points.
!                             Also: element type needs to be TRIANGLE for all
!                             zones of output triangulations, now that the
!                             triangulation_io package handles volume grids too.
!     03/03/15    "     "     Unintended usage on a left-half volume grid
!                             caused enough grief that it needs to be trapped.
!                             Trying to reflect the data in this case was
!                             attempted, but abandoned.  We stop early instead.
!     02/23/16    "     "     An unusually long wake hit the heuristic range
!                             limit for the line/boundary intersection.
!     03/10/16    "     "     Application to radio signal attenuation estimation
!                             suggested an option for uniform discretization of
!                             the radial lines.
!     11/10/16    "     "     Repeat of the 02/23/16 issue (Red Dragon case).
!                             Apologies to Josh Monk.  This location is worst-
!                             case because the reference length of the grid line
!                             closest to the body normal line is short yet the
!                             lines in the wake can be very long.  Raise the
!                             limit for body points where the unit normal points
!                             aft from 3x40 reference lengths to 4x40 lengths.
!     11/11/16    "     "     Same issue: now, use the volume grid data range to
!                             limit excessive intervals for the intersections.
!     11/14/16    "     "     New problem stumbled on: concavities in the outer
!                             boundary can lead to spurious local minima where
!                             the minimum squared intersection distance is not
!                             close to 0.  One retry with a revised interval
!                             further outboard seems to be a workaround.
!     12/16/16    "     "     Josh Monk requested changing 4x40 L to 8x40L.
!                             Try that in the hope of not hurting less extreme
!                             cases.
!     02/14/18    "     "     Started arranging to suppress unneeded lines for
!                             the case of radiometers with some view cone angle
!                             less than 90 deg.  Enter the needed angle on the
!                             body point x/y/z line.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!           Now with ERC, Inc. at  NASA ARC.
!           Now with AMC, Inc. at  NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Employed by the XYZQ_IO package
   use xyzq_io_module        ! PLOT3D-type I/O package
   use surface_patch_utilities  ! For the outer boundary data ranges
   use tri_header_structure  ! For multizone unstructured Tecplot files
   use tri_zone_structure    ! For one zone of an unstructured Tecplot file
   use triangulation_io      ! I/O package for multizone surface triangulations
   use adt_utilities         ! All variants of the ADT build & search utilities
   use trigd

   implicit none

!  Local constants:

   integer, parameter :: &
      lunpoints = 1,     &   ! Input list of surface points
      lunvolg   = 2,     &   ! Input volume grid
      lunout1   = 3,     &   ! Output triangulated unit hemisphere quadrant
      lunout2   = 4,     &   ! Output body-tangent [half-]hemisphere trianglns.
      lunout3   = 7,     &   ! Output triangulations of LOS/bndry. intersectns.
      lunout4   = 8,     &   ! Output lines (PLOT3D multiblock, 1 LOS/block)
      lunkbd    = 5,     &   ! Keyboard inputs
      luncrt    = 6,     &   ! Screen
      nl        = 2          ! 2-point lines are passed to INTSEC6

   real, parameter ::    &
      big       = 1.e+9, &
      eps       = 1.e-5, &   ! For assuming a target is on the symmetry plane
      half      = 0.5,   &
      ninety    = 90.,   &
      one       = 1.0,   &
      zero      = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character, parameter :: &
      lcs_method * 1 = 'L'    ! But LCSFIT isn't needed (hook for general case)

!  Local variables:

   integer :: &
      i, ib, ibp, ic, ier, in, ios, iquad, j, jc, jn, k, l, lbp, lid, m, &
      n, ni, nj, nblocks, ne, nelements, nf, ninside, nk, nlines, nn, nnodes, &
      noutside, npts, nquad, nquadrants, ntargets

   integer, allocatable :: &
      conn(:,:),           &  ! For (patch,i,j) of surface quads. to search
      nij(:,:)                ! For (patch,i,j) of results of searches

   real :: &
      d1, d2, dmax, dmean, dsqmin, data_range, dtolsq, &
      obxmax, obxmin, obymax, obymin, obzmax, obzmin, &
      p, q, s1, s2, t, total, xyz_interp(3), ymax, ymin

   real, dimension (nl) :: &
      tl, xl, yl, zl          ! For 2-point lines passed to INTSEC6

   real, dimension (3) ::  &
      un, v1, v2              ! For unit normals and rotation axis end pts.

   real, allocatable, dimension (:) :: &
      cone_angle,                      & ! 90 or missing => full hemisphere
      tline, xline, yline, zline,      & ! For one radial volume grid line
      tnorm

   real, allocatable, dimension (:,:) :: &
      pq, xyz                 ! (p,q)s & (x,y,z)s of list of surface points

   logical :: &
      cell_centered, formatted, full_90, indices, uniform, yeq0

   character :: &
      answer*1, identifier*64

   type (grid_type), pointer, dimension (:) :: &
      inner_surface, outer_surface, radial_lines, volume_grid

   type (tri_header_type) :: &
      tri_header              ! One header for all transformations suffices

   type (tri_type), pointer, dimension (:) :: &
      intersected_quadrant, transformed_quadrant, unit_hemi_quadrant

!  Execution:
!  ----------

   write (luncrt, '(/, a)', advance='no') &
      ' # hemisphere points, pole to equator, > 0: '
   read (lunkbd, *) ne

!  Read the target body point indices or coordinates [& cone angles?]:
!  -------------------------------------------------------------------

   call get_bp_data ()
   if (ios /= 0) go to 99

!  Read the volume grid:
!  ---------------------

   call file_prompt (lunvolg, 0, 'volume grid', 'old', true, formatted, ios)
   if (ios /= 0) go to 99

   call xyzq_read (lunvolg, -lunvolg, formatted, nblocks, nf, cell_centered,   &
                   volume_grid, ios)
   if (ios /= 0) go to 99

!  Trap unwitting input of a left-half volume grid:

   ymin = huge (ymin)
   ymax = -ymin
   do ib = 1, nblocks
      ymax = max (ymax, maxval (volume_grid(ib)%y))
      ymin = min (ymin, minval (volume_grid(ib)%y))
   end do

   if (ymax < (ymax - ymin)*0.1) then
      write (luncrt, '(/, (a))') &
         ' The volume grid appears to be the left half (y < 0).', &
         ' If all body points are on the centerline, use REFLECT_BLOCKS',  &
         ' and input the right half volume grid here.  Otherwise, use it', &
         ' to provide both halves here.', &
         ' Remember to reflect any associated volume flow data to match.'
      go to 99
   end if

   write (luncrt, '(/, a)', advance='no') &
      ' Identifier for the triangulation & line-per-block outputs: '
   read  (lunkbd, *) identifier  ! Here because it may relate to the volume grid

   lid = len_trim (identifier) + 1
   identifier(lid:lid) = '.'

!  Allow for user-specified discretizations of the lines of sight, or uniform:

   nk = volume_grid(1)%nk

   write (luncrt, '(a, i4, a)', advance='no') &
      ' # radial volume grid pts. found:', nk, &
      '   Desired # along lines of sight: '
   read  (lunkbd, *) npts
   write (luncrt, '(2a)', advance='no') ' Multiples of inner & outer relative',&
      ' grid spacings to discretize with [say 100 10, or 0 0 => uniform]: '
   read  (lunkbd, *) s1, s2
   uniform = s1 <= zero

!  Extract the inner and outer boundaries as multiblock surface grids:
!  -------------------------------------------------------------------

   allocate (inner_surface(nblocks), outer_surface(nblocks))

   do ib = 1, nblocks
      inner_surface(ib)%ni = volume_grid(ib)%ni
      outer_surface(ib)%ni = volume_grid(ib)%ni
      inner_surface(ib)%nj = volume_grid(ib)%nj
      outer_surface(ib)%nj = volume_grid(ib)%nj
      inner_surface(ib)%nk = 1
      outer_surface(ib)%nk = 1
   end do

   obxmax = -big;  obxmin = big
   obymax = -big;  obymin = big
   obzmax = -big;  obzmin = big

   do ib = 1, nblocks
      call xyz_allocate (inner_surface(ib), ios)

      inner_surface(ib)%x(:,:,1) = volume_grid(ib)%x(:,:,1)
      inner_surface(ib)%y(:,:,1) = volume_grid(ib)%y(:,:,1)
      inner_surface(ib)%z(:,:,1) = volume_grid(ib)%z(:,:,1)

      call xyz_allocate (outer_surface(ib), ios)

      nk = volume_grid(ib)%nk
      outer_surface(ib)%x(:,:,1) = volume_grid(ib)%x(:,:,nk)
      outer_surface(ib)%y(:,:,1) = volume_grid(ib)%y(:,:,nk)
      outer_surface(ib)%z(:,:,1) = volume_grid(ib)%z(:,:,nk)

      call patch_data_range (outer_surface(ib))

      obxmax = max (obxmax, outer_surface(ib)%xmax)
      obxmin = min (obxmin, outer_surface(ib)%xmin)
      obymax = max (obymax, outer_surface(ib)%ymax)
      obymin = min (obymin, outer_surface(ib)%ymin)
      obzmax = max (obzmax, outer_surface(ib)%zmax)
      obzmin = min (obzmin, outer_surface(ib)%zmin)
   end do

   data_range = sqrt ((obxmax - obxmin)**2 + &
                      (obymax - obymin)**2 + &
                      (obzmax - obzmin)**2)
   write (luncrt, '(/, a, es13.6)') ' Outer boundary data range:', data_range

!  Count the surface quads (same for inner and outer surfaces):
!  ------------------------------------------------------------

   nquad = 0;  ymin = 1.e+30;  ymax = -ymin

   do ib = 1, nblocks
      ni = inner_surface(ib)%ni
      nj = inner_surface(ib)%nj
      nquad = (ni - 1) * (nj - 1) + nquad

      do j = 1, nj
         do i = 1, ni
            ymin = min (ymin, inner_surface(ib)%y(i,j,1))
            ymax = max (ymax, inner_surface(ib)%y(i,j,1))
         end do
      end do
   end do

   dtolsq = ((ymax - ymin) * 0.00001) ** 2  ! Tolerance for search diagnostics

   allocate (conn(3,nquad)) ! For (patch,i,j) of each surface quad.

!  Search the inner surface only if necessary:
!  -------------------------------------------

   if (indices) then ! No need to search the inner surface

      do n = 1, ntargets
         ib = nij(1,n)
         i  = nij(2,n)
         j  = nij(3,n)
         xyz(1,n) = inner_surface(ib)%x(i,j,1)
         xyz(2,n) = inner_surface(ib)%y(i,j,1)
         xyz(3,n) = inner_surface(ib)%z(i,j,1)
         pq (1,n) = zero
         pq (2,n) = zero
      end do

   else ! Build a search tree for the inner surface and do the searching:

      call build_adt (nblocks, inner_surface, nquad, conn)

!     For all target points, search the inner boundary
!     and save the relevant patch number and cell indices:
!     ----------------------------------------------------

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      do n = 1, ntargets

         call search_adt (xyz(:,n), iquad, p, q, dsqmin, true, nblocks, &
                          inner_surface, nquad, conn, xyz_interp)

         if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
            ninside  = ninside + 1
         else
            noutside = noutside + 1
         end if

         dmax  = max (dmax, dsqmin)
         dmean = dmean + sqrt (dsqmin)

         nij(:,n) = conn(:,iquad)
         xyz(:,n) = xyz_interp(:) ! This is what FLOW_INTERP should also get
         pq (1,n) = p
         pq (2,n) = q

      end do ! Next target surface point

      dmax = sqrt (dmax);  dmean = dmean / real (ntargets)

      write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
         ' # surface points inside/outside tolerance:', ninside, noutside, &
         ' max & mean distance:', dmax, dmean

      call release_adt ()

   end if

!  Build a [new?] search tree from the outer boundary patches:
!  -----------------------------------------------------------

   call build_adt (nblocks, outer_surface, nquad, conn)

   ninside  = 0;  dmax  = zero
   noutside = 0;  dmean = zero

!  Generate the indicated triangulation of a quadrant of a unit hemisphere:
!  ------------------------------------------------------------------------

   allocate (unit_hemi_quadrant(1))  ! Can't be a scalar because it's a pointer

   call spherical_triangulation (ne, unit_hemi_quadrant(1))

   nelements = unit_hemi_quadrant(1)%nelements  ! (ne - 1)**2   = # triangles
   nnodes    = unit_hemi_quadrant(1)%nnodes     ! (ne + 1)*ne/2 = # vertices

   write (luncrt, '(/, (4x, a, i6))') &
      '# edge pts. per hemisphere quadrant:', ne,        &
      '# quadrant vertices:                ', nnodes,    &
      '# quadrant elements:                ', nelements, &
      'Primary line number:                ', 1

!  This unit quadrant has all of x, y, z in [0, 1], and the rotational trans-
!  formations that complete the hemisphere should be about its Oz axis, not
!  Ox, because the triangulation is produced via constant-z slices.  At the
!  nose it needs to start outside the body, not inside as initially, then be
!  shifted and rotated to align its original Oz axis with the body normal.
!  It is less confusing if we rotate the initial quadrant 90 degrees so its
!  original Oz axis aligns with -Ox.

   v1(:) = zero  ! End points of rotation axis
   v2(:) = zero
   v2(2) = -one

   call rotate_xyz (nnodes, unit_hemi_quadrant(1)%xyz, ninety, v1, v2)

!! unit_hemi_quadrant(1)%xyz(1,:)     = -unit_hemi_quadrant(1)%xyz(1,:)
   unit_hemi_quadrant(1)%zone_title   = 'Unit hemisphere quadrant'
   unit_hemi_quadrant(1)%element_type = 'TRIANGLE'

!  Save it for visualization as a Tecplot surface triangulation:

   tri_header%filename    = 'unit_hemi_quadrant_vertices.dat'
   tri_header%fileform    = 1        ! Vertex-centered
   tri_header%formatted   = true
   tri_header%nvertices   = 3        ! Triangles, not tets
   tri_header%numf        = 0
   tri_header%nzones      = 1
   tri_header%datapacking = 0        ! Point order
   tri_header%title       = '# edge points:'
   write (tri_header%title(15:18), '(i4)') ne
   allocate (tri_header%varname(3))
   tri_header%varname(1)  = 'x'
   tri_header%varname(2)  = 'y'
   tri_header%varname(3)  = 'z'

   call tri_write (lunout1, tri_header, unit_hemi_quadrant, ios)

!  ------------------------------------------------------------------------
!  Every target point produces three files of either 2 or 4 quadrants each:
!     (a) the body-tangent [half-]hemisphere (2|4-zone triangulation);
!     (b) the intersection surface (similar triangulation);
!     (c) lines of sight, all 1 x 1 x npts (PLOT3D multiblock form)
!  ------------------------------------------------------------------------

   allocate (transformed_quadrant(1))  ! Work-space for all real-space quadrants
   allocate (intersected_quadrant(1))  ! ... and for line-surface intersections

   transformed_quadrant(1)%nelements = nelements
   intersected_quadrant(1)%nelements = nelements
   transformed_quadrant(1)%nnodes    = nnodes
   intersected_quadrant(1)%nnodes    = nnodes

   call tri_zone_allocate (tri_header, transformed_quadrant(1), ios)
   if (ios /= 0) go to 99

   call tri_zone_allocate (tri_header, intersected_quadrant(1), ios)
   if (ios /= 0) go to 99

   transformed_quadrant(1)%conn = unit_hemi_quadrant(1)%conn
   intersected_quadrant(1)%conn = unit_hemi_quadrant(1)%conn

   deallocate (unit_hemi_quadrant(1)%conn)

   allocate (xline(nk), yline(nk), zline(nk), tline(nk), tnorm(npts))

   tnorm(1) = zero;  tnorm(npts) = one

   if (uniform) call xgrid (npts, 0, zero, one, tnorm)

   tri_header%filename = identifier(1:lid)

   do ibp = 1, ntargets  ! Probably only one target body point

      call int_to_char (ibp, tri_header%filename(lid+1:), lbp)

      tri_header%filename(lid+lbp+1:) = '.hemi.dat'

      call tri_header_write (lunout2, tri_header, transformed_quadrant, ios)
      if (ios /= 0) go to 99

      tri_header%filename(lid+lbp+1:) = '.boundary.dat'

      call tri_header_write (lunout3, tri_header, transformed_quadrant, ios)
      if (ios /= 0) go to 99

      tri_header%filename(lid+lbp+1:) = '.lines.g'

      open (lunout4, file=tri_header%filename, status='unknown')

      yeq0 = abs (xyz(2,ibp)) < eps

      if (yeq0) then
         nquadrants = 2
      else
         nquadrants = 4
      end if

      nlines = nnodes * nquadrants

      write (luncrt, '(/, a, i5, a, i4)') &
         ' # lines of sight indicated:', nlines, '  for body point', ibp

      allocate (radial_lines(nlines))

      radial_lines(:)%ni = 1
      radial_lines(:)%nj = 1
      radial_lines(:)%nk = npts

!     Construct a carefully-calculated unit normal at the current body point.

      n = nij(1,ibp)  ! Patch # at inner surface for this target surface point
      i = nij(2,ibp)  ! Corresponding cell indices (lower left)
      j = nij(3,ibp)
      p = pq (1,ibp); in = i; if (p > half) in = i + 1 ! Nearest radial line
      q = pq (2,ibp); jn = j; if (q > half) jn = j + 1 ! from indicated cell

      do k = 1, nk  ! Nearby radial grid line to use for the point distribution
         xline(k) = volume_grid(n)%x(in,jn,k)
         yline(k) = volume_grid(n)%y(in,jn,k)
         zline(k) = volume_grid(n)%z(in,jn,k)
      end do

      call chords3d (nk, xline, yline, zline, true, total, tline) ! Normalized

      write (luncrt, '(/, a, es13.6)') ' Length of primary line:', total

!     Generate the relative distribution specified by the input scale factors?

      if (.not. uniform) then
         d1 = s1*tline(2);  d2 = s2*(one - tline(nk-1))

         call vinokur (1, npts, d1, d2, tnorm, luncrt, ier)

         if (ier /= 0) then
            write (luncrt, '(/, a, i5, a, i4, 2es16.8)') &
               ' Vinokur distribution trouble at target pt. #', l, &
               '; npts, d1, d2:', npts, d1, d2
            go to 99
         end if
      end if

!     Unit normal at the body point:

      call surface_normal (inner_surface(n)%ni, inner_surface(n)%nj, &
                           inner_surface(n)%x,  inner_surface(n)%y,  &
                           inner_surface(n)%z,  i, j, p, q, un)

      xl(1) = xyz(1,ibp)  ! Body-point end of all 2-pt LOSs to intersect
      yl(1) = xyz(2,ibp)
      zl(1) = xyz(3,ibp)

!     Ensure a normal in the symmetry plane where appropriate:

      if (abs (yl(1)) < eps) then ! Presumably ...
         write (luncrt, '(a, 3f12.8)') ' Body normal at y ~ 0:', un(:)
         yl(1) = zero
         un(2) = sqrt (un(1)**2 + un(3)**2)
         un(1) = un(1) / un(2)
         un(3) = un(3) / un(2)
         un(2) = zero
         write (luncrt, '(a, 3f12.8)') ' Adjusted body normal:', un(:)
      end if

      transformed_quadrant(1)%zone_title   = 'Quadrant x'
      transformed_quadrant(1)%element_type = 'TRIANGLE'
      intersected_quadrant(1)%element_type = 'TRIANGLE'

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      l = 0
      do m = 1, nquadrants

         write (transformed_quadrant(1)%zone_title(10:10), '(i1)') m
         intersected_quadrant(1)%zone_title = transformed_quadrant(1)%zone_title

         call transform_quadrant (m)

         call tri_zone_write (lunout2, tri_header, transformed_quadrant(1), ios)
         if (ios /= 0) go to 99

         do nn = 1, nnodes

            l = l + 1  ! LOS count over all quadrants

!           2-pt. LOS of length ~ grid line nearest to primary line:

            xl(2) = xl(1) + total*(transformed_quadrant(1)%xyz(1,nn) - xl(1))
            yl(2) = yl(1) + total*(transformed_quadrant(1)%xyz(2,nn) - yl(1))
            zl(2) = zl(1) + total*(transformed_quadrant(1)%xyz(3,nn) - zl(1))

            call xyz_allocate (radial_lines(l), ios)
            if (ios /= 0) go to 99

            call intersect_line ()

            intersected_quadrant(1)%xyz(1,nn) = radial_lines(l)%x(1,1,npts)
            intersected_quadrant(1)%xyz(2,nn) = radial_lines(l)%y(1,1,npts)
            intersected_quadrant(1)%xyz(3,nn) = radial_lines(l)%z(1,1,npts)

         end do  ! Next line of this quadrant

         call tri_zone_write (lunout3, tri_header, intersected_quadrant(1), ios)
         if (ios /= 0) go to 99

      end do  ! Next quadrant

      call xyz_write (lunout4, true, nlines, radial_lines, ios)

      if (cone_angle(ibp) < 90.) call save_cone_only ()  ! Reduced set of lines

      do i = 1, nlines
         deallocate (radial_lines(i)%x, radial_lines(i)%y, radial_lines(i)%z)
      end do

      deallocate (radial_lines)

      dmax = sqrt (dmax);  dmean = dmean / real (nlines)

      write (luncrt, '(/, a, i4, a, 2i5, /, a, 2es12.5)') &
         ' Body pt.', ibp, &
         ':  # outer intersection points inside/outside tolerance:', &
         ninside, noutside, ' max & mean distance:', dmax, dmean

   end do  ! Next target body point

!  Let Fortran 90 do the other deallocates.

99 continue ! Avoid system dependencies where STOP is concerned

!  Internal procedures for program HEMISPHERES_OF_SIGHT:
!  -----------------------------------------------------

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bp_data ()  ! The body pt. file may include cone angles now
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: lbuf, nitems
      logical :: ascii
      character (3)  :: seps
      character (64) :: buffer

!     Execution:

      ascii = true  ! Can't pass a constant as an inout argument.

      call file_prompt (lunpoints, 0, 'body point x/y/z[/cone angle]', 'old', &
                        false, ascii, ios)
      if (ios /= 0) go to 99

      ntargets = 0
      do ! Count the target points till EOF
         read (lunpoints, '(a)', iostat=ios) buffer ! We can look for '.' below
         if (ios /= 0) exit
         ntargets = ntargets + 1
      end do

      write (luncrt, '(/, a, i4)') ' # target points specified: ', ntargets

      allocate (nij(3,ntargets), xyz(3,ntargets), pq(2,ntargets), &
                cone_angle(ntargets))

!     Assume the user doesn't mix indices with cone angles that include '.':

      lbuf = index (buffer, '!') - 1  ! Avoid counting comment tokens
      if (lbuf < 0) lbuf = len_trim (buffer)

      indices =  index (buffer(1:lbuf), '.') == 0
      seps    = ' ,' // char (9)  ! Space, comma, or tab

      call token_count (buffer(1:lbuf), seps, nitems)  ! Should be 3 or 4

      write (luncrt, '(a, l1)') ' Body point indices specified? ', indices
      if (nitems == 3) then
         write (luncrt, '(a)') ' Not expecting cone angle with BP inputs.'
      else
         write (luncrt, '(a)') ' Expecting cone angle with BP inputs.'
      end if

      rewind (lunpoints)
      cone_angle(:) = 90.
      if (indices) then
         do n = 1, ntargets ! Allow for trailing comments
            if (nitems == 3) then
               read (lunpoints, *) nij(:,n)
            else
               read (lunpoints, *) nij(:,n), cone_angle(n)
            end if
         end do
      else
         do n = 1, ntargets
            if (nitems == 3) then
               read (lunpoints, *) xyz(:,n)
            else
               read (lunpoints, *) xyz(:,n), cone_angle(n)
            end if
         end do
      end if
      close (lunpoints)
      ios = 0
99    return

      end subroutine get_bp_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine transform_quadrant (m)
!
!     A quadrant of a triangulated unit hemisphere is transformed 2 or 4 times
!     to define all the lines of sight for one body point.  Processing for each
!     quadrant is completed before the next.  For the first quadrant (m = 1),
!     the unit quadrant is transformed so that its primary axis along Ox is
!     coincident with the unit normal at the body point. Each next quadrant is
!     obtained by rotating the previous quadrant 90 degrees about the body
!     normal.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Argument:

      integer, intent (in) :: m  ! Quadrant number (1:nquadrants)

!     Local variables:

      real, parameter :: small_value = 0.001

!     Local variables:

      real, save    :: theta
      logical, save :: parallel

!     Execution:

      if (m == 1) then  ! Align the unit hemisphere quadrant w/ the body normal

         v1(1) = -one;  v1(2) = zero;  v1(3) = zero  ! Unit vector along -Ox

         write (luncrt, '(a, i4)')     ' Body pt:', ibp
         write (luncrt, '(a, 3i4)')    ' ib,i,j :', nij(:,ibp)
         write (luncrt, '(a, 3f12.8)') ' x,y,z  :', xyz(:,ibp)
         write (luncrt, '(a, 3f12.8)') ' Unit nm:', un(:)
         write (luncrt, '(a, 3f12.8)') ' v1(:)  :', v1(:)

         call cross (v1, un, v2)     ! v2 = v1 x un provides a second point on
                                     ! a rotation axis through the body point,
                                     ! and also the initial angle of rotation

         write (luncrt, '(a, 3f12.8)') ' v1 x un:', v2(:)

!        Avoid theta = 90, where tan (theta) is infinite:

         if (abs (un(1)) < small_value) then
            theta = sign (ninety, un(3))
         else                                                      ! |v1 x un| /
            theta = atand (sqrt (un(2)**2 + un(3)**2) / (-un(1)))  !  v1 . un
            if (un(1) > zero) theta = theta + ninety + ninety      ! Aft body
         end if

         write (luncrt, '(a, 3f14.8)') ' theta  :', theta

!        Avoid parallel vectors, with zero-magnitude cross-product:

         parallel = dot_product (v2, v2) < small_value  ! Normal is along Ox

         v1(:) = xyz(:,ibp)          ! First point on initial rotation axis

         if (parallel) then          ! Second  "   "   "   "   "   "   "
            v2(:) = v1(:)
            v2(2) = v2(2) + one      ! Axis is parallel to Oy
         else
            v2(:) = v1(:) + v2(:)
         end if

         write (luncrt, '(a, 3f12.8)') ' Axis p2:', v2(:)

!        The first transformation requires an initial shift:

         do i = 1, nnodes
            transformed_quadrant(1)%xyz(:,i) = &
               unit_hemi_quadrant(1)%xyz(:,i) + v1(:)
         end do

      else if (m == 2) then     ! The extra quadrants rotate the previous
                                ! quadrant about the body normal
         v2(:) = v1(:) + un(:)  ! Second point along normal
         theta = ninety

      end if

!     Rotate in-place:

      call rotate_xyz (nnodes, transformed_quadrant(1)%xyz, theta, v1, v2)

      end subroutine transform_quadrant

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine intersect_line ()
!
!     Intersect the line of sight defined by the current vertex of the
!     current transformed quadrant, and impose the specified distribution.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: search_a = 0.1  ! Conservative search interval limits
      real, parameter :: search_b = 40.  ! as multiples of the normalized
                                         ! primary line length
      real, parameter :: aftmultb = 8.0  ! Fudge factor for handling long wake

!     Local variables:

      real :: dx, dy, dz

!     Execution:

      tl(1) = search_a      ! Fractions of primary line length to search in
      tl(2) = search_b      ! Err on the high side

!     Aft-body points near the shoulder can produce a wide range of intersected
!     line lengths.  Don't compromise forebody cases:

      if (un(1) >= zero) then
         tl(1) = tl(1)*half
         tl(2) = tl(2)*aftmultb
      end if

!!    write (66, '(a, 4f10.4)') 'u1, tl(1:2), tl(2) adjusted:', &
!!       un(1), tl(:), 2.0*data_range/total

!     Avoid excessively large search intervals, but not large enough risks
!     a bad local minimum distance that isn't very close to zero, with t at
!     tl(1) + eps because of the way FMINRC in INTSEC6 works.  :(

      tl(2) = min (tl(2), 2.0*data_range/total)  ! tl(1), tl(2) are normalized

      call intsec6_2pt (nblocks, outer_surface, nquad, conn, nl, xl, yl, zl, &
                        tl, lcs_method, iquad, p, q, t, xyz_interp, dsqmin)

!!!   write (luncrt, '(a, 2i4, 4i10)') &
!!!      ' vertex, iquad, conn(:,iquad):', nn, iquad, conn(:,iquad)
!!!   write (luncrt, '(a, 3es19.11)') ' p, q, t:   ', p, q, t, &
!!!                                   ' xyz_interp:', xyz_interp

      if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
         ninside  = ninside + 1
!!!      write (luncrt, '(a, i6, 4es13.6)') &
!!!      ' No  intersection trouble.  Line #, ta, tb, dsq, t:', l, tl, dsqmin, t
      else
         noutside = noutside + 1
         write (luncrt, '(a, i6, 4es13.6)') &
         ' *** Intersection trouble.  Line #, ta, tb, dsq, t:', l, tl, dsqmin, t

         write (luncrt, '(a)') &
         '     Assume local min. due to boundary concavity.  Do one retry.'
         tl(1) = tl(1)*1.000001
         tl(2) = tl(2)*2.

         call intsec6 (nblocks, outer_surface, nquad, conn, nl, xl, yl, zl, &
                       tl, lcs_method, iquad, p, q, t, xyz_interp, dsqmin)

         if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
            ninside  = ninside + 1;  noutside = noutside - 1
            write (luncrt, '(a, i6, 4es13.6)') &
               '     Recovered.       Line #, ta, tb, dsq, t:', l, tl, dsqmin, t
         else
            write (luncrt, '(a, i6, 4es13.6)') &
         ' *** Intersection failure.  Line #, ta, tb, dsq, t:', l, tl, dsqmin, t

            ! Keep going
         end if
      end if

      dmax  = max (dmax, dsqmin)
      dmean = dmean + sqrt (dsqmin)

!     Transform the distribution of the primary line to the secondary:

      dx = (xl(2) - xl(1)) * t ! Length components
      dy = (yl(2) - yl(1)) * t
      dz = (zl(2) - zl(1)) * t

      do k = 1, npts
         radial_lines(l)%x(1,1,k) = xl(1) + dx * tnorm(k)
         radial_lines(l)%y(1,1,k) = yl(1) + dy * tnorm(k)
         radial_lines(l)%z(1,1,k) = zl(1) + dz * tnorm(k)
      end do

      end subroutine intersect_line

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_cone_only ()  ! Suppress lines outside the cone angle
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, ncone, nk
      logical, allocatable :: keep(:)
      real :: angle, vnormal(3), v2(3)
      type (grid_type), pointer, dimension (:) :: cone_lines
      character (64) :: filename

!     Execution:

      filename(1:11)  = 'cone.xx.xx.'
      write (filename(6:10), '(f5.2)') cone_angle(ibp)
      filename(12:)   = tri_header%filename
      open (lunout4, file=filename, status='unknown')

      nk = radial_lines(1)%nk
      vnormal(1) = radial_lines(1)%x(1,1,nk) - radial_lines(1)%x(1,1,1)
      vnormal(2) = radial_lines(1)%y(1,1,nk) - radial_lines(1)%y(1,1,1)
      vnormal(3) = radial_lines(1)%z(1,1,nk) - radial_lines(1)%z(1,1,1)
      allocate (keep(nlines))
      keep(1)  = true
      keep(2:) = false

      ncone = 1
      do i = 2, nlines
         v2(1) = radial_lines(i)%x(1,1,nk) - radial_lines(1)%x(1,1,1)
         v2(2) = radial_lines(i)%y(1,1,nk) - radial_lines(1)%y(1,1,1)
         v2(3) = radial_lines(i)%z(1,1,nk) - radial_lines(1)%z(1,1,1)
         call angle_between_vectors (vnormal, v2, angle)
         if (angle > cone_angle(ibp)) cycle
            ncone = ncone + 1
            keep(i) = true
      end do

      allocate (cone_lines(ncone))
      cone_lines(:)%ni = 1
      cone_lines(:)%nj = 1
      cone_lines(:)%nk = nk

      ncone = 0
      do i = 1, nlines
         if (.not. keep(i)) cycle
            ncone = ncone + 1
            allocate (cone_lines(ncone)%x(1,1,nk), &
                      cone_lines(ncone)%y(1,1,nk), &
                      cone_lines(ncone)%z(1,1,nk))
            cone_lines(ncone)%x(:,:,:) = radial_lines(i)%x(:,:,:)
            cone_lines(ncone)%y(:,:,:) = radial_lines(i)%y(:,:,:)
            cone_lines(ncone)%z(:,:,:) = radial_lines(i)%z(:,:,:)
      end do

      call xyz_write (lunout4, true, ncone, cone_lines, ios)

      do i = 1, ncone
         deallocate (cone_lines(i)%x, cone_lines(i)%y, cone_lines(i)%z)
      end do
      deallocate (keep, cone_lines)

      end subroutine save_cone_only

   end program hemispheres_of_sight
