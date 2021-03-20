!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program Stardust_Lines

!  This version of Stardust_Lines is generalized to handle axisymmetric grid
!  topologies that do not necessarily have just one layer of blocks around the
!  body.  It was prompted by meteor simulations with a 5-block topology that
!  allows arbitrarily long wakes.
!
!  Reasonable Stipulations:
!
!     o  The CFD coordinate system has Ox pointing downstream, with the nose
!        at (0,0) in 2D and (0,0,0) in 3D.  The revolved 2D grid has Oy > 0
!        on the starboard side and Oz points up.
!
!     o  The axisymmetry means the starboard half of the 3D grid suffices;
!        integration results can be doubled, and the flow solution looks the
!        same from any azimuthal angle (side view), so we work with the view
!        along the symmetry plane, rotated about Oy counterclockwise by the
!        input viewing angle.
!
!     o  The input viewing angle is relative to the CFD coordinate system's Ox.
!        If the flight patch angle were zero (parallel to the horizon) this
!        angle would an angle of elevation (positive looking from below/front).
!        In practice, the actual trajectory flight path angle has to be taken
!        into account along with the observing aircraft position. These details
!        are TBD at the time of writing.
!
!     o  Grid block wall and free-stream boundaries are either jmin or jmax,
!        and ibc >= 25 means wall, while ibc = 1 or 3 means outer boundary.
!
!     o  Block BC data are input as in these two examples from DPLR inputs, as
!        the file 'block-bc.dat.  If this is not present, a one-layer block
!        topology is assumed.
!
!        >  2 blocks:  Here, both blocks have wall and outer boundary faces.
!
!              Boundary condition type [ibc]
!              imin imax jmin jmax kmin kmax
!                14   20   25    1   19   19
!                20   14   25    3   19   19
!
!        >  5 blocks:  Here, only 2 blocks have a wall face and 3 have an outer
!                      boundary face.
!
!              Boundary condition type [ibc]
!              imin imax jmin jmax kmin kmax
!                14   20   25   20   19   19
!                14   20   20    1   19   19
!                20   14   25   20   19   19
!                20   20   20    1   19   19
!                20    3   14    1   19   19
!
!  Internally, the 2D blocks are reindexed so that their rotated forms have any
!  inner and outer boundaries at kmin and kmax.  Thus, the only real change is
!  that some (3D) blocks don't have a wall boundary at kmin and some don't have
!  a free-stream boundary at kmax.
!
!  The following description is carried over from the one-layer version:
!
!  Generate lines of sight suitable for estimating radiative energy from an
!  axisymmetric 2-D grid of the type used to calculate Stardust flow solutions.
!  The inner and outer grid boundaries are assumed to be closed with j = 1 at
!  the wall.  The rotated form produces closed multipatch surfaces with k = 1
!  at the wall.  (Actually, only half the volume grid is needed, as all lines
!  of sight can be considered to be in the vertical plane.)  If a flow solution
!  accompanies the grid, it should contain two temperatures followed by the
!  species number densities.  (Later:  Any flow solution is ignored here.  That
!  is needed after the lines of sight are produced here, for interpolation of
!  the temperatures and number densities on to them, by FLOW_INTERP or similar.)
!
!  Outline:
!
!  o  Read control info. such as the number of increments for the 3rd dimension,
!     the appropriate viewing angle, and the density of the lines-of-sight grid
!  o  Read the 2-D grid and associated flow solution if present (PLOT3D)
!  o  For each 2-D grid block
!     o  Rotate all grid points through the clock angle increments and store
!        as radial lines in the k direction forming a 3-space grid to be output
!        for FLOW_INTERP purposes; the inner and outer boundaries will be
!        processed here to identify lines of sight, and they are saved for
!        plotting purposes
!     o  If the flow solution is present, store its k lines similarly and save
!        it for FLOW_INTERP purposes; the flow solution isn't used here, but
!        the 3-D form of it is needed for the FLOW_INTERP step
!  o  Write the 3-D grid and optional flow solution (PLOT3D ASCII or not)
!  o  Establish a rectangular grid of points in a plane normal to the intended
!     lines of sight at the indicated viewing angle; save this as a PLOT3D
!     surface grid file for the integration of radiation solver results; it is
!     tangent to the CFD volume grid at the upstream and/or lower outer boundary
!  o  For each line of sight
!     o  Determine 2 or 0 intersections with the 3-D volume grid outer boundary
!  o  For each line of sight
!     o  Determine 2 or 0 intersections with the 3-D volume grid inner boundary
!     o  Combine the 4, 2, or 0 intersection points into 2, 1, or 0 segments
!     o  Discretize the segments using roughly the average radial grid spacing
!  o  Write results, with each segment a block of the output lines of sight file
!
!  Control File Format (Standard Input):
!
!     Stardust_Lines control file
!     ---------------------------
!     stardust.g           ! Input 2-D grid
!     T                    ! Formatted? [T|F]
!     none stardust.f      ! Input 2-D solution, or none (not used either way)
!     T                    ! Formatted? [T|F]
!     37                   ! # points in 3rd index direction via rotation
!     45.6                 ! Viewing angle above Ox, deg.
!     49    0.  ! d1 added ! # LOS grid pts.: vert. <-> 0 theta; x dirn. <-> 90
!     33    0.  ! d1 added ! # cross-stream LOS pts. & relative dy at y = 0.
!     stardust.3D.gu       ! Output volume grid
!     F                    ! Formatted? [T|F]
!     stardust.3D.fu       ! Output volume flow field, or none
!     F                    ! Formatted? [T|F]
!     line_segments.g      ! Output line segments, one segment per block
!     T                    ! Formatted? [T|F]
!     99.                  ! Wake cut-off (x body D from nose; big = no cutoff)
!     100  150             ! Min. & max # points on discretized lines
!
!  Control File Notes:
!
!     o   Theta = 90. gives a side view, and njLOS can be just 1.  All other
!         rotational slices through the CFD solution can be derived from that.
!         This possibility applies to meteor cases.
!     o   The rectangular grid defining the lines of sight was originally
!         uniform, but now it may contain one-sided stretching.  A second
!         input on the niLOS and njLOS lines controls the stretching as an
!         initial relative increment on [0, 1], with 0. used to indicate
!         uniform spacing.  If theta = 90., niLOS and its d1 input produce a
!         stretched grid in the X direction on [xmin, xmax] of the outer
!         CFD grid boundary, and d1 should be small (say 1.e-5).
!     o   The "pixel" array needed for integration of the radiation solver
!         outputs is oriented orthogonal to the viewing angle above Ox.
!         The pixel plane written here as pixels.g is tangent to the 3D volume
!         grid at the upstream side (or nearest to the ground if the viewing
!         angle is 90 degrees). In this case, enter something like 89.99999 deg
!         if you want to avoid the "meteor" option.
!     o   The lines formed by the rectangular pixel array that miss the 3D
!         volume grid are indicated in the output line segments file by nk = 2.
!         They should be ignored by the radiation solver.  The header of this
!         file can be used by the integration utility (Stardust_Integration by
!         the present author) to insert zero radiation at those pixels and to
!         identify the (i,j) of each meaningful line of sight.
!     
!  History:
!
!     03/21/06  D. Saunders  Initial design.
!     03/22/06   "       "   Initial code and test.
!     04/03/06   "       "   Introduced wake region cut-off option.
!     04/07/06   "       "   The optional flow-field I/O is functional now.
!     04/18/06   "       "   Tabulate the rectangular viewing grid details.
!     05/08/06   "       "   Additional outputs group the line segments that
!                            encounter the inner surface and those that pass
!                            through the outer boundary (only). The file names
!                            are hard-coded as lines_inner.g & lines_outer.g.
!     05/09/06   "       "   Yen Liu asked for the (i,j)s to be written to
!                            another pair of text files.
!     06/28/06   "       "   Yen asked for control over the number of points
!                            on the discretized lines.  Change nmin = 30 and
!                            nmax = 100 from constants to new inputs.
!     08/08/13   "       "   All ADT variants have been merged into one module
!                            with generic build_adt and search_adt interfaces.
!     10/03/14   "       "   Made use of xyq_io 2D I/O utilities (not available
!                            originally) with possible revival for asteroid
!                            studies in mind.
!     10/16/14   "       "   (After a hiatus:) Needed to add xyq_2d_to_3d to
!                            xyq_io to match original indexing.
!     05/15/15   "       "   Application to "meteors" with long wakes prompted
!                            several refinements: a 90-degree viewing angle
!                            with just 1 cross-stream point produces a single
!                            plane of lines applicable to every rotational
!                            station; 1-sided stretching of the stations in
!                            the axial direction is preferable to uniform
!                            spacing; wake cut-off should be a multiple of
!                            body diameter, not body length; separating nose-
!                            cap lines from wake lines when theta = 90. allows
!                            distinguishing wake radiation from nose-cap +
!                            body.
!     07/27/15   "       "   Generalization prompted by meteor cases with long
!                            wakes: allow for more than one layer of grid blocks
!                            in the input 2D grid, but assume that the wall and
!                            outer boundaries are still at jmin and jmax
!                            (no need to handle i boundaries).
!     02/20/20   "       "   Don't move the first column of pixels away from
!                            the symmetry plane.  The line-surface intersections
!                            should still be OK. (Revisited for Hayabusa 2.)
!                            Echo the control file to standard output.
!     02/24/20   "       "   Clarified viewing angle issues in the documentation
!                            above. Save the pixel array needed for the intended
!                            integrations, as pixels.g.
!     09/18/20   "       "   Argument linenum has been added to the line_surface
!                            utility.
!     09/19/20   "       "   Added printing of (i,j), pixel #, active line #.
!     10/15/20   "       "   Just a typo in a write statement.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!           Later with ERC, Inc. at ARC, then with AMA, Inc. at ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! 3-space grid block derived data type
   use xyq_io_module         ! 2-space PLOT2D-type file I/O package
   use xyzq_io_module        ! 3-space PLOT3D-type file I/O package
   use adt_utilities         ! All ADT variants
   use surface_patch_utilities

   implicit none

!  Local constants:

   integer, parameter :: &
      lun2Dg = 1,        &   ! Input 2-D grid ...
      lun2Df = 2,        &   ! ... and optional flow solution
      lunbcs = 3,        &   ! Optional file 'block-bc.dat'
      lunsrf = 4,        &   ! Inner & outer 3-D boundaries for visualization
      lunctl = 5,        &   ! Control file (standard input)
      lunlog = 6,        &   ! Log file (standard output)
      lun3Dg = 7,        &   ! Output 3-D grid ...
      lun3Df = 8,        &   ! ... and flow solution
      lunLOS = 9,        &   ! Output lines of sight (PLOT3D grid)
      lunLI  = 10,       &   ! Output lines meeting the inner surface ...
      lunLO  = 11,       &   ! ... and the outer surface only
      lunLN  = 12,       &   ! Added to distinguish nose-cap lines if theta = 90
      lunTI  = 13,       &   ! Corresponding lists of (i,j)s (text)
      lunTO  = 14,       &
      lunTN  = 15,       &
      lunscr = 16,       &   ! Scratch file avoids possible rewind issues
      lunpix = 17,       &   ! Save pixel array for later integrations
      nline  = 2             ! 2-point lines passed to line_surface

   real, parameter ::    &
      big    = 1.e10,    &
      ninety = 90.,      &
      one    = 1.,       &
      zero   = 0.

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character, parameter :: &
      bc_file*12   = 'block-bc.dat', &
      format*11    = 'unformatted',  &
      lcs_method*1 = 'L'   ! But LCSFIT isn't needed for 2-point lines (hook)

!  Local variables:

   integer :: &
      i, i1, ib, ios, j, k, lactive, line, luno, lunq, n, n3D, nblocks, &
      nb_inner, nb_outer, nf, ni, niLOS, nj, njLOS, nk, nlines_inner, &
      nlines_outer, nlines_nosecap, nmax, nmin, npts, nquad, ns, nspaces

   integer, allocatable, dimension (:,:) :: &
      conn, &                ! For (patch,i,j) of surface quads. to search
      ibc                    ! For (block,bc) data

   real :: &
      arc, d1_i, d1_j, dphi, ds, ds_nominal, dx, dy, dz, phi, r, s, s1, s2, &
      sum, theta, theta2, x, x1, x2, xcutoff, xless, xmore, xmax, xmin, xnose, &
      y, ymax, z, zmax, zmin

   real, dimension (nline) :: &
      xline, yline, zline, sline

   real, allocatable, dimension (:) :: &
      relative_i, relative_j

   real, allocatable, dimension (:,:) :: &
      yrotate, zrotate

   logical :: &
      cell_centered, cr, eof, flow_present, formatted, formatted_lines, &
      side_view

   logical, allocatable :: &
      has_inner(:), has_outer(:)

   character (80) :: &
      filename

   type (grid_type) :: &
      bounding_box           ! For rectangular grids defining lines of sight

   type (grid_type), pointer, dimension (:) :: &
      xyzf2D, xyzf3D, inner_surface, outer_surface

   type (grid_type), allocatable, dimension (:,:) :: &
      line_segments

!  Execution:
!  ----------

!  Echo the control file to standard output:

   call echo (lunctl, lunlog, lunscr)

   read (lunscr, *) ! Standard input
   read (lunscr, *)
   read (lunscr, *) filename    ! Input 2-D grid
   read (lunscr, *) formatted;  i1 = 1;  if (formatted) i1 = 3

   open (lun2Dg, file=filename, status='old', form=format(i1:11), iostat=ios)
   if (ios /= 0) then
      write (lunlog, '(/, 2a)') ' Unable to open the 2-D grid file: ',         &
         trim (filename)
      go to 99
   end if

   read (lunscr, *) filename    ! Input flow (optional)
   read (lunscr, *) formatted;  i1 = 1;  if (formatted) i1 = 3

   flow_present = filename(1:4) /= 'none'
   lunq = -lun2Df

   if (flow_present) then
      open (lun2Df, file=filename, status='old', form=format(i1:11), iostat=ios)
      if (ios /= 0) then
         write (lunlog, '(/, 2a)') ' Unable to open the 2-D flow file: ',      &
            trim (filename)
         go to 99
      end if
      lunq = lun2Df
   end if

   call xyq_read (lun2Dg, lunq, formatted, nblocks, nf, cell_centered, &
                  xyzf2D, ios)
   if (ios /= 0) go to 99

!  Before xyq_io was installed, the 2D grid and flow were dimensioned (ni,1,nk)
!  to help the 3D rotation, so now we need to convert (ni,nj,1) to (ni,1,nj):

   call xyq_2d_to_3d (nblocks, nf, xyzf2D, ios)  ! Added to xyq_io to help here
   if (ios /= 0) go to 99

   read (lunscr, *) n3D          ! # pts. in 3rd index dir. via 180 deg rotation
   read (lunscr, *) theta        ! Viewing angle, degrees
   side_view = theta == 90.      ! Meteor case; avoid messing up original scheme
   read (lunscr, *) niLOS, d1_i  ! # vertical pts. in LOS grid (if theta = 0)
   read (lunscr, *) njLOS, d1_j  ! # cross-stream pts.  "   "; see Notes for d1s
   read (lunscr, *) filename     ! Output volume grid
   read (lunscr, *) formatted;  i1 = 1;  if (formatted) i1 = 3

   open (lun3Dg, file=filename, status='unknown', form=format(i1:11),iostat=ios)
   if (ios /= 0) then
      write (lunlog, '(/, 2a)') ' Unable to open the output volume grid: ',    &
         trim (filename)
      go to 99
   end if

   read (lunscr, *) filename      ! or none (output flow file)
   read (lunscr, *) ! Same format as grid

   if (flow_present) then
      open (lun3Df, file=filename, status='unknown', form=format(i1:11),       &
            iostat=ios)
      if (ios /= 0) then
         write (lunlog, '(/, 2a)') ' Unable to open the output flow file: ',   &
            trim (filename)
         go to 99
      end if
   end if

   read (lunscr, *) filename    ! Output line segments
   read (lunscr, *) formatted_lines;   i1 = 1;  if (formatted_lines) i1 = 3
   read (lunscr, *) xcutoff     ! Multiple of D aft of xmin; converted to an X
   read (lunscr, *) j, k        ! Min. & max. # points on discretized lines

   nmin = min (j, k)
   nmax = max (j, k)

   close (lunscr, status='delete')

   open (lunLOS, file=filename, status='unknown', form=format(i1:11),iostat=ios)

   if (ios /= 0) then
      write (lunlog, '(/, 2a)') ' Unable to open the output lines grid: ',     &
         trim (filename)
      go to 99
   end if

!  Optional BC file for input grids that are not necessarily single-layer:

   allocate (has_inner(nblocks), has_outer(nblocks))

   open (lunbcs, file=bc_file, status='old', iostat=ios)

   if (ios /= 0) then      ! Assume single-layer topology
      has_inner(:) = true  ! Wall is at jmin for all blocks
      has_outer(:) = true  ! Outer boundary is at jmax for all blocks
   else
      allocate (ibc(6,nblocks))
      read (lunbcs, *, iostat=ios)
      read (lunbcs, *, iostat=ios)
      read (lunbcs, *, iostat=ios) ibc
      if (ios /= 0) then
         write (lunlog, '(/, 2a)') ' Trouble reading BC data from file ', &
         bc_file
         go to 99
      end if
      close (lunbcs)
      do ib = 1, nblocks
         has_inner(ib) = ibc(3,ib) >= 25
         has_outer(ib) = ibc(4,ib) <= 3
      end do
      deallocate (ibc)
   end if

   nb_inner = 0
   nb_outer = 0
   do ib = 1, nblocks
      if (has_inner(ib)) nb_inner = nb_inner + 1
      if (has_outer(ib)) nb_outer = nb_outer + 1
   end do

!  Rotate the 2-D grid to produce half a volume grid; find an average spacing:
!  ---------------------------------------------------------------------------

   allocate (xyzf3D(nblocks))

   dphi = 180. / real (n3D - 1);  nj = n3D;  sum = zero;  nspaces = 0

   do ib = 1, nblocks

      ni = xyzf2D(ib)%ni;  xyzf3D(ib)%ni = ni
                           xyzf3D(ib)%nj = nj
      nk = xyzf2D(ib)%nk;  xyzf3D(ib)%nk = nk

      npts = ni * nk;  nspaces = ni * (nk - 1) + nspaces

      do i = 1, ni
         do k = 2, nk
            sum = sqrt ((xyzf2D(ib)%x(i,1,k) - xyzf2D(ib)%x(i,1,k-1))**2  +    &
                        (xyzf2D(ib)%z(i,1,k) - xyzf2D(ib)%z(i,1,k-1))**2) + sum
         end do
      end do

      allocate (yrotate(ni,nk), zrotate(ni,nk)) ! To get a j-plane contiguous

      allocate (xyzf3D(ib)%x(ni,nj,nk), xyzf3D(ib)%y(ni,nj,nk),                &
                xyzf3D(ib)%z(ni,nj,nk))

      if (flow_present) then ! %mi/j/k are needed by q_write

         xyzf3D(ib)%mi = ni;  xyzf3D(ib)%mj = nj;  xyzf3D(ib)%mk = nk

         allocate (xyzf3D(ib)%q(nf,ni,nj,nk))

      end if

      do j = 1, n3D

         phi = -dphi * real (j - 1)  ! < 0 gives the half for which y > 0

         do k = 1, nk
            do i = 1, ni
               yrotate(i,k) = zero
               zrotate(i,k) = xyzf2D(ib)%z(i,1,k)
            end do
         end do

!        Counterclockwise, in place; assume the x-axis is the line of symmetry

         call rotate2d (npts, yrotate, zrotate, phi, zero, zero)

         do k = 1, nk
            do i = 1, ni
               xyzf3D(ib)%x(i,j,k) = xyzf2D(ib)%x(i,1,k)
               xyzf3D(ib)%y(i,j,k) = yrotate(i,k)
               xyzf3D(ib)%z(i,j,k) = zrotate(i,k)
            end do
         end do

         if (flow_present) xyzf3D(ib)%q(:,:,j,:) = xyzf2D(ib)%q(:,:,1,:)

      end do

      deallocate (yrotate, zrotate)

   end do ! Next block to rotate

   ds_nominal = 0.05 * sum / real (nspaces)

   write (lunlog, '(/, a, es14.6)') &
      ' Nominal spacing to be used for discretizing line segments:',    &
      ds_nominal
   write (lunlog, '(a, 2i5)') ' Min. & max. # pts. per segment:', nmin, nmax

!  Save the volume grid and optional flow-field:

   call xyz_write (lun3Dg, formatted, nblocks, xyzf3D, ios)

   if (ios /= 0) then
      write (lunlog, '(/, a)') ' Trouble writing the volume grid.'
      go to 99
   end if

   if (flow_present) then

      call q_write (lun3Df, formatted, nblocks, nf, xyzf3D, ios)

      if (ios /= 0) then
         write (lunlog, '(/, a)') ' Trouble writing the volume flow field.'
         go to 99
      end if

      do ib = 1, nblocks
         deallocate (xyzf3D(ib)%q)  ! The flow field will be used by FLOW_INTERP
      end do

   end if

!  Extract the inner and outer boundaries as surface grids for intersecting:
!  -------------------------------------------------------------------------

   allocate (inner_surface(nb_inner))

   i = 0
   do ib = 1, nblocks
      if (.not. has_inner(ib)) cycle
      i = i + 1
      ni = xyzf3D(ib)%ni
      inner_surface(i)%ni=ni;  inner_surface(i)%nj=nj;  inner_surface(i)%nk=1
      call xyz_allocate (inner_surface(i), ios)
      inner_surface(i)%x(:,:,1) = xyzf3D(ib)%x(:,:,1)
      inner_surface(i)%y(:,:,1) = xyzf3D(ib)%y(:,:,1)
      inner_surface(i)%z(:,:,1) = xyzf3D(ib)%z(:,:,1)
   end do

   allocate (outer_surface(nb_outer))

   i = 0
   do ib = 1, nblocks
      if (.not. has_outer(ib)) cycle
      i = i + 1
      ni = xyzf3D(ib)%ni
      nk = xyzf3D(ib)%nk
      outer_surface(i)%ni=ni;  outer_surface(i)%nj=nj;  outer_surface(i)%nk=1
      call xyz_allocate (outer_surface(i), ios)
      outer_surface(i)%x(:,:,1) = xyzf3D(ib)%x(:,:,nk)
      outer_surface(i)%y(:,:,1) = xyzf3D(ib)%y(:,:,nk)
      outer_surface(i)%z(:,:,1) = xyzf3D(ib)%z(:,:,nk)
   end do

!  Turn the wake cut-off input into an X.  We need the OML data ranges:

   xmax = -big;  xmin = big;  zmax = -big

   do ib = 1, nb_inner
      call patch_data_range (inner_surface(ib))

      xmax = max (xmax, inner_surface(ib)%xmax)
      xmin = min (xmin, inner_surface(ib)%xmin)
      zmax = max (zmax, inner_surface(ib)%zmax)
   end do

   xcutoff = xmin + (zmax + zmax) * xcutoff  ! Body diam. multiple aft of nose
   xnose   = xmin

!  Save the boundaries for graphics purposes:

   open (lunsrf, file='inner_surface.gu', form='unformatted', status='unknown')

   call xyz_write (lunsrf, false, nb_inner, inner_surface, ios)

   open (lunsrf, file='outer_surface.gu', form='unformatted', status='unknown')

   call xyz_write (lunsrf, false, nb_outer, outer_surface, ios)

!  Establish a rectangular grid of pts. in a plane normal to the lines of sight:
!  -----------------------------------------------------------------------------

!  First, to locate the grid, temporarily rotate the outer boundary to find the
!  the data range easily:

   theta2 = ninety - theta

   xmax = -big;  xmin = big
   ymax = -big
   zmax = -big;  zmin = big

   do ib = 1, nb_outer

      npts = outer_surface(ib)%ni * outer_surface(ib)%nj

      call rotate2d (npts, outer_surface(ib)%x, outer_surface(ib)%z, theta2,   &
                     zero, zero)

      call patch_data_range (outer_surface(ib))

      xmax = max (xmax, outer_surface(ib)%xmax)
      xmin = min (xmin, outer_surface(ib)%xmin)
      ymax = max (ymax, outer_surface(ib)%ymax)
      zmax = max (zmax, outer_surface(ib)%zmax)
      zmin = min (zmin, outer_surface(ib)%zmin)

   end do

!  Recopy the outer boundary rather than rotating it back:

   i = 0
   do ib = 1, nblocks
      if (.not. has_outer(ib)) cycle
      i = i + 1
      nk = xyzf3D(ib)%nk
      outer_surface(i)%x(:,:,1) = xyzf3D(ib)%x(:,:,nk)
      outer_surface(i)%y(:,:,1) = xyzf3D(ib)%y(:,:,nk)
      outer_surface(i)%z(:,:,1) = xyzf3D(ib)%z(:,:,nk)
   end do

!  Generate a rectangular bounding box grid for the rotated grid.
!  It may or may not be a uniform grid.

   dx = (xmax - xmin) / real (niLOS - 1)
   dy =  ymax         / real (max (njLOS - 1, 1))
   dz =  zmax - zmin

   bounding_box%ni = niLOS;  bounding_box%nj = njLOS;  bounding_box%nk = 2

   call xyz_allocate (bounding_box, ios)

   allocate (relative_i(niLOS))

!  One-sided stretching, with handling of [near-]uniform cases:

   call expdis5 (1, zero, one, d1_i, niLOS, relative_i, -lunlog)

   allocate (relative_j(njLOS))

   if (njLOS > 1) then  ! It should be 1 if theta = 90 deg
      call expdis5 (1, zero, one, d1_j, njLOS, relative_j, -lunlog)
!!!   relative_j(1) = 0.5 * relative_j(2)  ! Avoid intersection trouble at y = 0
!     Don't like this; shouldn't be a problem.
   else
      relative_j(1) = 1.e-6
   end if

   do j = 1, njLOS
      y = ymax * relative_j(j)
      do i = 1, niLOS
         bounding_box%x(i,j,:) = xmin + (xmax - xmin) * relative_i(i)
         bounding_box%y(i,j,:) = y
         bounding_box%z(i,j,1) = zmax
         bounding_box%z(i,j,2) = zmin
      end do
   end do

!  The viewing grid details are needed for some quadrature:

   write (lunlog, &
      '(/, a, //, a, f8.2, a, /, a, i4, a, i4, /, (a, 3es15.6))') &
      ' Bounding box grid details (rotated):', &
      ' Viewing angle: ', theta, ' degrees',   &
      ' Discretization:', niLOS, '  x', njLOS, &
      ' xmin, xmax, dx:', xmin, xmax, dx, &
      ' ymin, ymax, dy:', zero, ymax, dy, &
      ' zmin, zmax, dz:', zmin, zmax, dz

!  Rotate this bounding box back to the original coordinate system:

   npts = niLOS * njLOS * 2

   call rotate2d (npts, bounding_box%x, bounding_box%z, -theta2, zero, zero)

!  The pixel array needed for integration of radiation solver results is
!  the k = 2 plane of this bounding box.  Save it in PLOT3D surface form.

   open (lunpix, file='pixels.g', status='unknown')

   write (lunpix, '(i1)') 1
   write (lunpix, '(3i5)') niLOS, njLOS, 1
   write (lunpix, '(6es19.11)') bounding_box%x(:,:,2), bounding_box%y(:,:,2), &
                                bounding_box%z(:,:,2)
   close (lunpix)

!  Set up for intersecting the outer boundary with the lines of sight:
!  -------------------------------------------------------------------

   allocate (bounding_box%q(3,niLOS,njLOS,2))  ! q(1:3,:,:,1) = outer ns, s1, s2
                                               ! q(1:3,:,:,2) = inner ns, s1, s2
   nquad = 0
   do ib = 1, nb_outer
      nquad = (outer_surface(ib)%ni - 1) * (outer_surface(ib)%nj - 1) + nquad
   end do

   allocate (conn(3,nquad))  ! for patch # and (i,j)

   call build_adt (nb_outer, outer_surface, nquad, conn)

   write (lunlog, '(/, a, /)') ' Outer surface intersections:'

   sline(1) = zero;  sline(2) = dz;  ds = dz;  line = 0

   do j = 1, njLOS

      do i = 1, niLOS

         line = line + 1
         xline(1) = bounding_box%x(i,j,1);  xline(2) = bounding_box%x(i,j,2)
         yline(1) = bounding_box%y(i,j,1);  yline(2) = bounding_box%y(i,j,2)
         zline(1) = bounding_box%z(i,j,1);  zline(2) = bounding_box%z(i,j,2)

         call line_surface (line, nb_outer, outer_surface, nquad, conn,        &
                            nline, xline, yline, zline, sline, lcs_method,     &
                            ns, s1, s2)

!!!      write (52, '(a, 3i5)') 'outer i,j,ns:', i, j, ns

!        If intersections were found, they may still be beyond xcutoff:

         if (ns > 0) then  ! ns = 2;  s1 <= s2;  lower s is at higher altitude

            dx = xline(2) - xline(1)
            if (abs (dx) < 1.e-10) then  ! Theta = 90 deg possibility
               dx = sign (1.e-10, dx)
            end if
            x1 = xline(1) + (s1 / ds) * dx
            x2 = xline(1) + (s2 / ds) * dx
            xless = min (x1, x2)
            xmore = max (x1, x2)

            if (xless >= xcutoff) then
               ns = 0
            else
               if (xmore >= xcutoff) s1 = (xcutoff - xline(1)) * (ds / dx)
            end if

            bounding_box%q(2,i,j,1) = s1
            bounding_box%q(3,i,j,1) = s2

!!          write (lunlog, '(a, 2i4, a, 5es15.7)') &
!!             ' i,j:', i,j, '  x1, x2, xc, s1, s2:', x1, x2, xcutoff, s1, s2
         end if

         bounding_box%q(1,i,j,1) = real (ns)

      end do

   end do

!  Repeat for the inner boundary, skipping if ns = 0 for the outer boundary:

   call release_adt ()              ! Deallocate work-space from the first ADT

   call build_adt (nb_inner, inner_surface, nquad, conn)

   write (lunlog, '(/, a, /)') ' Inner surface intersections:'

   line = 0

   do j = 1, njLOS

      do i = 1, niLOS

         line = line + 1
         xline(1) = bounding_box%x(i,j,1);  xline(2) = bounding_box%x(i,j,2)
         yline(1) = bounding_box%y(i,j,1);  yline(2) = bounding_box%y(i,j,2)
         zline(1) = bounding_box%z(i,j,1);  zline(2) = bounding_box%z(i,j,2)

         if (bounding_box%q(1,i,j,1) > zero) then

            call line_surface (line, nb_inner, inner_surface, nquad, conn,     &
                               nline, xline, yline, zline, sline, lcs_method,  &
                               ns, s1, s2)

!!!         write (52, '(a, 3i5)') 'inner i,j,ns:', i, j, ns

            if (ns > 0) then  ! Just check for a bad X cut-off value

               dx = xline(2) - xline(1)
               x1 = xline(1) + (s1 / ds) * dx
               x2 = xline(1) + (s2 / ds) * dx
               xmore = max (x1, x2)

               if (xmore >= xcutoff) then
                  write (lunlog, '(/, a, /, a, 3es15.6)') &
                  ' Cutoff X is less than inner surface intersection X(s): ',  &
                  ' xcutoff, x1, x2: ', xcutoff, x1, x2
                  go to 99
               end if

               bounding_box%q(2,i,j,2) = s1
               bounding_box%q(3,i,j,2) = s2

!!             write (lunlog, '(a, 2i4, a, 5es15.7)') &
!!             ' i,j:', i,j, '  x1, x2, xc, s1, s2:', x1, x2, xcutoff, s1, s2

            end if

            bounding_box%q(1,i,j,2) = real (ns)

         else ! No inner surface intersections

            bounding_box%q(:,i,j,2) = zero

         end if

      end do

   end do

!  Discretize the meaningful line segments:
!  ----------------------------------------

   nlines_inner = 0;  nlines_outer = 0

   write (lunlog, '(/, a, /)') &
      ' Arc lengths for lines intersecting both boundaries:'

   allocate (line_segments(niLOS,njLOS))

   do j = 1, njLOS

      do i = 1, niLOS

         line_segments(i,j)%ni = 1;  line_segments(i,j)%nj = 1

         if (bounding_box%q(1,i,j,1) == zero) then ! No intersections
            s1 = zero
            s2 = zero  ! It's not clear yet what to do here
            nk = 2
         else          ! s1 < s2; s1 is at upper boundary; s2 is at lower bndry.
            if (bounding_box%q(1,i,j,2) == zero) then ! No inner surface contact
               nlines_outer = nlines_outer + 1
               s1 = bounding_box%q(2,i,j,1) ! for outer surface
               s2 = bounding_box%q(3,i,j,1)
               nk = nint ((s2 - s1) / ds_nominal)
               nk = min (nmax, max (nk, nmin))
            else ! Inner & outer bndry. intersections; select lower line segment
               nlines_inner = nlines_inner + 1
               s1 = bounding_box%q(3,i,j,2) ! s2 for inner surface
               s2 = bounding_box%q(3,i,j,1) ! s2 for outer surface
               nk = nint ((s2 - s1) / ds_nominal)
               nk = min (nmax, max (nk, nmin))
               write (lunlog, '(a, 2i4, a, 4es14.6)')  &
                  ' i,j:', i, j, '  s1o, s1i, s2i, s2o:', &
                  bounding_box%q(2,i,j,1), bounding_box%q(2,i,j,2), s1, s2
            end if
         end if

         line_segments(i,j)%nk = nk

         call xyz_allocate (line_segments(i,j), ios)

         dx = bounding_box%x(i,j,2) - bounding_box%x(i,j,1)
         dy = bounding_box%y(i,j,2) - bounding_box%y(i,j,1)
         dz = bounding_box%z(i,j,2) - bounding_box%z(i,j,1)
         ds = (s2 - s1) / real (nk - 1)

         do k = 1, nk
            s = s1 + ds * real (k - 1)
            r = s / sline(nline)
            line_segments(i,j)%x(1,1,k) = bounding_box%x(i,j,1) + r * dx
            line_segments(i,j)%y(1,1,k) = bounding_box%y(i,j,1) + r * dy
            line_segments(i,j)%z(1,1,k) = bounding_box%z(i,j,1) + r * dz
         end do

      end do

   end do

!  Save all lines (including those that miss), 1 segment per grid block:

   nblocks = niLOS * njLOS

   if (formatted_lines) then
      write (lunLOS, '(i6)') nblocks
      write (lunLOS, '(i1, i2, i4)') &
         ((1, 1, line_segments(i,j)%nk, i = 1, niLOS), j = 1, njLOS)
      do j = 1, njLOS
         do i = 1, niLOS
            write (lunLOS, '(6es18.10)') line_segments(i,j)%x
            write (lunLOS, '(6es18.10)') line_segments(i,j)%y
            write (lunLOS, '(6es18.10)') line_segments(i,j)%z
         end do
      end do
   else
      write (lunLOS) nblocks
      write (lunLOS) &
         ((1, 1, line_segments(i,j)%nk, i = 1, niLOS), j = 1, njLOS)
      do j = 1, njLOS
         do i = 1, niLOS
            write (lunLOS) line_segments(i,j)%x, &
                           line_segments(i,j)%y, &
                           line_segments(i,j)%z
         end do
      end do
   end if

   close (lunLOS)

!  Print the association of each active line with an (i,j) and a pixel #:

   write (lunlog, '(/, a)') '    i    j pixel active'
   line = 0;  lactive = 0
   do j = 1, njLOS
      do i = 1, niLOS
         line = line + 1
         if (line_segments(i,j)%nk > 2) lactive = lactive + 1
         write (lunlog, '(2i5, 2i6)') i, j, line, lactive
      end do
   end do

   open (lunLI, file='lines_inner.g',   status='unknown')
   open (lunLO, file='lines_outer.g',   status='unknown')
   open (lunTI, file='lines_inner.txt', status='unknown')
   open (lunTO, file='lines_outer.txt', status='unknown')

!  Adjust nlines_outer?  Some lines may belong to the nose cap:

   if (side_view) then
      nlines_nosecap = 0
      do i = 1, niLOS
         if (bounding_box%q(1,i,1,1) > zero) then  ! Outer boundary intersection
            if (bounding_box%q(1,i,1,2) == zero) then  ! ... but no inner
               if (line_segments(i,1)%x(1,1,1) < xnose) then
                  nlines_nosecap = nlines_nosecap + 1
               end if
            end if
         end if
      end do
      nlines_outer = nlines_outer - nlines_nosecap
      if (nlines_nosecap > 0) then
         open (lunLN, file='lines_nosecap.g',   status='unknown')
         open (lunTN, file='lines_nosecap.txt', status='unknown')
         write (lunLN, '(i5)') nlines_nosecap
         write (lunTN, '(i5)') nlines_nosecap
      else
         side_view = false
      end if
   end if

   write (lunLI, '(i5)') nlines_inner
   write (lunLO, '(i5)') nlines_outer
   write (lunTI, '(i5)') nlines_inner
   write (lunTO, '(i5)') nlines_outer

   do j = 1, njLOS
      do i = 1, niLOS
         if (bounding_box%q(1,i,j,2) > zero) then ! Inner surface intersection
            write (lunLI, '(i1, i2, i4)') 1, 1, line_segments(i,j)%nk
            write (lunTI, '(i3, i4)') i, j
         else if (bounding_box%q(1,i,j,1) > zero) then ! Outer boundary only
            if (side_view .and. line_segments(i,1)%x(1,1,1) < xnose) then
               write (lunLN, '(i1, i2, i4)') 1, 1, line_segments(i,j)%nk
               write (lunTN, '(i3, i4)') i, j
            else
               write (lunLO, '(i1, i2, i4)') 1, 1, line_segments(i,j)%nk
               write (lunTO, '(i3, i4)') i, j
            end if
         end if
      end do
   end do

   close (lunTI)
   close (lunTO)
   if (side_view) close (lunTN)

   do j = 1, njLOS
      do i = 1, niLOS
         if (bounding_box%q(1,i,j,2) > zero) then ! Inner surface intersection
            write (lunLI, '(6es18.10)') line_segments(i,j)%x
            write (lunLI, '(6es18.10)') line_segments(i,j)%y
            write (lunLI, '(6es18.10)') line_segments(i,j)%z
         else if (bounding_box%q(1,i,j,1) > zero) then ! Outer boundary only
            luno = lunLO
            if (side_view .and. line_segments(i,1)%x(1,1,1) < xnose) then
               luno = lunLN
            end if
            write (luno, '(6es18.10)') line_segments(i,j)%x
            write (luno, '(6es18.10)') line_segments(i,j)%y
            write (luno, '(6es18.10)') line_segments(i,j)%z
         end if
      end do
   end do

   close (lunLI)
   close (lunLO)
   if (side_view) close (lunLN)

99 continue

   end program Stardust_Lines
