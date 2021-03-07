!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program lines_of_sight_2D

!  Description:
!
!     This is the 2D analogue of the earlier LINES_OF_SIGHT 3-space utility.
!     It can handle either an (x,y) input dataset with z = 0 assumed or an
!     (x,y,z) dataset with z assumed constant (and ignored).  See paragraph 4
!     below for the most recent enhancement in this version.
!
!     For a list of surface points (grid indices or (x,y) coordinates) and the
!     associated 2D volume grid, generate lines of sight - i.e., straight lines
!     normal to the surface and extending to the outer boundary with point dis-
!     tributions close to those of the local radial grid lines.  The results are
!     saved in PLOT2D multiblock form (one line of sight per block, ASCII)
!     compatible with FLOW_INTERP_2D, which can perform the flow interpolations
!     and tabulations that are normally what are really desired for application
!     to hypersonic flows.
!
!     This version has the option to produce lines parallel to Ox rather than
!     normal to the wall.  More recently, it also has the option to produce
!     lines normal to the outer shock boundary, which may or may not be the best
!     choice for tangent-slab radiation calculations.
!
!     Most recently, radio signal attenuation analyses have called for the
!     option to generate lines at arbitrary angles off the body (pointing from
!     transmitter to receiving satellite at some instant in time during a
!     spacecraft entry into an atmosphere).  In this case, the angle is taken
!     to be with Ox (positive meaning pointing above Ox; negative meaning
!     below Ox). The angle(s) should be entered in a third column: x, y, angle.
!
!     Always plot the results from LINES_OF_SIGHT[_2D] before proceeding with
!     related computations, especially if the geometry includes an aft body.
!
!  Assumptions (probably generalizable, but they may never need to be):
!
!     o  The structured 2D volume grid contains one layer of blocks, with j = 1
!        at the wall.  This simplifies determination of the inner and outer
!        boundary curves.  (To overcome these restrictions, one could use the
!        boundary condition data employed by the relevant flow solver.)
!
!     o  Normally the upper half of the geometry and 2D volume grid suffices
!        for axisymmetric bodies at zero angle of attack.  However, a body such
!        as a sample return capsule may be so non-convex that some lines of
!        sight may cross the X axis.  In such a case, use companion utility
!        REFLECT_BLOCKS_2D to reflect the grid and flow data before running
!        LINES_OF_SIGHT_2D and FLOW_INTERP_2D.
!
!  Strategy:
!
!     o  Prompt for all inputs (no control file).
!
!     o  Read the entire 2D volume grid (ASCII|binary, 2D|3D; any zs ignored
!        unless they are really angles for the radio attenuation case).
!
!     o  For all lines of sight, search the inner boundary and save the relevant
!        block number and cell indices if the body points are defined by (x,y)
!        coordinates.  (They may also be defined by a list of (block, i) pairs.)
!        This is now done with an ADT build of the body surface cells and an ADT
!        search for each target body point.
!
!     o  Release the inner surface search tree and build one for the outer
!        boundary, as needed for INTSEC2D (body-normal or Ox-parallel cases)
!        and also for the shock-normal case.
!
!     o  For each line of sight:
!
!          If body-normal:
!
!            > Construct a 2-point line normal to the wall with length that of a
!              local radial grid line.  This should be at least as long as the
!              straight line distance to the outer boundary.
!
!            > Intersect the line with the outer boundary and transform the
!              point distribution of the radial grid line to the relevant
!              portion of the straight line.
!
!          If shock-normal:
!
!            > Simply apply the SEARCH_ADT curve utility to each body point
!              and the outer grid boundary: this finds the closest point on the
!              shock boundary, and the associated line is orthogonal to that
!              boundary.  BEWARE OF THIS OPTION IF AN AFT BODY IS PRESENT.
!
!            > Discretize the 2-point line very simply.
!
!          If parallel to Ox:
!
!            > Adjust the body-normal method to work with unit vector (-1, 0)'
!              instead of the unit normal at the body point, and perform the
!              same intersection calculation and discretization.
!
!          If radio-attenuation case:
!
!            > The angle in column three defines a direction vector, and the
!              same intersection calculation and discretization is performed.
!
!  Input surface point format (ASCII, read to EOF):
!
!     Either     or        or (attenuation case)
!
!     n   i      x   y     x   y   angle
!     n   i      x   y     x   y   angle
!     n   i      x   y     x   y   angle
!     :   :      :   :     :   :   :
!
!     where n = block number and j = 1 is implied.
!
!     Note that trailing comments may be safely added to a body point list of
!     indices or coordinates.  The presence of a decimal point in the first
!     token of the LAST list line is used to distinguish the list as real,
!     with the exception of a lone '0', which is also interpreted as real.
!
!  XYZ Conventions:
!
!     Since the DPLR postprocessor extracts only x & y for a 2D grid, the input
!     volume grid to be interpolated may be either 2D/xy or 3D/xyz with z all 0.
!     Thus y is "up" for input and output files here.
!
!  Control:
!
!     A handful of prompts suffice.
!
!  History:
!
!     10/07/05  D.A.Saunders  Initial implementation of 3-space LINES_OF_SIGHT.
!     02/14/12    "     "     LINES_OF_SIGHT_2D adapted from the 3D form.
!     07/09/13    "     "     Dinesh Prabhu proposed making the lines of sight
!                             orthogonal to the shock as the proper thing to do
!                             for tangent-slab radiation calculations.  The
!                             earlier body-normal and Ox-parallel options have
!                             been retained, and the starting guesses for
!                             line-line intersections have been improved.
!     08/16/14    "     "     A trailing comment in the body-point file allowed
!                             indices to be interpreted as reals if the last
!                             line comment contained a period.  This has been
!                             remedied.
!     11/12/14    "     "     Tabulate x,y with body point indices to help
!                             plotting of radiative heating along surfaces.
!     04/08/15    "     "     Logic failure on a reflected axisymmetric 2D
!                             volume grid precipitated 2-space curve forms of
!                             the Alternating Digital Tree utilities and also
!                             implementation of INTSEC2D which uses those.
!                             Now, the 2-space ADT method applied to the outer
!                             grid boundary works for both the body-normal case
!                             (INTSEC2D) and the shock-normal case (SEARCH_ADT).
!                             The 2-space ADT method applied to the inner
!                             boundary also replaces NEAREST_CURVE_POINT.
!     05/07/19    "     "     Added the option to generate lines at any angle
!                             w.r.t. Ox for radio attenuation applications.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  --------

   use grid_block_structure  ! Employed by the XYZQ_IO package
   use xyzq_io_module        ! PLOT3D-type I/O package
   use adt_utilities         ! Search utilities, incl. 2-space curve variant
   use trigd

   implicit none

!  Global constants:
!  -----------------

   integer, parameter :: &
      lunpoints = 1,     &   ! Input list of surface points
      lunvolg   = 2,     &   ! Input volume grid
      lunout    = 3,     &   ! Output lines of sight
      lunkbd    = 5,     &   ! Keyboard inputs
      luncrt    = 6          ! Screen outputs

   real, parameter ::    &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Global variables:
!  -----------------

   integer :: &
      ios, nblocks, nlines

   logical :: &
      any_angle, body_normal, indices, parallel_to_Ox, shock_normal

   integer, allocatable, dimension (:,:) :: &
      n_and_i_target        ! For (block,i) of results of surface searches

   real, allocatable, dimension (:) :: &
      direction             ! Direction towards receiving antenna, off Ox, deg.

   real, allocatable, dimension (:,:) :: &
      pq, xy_target         ! (p,q)s & (x,y)s of target surface points

   type (grid_type), pointer, dimension (:) :: &
      radial_lines, volume_grid

!  Execution:
!  ----------

   call control ()          ! Prompt for and read input files, etc.
   if (ios /= 0) go to 99

   call construct_lines ()  ! Intersect 2-point normals at each target point
                            ! with the outer boundary and discretize them

   call save_results ()     ! Write the lines of sight and wrap up

99 continue

!  Internal procedures for program LINES_OF_SIGHT_2D:
!  --------------------------------------------------

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine control ()  ! Prompt for input/output files and other controls;
!                            ! read input files and set up output file.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      character, parameter :: &
         form*11 = 'unformatted'

!     Local variables:

      integer :: &
         l, n, npts

      logical :: &
         ascii, formatted, twod

      character :: &
         answer*1, buffer*80

!     Execution:

!     Lines of sight definition:

      write (luncrt, '(/, 2a)', advance='no') &
         ' Make lines body-normal (b), shock-normal (s)', &
         ' parallel to -Ox (x), or at given angles with Ox (a)?: '
      read (lunkbd, *) answer
      body_normal    = answer == 'b' .or. answer == 'B'
      shock_normal   = answer == 's' .or. answer == 'S'
      parallel_to_Ox = answer == 'x' .or. answer == 'X'
      any_angle      = answer == 'a' .or. answer == 'A'

!     Distinguish between a list of surface indices and a list of coordinates:

      ascii = true ! Can't pass a constant as an inout argument.

      call file_prompt (lunpoints, 0, 'list of surface points', 'old', false,  &
                        ascii, ios)
      if (ios /= 0) go to 99

      nlines = 0
      do ! Count the lines till EOF
         read (lunpoints, '(a)', iostat=ios) buffer ! We can look for '.' below
         if (ios /= 0) exit
         nlines = nlines + 1
      end do
      rewind (lunpoints)

!     Distinguish (n,i) block/point indices from (x,y) coordinates:

      call index_or_not (buffer, indices)

      write (luncrt, '(/, a, i4)') ' # lines of sight specified: ', nlines
      write (luncrt, '(a, l1)') ' indices specified? ', indices

      allocate (n_and_i_target(2,nlines), xy_target(2,nlines), pq(2,nlines))

      if (indices) then
         do l = 1, nlines ! Allow for trailing comments
            read (lunpoints, *) n_and_i_target(:,l)
         end do
      else if (.not. any_angle) then
         do l = 1, nlines  ! Avoid any z coordinates
            read (lunpoints, *) xy_target(:,l)
         end do
      else  ! Arbitrary angles from Ox
         allocate (direction(nlines))
         do l = 1, nlines  ! Avoid any trailing comments
            read (lunpoints, *) xy_target(:,l), direction(l)
         end do
      end if
      close (lunpoints)

!     Prompt for and open the input 2D volume grid and the output lines file:

      call file_prompt (lunvolg, 0, '2D volume grid', 'old', true, formatted,  &
                        ios)
      if (ios /= 0) go to 99

      call file_prompt (lunout, 0, 'output lines of sight', 'unknown', false,  &
                        ascii, ios)
      if (ios /= 0) go to 99

!     Read the 2D volume grid.  We eschew the packaged I/O utilities to allow
!     for the presence of [constant] Z or not.

      if (formatted) then
         read (lunvolg, *, iostat=ios) nblocks
      else
         read (lunvolg,    iostat=ios) nblocks
      end if
      if (ios /= 0) then
         write (luncrt, '(2a)') 'Trouble reading # input 2D volume grid blocks.'
         go to 99
      end if

      allocate (volume_grid(nblocks))

      twod = false
      if (formatted) then
         read (lunvolg, *, iostat=ios) &  ! Try for 3D first
            (volume_grid(n)%ni, volume_grid(n)%nj, volume_grid(n)%nk, &
            n = 1, nblocks)
         do n = 1, nblocks
            if (volume_grid(n)%nk /= 1) ios = 1  ! Shouldn't be necessary
         end do
         if (ios /= 0) then  ! Maybe it's 2D
            twod = true
            rewind (lunvolg);  read (lunvolg, *) ! Skip # blocks
            read (lunvolg, *, iostat=ios) &
               (volume_grid(n)%ni, volume_grid(n)%nj, n = 1, nblocks)
         end if
      else
         read (lunvolg, iostat=ios) &
            (volume_grid(n)%ni, volume_grid(n)%nj, volume_grid(n)%nk, &
            n = 1, nblocks)
         do n = 1, nblocks
            if (volume_grid(n)%nk /= 1) ios = 1  ! Shouldn't be necessary
         end do
         if (ios /= 0) then
            twod = true
            rewind (lunvolg);  read (lunvolg)
            read (lunvolg, iostat=ios) &
               (volume_grid(n)%ni, volume_grid(n)%nj, n = 1, nblocks)
         end if
      end if
      if (ios /= 0) then
         write (luncrt, '(2a)') &
            'Trouble reading input 2D volume grid block dimensions.'
         go to 99
      end if

      volume_grid(:)%nk = 1
      write (luncrt, '(a, /, (3i5))') &
         ' Input 2D volume grid block dimensions found:', &
         (volume_grid(n)%ni, volume_grid(n)%nj, 1, n = 1, nblocks)

!     Read the input 2D volume grid coordinates.
!     The intersection scheme now uses 3-space ADT techniques, so allocate %z
!     and set it to zero for BUILD_ADT:

      do n = 1, nblocks
         call xyz_allocate (volume_grid(n), ios)

         if (twod) then
            if (formatted) then
               read (lunvolg, *, iostat=ios) volume_grid(n)%x, volume_grid(n)%y
            else
               read (lunvolg,    iostat=ios) volume_grid(n)%x, volume_grid(n)%y
            end if
            if (ios /= 0) then
               write (luncrt, '(a, i4)') &
                  'Trouble reading (x,y) from 2D input volume; block #:', n
               go to 99
            end if
            volume_grid(n)%z = zero
         else
            npts = volume_grid(n)%ni * volume_grid(n)%nj
            call xyz_block_io (1, lunvolg, formatted, npts, volume_grid(n)%x,  &
                               volume_grid(n)%y, volume_grid(n)%z, ios)
            if (ios /= 0) then
               write (luncrt, '(a, i4)') &
                  'Trouble reading (x,y,z) from 3D input volume; block #:', n
               go to 99
            end if
         end if
      end do

      close (lunvolg)

99    continue

      end subroutine control

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine construct_lines ()
!
!     Body-normal case:
!
!     Originally, INTSEC2T and NEAREST_POINT_POINT were employed here, but
!     having to work with a reflected axisymmetric grid prompted introduction
!     of a 2-space curve variant of the Alternating Digital Tree rapid search
!     utilities, now used at both the inner boundary (if body points are given
!     as (x,y) coordinates) and the outer volume grid boundary.  The boundaries
!     no longer have to be assembled as continuous curves.
!
!     For each line of sight, we construct a normal to the inner surface
!     (j = 1 for all 2D volume grid blocks) and a straight line of length
!     equal to that of the nearest radial grid line.  INTSEC2D performs the
!     intersection with the outer boundary.  The point distribution of that
!     nearest grid line is then transformed to the straight line, forming one
!     block of the PLOT2D formatted output lines of sight file.
!
!     Shock-normal case:
!
!     Simply calculate the shortest distance from the body point to the outer
!     boundary via SEARCH_ADT, and discretize the two-point line.
!
!     Ox-parallel case:
!
!     Simply adjust the unit vector of the body-normal method then proceed.
!
!     Arbitrary angle case:
!
!     The angle defines a different unit direction vector.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      integer, parameter :: &
         nl = 2             ! INTSEC2T needed 3 line pts. as starting guess was
                            ! t at point (nl + 1)/2; now 2 points suffice
      real, parameter :: &
         big = 1.e+20, eps = 1.e-6, half = 0.5, one = 1.0

      character (1), parameter :: &
         method = 'L'       ! No need for monotonic spline fits along 2-pt. line

!     Local variables:

      integer :: &
         i, ib, icell, in, l, n, ncell, ni, ninside, nj, noutside

      integer, allocatable :: &
         conn(:,:)

      real :: &
         degrees, dmax, dmean, dsqmin, dtol, dtolsq, p, q, t, teval, tint, &
         total, xmin, xsmax, xy_interp(2), ymax, ymin

      real, dimension (nl) :: &
         tl, xl, yl, zl          ! For 2-point lines passed to INTSEC2D

      real, dimension (2) ::  &
         un                      ! For unit normals

      real, allocatable, dimension (:) :: &
         tline, xline, yline     ! For one radial volume grid line

!     Execution:

!     Count the number of 2-point boundary cells; determine a search tolerance:
!     -------------------------------------------------------------------------

      ymin  =  big
      ymax  = -big
      ncell =  0

      do ib = 1, nblocks
         ni = volume_grid(ib)%ni
         nj = volume_grid(ib)%nj
         ncell = (ni - 1) + ncell
         do i = 1, ni
            ymin = min (ymin, volume_grid(ib)%y(i,nj,1))
            ymax = max (ymax, volume_grid(ib)%y(i,nj,1))
         end do
      end do

      dtol   = (ymax - ymin) * eps
      dtolsq = dtol**2              ! For search diagnostics
      dtol   = dtol * 0.01          ! For axisymmetric body point on X axis

      allocate (conn(2,ncell))  ! BUILD_ADT sets 2-pt. cell block #s & i indices

!     Search the inner surface only if necessary:
!     -------------------------------------------

      if (indices) then ! No need to search the inner surface

!        However, the unit normal procedure expects cell indices to point to
!        the left end.

         write (luncrt, '(/, a)') '#  BP   ib    i           x,y'
         do l = 1, nlines
            ib = n_and_i_target(1,l)
            i  = n_and_i_target(2,l)
            xy_target(1,l) = volume_grid(ib)%x(i,1,1)
            xy_target(2,l) = volume_grid(ib)%y(i,1,1)

            write (luncrt, '(3i5, 2es16.8)') l, ib, i, xy_target(:,l)

            pq(1,l) = one
            if (i == volume_grid(ib)%ni) then
               n_and_i_target(2,l) = i - 1;  pq(1,l) = zero
            end if
            pq(2,l) = one - pq(1,l)
!!          write (luncrt, '(a, 3i5)') &
!!             'l, n_and_i_target(:,l):', l, n_and_i_target(:,l)
!!          write (luncrt, '(a, 2f15.8)') '    x/y_target(:,l):', xy_target(:,l)
!!          write (luncrt, '(a, 2f15.8)') '           p,q(:,l):', pq(:,l)
         end do

      else  ! Search the inner surface for the nearest pt. to each target pt.

         call build_adt (nblocks, volume_grid, ncell, 1, conn)

!        For all lines of sight, search the inner boundary
!        and save the relevant block number and cell index:
!        --------------------------------------------------

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         do l = 1, nlines  ! For each target surface point

            call search_adt (xy_target(:,l), icell, p, q, dsqmin, true, &
                             nblocks, volume_grid, ncell, 1, conn, xy_interp)

            if (dsqmin < dtolsq) then ! The nearest cell was within tolerance
               ninside  = ninside + 1
            else
               noutside = noutside + 1
            end if

            dmax  = max (dmax, dsqmin)
            dmean = dmean + sqrt (dsqmin)

            n_and_i_target(:,l) = conn(:,icell)
                 xy_target(:,l) = xy_interp(:)
                        pq(1,l) = p
                        pq(2,l) = q

         end do ! Next target surface point

         dmax = sqrt (dmax);  dmean = dmean / real (nlines)

         write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
            ' # surface points inside/outside tolerance:', ninside, noutside, &
            ' max & mean distance:', dmax, dmean

         call release_adt ()  ! An outer boundary search tree is needed next

      end if

      write (luncrt, '(/, 2a)') &
         'Line              x              y   ib    i', &
         '              p              q'
      write (luncrt,'(i4, 2f15.8, 2i5, 2f15.8)') &
         (l, xy_target(:,l), n_and_i_target(:,l), pq(:,l), l = 1, nlines)

!     Set up storage for the output lines of sight.
!     They are saved in (1,nj) multiblock PLOT2D form.
!     ------------------------------------------------

      allocate (radial_lines(nlines))

      nj = 1
      do l = 1, nlines
         radial_lines(l)%ni = 1
         ib = n_and_i_target(1,l)
         nj = max (nj, volume_grid(ib)%nj)
         radial_lines(l)%nj = nj         ! Expected to be all the same

         allocate (radial_lines(l)%x(1,nj,1), radial_lines(l)%y(1,nj,1))
      end do

      allocate (xline(nj), yline(nj), tline(nj))

!     Build a [new?] search tree for the outer boundary:
!     --------------------------------------------------

      call build_adt (nblocks, volume_grid, ncell, nj, conn)

!     Body pt. normals require surface arc lengths.  Store them in %z(:,1,1).
!     (The ADT scheme no longer needs %z.)

      do ib = 1, nblocks
         ni = volume_grid(ib)%ni
         nj = volume_grid(ib)%nj
         call chords2d (ni, volume_grid(ib)%x(:,1,1), &
                            volume_grid(ib)%y(:,1,1), false, total, &
                            volume_grid(ib)%z(:,1,1))        ! Unnormalized
      end do

!     Calculate each line of sight of the indicated type:
!     ---------------------------------------------------

      tl(1) = zero

      do l = 1, nlines

         n = n_and_i_target(1,l)  ! Block # for this target surface point
         i = n_and_i_target(2,l)  ! Corresponding cell left index
         p = pq (1,l)
         q = pq (2,l)
         in = i
         if (p > half) in = i + 1 ! Nearest radial line
         xl(1) = p * volume_grid(n)%x(i,1,1) + q * volume_grid(n)%x(i+1,1,1)
         yl(1) = p * volume_grid(n)%y(i,1,1) + q * volume_grid(n)%y(i+1,1,1)

         ni = volume_grid(n)%ni
         nj = volume_grid(n)%nj

         xline(:) = volume_grid(n)%x(in,:,1)  ! Nearby radial grid line to use
         yline(:) = volume_grid(n)%y(in,:,1)  ! for the discretization

         call chords2d (nj, xline, yline, true, total, tline) ! Normalized

         if (parallel_to_Ox) then  ! Retrofitted option most like body-normal

            un(1) = -one  ! Presumably not an aft body point
            un(2) = zero

         else if (any_angle) then  ! Retrofitted option for angle from Ox

            degrees = direction(l)
            un(1)   = cosd (degrees)
            un(2)   = sind (degrees)

         else if (body_normal) then

!           Construct a unit normal to the wall at the body point.
!           Wall segment volume_grid(n)%z(:,1,1) contains wall arc lengths.

            teval = p * volume_grid(n)%z(i,1,1) + q * volume_grid(n)%z(i+1,1,1)

            call unitnorm2d (method, ni, &
                             volume_grid(n)%x(:,1,1), volume_grid(n)%y(:,1,1), &
                             volume_grid(n)%z(:,1,1), teval, un)

            if (abs (xy_target(2,l)) < dtol) then
               un(1) = sign (one, un(1))
               un(2) = zero
            end if

            write (luncrt, '(/, a, i4, a, 2f15.8, a, es19.11)') &
               'Unit normal, BP', l, ':', un(:), '  teval along wall:', teval
         end if

         if (shock_normal) then  ! Simply find the nearest outer boundary point

            call search_adt (xy_target(:,l), icell, p, q, dsqmin, true, &
                             nblocks, volume_grid, ncell, nj, conn, xy_interp)

            xl(1) = xy_target(1,l);  yl(1) = xy_target(2,l)

            radial_lines(l)%x(1,:,1) = xl(1) + (xy_interp(1) - xl(1)) * tline(:)
            radial_lines(l)%y(1,:,1) = yl(1) + (xy_interp(2) - yl(1)) * tline(:)

         else  ! Body-normal, angle off Ox, or Ox-parallel method

!           Set up a 2-point straight line that should cross the outer bndry.:

            total = total + total
            if (un(1) >= zero) total = total * 5.     ! Wake? Reflected grid?

            xl(2) = xl(1) + total * un(1)
            yl(2) = yl(1) + total * un(2)
            tl(2) = total

!           Intersect the 2-point straight line with the outer boundary:

            call intsec2d (nblocks, volume_grid, ncell, nj, conn, &
                           nl, xl, yl, tl, method, icell, p, q, tint, &
                           xy_interp, dsqmin)

            write (luncrt, '(a, 3i10, a, es19.11)') &
               ' l, ib, i:', l, conn(:,icell), ' tint:', tint

!           Transform the point distribution of the radial grid line to the
!           relevant portion of the straight line of sight.

            t = tint / total
            radial_lines(l)%x(1,:,1) = xl(1) + ((xl(nl) - xl(1))*t) * tline(:)
            radial_lines(l)%y(1,:,1) = yl(1) + ((yl(nl) - yl(1))*t) * tline(:)

         end if

      end do ! Next line of sight

      end subroutine construct_lines

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine save_results ()
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: l

      write (lunout, '(i5)') nlines
      write (lunout, '(i1, i4)') (1, volume_grid(1)%nj, l = 1, nlines)

      do l = 1, nlines
         write (lunout, '(6es19.11)') radial_lines(l)%x, radial_lines(l)%y
      end do

      write (luncrt, '(/,a)') '*** Be sure to plot the computed lines of sight.'

      end subroutine save_results

   end program lines_of_sight_2D
