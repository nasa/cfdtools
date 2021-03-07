!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program cones_of_sight

!     CONES_OF_SIGHT is an adaptation of LINES_OF_SIGHT, q.v.  It retains the
!     functionality of its predecessor if the indicated cone angle is zero.
!     At the time of writing, a third variation working with lines defined by
!     points on a hemisphere is planned - see HEMISPHERES_OF SIGHT.
!
!     Most of the following is carried over from LINES_OF_SIGHT.
!
!  Description:
!
!     For a list of surface points (grid indices or (x,y,z) coordinates) and the
!     associated volume grid, generate lines of sight - i.e., straight lines
!     normal to the surface and extending to the outer boundary with point dis-
!     tributions close to those of the local radial grid lines.  The results are
!     saved in PLOT3D multiblock form (one line of sight per block) compatible
!     with the earlier FLOW_INTERP, which can perform the flow interpolations
!     and tabulations that are normally what are really desired for application
!     to hypersonic flows.
!
!     (Later:) If a cone angle > 0 is specified, the output blocks contain the
!     9 lines defined by the 9-point O(h^6) formula for integration on a circle
!     of radius h on pp. 891-892 of Abramowitz and Stegun.  The lines are in
!     clockwise order from the center out, so the corresponding quadrature
!     coefficients are 1/6, 1/24, 1/24, 1/24, 1/24, 1/6, 1/6, 1/6, 1/6.
!
!     This version can also read a structured surface grid for the target points
!     to simplify using a thinned form of the relevant volume grid's surface.
!
!  Rationale:
!
!     Generating the lines of sight directly within FLOW_INTERP is more awkward
!     than it might seem.  FLOW_INTERP's use of the VOLUME grid form of the ADT
!     search package precludes use of the SURFACE grid form needed here for both
!     the wall and the outer boundary.  (The subroutine names are the same for
!     all forms of the ADT search packages.)
!
!  Initial assumptions (probably generalizable, but they may never need to be):
!
!     o  The structured volume grid contains one layer of blocks, with k = 1
!        at the wall.  This simplifies determination of the inner and outer
!        boundary patches.  (To overcome these restrictions, one could use the
!        boundary condition data employed by the relevant flow solver.)
!
!  Strategy:
!
!     o  Read the entire volume grid and extract the inner and outer boundaries
!        as multiblock surface grids.
!
!     o  For all points defining lines of sight, search the inner boundary and
!        save the relevant patch number and cell indices.
!
!     o  Build a new search tree from the outer boundary.
!
!     o  For each inner surface point defining a line or lines of sight:
!
!        > Construct a two-point line normal to the wall with length that of
!          a local radial line.  This should be at least as long as the straight
!          line distance to the outer boundary.
!
!        > Intersect the line with the outer boundary and transform the point
!          distribution of the radial line to the relevant portion of the line.
!
!        > If the specified cone angle is not zero, derive 8 more lines of
!          sight from the first and impose the same relative distributions.
!
!  Input surface point format (read to EOF):
!
!     Either                                      or
!
!     n   i   j                                   x   y   z
!     n   i   j                                   x   y   z
!     n   i   j                                   x   y   z
!     :   :   :                                   :   :   :
!
!     where n = block number and k = 1 is implied.
!
!  History:
!
!     10/07/05  D.A.Saunders  Initial design of LINES_OF_SIGHT.
!     11/14/05    "     "     Fixed a bug in allocating radial_lines(ib).
!     11/22/05    "     "     Expanded the range of t for INTSEC6 to 2, not 1.1;
!                             Added the option to read a structured surface grid
!                             rather than a list of indices or coordinates.
!     03/06/06    "     "     Added the cone angle option.
!     03/07/06    "     "     Renamed it as CONES_OF_SIGHT.
!     01/31/14    "     "     All variants of the ADT routines have been merged
!                             into a module (generic build & search calls).
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Employed by the XYZQ_IO package
   use xyzq_io_module        ! PLOT3D-type I/O package
   use adt_utilities         ! All variants of the ADT build & search utilities
   use trigd

   implicit none

!  Local constants:

   integer, parameter :: &
      lunpoints = 1,     &   ! Input list of surface points
      lunsurfg  = 2,     &   ! Alternative input structured surface grid
      lunvolg   = 3,     &   ! Input volume grid
      lunout    = 4,     &   ! Output lines of sight (multiblock PLOT3D ASCII)
      lunkbd    = 5,     &   ! Keyboard inputs
      luncrt    = 6,     &   ! Screen
      ncone     = 9,     &   ! # lines per cone if cone_angle > 0
      nl        = 2          ! 2-point lines are passed to INTSEC6

   real, parameter ::    &
      half      = 0.5,   &
      one       = 1.0,   &
      zero      = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character, parameter :: &
      lcs_method * 1 = 'L'    ! But LCSFIT isn't needed (hook for general case)

!  Local variables:

   integer :: &
      i, ib, ic, in, ios, iquad, j, jc, jn, k, l, n, ni, nj, nb_target,        &
      nblocks, nf, ninside, nk, nlines, noutside, nquad, ntargets

   integer, allocatable :: &
      conn(:,:),           &  ! For (patch,i,j) of surface quads. to search
      nij(:,:)                ! For (patch,i,j) of results of searches

   real :: &
      cone_angle, dmax, dmean, dsqmin, dtolsq, p, q, t, tan_cone_angle, total, &
      xyz_interp(3), ymax, ymin

   real, dimension (nl) :: &
      tl, xl, yl, zl          ! For 2-point lines passed to INTSEC6

   real, dimension (3) ::  &
      un, vec1, vec2          ! For unit normals

   real, allocatable, dimension (:) :: &
      tline, xline, yline, zline       ! For one radial volume grid line

   real, allocatable, dimension (:,:) :: &
      pq, xyz                 ! (p,q)s & (x,y,z)s of list of surface points

   logical :: &
      ascii, cell_centered, cone, formatted, indices, structured_targets

   character :: &
      answer * 1, buffer * 64

   type (grid_type), pointer, dimension (:) :: &
      inner_surface, outer_surface, radial_lines, surface_grid, volume_grid

!  Execution:

!  Structured target points, or just a list of indices or coordinates?
!  -------------------------------------------------------------------

   write (luncrt, '(/, a)', advance='no') &
      ' Structured target points, or just a list? [s|l]: '
   read (lunkbd, *) answer
   structured_targets = answer == 's' .or. answer == 'S'

   write (luncrt, '(/, a)', advance='no') ' Cone angle >= 0 degrees: '
   read (lunkbd, *) cone_angle
   cone = cone_angle > zero

   if (structured_targets) then

!     Deal with a structured surface grid:
!     ------------------------------------

      ascii = true
      call file_prompt (lunsurfg, 0, '[formatted] target surface', 'old',      &
                        false, ascii, ios)
      if (ios /= 0) go to 99

!     Read all patches:

      call xyzq_read (lunsurfg, -lunsurfg, ascii, nb_target, nf,               &
                      cell_centered, surface_grid, ios)
      if (ios /= 0) go to 99

!     Convert the structured surface points to a list.
!     Should we try to suppress duplicates??  Not now.

      ntargets = 0
      do ib = 1, nb_target
         ntargets = ntargets + surface_grid(ib)%ni * surface_grid(ib)%nj
      end do

      allocate (nij(3,ntargets), xyz(3,ntargets), pq(2,ntargets))

      l = 0
      do ib = 1, nb_target
         do j = 1, surface_grid(ib)%nj
            do i = 1, surface_grid(ib)%ni
               l = l + 1
               xyz(1,l) = surface_grid(ib)%x(i,j,1)
               xyz(2,l) = surface_grid(ib)%y(i,j,1)
               xyz(3,l) = surface_grid(ib)%z(i,j,1)
            end do
         end do
         deallocate (surface_grid(ib)%x, surface_grid(ib)%y, surface_grid(ib)%z)
      end do

      indices = false ! To reuse original code below

   else

!     Deal with a list of surface indices or coordinates:
!     ---------------------------------------------------

      ascii = true ! Can't pass a constant as an inout argument.

      call file_prompt (lunpoints, 0, 'list of surface points', 'old', false,  &
                        ascii, ios)
      if (ios /= 0) go to 99

      ntargets = 0
      do ! Count the target points till EOF
         read (lunpoints, '(a)', iostat=ios) buffer ! We can look for '.' below
         if (ios /= 0) exit
         ntargets = ntargets + 1
      end do
      rewind (lunpoints)

!     Distinguish (n,i,j) indices from (x,y,z) coordinates:

      indices = index (buffer, '.') == 0

      write (luncrt, '(/, a, i4)') ' # target points specified: ', ntargets
      write (luncrt, '(a, l1)') ' indices specified? ', indices

      allocate (nij(3,ntargets), xyz(3,ntargets), pq(2,ntargets))

      if (indices) then
         do l = 1, ntargets ! Allow for trailing comments
            read (lunpoints, *) nij(:,l)
         end do
      else
         do l = 1, ntargets
            read (lunpoints, *) xyz(:,l)
         end do
      end if
      close (lunpoints)

   end if

   if (cone) then
      nlines = ntargets * ncone
      tan_cone_angle = tand (cone_angle)
   else
      nlines = ntargets
   end if

   write (luncrt, '(/, a, i4)') ' # lines of sight indicated: ', nlines

!  Deal with the volume grid (and initializing the output results file):
!  ---------------------------------------------------------------------

   call file_prompt (lunvolg, 0, 'volume grid', 'old', true, formatted, ios)

   if (ios /= 0) go to 99

   call file_prompt (lunout, 0, 'output lines of sight', 'unknown', false,     &
                     ascii, ios)

!  Read the entire volume grid:

   call xyzq_read (lunvolg, -lunvolg, formatted, nblocks, nf, cell_centered,   &
                   volume_grid, ios)
   if (ios /= 0) go to 99

!  Extract the inner and outer boundaries as multiblock surface grids.

   allocate (inner_surface(nblocks), outer_surface(nblocks))

   do ib = 1, nblocks
      inner_surface(ib)%ni = volume_grid(ib)%ni
      outer_surface(ib)%ni = volume_grid(ib)%ni
      inner_surface(ib)%nj = volume_grid(ib)%nj
      outer_surface(ib)%nj = volume_grid(ib)%nj
      inner_surface(ib)%nk = 1
      outer_surface(ib)%nk = 1
   end do

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
   end do

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

      do l = 1, ntargets
         ib = nij(1,l)
         i  = nij(2,l)
         j  = nij(3,l)
         xyz(1,l) = inner_surface(ib)%x(i,j,1)
         xyz(2,l) = inner_surface(ib)%y(i,j,1)
         xyz(3,l) = inner_surface(ib)%z(i,j,1)
         pq (1,l) = zero
         pq (2,l) = zero
      end do

   else ! Build a search tree for the inner surface and do the searching:

      call build_adt (nblocks, inner_surface, nquad, conn)

!     For all target points, search the inner boundary
!     and save the relevant patch number and cell indices:
!     ----------------------------------------------------

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      do l = 1, ntargets

         call search_adt (xyz(:,l), iquad, p, q, dsqmin, true, nblocks,        &
                          inner_surface, nquad, conn, xyz_interp)

         if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
            ninside  = ninside + 1
         else
            noutside = noutside + 1
         end if

         dmax  = max (dmax, dsqmin)
         dmean = dmean + sqrt (dsqmin)

         nij(:,l) = conn(:,iquad)
         xyz(:,l) = xyz_interp(:) ! This is what FLOW_INTERP should also get
         pq (1,l) = p
         pq (2,l) = q

      end do ! Next target surface point

      dmax = sqrt (dmax);  dmean = dmean / real (ntargets)

      write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
         ' # surface points inside/outside tolerance:', ninside, noutside,     &
         ' max & mean distance:', dmax, dmean

   end if

!  Set up storage for the output lines of sight.
!  They are saved in (1,1,nk) or (ncone,1,nk) multiblock PLOT3D form.
!  ------------------------------------------------------------------

   allocate (radial_lines(ntargets), stat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a, i10)') &
         ' Allocation trouble with radial_lines. ntargets:', ntargets
      go to 99
   end if

   nk = 1
   do l = 1, ntargets
      radial_lines(l)%ni = 1;  if (cone) radial_lines(l)%ni = ncone
      radial_lines(l)%nj = 1
      ib = nij(1,l)
      nk = max (nk, volume_grid(ib)%nk)
      radial_lines(l)%nk = nk        ! Probably all the same

      call xyz_allocate (radial_lines(l), ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i10)') &
            ' Allocation trouble with radial_lines(l). l:', l
         write (luncrt, '(/, a, 3i10)') ' ni, nj, nk: ', &
            radial_lines(l)%ni, radial_lines(l)%nj, nk
         go to 99
      end if

   end do

   allocate (xline(nk), yline(nk), zline(nk), tline(nk), stat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a, i10)') &
         ' Allocation trouble with x/y/z/tline. nk:', nk
      go to 99
   end if

!  Build a [new?] search tree from the outer boundary patches
!  and intersect each line of sight with the outer boundary:
!  ----------------------------------------------------------

   if (.not. indices) call release_adt ()

   call build_adt (nblocks, outer_surface, nquad, conn)

   ninside  = 0;  dmax  = zero
   noutside = 0;  dmean = zero

   do l = 1, ntargets

      n = nij(1,l)  ! Patch # at inner surface for this target surface point
      i = nij(2,l)  ! Corresponding cell indices (lower left)
      j = nij(3,l)
      p = pq (1,l); in = i; if (p > half) in = i + 1 ! Nearest radial line
      q = pq (2,l); jn = j; if (q > half) jn = j + 1 ! from indicated cell

!!!   write (luncrt, '(a, 3i5, 2f20.15)') ' n, i, j, p, q: ', n, i, j, p, q

      nk = volume_grid(n)%nk

      do k = 1, nk  ! Nearby radial grid line to use for the point distribution
         xline(k) = volume_grid(n)%x(in,jn,k)
         yline(k) = volume_grid(n)%y(in,jn,k)
         zline(k) = volume_grid(n)%z(in,jn,k)
      end do

      call chords3d (nk, xline, yline, zline, true, total, tline) ! Normalized

!     Construct a two-point line normal to the wall with length that of a
!     local radial line.  This should be at least as long as the straight
!     line distance to the outer boundary.
!     First, a carefully calculated unit normal:

      call surface_normal (inner_surface(n)%ni, inner_surface(n)%nj,           &
                           inner_surface(n)%x,  inner_surface(n)%y,            &
                           inner_surface(n)%z,  i, j, p, q, un)

      xl(1) = xyz(1,l);  xl(2) = xl(1) + total * un(1)
      yl(1) = xyz(2,l);  yl(2) = yl(1) + total * un(2)
      zl(1) = xyz(3,l);  zl(2) = zl(1) + total * un(3)

!     Intersect the 2-point line with the outer boundary:

      tl(1)  = 0.5  ! Tell INTSEC6 the range of normalized t to search in
      tl(nl) = 3.0

      call intsec6 (nblocks, outer_surface, nquad, conn, nl, xl, yl, zl, tl,   &
                    lcs_method, iquad, p, q, t, xyz_interp, dsqmin)

!!!   write (luncrt, '(a, 4i10)') ' iquad, conn(:,iquad):', iquad, conn(:,iquad)
!!!   write (luncrt, '(a, 3es19.11)') ' p, q, t:', p, q, t

      if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
         ninside  = ninside + 1
      else
         noutside = noutside + 1
      end if

      dmax  = max (dmax, dsqmin)
      dmean = dmean + sqrt (dsqmin)

!     Transform the point distribution of the radial grid line to the relevant
!     portion of the straight line of sight.

      xl(2) = (xl(2) - xl(1)) * t ! These are now the components of the length
      yl(2) = (yl(2) - yl(1)) * t
      zl(2) = (zl(2) - zl(1)) * t

      do k = 1, nk
         radial_lines(l)%x(1,1,k) = xl(1) + xl(2) * tline(k)
         radial_lines(l)%y(1,1,k) = yl(1) + yl(2) * tline(k)
         radial_lines(l)%z(1,1,k) = zl(1) + zl(2) * tline(k)
      end do

!     Derive further lines of sight for the specified cone angle?

      if (cone) then

         call cone_lines ()     ! Internal procedure below

      end if

   end do ! Next target point

   dmax = sqrt (dmax);  dmean = dmean / real (nlines)

   write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
      ' # outer intersection points inside/outside tolerance:', &
      ninside, noutside, ' max & mean distance:', dmax, dmean

!  Save the lines of sight as a PLOT3D formatted file:
!  ---------------------------------------------------

   call xyz_write (lunout, true, ntargets, radial_lines, ios)

!  Let Fortran 90 do all the deallocates.

99 continue ! Avoid system dependencies where STOP is concerned

!  Internal procedure for program CONES_OF_SIGHT:
!  ----------------------------------------------

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine cone_lines ()

!     Derive further lines of sight from the current normal line.  See the
!     program header for the chosen scheme (which may change).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

!     Local variables:

      integer :: i
      real    :: alpha, axis(3), d, radius
      real, dimension (2)     :: xs, ys, zs, ts
      real, dimension (ncone) :: xc, yc, zc, ycircle, zcircle

!     Storage:

!     Quadrature points in the yz plane  (xcircle will be primary length t):

      data ycircle / 0., 0.5,  0.5, -0.5, -0.5, 0., 1.,  0., -1. /
      data zcircle / 0., 0.5, -0.5, -0.5,  0.5, 1., 0., -1.,  0. /

!     Execution:

!!!   write (6, '(/, a, 3es16.8)') ' Principal inner xyz: ', &
!!!radial_lines(l)%x(1,1,1), radial_lines(l)%y(1,1,1), radial_lines(l)%z(1,1,1)
!!!   write (6, '(   a, 3es16.8)') ' Principal outer xyz: ', &
!!!radial_lines(l)%x(1,1,nk),radial_lines(l)%y(1,1,nk),radial_lines(l)%z(1,1,nk)

!     Direction angle of (shifted) primary line with x-axis:

      d = sqrt (xl(2)**2 + yl(2)**2 + zl(2)**2)
      alpha = acosd (xl(2) / d)

!     An axis of rotation for aligning Ox with the primary line is given
!     by the cross product of Ox with that line shifted to the origin:
!     This is degenerate, though if the primary line is already along Ox.

      axis(1) = zero
      axis(2) = -d * zl(2)
      axis(3) =  d * yl(2)

      radius = d * tan_cone_angle

!!!   write (6, '(/, a, 3es16.8)') ' d, alpha, radius:', d, alpha, radius
!!!   write (6, '(/, a, 3es16.8)') ' axis(1:3):       ', axis

      do i = 1, ncone
         xc(i) = d
         yc(i) = ycircle(i) * radius
         zc(i) = zcircle(i) * radius
      end do

!!!   write (6, '(/, a, /, (i2, 3es16.8))') ' Unrotated i, xc, yc, zc', &
!!!      (i, xc(i), yc(i), zc(i), i = 1, ncone)

      if (abs (alpha) > 0.1) then  ! Avoid a degenerate rotation axis

         call rotate3d (ncone, xc, yc, zc, alpha, zero, zero, zero, & ! In place
                        axis(1), axis(2), axis(3))
      end if

!!!   write (6, '(/, a, /, (i2, 3es16.8))') '   Rotated i, xc, yc, zc', &
!!!      (i, xc(i), yc(i), zc(i), i = 1, ncone)

!     For each secondary line in the cone, repeat the intersection with the
!     outer boundary, etc., as for the primary line:

      ts(1) = 0.5  ! Fractions of secondary line lengths to search in
      ts(2) = 10.0 ! Err on the high side

      do i = 2, ncone
         xs(1) = xl(1);  xs(2) = xl(1) + xc(i)
         ys(1) = yl(1);  ys(2) = yl(1) + yc(i)
         zs(1) = zl(1);  zs(2) = zl(1) + zc(i)

!!!      write (6, '(/, a, i2, 3es16.8)') ' Secondary inner i, xyz:', i, &
!!!         xs(1), ys(1), zs(1)
!!!      write (6, '(   a, i2, 3es16.8)') ' Secondary outer i, xyz:', i, &
!!!         xs(2), ys(2), zs(2)

         call intsec6 (nblocks, outer_surface, nquad, conn, nl, xs, ys, zs, ts,&
                       lcs_method, iquad, p, q, t, xyz_interp, dsqmin)

!!!      write (luncrt, '(a, i2, 4i10)') &
!!!         ' secondary i, iquad, conn(:,iquad):', i, iquad, conn(:,iquad)
!!!      write (luncrt, '(a, 3es19.11)') ' p, q, t:   ', p, q, t, &
!!!                                      ' xyz_interp:', xyz_interp

         if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
            ninside  = ninside + 1
         else
            noutside = noutside + 1
         end if

         dmax  = max (dmax, dsqmin)
         dmean = dmean + sqrt (dsqmin)

!        Transform the distribution of the primary line to the secondary line:

         xs(2) = (xs(2) - xs(1)) * t ! Length components
         ys(2) = (ys(2) - ys(1)) * t
         zs(2) = (zs(2) - zs(1)) * t

         do k = 1, nk
            radial_lines(l)%x(i,1,k) = xs(1) + xs(2) * tline(k)
            radial_lines(l)%y(i,1,k) = ys(1) + ys(2) * tline(k)
            radial_lines(l)%z(i,1,k) = zs(1) + zs(2) * tline(k)
         end do

      end do

      end subroutine cone_lines

   end program cones_of_sight
