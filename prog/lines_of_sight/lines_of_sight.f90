!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program lines_of_sight

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
!     This version can also read a structured surface grid for the target points
!     to simplify using a thinned form of the relevant volume grid's surface.
!
!     This version has the option to produce lines parallel to Ox rather than
!     normal to the wall.  Most recently, it also has the option to produce
!     lines normal to the outer shock boundary, which may be the best choice
!     for tangent-slab radiation calculations.
!
!     For the Ox-parallel option, the targets are assumed to be in list form.
!
!     It is understood that body points will normally be confined to forebodies,
!     since aft-body points are unlikely to have solutions in the usual sense.
!
!  Assumptions (probably generalizable, but they may never need to be):
!
!     o  The structured volume grid contains one layer of blocks, with k = 1
!        at the wall.  This simplifies determination of the inner and outer
!        boundary patches.  (To overcome these restrictions, one could use the
!        boundary condition data employed by the relevant flow solver.)
!
!  Strategy:
!
!     o  Prompt for all inputs (no control file).
!
!     o  Read the entire volume grid and extract the inner and outer boundaries
!        as multiblock surface grids.
!
!     o  For all lines of sight, search the inner boundary and save the relevant
!        patch number and cell indices.
!
!     o  Build a new search tree from the outer boundary.
!
!     o  For each line of sight:
!
!          If body-normal:
!
!            > Construct a 2-point line normal to the wall with length that of a
!              local radial line.  This should be at least as long as the
!              straight line distance to the outer boundary.  (Actually, it may
!              need to be longer off the shoulder of a capsule.)
!
!            > Intersect the line with the outer boundary and transform the
!              point distribution of the radial grid line to the relevant
!              portion of the straight line.
!
!          If shock-normal:
!
!            > Simply apply the ADT search utility to each body point and the
!              outer grid boundary: this finds the closest point on the shock
!              boundary, and the associated line is orthogonal to that boundary.
!
!            > Discretize the 2-point line very simply.
!
!          If parallel to Ox:
!
!            > Adjust the body-normal method to work with unit vector (-1,0,0)'
!              instead of the unit normal at the body point, and perform the
!              same intersection calculation and discretization.
!
!  Input surface point format (read to EOF):
!
!     Either                            or                          or
!
!     n   i   j                         x   y   z             nb
!     n   i   j                         x   y   z             ni  nj  1
!     n   i   j                         x   y   z             ni  nj  1
!     :   :   :                         :   :   :             :   :   :
!                                                             x11 x12 x13 ...
!     where n = block number and k = 1 is implied.            ...............
!
!     Note that trailing comments may be safely added to a body point list of
!     indices or coordinates.  The presence of a decimal point in the first
!     token of the LAST list line is used to distinguish the list as real,
!     with the exception of a lone '0', which is also interpreted as real.
!
!  History:
!
!     10/07/05  D.A.Saunders  Initial design.
!     11/14/05    "     "     Fixed a bug in allocating radial_lines(ib).
!     11/22/05    "     "     Expanded the range of t for INTSEC6 to 2, not 1.1;
!                             Added the option to read a structured surface grid
!                             rather than a list of indices or coordinates.
!     08/21/09    "     "     If a point entered via patch indices were at an
!                             upper index boundary, the surface normal utility
!                             was not being given the lower-left cell indices.
!     07/11/13    "     "     Dinesh Prabhu proposed making the lines of sight
!                             orthogonal to the shock as the proper thing to do
!                             for tangent-slab radiation calculations.  The
!                             earlier body-normal and Ox-parallel options have
!                             been retained.  (Later: Shock-normal is NOT the
!                             recommended usage.)
!     08/06/13    "     "     All ADT variants have been combined into a module
!                             for distribution reasons (generic build & search).
!     08/16/14    "     "     A trailing comment in the body-point file allowed
!                             indices to be interpreted as reals if the last
!                             line comment contained a period.  This has been
!                             remedied.
!     11/11/14    "     "     Tabulate x,y,z with body point indices to help
!                             plotting of radiative heating along surfaces.
!     08/25/17    "     "     A 20-meter search range was occasionally not
!                             enough.  Make it 20 diameters as it should have
!                             been from the start for aft body points.
!     09/05/17    "     "     The calculation of diameter was flawed for left
!                             half body.
!     10/22/18    "     "     The unstructured analogue, USLOS, needs the
!                             release_adt call that has been missing here all
!                             along.  It surely belongs.  Also: use the outer
!                             bounding box diagonal to set the upper limit for
!                             the intersections instead of heuristics involving
!                             the body diameter.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!                           ERC, Inc. at ARC (08/2010 through 06/2015).
!                           AMA, Inc. at ARC (from 07/2015).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Employed by the XYZQ_IO & ADT packages
   use adt_utilities         ! All variants of the ADT build & search routines
   use xyzq_io_module        ! PLOT3D-type I/O package

   implicit none

!  Local constants:
!  ----------------

   integer, parameter :: &
      luntarget = 1,     &   ! Target body points in one of three forms
      lunvolg   = 2,     &   ! Input volume grid
      lunout    = 3,     &   ! Output lines of sight (multiblock PLOT3D ASCII)
      lunkbd    = 5,     &   ! Keyboard inputs
      luncrt    = 6,     &   ! Screen
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
!  ----------------

   integer :: &
      i, ib, ic, in, ios, iquad, j, jc, jn, k, l, n, ni, nj, nb_target,        &
      nblocks, nf, ninside, nk, nlines, noutside, nquad

   integer, allocatable :: &
      conn(:,:),           &  ! For (patch,i,j) of surface quads. to search
      nij(:,:)                ! For (patch,i,j) of results of searches

   real :: &
      bbox_diagonal, dmax, dmean, dsqmin, dtolsq, p, q, t, total, &
      xyz_interp(3), xmax, xmin, ymax, ymin, zmax, zmin

   real, dimension (nl) :: &
      tl, xl, yl, zl          ! For 2-point lines passed to INTSEC6

   real, dimension (3) ::  &
      un, vec1, vec2          ! For unit normals

   real, allocatable, dimension (:) :: &
      tline, xline, yline, zline       ! For one radial volume grid line

   real, allocatable, dimension (:,:) :: &
      pq, xyz                 ! (p,q)s & (x,y,z)s of list of surface points

   logical :: &
      ascii, body_normal, cell_centered, formatted, indices, parallel_to_Ox,   &
      shock_normal, structured_targets

   character :: &
      answer*1, buffer*80

   type (grid_type), pointer, dimension (:) :: &
      inner_surface, outer_surface, radial_lines, surface_grid, volume_grid

!  Execution:
!  ----------

   write (luncrt, '(/, 2a)', advance='no') &
      ' Make lines body-normal (b), shock-normal (s), or parallel to -Ox (x)?: '
   read (lunkbd, *) answer
   body_normal    = answer == 'b' .or. answer == 'B'
   shock_normal   = answer == 's' .or. answer == 'S'
   parallel_to_Ox = answer == 'x' .or. answer == 'X'

!  Avoid a prompt for distinguishing the target point file type.
!  We also avoid the available token_count utility by looking for just 1 token.

   ascii = true
   call file_prompt (luntarget, 0, 'target body point', 'old', &
                     false, ascii, ios)
   if (ios /= 0) go to 99

   read (luntarget, '(a)') buffer
   buffer = adjustl (buffer);  l = len_trim (buffer)
   structured_targets = index (buffer(1:l), ' ')      == 0 .and. & ! No embedded
                        index (buffer(1:l), char (9)) == 0         ! blank | tab
   rewind (luntarget)

   if (structured_targets) then  ! Must have been only 1 item on line 1: # blks.

      call xyzq_read (luntarget, -luntarget, ascii, nb_target, nf,             &
                      cell_centered, surface_grid, ios)
      if (ios /= 0) go to 99

!     Convert the structured surface points to a list.
!     Should we try to suppress duplicates??  Not now.

      nlines = 0
      do ib = 1, nb_target
         nlines = nlines + surface_grid(ib)%ni * surface_grid(ib)%nj
      end do
      write (luncrt, '(/, a, i4)') ' # lines of sight indicated: ', nlines

      allocate (nij(3,nlines), xyz(3,nlines), pq(2,nlines))

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

   else  ! Deal with a list of surface indices or coordinates

      nlines = 0
      do ! Count the lines till EOF
         read (luntarget, '(a)', iostat=ios) buffer ! We can look for '.' below
         if (ios /= 0) exit
         nlines = nlines + 1
      end do
      rewind (luntarget)

!     Distinguish (n,i,j) block/point indices from (x,y,z) coordinates:

      call index_or_not (buffer, indices)

      write (luncrt, '(/, a, i4)') ' # lines of sight specified: ', nlines
      write (luncrt, '(a, l1)') ' indices specified? ', indices

      allocate (nij(3,nlines), xyz(3,nlines), pq(2,nlines))

      if (indices) then
         do l = 1, nlines ! Allow for trailing comments
            read (luntarget, *) nij(:,l)
         end do
      else
         do l = 1, nlines
            read (luntarget, *) xyz(:,l)
         end do
      end if
      close (luntarget)

   end if

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

   call get_data_range ()  ! Bounding box diagonal, for intersection purposes

!  Count the surface quads (same for inner and outer surfaces):
!  ------------------------------------------------------------

   nquad = 0
   ymin = 1.e+30;  ymax = -ymin
   zmin = ymin;    zmax = -ymin

   do ib = 1, nblocks
      ni = inner_surface(ib)%ni
      nj = inner_surface(ib)%nj
      nquad = (ni - 1) * (nj - 1) + nquad

      do j = 1, nj
         do i = 1, ni
            ymin = min (ymin, inner_surface(ib)%y(i,j,1))
            ymax = max (ymax, inner_surface(ib)%y(i,j,1))
            zmin = min (zmin, inner_surface(ib)%z(i,j,1))
            zmax = max (zmax, inner_surface(ib)%z(i,j,1))
         end do
      end do
   end do

   dtolsq = (bbox_diagonal*1.e-6)**2          ! Tolerance for search diagnostics

   allocate (conn(3,nquad)) ! For (patch,i,j) of each surface quad.

!  Search the inner surface only if necessary:
!  -------------------------------------------

   if (indices) then ! No need to search the inner surface

!     However, the unit normal routine expects cell indices to point to
!     the "lower left" corner.

      write (luncrt, '(/, a)') '#  BP   ib    i    j           x,y,z'
      do l = 1, nlines
         ib = nij(1,l)
         i  = nij(2,l)
         j  = nij(3,l)
         xyz(1,l) = inner_surface(ib)%x(i,j,1)
         xyz(2,l) = inner_surface(ib)%y(i,j,1)
         xyz(3,l) = inner_surface(ib)%z(i,j,1)

         write (luncrt, '(4i5, 3es16.8)') l, nij(:,l), xyz(:,l)

         pq (1,l) = zero
         pq (2,l) = zero
         if (i == inner_surface(ib)%ni) then
            nij(2,l) = i - 1;  pq(1,l) = one
         end if
         if (j == inner_surface(ib)%nj) then
            nij(3,l) = j - 1;  pq(2,l) = one
         end if
      end do

   else ! Build a search tree for the inner surface and do the searching:

      call build_adt (nblocks, inner_surface, nquad, conn)

!     For all lines of sight, search the inner boundary
!     and save the relevant patch number and cell indices:
!     ----------------------------------------------------

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      do l = 1, nlines

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

      dmax = sqrt (dmax);  dmean = dmean / real (nlines)

      write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
         ' # surface points inside/outside tolerance:', ninside, noutside,     &
         ' max & mean distance:', dmax, dmean

      call release_adt ()  ! This has long been overlooked
   end if

!! write (luncrt, '(a, 3i5)') 'nij(1:3,1):', nij(1:3,1)
!! write (luncrt, '(a, 3f15.8)') 'x,y,z(1):', xyz(1:3,1)
!! write (luncrt, '(a, 2f15.8)') '  p,q(1):',  pq(1:2,1)

!  Set up storage for the output lines of sight.
!  They are saved in (1,1,nk) multiblock PLOT3D form.
!  --------------------------------------------------

   allocate (radial_lines(nlines), stat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a, i10)') &
         ' Allocation trouble with radial_lines. nlines:', nlines
      go to 99
   end if

   nk = 1
   do l = 1, nlines
      radial_lines(l)%ni = 1
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

!  Build a [new?] search tree from the outer boundary patches:
!  ----------------------------------------------------------

   call build_adt (nblocks, outer_surface, nquad, conn)

!  Calculate each line of sight of the indicated type:
!  ---------------------------------------------------

   ninside  = 0;  dmax  = zero
   noutside = 0;  dmean = zero

   do l = 1, nlines

      n = nij(1,l)  ! Patch # at inner surface for this target surface point
      i = nij(2,l)  ! Corresponding cell indices (lower left)
      j = nij(3,l)
      p = pq (1,l); in = i; if (p > half) in = i + 1 ! Nearest radial line
      q = pq (2,l); jn = j; if (q > half) jn = j + 1 ! to target point

!!    write (luncrt, '(a, 3i5, 2f20.15)') ' n, i, j, p, q: ', n, i, j, p, q

      nk = volume_grid(n)%nk

      do k = 1, nk  ! Nearby radial grid line to use for the point distribution
         xline(k) = volume_grid(n)%x(in,jn,k)
         yline(k) = volume_grid(n)%y(in,jn,k)
         zline(k) = volume_grid(n)%z(in,jn,k)
      end do

      call chords3d (nk, xline, yline, zline, true, total, tline) ! Normalized

      if (parallel_to_Ox) then  ! Retrofitted kludge

         un(1) = -one  ! ? What about an aft-body point?
         un(2) = zero
         un(3) = zero

      else if (body_normal) then

!        Construct a two-point line normal to the wall with length that of a
!        local radial line.  This should be at least as long as the straight
!        line distance to the outer boundary.
!        First, a carefully calculated unit normal:

         call surface_normal (inner_surface(n)%ni, inner_surface(n)%nj,        &
                              inner_surface(n)%x,  inner_surface(n)%y,         &
                              inner_surface(n)%z,  i, j, p, q, un)
      end if

      if (shock_normal) then  ! Simply find the nearest outer boundary point

         call search_adt (xyz(:,l), iquad, p, q, dsqmin, true, nblocks,        &
                          outer_surface, nquad, conn, xyz_interp)

          xl(1) = xyz(1,l);  yl(1) = xyz(2,l);  zl(1) = xyz(3,l)

         do k = 1, nk
            radial_lines(l)%x(1,1,k) = xl(1) + (xyz_interp(1) - xl(1))*tline(k)
            radial_lines(l)%y(1,1,k) = yl(1) + (xyz_interp(2) - yl(1))*tline(k)
            radial_lines(l)%z(1,1,k) = zl(1) + (xyz_interp(3) - zl(1))*tline(k)
         end do

      else  ! ! Body-normal or Ox-parallel method

         xl(1) = xyz(1,l);  xl(2) = xl(1) + bbox_diagonal*un(1)
         yl(1) = xyz(2,l);  yl(2) = yl(1) + bbox_diagonal*un(2)
         zl(1) = xyz(3,l);  zl(2) = zl(1) + bbox_diagonal*un(3)

!        Intersect the 2-point line with the outer boundary:

         tl(1)  = zero           ! Tell INTSEC6 the range of t to search in
         tl(nl) = bbox_diagonal  ! These are NOT redundant if nl = 2

         call intsec6_2pt (nblocks, outer_surface, nquad, conn, nl, xl, yl, zl,    &
                           tl, lcs_method, iquad, p, q, t, xyz_interp, dsqmin)

!        The output t is normalized.
!!       write (luncrt, '(a,4i10)') 'iquad,conn(:,iquad):', iquad,conn(:,iquad)
!!       write (luncrt, '(a,4es19.11)') 'p, q, t, bbox:', p, q, t, bbox_diagonal
!!       write (luncrt, '(a,4es19.11)') 'xyzint, dsqmin:', xyz_interp, dsqmin

         if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
            ninside  = ninside + 1
         else
            noutside = noutside + 1
         end if

         dmax  = max (dmax, dsqmin)
         dmean = dmean + sqrt (dsqmin)

!        Transform the point distribution of the radial grid line to the
!        relevant portion of the straight line of sight:

         do k = 1, nk
            radial_lines(l)%x(1,1,k) = xl(1) + ((xl(2) - xl(1))*t)*tline(k)
            radial_lines(l)%y(1,1,k) = yl(1) + ((yl(2) - yl(1))*t)*tline(k)
            radial_lines(l)%z(1,1,k) = zl(1) + ((zl(2) - zl(1))*t)*tline(k)
         end do

      end if

   end do ! Next line of sight

   dmax = sqrt (dmax);  dmean = dmean / real (nlines)

   write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
      ' # outer intersection points inside/outside tolerance:',                &
      ninside, noutside, ' max & mean distance:', dmax, dmean

!  Save the lines of sight as a PLOT3D formatted file:
!  ---------------------------------------------------

   call xyz_write (lunout, true, nlines, radial_lines, ios)

!  Let F90 do all the deallocates.

99 continue

!  Internal procedures for LINES_OF_SIGHT:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_data_range ()  ! .. of the shock boundary, for intersectns.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      xmin = 1.e+10
      xmax = -xmin
      ymin =  xmin
      ymax =  xmax
      zmin =  xmin
      zmax =  xmax

      do ib = 1, nblocks
         do j = 1, outer_surface(ib)%nj
            do i = 1, outer_surface(ib)%ni
               xmax = max (xmax, outer_surface(ib)%x(i,j,1))
               xmin = min (xmin, outer_surface(ib)%x(i,j,1))
               ymax = max (ymax, outer_surface(ib)%y(i,j,1))
               ymin = min (ymin, outer_surface(ib)%y(i,j,1))
               zmax = max (zmax, outer_surface(ib)%z(i,j,1))
               zmin = min (zmin, outer_surface(ib)%z(i,j,1))
            end do
         end do
      end do

      bbox_diagonal = sqrt ((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)

      end subroutine get_data_range

   end program lines_of_sight
