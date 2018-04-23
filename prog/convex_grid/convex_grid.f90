!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program convex_grid
!
!  Description:
!
!     CONVEX_GRID implements a simplistic approach to volume gridding for a
!  convex body in 3-space defined by a structured multiblock surface grid.  It
!  requires a structured outer boundary surface that may or may not be aligned
!  with the bow shock for an intended hypersonic flow calculation.  It simply
!  constructs straight radial volume grid lines from the inner surface to the
!  outer surface using an existing body normal vector utility that ensures
!  unique results even at surface patch edges and corners.  (SURFACE_VECTORS
!  constructs the average unit vector if there is more than one possibility,
!  and ensures that matching boundary points are assigned the same unit normal.)
!
!     The quality of the volume grid will depend closely on the quality of the
!  input body surface grid, which should resolve any high-curvature regions such
!  as the shoulder of a typical capsule more thoroughly than is usually the case
!  (meaning plenty of points inboard of the shoulder tangency point as well as
!  on the shoulder itself).  The application that prompted this fairly general
!  capability is a spherical quadrant surface grid with many small circular
!  features along the centerline edge that had been carefully resolved.  The
!  associated initial volume grid produced by GridPro was less than satisfactory
!  as it distorted the initial circular surface grid features.  Orthogonality at
!  the wall also suffered strangely.  (Talk to Dinesh Prabhu for the details.)
!
!     A related issue is that tangent-slab radiation estimations should be
!  performed along body-normal grid lines, so this kind of grid may facilitate
!  loose coupling of flow and radiation computations.
!
!  Method:
!
!     For each point of each block of the inner surface grid, construct the unit
!  normal vector (carefully).  Define a second point along it suitably beyond
!  the outer boundary grid, and perform a line-surface intersection.  (Such an
!  intersection is done robustly as a 1-D minimization problem: minimize the
!  the distance from a point on the straight line to the target surface.)
!  Discretize the radial grid line in the usual way: specified number of points
!  along with a (constant) wall spacing and a fraction of the outermost spacing
!  from one-sided stretching to define a two-sided stretching.
!
!  History:
!
!     05/17/2016  D.A.Saunders  Initial implementation, at Dinesh's request.
!     05/18/2016    "     "     It works with surface_normal (initial test).
!                               Now it uses surface_vectors to treat common
!                               edge points properly.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure     ! Employed by the XYZQ_IO package
   use xyzq_io_module           ! PLOT3D-type I/O package
   use adt_utilities            ! All variants of the ADT build/search utilities
   use surface_patch_utilities  ! Utilities for structured surface grids

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2,        &
      lunkbd = 5,        &
      luncrt = 6

   character (11), parameter :: &
      format = 'unformatted'

!  Variables:

   integer :: &
      ios, nk, numf, nbinner, nbouter

   real :: &
      ds1, ds2, ds2_fraction

   logical :: &
      cell_centered, cr, eof, formatted

   character :: &
      answer*1, filename*64

   type (grid_type), pointer, dimension (:) :: &
      surface_inner, surface_outer, volume_grid

!  Execution:

   call initialize ()  ! Prompt for a few options and read the surface grids
   if (ios /= 0) go to 99

   call construct_volume_grid ()  ! As a bunch of straight lines

99 continue

!  Internal procedures for program convex_grid:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialize ()  ! Prompt for controls and file names; read data
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: &
         i1, ib

!     Execution:

!     Read the inner surface mesh:

      ios = 1
      do while (ios /= 0)
         call reads (luncrt, 'Input surface grid file (PLOT3D /mgrid): ', &
                     lunkbd, filename, cr, eof)
         if (eof) go to 99  ! Quit
         if (cr) cycle

         call determine_grid_form (filename, lunin, formatted, ios)
      end do

      i1 = 1;  if (formatted) i1 = 3

      open (lunin, file=filename, form=format(i1:11), status='OLD')

      call xyzq_read (lunin, -1, formatted, nbinner, numf, cell_centered, &
                      surface_inner, ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Trouble reading ', trim (filename)
         go to 99
      end if

      write (luncrt, '(/, (i6, 2x, 3i5))') &
         (ib, surface_inner(ib)%ni, surface_inner(ib)%nj, 1, ib = 1, nbinner)

!     Read the outer boundary mesh:

      ios = 1
      do while (ios /= 0)
         call reads (luncrt, 'Outer boundary grid file (PLOT3D /mgrid): ', &
                     lunkbd, filename, cr, eof)
         if (eof) go to 99  ! Quit
         if (cr) cycle

         call determine_grid_form (filename, lunin, formatted, ios)
      end do

      i1 = 1;  if (formatted) i1 = 3

      open (lunin, file=filename, form=format(i1:11), status='OLD')

      call xyzq_read (lunin, -1, formatted, nbouter, numf, cell_centered, &
                      surface_outer, ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Trouble reading ', trim (filename)
         go to 99
      end if

      write (luncrt, '(/, (i6, 2x, 3i5))') &
         (ib, surface_outer(ib)%ni, surface_outer(ib)%nj, 1, ib = 1, nbouter)

!     Set up the output volume grid:

      call reads (luncrt, 'Output volume grid file (PLOT3D /mgrid): ', &
                  lunkbd, filename, cr, eof)
      formatted = .false.
      call ready (luncrt, 'Formatted [y] or unformatted [n|cr]: ', &
                  lunkbd, formatted, cr, eof)

      nk = 129
      call readi (luncrt, 'Number of radial volume grid points, nk [cr=129]: ',&
                  lunkbd, nk, cr, eof)
      ios = 1
      if (eof) go to 99

      ds1 = 1.e-6
      call readr (luncrt, 'Wall spacing (constant everywhere) [cr=1.e-6]: ',&
                  lunkbd, ds1, cr, eof)
      if (eof) go to 99

      ds2_fraction = 0.1
      call readr (luncrt, 'Fraction applied to 1-sided outer ds [cr=0.1]: ',&
                  lunkbd, ds2_fraction, cr, eof)
      if (eof) go to 99

      i1 = 1;  if (formatted) i1 = 3
      open (lunout, file=filename, form=format(i1:11), status='UNKNOWN')

      allocate (volume_grid(nbinner))

      do ib = 1, nbinner
         volume_grid(ib)%ni = surface_inner(ib)%ni
         volume_grid(ib)%nj = surface_inner(ib)%nj
         volume_grid(ib)%nk = nk
      end do

      call xyz_header_io (2, lunout, formatted, nbinner, volume_grid, ios)

 99   return

      end subroutine initialize

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine construct_volume_grid ()  ! As a bunch of straight lines
!
!     SHOCK_STANDOFF is the nearest similar application.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use surface_patch_utilities  ! For structured surface grids

!     Local constants:

      integer,       parameter :: nl = 2  ! 2-point lines are passed to INTSEC6
      real,          parameter :: big = 1.e+32, zero = 0.
      logical,       parameter :: false = .false., true  = .true.
      character (1), parameter :: lcs_method = 'L'  ! 2-pt. straight lines

!     Local variables:

      integer :: i, j, k, ib, ic, jc, iquad, ni, nj
      integer :: ninside, noutside, np, npts, nquad
      real    :: dmax, dmean, dsqmin, dtolsq, dx, dy, dz, p, q
      real    :: rmaxo, xmaxo, ymaxo, zmaxo
      real    :: rmino, xmino, ymino, zmino
      real    :: t, tani(3), tanj(3), un(3), xyz_interp(3)
      real, dimension (nl) :: tl, xl, yl, zl ! For 2-pt. lines passed to INTSEC6

      integer, allocatable :: &
         conn(:,:)               ! For (patch,i,j) of boundary quads. to search

      real, allocatable, dimension (:) :: &
         tline   ! For one discretized radial line distribution

!     Execution:
!     ----------

!     We need the data ranges to set the intersection tolerance appropriately:

      do ib = 1, nbinner
         call patch_data_range (surface_inner(ib))
      end do

      rmino =  big;  xmino =  big;  ymino =  big;  zmino =  big
      rmaxo = -big;  xmaxo = -big;  ymaxo = -big;  zmaxo = -big
      nquad = 0

      do ib = 1, nbouter
         nquad = (surface_outer(ib)%ni - 1)*(surface_outer(ib)%nj - 1) + nquad
         call patch_data_range (surface_outer(ib))
         xmino = min (surface_outer(ib)%xmin, xmino)
         xmaxo = max (surface_outer(ib)%xmax, xmaxo)
         ymino = min (surface_outer(ib)%ymin, ymino)
         ymaxo = max (surface_outer(ib)%ymax, ymaxo)
         zmino = min (surface_outer(ib)%zmin, zmino)
         zmaxo = max (surface_outer(ib)%zmax, zmaxo)
      end do
      rmaxo = max (xmaxo - xmino, ymaxo - ymino, zmaxo - zmino)

      dtolsq = (rmaxo * 0.00001)**2  ! Tolerance for search diagnostics

      allocate (tline(nk))

      allocate (conn(3,nquad)) ! For (patch,i,j) of each boundary quad.

!     Build a search tree from the outer boundary patches:
!     ----------------------------------------------------

      call build_adt (nbouter, surface_outer, nquad, conn)

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero
      npts = 0
      first_unit_normal = true  ! See surface_patch_utilities

!     For each inner surface point ...
!     --------------------------------

      do ib = 1, nbinner

         call xyz_allocate (volume_grid(ib), ios)

         ni   = surface_inner(ib)%ni
         nj   = surface_inner(ib)%nj
         npts = ni*nj + npts

         do j = 1, nj

            do i = 1, ni

!              Carefully calculated unit normal:

!!!            ic = i;  jc = j;  p = zero;  q = zero  ! ni, nj points are now
!!!                                                   ! handled internally
!!!                                                   ! but not common edge pts.

!!!            call surface_normal (ni, nj, surface_inner(ib)%x, &
!!!                                 surface_inner(ib)%y, surface_inner(ib)%z, &
!!!                                 ic, jc, p, q, un)

               call surface_vectors (nbinner, surface_inner, ib, i, j, &
                                     tani, tanj, un)
               first_unit_normal = false

!              Two-point line off the wall:

               xl(1) = surface_inner(ib)%x(i,j,1);  xl(2) = xl(1) + rmaxo*un(1)
               yl(1) = surface_inner(ib)%y(i,j,1);  yl(2) = yl(1) + rmaxo*un(2)
               zl(1) = surface_inner(ib)%z(i,j,1);  zl(2) = zl(1) + rmaxo*un(3)

               tl(1) = 0.001   ! Tell INTSEC6 the normalized t search range;
               tl(2) = 1.1     ! use of rmaxo should make this more than enough

               call intsec6 (nbouter, surface_outer, nquad, conn, &
                             nl, xl, yl, zl, tl, lcs_method, iquad, p, q, &
                             t, xyz_interp, dsqmin)

!!!            write (luncrt, '(a, 3i4, 3es16.8, es16.6)') &
!!!               'ib, i, j, un, dsqmin:', ib, i, j, un(:), dsqmin

               if (t > 0.99 * tl(nl)) write (luncrt, '(a, 3i5, 2f11.7)') &
                  ' Intersection warning. ib, i, j:', ib, i, j, t, tl(nl)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if

               dmax  = max (dmax, dsqmin)
               dmean = dmean + sqrt (dsqmin)

!              1-sided stretching first, then derive 2-sided d2 from it:

!              EXPDIS5 is equivalent to Vinokur's 1-sided tanh-based method.
!              MODE = 1 bunches towards the low end.
!              If increasing spacing (or close to it) is implied, it switches
!              to geometric spacing, with blending for d1 in [0.94 du, 0.97 du].

               call expdis5 (1, zero, t, ds1, nk, tline, -luncrt)

               ds2 = (t - tline(nk-1)) * ds2_fraction

!              Pure 2-sided Vinokur if ngeometric = 2 else geometric in b.layer.
!!!            tline(1)  = zero  ! Always
!!!            tline(nk) = t     ! BLGRID expects the end-points in place:

               call blgrid (nk, ds1, ds2, 2, 1.1, tline, luncrt, ios)

!              Take advantage of the fact that the radial line is straight:

               dx = xyz_interp(1) - xl(1)
               dy = xyz_interp(2) - yl(1)
               dz = xyz_interp(3) - zl(1)

               do k = 1, nk
                  volume_grid(ib)%x(i,j,k) = xl(1) + tline(k)*dx
                  volume_grid(ib)%y(i,j,k) = yl(1) + tline(k)*dy
                  volume_grid(ib)%z(i,j,k) = zl(1) + tline(k)*dz
               end do

            end do  ! Next i

         end do  ! Next j

         np = ni*nj*nk
         call xyz_block_io (2, lunout, formatted, np, volume_grid(ib)%x, &
                            volume_grid(ib)%y, volume_grid(ib)%z, ios)
         if (ios /= 0) then
            write (*, '(/, a, i4, a, i4)') ' Trouble writing volume block #', &
               ib, '.  I/O status:', ios
            write (*, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 99
         end if

         deallocate (volume_grid(ib)%x, volume_grid(ib)%y, volume_grid(ib)%z)

      end do  ! Next block

      dmax = sqrt (dmax);  dmean = dmean / real (npts)

      write (luncrt, '(/, a, 2i6, /, a, 1p, 2e12.5)') &
         ' # outer intersection points inside/outside tolerance:', &
         ninside, noutside, ' max & mean distance:', dmax, dmean

      deallocate (tline)

 99   return

      end subroutine construct_volume_grid

   end program convex_grid
