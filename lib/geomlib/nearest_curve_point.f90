!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_curve_point (nc, xc, yc, zc, xt, yt, zt, &
                                   ic, xb, yb, zb, p, q, dsqmin)

!  One-liner:  Locate the nearest point on a finite discretized curve in 3-space
!
!  Description:
!
!     Subroutine NEAREST_CURVE_POINT applies the earlier NEAREST_EDGE_POINT to
!  each cell (line segment) of a 3-space curve defined by xc/yc/zc(1:nc).  In
!  the spirit of the ADT (Alternating Digital Trees) search utilities available
!  for surface and volume grids, it always returns the best point on the
!  piecewise linear curve "grid" that it can, meaning a curve end point if the
!  target point is beyond one end or the other, else the foot of the perpen-
!  dicular to the nearest straight line cell segment of the curve.
!
!     No attempt is made to use a search tree, as the arithmetic for visiting
!  all cells is minimal for any likely application.  At a higher level, the
!  intended application is to all surface (wall) edges of a 2-space CFD grid,
!  where the present routine is applied to each block of a multiblock grid.
!  At an even higher level, the intent is to facilitate turning a given 2D
!  volume grid into the equivalent one for a different surface grid, which
!  may or may not be mismatched in some way.
!
!     The third coordinate is retained in case it ever proves handy (and the
!  lower level utility requires it anyway).
!
!  History:
!
!     02/03/12  D.A.Saunders  Initial adaptation of NEAREST_EDGE_POINT, as
!                             needed by RADIAL_INTERP_2D.
!
!  Author:
!
!     David Saunders, ERC, Inc. at NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nc                     ! # pts. on searched curve
   real,    intent (in)  :: xc(nc), yc(nc), zc(nc) ! Curve point coordinates
   real,    intent (in)  :: xt, yt, zt             ! Target point coordinates
   integer, intent (out) :: ic                     ! Lower end index of the
                                                   ! best curve segment found;
                                                   ! 1 <= ic < nc
   real,    intent (out) :: xb, yb, zb             ! Point on the piecewise
                                                   ! linear curve found to be
                                                   ! nearest to (xt,yt,zt)
   real,    intent (out) :: p, q                   ! Interpolation coefficients
                                                   ! applicable to the ic:ic+1
                                                   ! curve segment or cell;
                                                   ! xb = p xc(ic) + q xc(ic+1)
                                                   ! and likewise for yb, zb;
                                                   ! p + q = 1.; 0 <= p, q <= 1.
   real,    intent (out) :: dsqmin                 ! The corresponding smallest
                                                   ! squared distance from the
                                                   ! target point to (xb,yb,zb)
!  Local constants:

   real, parameter :: one = 1., zero = 0.

!  Local variables:

   integer :: i
   real    :: dsq, pi, qi
   real, dimension (3) :: x1, x2, xtarget, xnear

!  Execution:

   dsqmin = 1.e+20
   xtarget(1) = xt;  xtarget(2) = yt;  xtarget(3) = zt

   do i = 1, nc - 1

      x1(1) = xc(i);    x1(2) = yc(i);    x1(3) = zc(i)
      x2(1) = xc(i+1);  x2(2) = yc(i+1);  x2(3) = zc(i+1)

      call nearest_edge_point (x1, x2, xtarget, xnear, pi, qi, dsq)

      if (dsq < dsqmin) then
          dsqmin = dsq
          ic = i
          p = pi
          q = qi
          xb = xnear(1);  yb = xnear(2);  zb = xnear(3)
      end if

   end do

   end subroutine nearest_curve_point
