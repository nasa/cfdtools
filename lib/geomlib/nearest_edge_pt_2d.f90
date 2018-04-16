!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_edge_pt_2d (x1, x2, xtarget, xnear, p, q, dsq)

!  One-liner:  Locate the nearest point on a finite straight edge in 2-space
!
!  Overview:
!
!     Routines NEAREST_EDGE_PT_2D and NEAREST_TRI_PT_2D are 2-space analogues
!  of earlier 3-space utilities, designed to facilitate searching of meshes
!  consisting of triangles defined by integer triples pointing into a list of
!  (x,y) vertices (possibly from a Delaunay triangulation of a set of scattered
!  points in 2-space).
!
!     Discussion of rapid searching of 3-space meshes has been removed here,
!  but the same strategy of always returning a result that is NOT OUTSIDE the
!  target cell (edge or triangle in this case) still applies.  Then, a search of
!  a mesh boils down to locating the point with the smallest squared distance
!  from the target point, which may or may not be zero.  In the following, the
!  term "cell" can refer to a triangle (companion utility) or to an edge line
!  segment (this utility).
!
!     If a target point is not outside the current cell, then the distance to
!  the interpolated point is zero and the cell has been found.  More often, the
!  target point is outside a cell being tested.  For geometric applications at
!  least, perpendicular distances to cell boundaries are then clearly the
!  correct measure of closeness.  Sometimes the foot of the perpendicular to a
!  cell boundary lies within the finite edge or face.  More often, it is not and
!  it needs to be adjusted to provide the boundary point closest to the target.
!
!     Each of these two utilities combines an initial calculation of linear
!  interpolation coefficients with the adjustment that is most commonly needed
!  to produce an interpolated point not outside the given cell.
!
!  NEAREST_EDGE_PT_2D:
!                      xnear = p x1 + q x2                     p + q = 1
!  NEAREST_TRI_PT_2D:
!                      xnear = p x1 + q x2 + r x3              p + q + r = 1
!
!     Note:  In the interest of efficiency, degenerate cells are not handled at
!  this low level.  They should be treated by the application.
!
!  History:
!
!     08/26/04  D.A.Saunders  Initial implementation of NEAREST_EDGE_POINT (3D).
!     07/25/16    "     "     NEAREST_EDGE_PT_2D adapted from the 3D original.
!
!  Author:
!
!     David Saunders, AMA, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in),  dimension (2) :: &
      x1, x2, xtarget                   ! Two points defining the line, and
                                        ! the target point
   real, intent (out), dimension (2) :: &
      xnear                             ! The point on the straight line
                                        ! between x1 and x2 that is nearest to
                                        ! xtarget and not "off an end"
   real, intent (out) :: &
      p, q                              ! Interpolation coefficients such that
                                        ! p + q = 1 and p x1 + q x2 = xnear
   real, intent (out) :: &
      dsq                               ! The corresponding shortest distance
                                        ! between xtarget and xnear, squared

!  Local constants:

   real, parameter :: one = 1., zero = 0.

!  Local variables:

   real, dimension (2) :: a, b

!  Execution:

!  (Note that the coding is identical to the earlier 3-space variant except for
!  the vector dimensions, which are 2 here instead of 3.)
!
!  The square of the distance between (1 - q) x1 + q x2 and xtarget can be
!  represented as a scalar product and differentiated w.r.t. q to find where
!  it is smallest.  Adjusting the optimized q to the range [0, 1] ensures that
!  we don't return a point off an end.  Then the squared distance to the
!  interpolated point is a similar squared vector norm or scalar product.

   a = x2 - x1;                                   b = xtarget - x1

   q = dot_product (a, b) / dot_product (a, a);   q = max (zero, min (one, q))

   p = one - q;                                   xnear = p * x1 + q * x2

   a = q * a - b;                                 dsq = dot_product (a, a)

   end subroutine nearest_edge_pt_2d
