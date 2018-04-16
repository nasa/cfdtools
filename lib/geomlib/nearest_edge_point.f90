!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_edge_point (x1, x2, xtarget, xnear, p, q, dsq)

!  One-liner:  Locate the nearest point on a finite straight edge in 3-space
!
!  Overview:
!
!     Routines NEAREST_EDGE_POINT, NEAREST_TRI_POINT, and NEAREST_TET_POINT are
!  low-level utilities designed to facilitate rapid searching of unstructured
!  surface and volume meshes (surface triangulations and tetrahedral meshes).
!  The rapid part is done well by ADT techniques (Alternating Digital Trees).
!  An ADT search package is typically not affected by whether a target point is
!  within some mesh cell or not.  Apart from avoiding expensive testing of the
!  bulk of the cells via tree structures, the strategy for each cell that IS
!  tested is to return interpolation coefficients corresponding to the nearest
!  point not outside that cell.  Then, even if no cell is found to contain a
!  given target point, the search returns the best possible cell along with
!  interpolation coefficients to go with it that do not imply extrapolation.
!  This is typically what is really wanted in the presence of unintended
!  geometric discrepancies, or where best mesh boundary values are appropriate.
!
!     If a target point is not outside the current cell, then the distance to
!  the interpolated point is zero and the cell has been found.  More often, the
!  target point is outside a cell being tested.  For geometric applications at
!  least, perpendicular distances to cell boundaries are then clearly the
!  correct measure of closeness.  Sometimes the foot of the perpendicular to a
!  cell boundary lies within the finite edge or face.  More often, it is not and
!  it needs to be adjusted to provide the boundary point closest to the target.
!
!     Each of these three utilities combines an initial calculation of linear
!  interpolation coefficients with the adjustment that is most commonly needed
!  to produce an interpolated point not outside the given cell.
!
!  NEAREST_EDGE_POINT:
!                      xnear = p x1 + q x2                     p + q = 1
!  NEAREST_TRI_POINT:
!                      xnear = p x1 + q x2 + r x3              p + q + r = 1
!  NEAREST_TET_POINT:
!                      xnear = p x1 + q x2 + r x3 + s x4       p + q + r + s = 1
!
!     Note:  In the interest of efficiency, degenerate cells are not handled at
!  this low level.  They should be treated by the application.
!
!  History:
!
!     08/26/04  D.A.Saunders  Initial implementation.
!
!  Author:
!
!     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in),  dimension (3) :: &
      x1, x2, xtarget                   ! Two points defining the line, and
                                        ! the target point
   real, intent (out), dimension (3) :: &
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

   real, dimension (3) :: a, b

!  Execution:

!  The square of the distance between (1 - q) x1 + q x2 and xtarget can be
!  represented as a scalar product and differentiated w.r.t. q.  Adjusting
!  the optimized q to the range [0, 1] ensures we don't return a point off
!  an end.  Then the squared distance to the interpolated point is a similar
!  squared vector norm or scalar product.

   a = x2 - x1;                                   b = xtarget - x1

   q = dot_product (a, b) / dot_product (a, a);   q = max (zero, min (one, q))

   p = one - q;                                   xnear = p * x1 + q * x2

   a = q * a - b;                                 dsq = dot_product (a, a)

   end subroutine nearest_edge_point
