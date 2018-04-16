!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_tri_pt_2d (x1, x2, x3, xtarget, xnear, p, q, r, dsq)

!  One-liner:  Locate the nearest point on a triangle in 2-space
!
!  Overview:
!
!     Routines NEAREST_EDGE_PT_2D and NEAREST_TRI_PT_2D are the 2-space analogs
!  of the earlier NEAREST_EDGE_POINT and NEAREST_TRI_POINT 3-space utilities.
!  As explained in NEAREST_EDGE_PT_2D, the strategy is to return interpolation
!  coefficients corresponding to the nearest point NOT OUTSIDE the current cell.
!  Then, even if no cell is found to contain a given target point, a search
!  returns the best possible cell along with interpolation coefficients to go
!  with it that do not imply extrapolation.  These variants were prompted by a
!  need to search a 2-space Delaunay triangulation of scattered points as part
!  of avoiding (excessive) extrapolation when the target is (well) outside the
!  data range.  In the following, "cell" refers to an edge seqment or to a
!  triangle.
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
!  NEAREST_EDGE_POINT:
!                      xnear = p x1 + q x2                     p + q = 1
!  NEAREST_TRI_POINT:
!                      xnear = p x1 + q x2 + r x3              p + q + r = 1
!
!     Note:  In the interest of efficiency, degenerate cells are not handled at
!  this low level.  They should be treated by the application.
!
!  History:
!
!     08/27/04  D.A.Saunders  Initial 3-space implementation of NEAREST_TRI-
!                             POINT after many years of flawed adjustment of
!                             PROJECT3 results in various applications.
!     07/25/16    "     "     2-space analog of NEAREST_TRI_POINT, prompted by
!                             the need to search Delaunay triangulations.
!  Author:
!
!     David Saunders, AMA, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in),  dimension (2) :: &
      x1, x2, x3, xtarget               ! Three points defining the triangle,
                                        ! and the target point
   real, intent (out), dimension (2) :: &
      xnear                             ! The point in the plane of the
                                        ! triangle that is nearest to
                                        ! xtarget and not outside the triangle
   real, intent (out) :: &
      p, q, r                           ! Interpolation coefficients such that
                                        ! p x1 + q x2 + r x3 = xnear      and
                                        ! p + q + r = 1
   real, intent (out) :: &
      dsq                               ! The corresponding shortest distance
                                        ! between xtarget and xnear, squared
!  Local constants:

   real, parameter :: one = 1., zero = 0.

!  Local variables:

   integer :: i, ier, j, k
   real    :: a(2,2), b(2)

!  Execution:

!  There is more than one way to interpret the calculation.  PROJECT3 took
!  the approach suggested by Strang in Linear Algebra and its Applications,
!  where the target vector is projected into the subspace or plane that the
!  vectors representing two sides of the triangle define.  However, the
!  following interpretation is simpler.
!
!  In 2-space, we seek p, q and r such that:
!
!                            x1 p + x2 q + x3 r = xt
!                            y1 p + y2 q + y3 r = yt
!                               p +    q +    r = 1
!
!  This 3x3 is reducible to a 2x2 and the solution uniquely defines p, q and r.
!
!  (In 3-space, we have a fourth equation for z, but still just 3 coefficients.
!  See NEAREST_TRI_POINT for the details.)

   do i = 1, 2
      a(i,1) = x1(i)      - x3(i)  ! LHS column 1
      a(i,2) = x2(i)      - x3(i)  !  "    "    2
      b(i)   = xtarget(i) - x3(i)  ! RHS
   end do

!  Solve for p and q:

   call lusolve (2, 2, a, b, ier)  ! Assume degenerate cases are handled at
                                   ! a higher level

   p = b(1)         ! Initial coefficients, possibly needing adjustment
   q = b(2)
   r = one - p - q

!  Each of p, q and r must be non-negative.  Since they can't all be zero if
!  their sum is 1, we have 2**3 - 1 cases treated explicitly here to avoid
!  projecting to more than one edge unnecessarily.

!  The following kludge avoids some if ... then ... elses.
!  Is there a smarter way?

   if (p < zero) then
      i = 0
   else
      i = 1
   end if

   if (q < zero) then
      j = 0
   else
      j = 1
   end if

   if (r < zero) then
      k = 0
   else
      k = 1
   end if

   k = 4*i + 2*j + k  ! Binary --> decimal
   dsq = -one

   select case (k)

      case (1) ! 001:  p < 0, q < 0, r > 0

         p = zero;  q = zero;  r = one;  xnear = x3

      case (2) ! 010:  p < 0, q > 0, r < 0

         p = zero;  q = one;  r = zero;  xnear = x2

      case (3) ! 011:  p < 0, q > 0, r > 0

         p = zero;  call nearest_edge_pt_2d (x2, x3, xtarget, xnear, q, r, dsq)

      case (4) ! 100:  p > 0, q < 0, r < 0

         p = one;  q = zero;  r = zero;  xnear = x1

      case (5) ! 101:  p > 0, q < 0, r > 0

         q = zero;  call nearest_edge_pt_2d (x3, x1, xtarget, xnear, r, p, dsq)

      case (6) ! 110:  p > 0, q > 0, r < 0

         r = zero;  call nearest_edge_pt_2d (x1, x2, xtarget, xnear, p, q, dsq)

      case (7) ! 111:  p > 0, q > 0, r > 0; no adjustment

!!!      xnear = p * x1 + q * x2 + r * x3  ! 3-space case
         xnear = xtarget                   ! Surely

   end select

   if (dsq < zero) then  ! Carried over from the 3-space case, but probably
      b = xnear - xtarget;  dsq = dot_product (b, b)             ! redundant
   end if

   end subroutine nearest_tri_pt_2d
