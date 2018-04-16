!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_tet_point (x1, x2, x3, x4, xtarget, xnear, p, q, r, s, &
                                 dsq)

!  One-liner:  Locate the nearest point in or on a tetrahedron
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
!     08/27/04  D.A.Saunders  Initial implementation, by analogy with the
!                             newly implemented triangle case.
!  Author:
!
!     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in),  dimension (3) :: &
      x1, x2, x3, x4, xtarget           ! Four points defining the tetrahedron,
                                        ! and the target point
   real, intent (out), dimension (3) :: &
      xnear                             ! The point nearest to xtarget that is
                                        ! not outside the tetrahedron
   real, intent (out) :: &
      p, q, r, s                        ! Interpolation coefficients such that
                                        ! p x1 + q x2 + r x3 + s x4 = xnear  and
                                        ! p + q + r + s = 1
   real, intent (out) :: &
      dsq                               ! The corresponding shortest distance
                                        ! between xtarget and xnear, squared
!  Local constants:

   real, parameter :: one = 1., zero = 0.

!  Local variables:

   integer :: i, j, k, l
   real    :: a(3,3), b(3)

!  Execution:

!  Calculate initial interpolation coefficients.  See TET_COEFS if this is all
!  that is needed.  Do it in-line here.

   do i = 1, 3
      a(i,1) = x1(i)      - x4(i)
      a(i,2) = x2(i)      - x4(i)
      a(i,3) = x3(i)      - x4(i)
      b(i)   = xtarget(i) - x4(i)
   end do

!  Solution of Ax = b by LU decomposition.  RHS is overwritten with solution.

   call lusolve (3, 3, a, b, i)

   if (i /= 0) THEN
      write (*, '(/, a)')           &
         ' NEAREST_TET_POINT:  Singular matrix.  Collapsed cell?  Aborting.'
      write (*, '(a, 1p, 3e19.11)') &
         ' x1:', x1, ' x2:', x2, ' x3:', x3, ' x4:', x4, 'xt:', xtarget
      stop
   end if

   p = b(1);  q = b(2);  r = b(3);  s = one - p - q - r

!  Each of p, q, r and s must be non-negative.  Since they can't all be zero
!  if their sum is 1, we have 2**4 - 1 cases treated explicitly here to avoid
!  projecting to more than one face or edge unnecessarily.

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

   if (s < zero) then
      l = 0
   else
      l = 1
   end if

   k = 8*i + 4*j + 2*k + l  ! Binary --> decimal
   dsq = -one

   select case (k)

      case (1)  ! 0001:  p < 0, q < 0, r < 0, s > 0

         p = zero;  q = zero;  r = zero;  s = one;  xnear = x4

      case (2)  ! 0010:  p < 0, q < 0, r > 0, s < 0

         p = zero;  q = zero;  r = one;  s = zero;  xnear = x3

      case (3)  ! 0011:  p < 0, q < 0, r > 0, s > 0

         p = zero;  q = zero

         call nearest_edge_point (x3, x4, xtarget, xnear, r, s, dsq)

      case (4)  ! 0100:  p < 0, q > 0, r < 0, s < 0

         p = zero;  q = one;  r = zero;  s = zero;  xnear = x2

      case (5)  ! 0101:  p < 0, q > 0, r < 0, s > 0

         p = zero;  r = zero

         call nearest_edge_point (x2, x4, xtarget, xnear, q, s, dsq)

      case (6)  ! 0110:  p < 0, q > 0, r > 0, s < 0

         p = zero;  s = zero

         call nearest_edge_point (x2, x3, xtarget, xnear, q, r, dsq)

      case (7)  ! 0111:  p < 0, q > 0, r > 0, s > 0

         p = zero

         call nearest_tri_point (x2, x3, x4, xtarget, xnear, q, r, s, dsq)

      case (8)  ! 1000:  p > 0, q < 0, r < 0, s < 0

         p = one;  q = zero;  r = zero;  s = zero;  xnear = x1

      case (9)  ! 1001:  p > 0, q < 0, r < 0, s > 0

         q = zero;  r = zero

         call nearest_edge_point (x4, x1, xtarget, xnear, s, p, dsq)

      case (10) ! 1010:  p > 0, q < 0, r > 0, s < 0

         q = zero;  s = zero

         call nearest_edge_point (x1, x3, xtarget, xnear, p, r, dsq)

      case (11) ! 1011:  p > 0, q < 0, r > 0, s > 0

         q = zero

         call nearest_tri_point (x3, x4, x1, xtarget, xnear, r, s, p, dsq)

      case (12) ! 1100:  p > 0, q > 0, r < 0, s < 0

         r = zero;  s = zero

         call nearest_edge_point (x1, x2, xtarget, xnear, p, q, dsq)

      case (13) ! 1101:  p > 0, q > 0, r < 0, s > 0

         r = zero

         call nearest_tri_point (x1, x2, x4, xtarget, xnear, p, q, s, dsq)

      case (14) ! 1110:  p > 0, q > 0, r > 0, s < 0

         s = zero

         call nearest_tri_point (x1, x2, x3, xtarget, xnear, p, q, r, dsq)

      case (15) ! 1111:  p > 0, q > 0, r > 0, s > 0

         xnear = xtarget;  dsq = zero  ! No adjustment; target is in this cell

   end select

   if (dsq < zero) then
      b = xnear - xtarget;  dsq = dot_product (b, b)
   end if

   end subroutine nearest_tet_point
