!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine unitnorm2d (method, n, x, y, s, seval, unormal)
!
!     UNITNORM2D determines the x and y components of the unit normal to a
!  2-space curve at the given arc length along the curve.  It makes use of
!  the derivative available from the local cubic spline routine, LCSFIT,
!  applied to x vs. s and y vs. s.  This avoids recoding derivatives w.r.t.
!  arc-length (possibly non-central or monotonic or for smoothly-closed
!  curves) in the presence of non-uniform data.
!
!     The direction of the normal is such that, using a right-hand-rule,
!  t x n is "out of the page" (parallel to the implied right-hand z axis),
!  where t is the unit tangent in the direction of increasing arc length.
!
!  History:
!
!  28 Sep. 1996  D.A.Saunders  Initial implementation of UNITNORM (results
!                              at the data points only).
!  13 Feb. 2012    "    "      Same idea but for just one normal per call, at
!                              a general point s = seval.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit  none

!  Arguments:

   character, intent (in)  :: method * 1  ! 'B' (loose), 'M' (monotonic) or
                                          ! 'C' (cyclic); cyclic means the
                                          ! curve is closed smoothly and
                                          ! point n must match point 1
   integer,   intent (in)  :: n           ! Number of data points
   real,      intent (in)  :: x(n), y(n)  ! The data points
   real,      intent (in)  :: s(n)        ! Cumulative arc lengths,
                                          ! normalized or not
   real,      intent (in)  :: seval       ! Target arc length
   real,      intent (out) :: unormal(2)  ! Unit normal x & y components

!  Local constants:

   logical, parameter :: true = .true.

!  Local variables:

   real :: size, utangent(2), xyeval

!  Procedures:

   external :: lcsfit

!  Execution:

!  dx/dt at teval:

   call lcsfit (n, s, x, true, method, 1, seval, xyeval, utangent(1))

!  Likewise for dy/dt:

   call lcsfit (n, s, y, true, method, 1, seval, xyeval, utangent(2))

   size = sqrt (dot_product (utangent, utangent))

!  Turn the partial derivatives into components of the unit normal:

   unormal(1) = -utangent(2) / size
   unormal(2) =  utangent(1) / size

   end subroutine unitnorm2d
