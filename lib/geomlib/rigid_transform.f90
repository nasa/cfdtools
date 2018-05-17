!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine rigid_transform (n, x1, y1, z1, x2, y2, z2)
!
!  For a 3-space curve defined by x1/y1/z1(1:n), and end points of the same
!  curve moved to some arbitrary location defined by points 1 and n of x2/y2/z2,
!  transform the interior points by the equivalent rigid body motion.  The new
!  end points should be the same distance apart as the starting end point, but
!  this is not checked for.
!
!  10/07/2011  D.A.Saunders  Initial implementation, prompted by a need to
!                            impose catenary-type deflections between the spokes
!                            of a sphere-cone-type heat shield that opens out
!                            like an umbrella from a folded configuration.
!  10/08/2011    "     "     Note that the two new end points are not enough to
!                            define a unique transformation except if the curve
!                            is a straight line.  However, whatever result we
!                            get should be either what we want or within some
!                            rotation about the line connecting the new end
!                            points, and the application can presumably make
!                            such an adjustment.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use trigd
   implicit none

!  Arguments:

   integer, intent (in)    :: n  ! Number of points
   real,    intent (in)    :: x1(n), y1(n), z1(n)  ! Curve to be transformed
   real,    intent (inout) :: x2(n), y2(n), z2(n)  ! Input with desired new
                                                   ! end points space as for
                                                   ! the first curve;  output
                                                   ! with all points transformed
!  Local variables:

   integer :: i
   real    :: shift(3), theta, xa(3), v1(3), v2(3), v3(3)

!  Execution:

   v1(1) = x1(n) - x1(1)  ! Vector in the direction of the straight line
   v1(2) = y1(n) - y1(1)  ! connectin the end points before ...
   v1(3) = z1(n) - z1(1)

   v2(1) = x2(n) - x2(1)  ! ... and after
   v2(2) = y2(n) - y2(1)
   v2(3) = z2(n) - z2(1)

   call cross (v1, v2, v3)

   xa(1) = x2(1) + v3(1)  ! Second point defining rotation axis
   xa(2) = y2(1) + v3(2)
   xa(3) = z2(1) + v3(3)

!  tan (angle between v1 and v2) = |v1 x v2| / v1.v2:

   theta = atand (sqrt (dot_product (v3, v3)) / dot_product (v1, v2))

!  This angle is in [-p1/2, +p1/2]; it's not clear if that's an issue.

   shift(1) = x2(1) - x1(1)
   shift(2) = y2(1) - y1(1)
   shift(3) = z2(1) - z1(1)

!! write (6, '(a, 3f25.16)') &
!!    'v1:   ', v1, &
!!    'v2:   ', v2, &
!!    'v3:   ', v3, &
!!    'xa:   ', xa, &
!!    'shift:', shift, &
!!    'theta:', theta

!  Rotate the interior points after translating them first:

   do i = 2, n - 1
      x2(i) = x1(i) + shift(1)
      y2(i) = y1(i) + shift(2)
      z2(i) = z1(i) + shift(3)
   end do

   call rotate3d (n-2, x2(2), y2(2), z2(2), theta, x2(1), y2(1), z2(1), &
                  xa(1), xa(2), xa(3))

   end subroutine rigid_transform
