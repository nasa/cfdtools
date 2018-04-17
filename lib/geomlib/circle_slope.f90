!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine circle_slope (xc, yc, rc, slope, theta, xt, yt)
!
!  For the circle centered at (xc, yc) with radius rc, determine the point on it
!  that has a given slope.  In general, there are two such points, so only the
!  upper semicircle is treated here.  Treating the circle as
!                           x  =  xc  +  r cos (theta)
!                           y  =  yc  +  r sin (theta)
!  means the second solution could be determined using 180 + theta if needed.
!
!  10/24/12  D.A.Saunders  Initial implementation, prompted by bent ribs of an
!                          "umbrella"-type capsule forebody (which need to be
!                          blended with the circular shoulder somehow).
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use trigd
   implicit none

!  Arguments:

   real, intent (in)  :: xc, yc   ! Circle center
   real, intent (in)  :: rc       ! Circle radius
   real, intent (in)  :: slope    ! Circle slope (dy/dx) at the (xt, yt) sought
   real, intent (out) :: theta    ! Parametric variable <-> desired pt., deg.
   real, intent (out) :: xt, yt   ! Desired coordinates on the circle where its
                                  ! gradient matches the given slope
!  Local constants:

   real, parameter :: eps = 1.e-5, one = 1., zero = 0.

!  Execution:

   if (abs (slope) < eps) then
      theta = 90.
   else
      theta = -atan2d (one, slope)  ! In [0, 180] for the upper semicircle
      if (theta < zero) theta = theta + 180.
   end if

   xt = xc + rc * cosd (theta)
   yt = yc + rc * sind (theta)

   end subroutine circle_slope

