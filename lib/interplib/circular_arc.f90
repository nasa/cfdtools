!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine circular_arc (x1, y1, x2, y2, x3, y3, r, n, xc, yc)

!  Evaluate the 2-space coordinates of n points distributed uniformly along a
!  circular arc as defined by two tangential points 1 and 2 and the point 3
!  where those tangents meet.  The associated radius is unique (and returned).
!  Point 3 must be equidistant from points 1 and 2 (assumed).
!
!  ALTERNATIVELY:  If the radius argument is input as a positive value, make
!                  that the circle radius instead of deriving the radius from
!  the implied tangent length.  This is checked for first, and if so, derive
!  the tangent length then proceed as for the original option.  The input
!  coordinates for points 1 and 2 should still be equidistant from the vertex
!  (because of tests for which side of the vertex to put the circle) but they
!  can be anywhere on the underlying tangents.  The radius argument is returned
!  as the implied tangent length in this case, for the convenience of the
!  calling program.
!
!  This functionality was prompted by a requirement to round angular vertices in
!  generatrices for the aft bodies of axisymmetric capsules (entry vehicles).
!  It may find other applications.
!
!  11/04/2011  D.A.Saunders  Initial implementation.
!  11/08/2011    "    "      Arranged for specifying the circle radius as an
!                            alternative to the original implied tangent length.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)    :: x1, y1, x2, y2  ! Tangential point coordinates
   real,    intent (in)    :: x3, y3          ! Underlying vertex coordinates
   real,    intent (inout) :: r               ! Radius of the defined circle,
                                              ! used if input > 0 then returned
                                              ! as the implied tangent length,
                                              ! else output with implied radius
   integer, intent (in)    :: n               ! Specified number of uniform pts.
   real,    intent (out)   :: xc(n), yc(n)    ! Requested discretization

!  Local constants:

   real, parameter :: half = 0.5, zero = 0.

!  Local variables:

   integer :: i
   real    :: costheta, da, db, dtheta, l1, l2, t, theta, ux, uy, &
              xa, xb, xcenter, xm, xl, xr, ya, yb, ycenter, ym, yl, yr, &
              v1(2), v2(2)
   logical :: radius_in_tangent_out

!  Execution:

   radius_in_tangent_out = r > zero

   if (radius_in_tangent_out) then  ! Use it rather than derive it
      v1(1) = x3 - x1;  v1(2) = y3 - y1;  l1 = sqrt (dot_product (v1, v1))
      v2(1) = x3 - x2;  v2(2) = y3 - y2;  l2 = sqrt (dot_product (v2, v2))
      theta = 90. - half * acosd (dot_product (v1, v2) / (l1 * l2))
      t     = r * tand (theta)
      ux = v1(1) / l1   ! Unit vector along the tangent
      uy = v1(2) / l1
      xl = x3 - t*ux    ! True tangent point
      yl = y3 - t*uy
      ux = v2(1) / l2   ! Unit vector along the tangent
      uy = v2(2) / l2
      xr = x3 - t*ux    ! Second true tangent point
      yr = y3 - t*uy
   else  ! Calculate the implied radius
      xl = x1  ! So that the revised code below works for both cases
      yl = y1
      xr = x2
      yr = y2
      xm = (xl + xr)*half  ! Mid point of chord between points 1-2
      ym = (yl + yr)*half
      t  = distance (xl, yl, x3, y3)       ! Tangent length, as for pts. 2-3
      r  = distance (xl, yl, xm, ym) * t / &
           distance (x3, y3, xm, ym)
      theta  = atand (t/r)                 ! Half the sector angle
   end if

!  The circle center along the radius through point 1 will be closer to
!  point 2 for one of two possibilities, so we compute both and pick one:

   ux = (x3 - xl) / t  ! Unit vector along the tangent
   uy = (y3 - yl) / t

   xa = xl + r*uy;  ya = yl - r*ux;  da = distance (xa, ya, xr, yr)
   xb = xl - r*uy;  yb = yl + r*ux;  db = distance (xb, yb, xr, yr)

   if (radius_in_tangent_out) r = t  ! Can be used by the calling program

   if (da < db) then
       xcenter = xa;  ycenter = ya
   else
       xcenter = xb;  ycenter = yb
   end if

   dtheta = (theta + theta) / real (n-1)

!  Rotating the radius through point 1 about the center clockwise and
!  anticlockwise will come closer to point 2 one way or the other:

   xa = xl;  ya = yl
   call rotate2d (1, xa, ya,  theta, xcenter, ycenter)
   da = distance (xa, ya, xr, yr)

   xb = xl;  yb = yl
   call rotate2d (1, xb, yb, -theta, xcenter, ycenter)
   db = distance (xb, yb, xr, yr)

   if (da > db) dtheta = -dtheta

!  Discretize the arc by rotating point 1 about the circle center:

   xc(:) = xl;  yc(:) = yl

   do i = 2, n - 1
      theta = dtheta * real (i-1)
      call rotate2d (1, xc(i), yc(i), theta, xcenter, ycenter)
   end do

   xc(n) = xr;  yc(n) = yr

   contains

      real function distance (xl, yl, xr, yr)  ! Obvious purpose and arguments

      real, intent (in) :: xl, yl, xr, yr

      distance = sqrt ((xl - xr)**2 + (yl - yr)**2)

      end function distance

   end subroutine circular_arc
