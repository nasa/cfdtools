!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine catenary_grid (semispan, deflection, ymax, mode, n, x, y, a, ier)
!
!  Generate a 2-space catenary curve (left, right or both halves) discretized as
!  n points uniformly spaced in x over the indicated span or semispan.  The peak
!  deflection and the semispan b define the parameter a of  y = a cosh (x/a) - a
!  passing through (0, 0), but the ymax argument allows placing the end point(s)
!  at y = 0 (or any other vertical location).
!
!  Such a grid line can be transformed as a 3-space surface grid line using the
!  rigid_transform subroutine by the same author.  Its discretization can then
!  be adjusted as necessary via spline interpolation versus arc length.
!
!  10/13/2011  D.A.Saunders  Initial implementation for program forebody_regrid.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)  :: semispan   ! Half width of catenary centered on Oy
   real,    intent (in)  :: deflection ! Maximum y droop below the end point(s)
   real,    intent (in)  :: ymax       ! Desired y of end point(s)
   integer, intent (in)  :: mode       ! mode = 1 => left half   on [-b, 0];
                                       ! mode = 2 => right half  on [ 0, b];
                                       ! mode = 3 => both halves on [-b, b]
   integer, intent (in)  :: n          ! Number of points in discretization
   real,    intent (out) :: x(n), y(n) ! Desired 2-space grid, uniform in x
   real,    intent (out) :: a          ! Catenary parameter calculated to give
                                       ! the desired shape, in case it's reqd.
   integer, intent (out) :: ier        ! ier /= 0 means the iterative
                                       ! calculation of parameter a failed;
                                       ! see subroutine catenary_parameters
!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i
   real    :: b, dx, ymin, yshift

!  Execution:

   b = semispan;  ymin = ymax - deflection

   call catenary_parameters (semispan, deflection, a, ier)

   if (ier /= 0) go to 99  ! Single return philosophy

   select case (mode)
      case (1)
         x(1) = -b;    x(n) = zero
         y(1) = ymax;  y(n) = ymin
      case (2)
         x(1) = zero;  x(n) = b
         y(1) = ymin;  y(n) = ymax
      case (3)
         x(1) = -b;    x(n) = b
         y(1) = ymax;  y(n) = ymax
   end select

   yshift = ymax - deflection - a
   dx = (x(n) - x(1)) / real (n - 1)

   do i = 2, n - 1
      x(i) = x(1) + dx*real (i - 1)
      y(i) = a*cosh (x(i)/a) + yshift
   end do

99 return

   end subroutine catenary_grid
