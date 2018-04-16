!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine derivs_123 (n, x, f, i, fp, fpp, fppp)

!  Calculate the first three derivatives of f(x) at x(i) for 3 <= i <= n-2.
!
!  Strategy:  Since finite differencing gets messy for nonuniform spacing and
!  higher-order polynomials, the easy way out is taken by calculating the unique
!  quartic through the five points i-2, i-1, i, i+1, i+2 and evaluating its
!  derivatives at point i.
!
!  12/11/2007  D.A.Saunders  Initial implementation as part of determining the
!                            quintic defined by f, f' and f, f', f'', f''' at
!                            two points respectively.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n              ! Length of discrete curve in 2-space
   real,    intent (in)  :: x(n), f(n)     ! Discrete curve
   integer, intent (in)  :: i              ! Central target point (see above)
   real,    intent (out) :: fp, fpp, fppp  ! Desired derivatives

!  Local constants:

   integer, parameter :: nw = 35           ! PNFIT needs 5 pts. * (4th deg + 3)
   logical, parameter :: thru00 = .false.  ! Curve doesn't pass through (0,0)

!  Local variables:

   integer :: ier  ! Needed by PNFIT; assumed never needed here
   real :: c(0:4)  ! Coefficients of the underlying quartic
   real :: w(nw)   ! Workspace for PNFIT
   real :: rmsdev  ! Should be zero for this polynomial fit
   real :: t       ! Target abscissa

!  Execution:

!  Fit the unique quartic:

   call pnfit (5, x(i-2), f(i-2), 4, thru00, nw, w, c, rmsdev, ier)

!  Evaluate the derivatives:

   t    = x(i)
   fp   = ((4.*c(4)*t + 3.*c(3))*t + 2.*c(2))*t + c(1)
   fpp  = (12.*c(4)*t + 6.*c(3))*t + 2.*c(2)
   fppp =  24.*c(4)*t + 6.*c(3)

   end subroutine derivs_123
