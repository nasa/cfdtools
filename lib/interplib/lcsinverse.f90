!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine lcsinverse (x, y, ytarget, xtarget)

!     LCSINVERSE performs a specialized inverse interpolation using a monotonic
!     local cubic spline on 4 points, where the solution is known to lie in the
!     middle interval.  "Inverse" here implies finding the x corresponding to a
!     given y, which is iterative if linear interpolation does not suffice.
!
!     Confining the iteration to the middle interval avoids potential multiple
!     solution difficulties as long as the monotonic form of the 4-point LCSFIT
!     method is employed.
!
!     Locating a boundary layer edge prompted this utility.  A sequential search
!     is used at the higher level to find the first data interval containing the
!     the target ordinate, and the associated 4 surrounding points should be
!     passed as x(*) and y(*).
!
!     No check is made for a bad ytarget outside the middle interval.  This
!     should be ensured at the higher level, including no chance that y2 = y3.
!
!     A safeguarded Newton iteration is used to solve the nonlinear equation
!     f(x) = y(x) - ytarget = 0.
!
!     07/21/2006  D.A.Saunders  Initial adaptation of PLXCUT.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      real, intent (in)  :: x(4), &  ! Four data points with monotonic x, and
                            y(4)     ! ytarget known to be between y(2) & y(3)

      real, intent (in)  :: ytarget  ! Ordinate for which the abscissa is sought

      real, intent (out) :: xtarget  ! Estimate of this abscissa, to within a
                                     ! tolerance suitable for 32-bit precision
                                     ! (to avoid determining machine epsilon
                                     ! on every call)
!     Procedures:

      External  lcsfit               ! Local cubic spline utility with a
                                     ! monotonic option, meaning y cannot
                                     ! exceed the data range in any interval
!     Local constants:

      integer,   parameter :: itmax = 15   ! 4 is typical for boundary layers
      integer,   parameter :: npoints = 4  ! Always
      real,      parameter :: eps = 5.e-7, half = 0.5, one = 1., zero = 0.
      character, parameter :: method*1 = 'M'

!     Local variables:

      integer :: iter
      real    :: alpha, dx, fi, fnorm, fnorm0, toler, xa, xb, xi, xlast, yi, yp
      logical :: new

!     Execution:

!!    ier = 0

      xa = x(2);  xb = x(3);  dx = xb - xa

!     Perform a safeguarded Newton iteration:

      iter   = 0
      new    = .true.
      alpha  = zero           ! For printout purposes only
      fnorm0 = 1.e+20         ! I.e., big to avoid iteration 0 test
      toler  = dx * eps
      xi     = dx * half + xa

  300 continue

!        Evaluate y and dy/dx = df/dx at the current x.

         call lcsfit (npoints, x, y, new, method, 1, xi, yi, yp)

         new   = .false.
         fi    = yi - ytarget
         fnorm = abs (fi)

!        Halve the step until | f | is reduced (except first time through).

         if (fnorm > fnorm0) then
            if (alpha > toler) then   ! Presume toler > 2 * machine epsilon
               alpha = half * alpha
               xi = max (xa, min (xb, xlast - alpha * dx))
               go to 300
            end if
            go to 820
         end if

!!       write (*, '(a, i3, a, 1p, e9.2, a, e9.2, a, e14.6)') &
!!             ' LCSINVERSE:', iter, '  |f|:', fnorm, '  step:', alpha, &
!!             '  x:', xi

         if (fnorm < toler) go to 900  ! Success


         iter = iter + 1
         if (iter == itmax) go to 810

         fnorm0 = fnorm

!        Newton step:

         dx    = fi / yp   ! Assume derivative cannot be zero
         xlast = xi
         xi    = max (xa, min (xb, xi - dx))

         alpha = one

      go to 300                           ! Do another iteration

!     Error handling:

  810 continue
!!    ier = 1
!!    write (*, '(a)') 'LCSINVERSE iteration limit reached.'
      go to 900

  820 continue
!!    ier = 2
!!    write (*, '(a)') 'LCSINVERSE step-halving iteration failed.'

  900 continue  ! Wrap up

!!    if (iter > 6) write (*, '(a, i3)') ' Iters:', iter

      xtarget = xi

      end subroutine lcsinverse
