!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine superfit2 (n, x, y, pmin, pmax, p, a, b, yc, rsq)

!     SUPERFIT2 performs the nonlinear optimization portion of fitting a super-
!  ellipse to data points where the major axis a is known, the center in the x
!  direction is assumed to be zero, and the unknowns are the minor axis b, the
!  center in the y direction yc, and the exponent p.  SUPERFIT1 calculates the
!  best values of b and yc for given exponent p, via linear least squares.
!  The nonderivative 1-D minimizer FMINRC calculates p here by calling SUPERFIT1
!  and minimizing the squared residual within the indicated exponent range.
!
!     This capability was prompted by a need to represent the cross-section of a
!  recessed arc-jet test article as a smooth discretized curve.  The data points
!  (x,y) here for the initial application are an averaged section derived from
!  laser scan measurements by other numerical techniques.
!
!     The data points may cover either one or two quadrants, but they must
!  represent either the top half of the [super]ellipse or the bottom half, but
!  not both.  If they are on the bottom half, the output b will be negative.
!  |x(i)| <= a is assumed for all i = 1:n.
!
!     Use SUPEREVAL to discretize the fitted curve.
!
!  07/29/2008  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n        ! Number of data points

   real,    intent (in)  :: x(n), &  ! Data coordinates on either the upper or
                            y(n)     ! lower half (or quarter) of the expected
                                     ! [super]ellise, but not both

   real,    intent (in)  :: pmin, &  ! Exponent range allowed, e.g., [2., 6.]
                            pmax

   real,    intent (out) :: p        ! Exponent of the superellipse

   real,    intent (in)  :: a        ! Major axis, with center x = 0 assumed

   real,    intent (out) :: b        ! Minor axis calculated here; b < 0
                                     ! if the lower portion of the curve is
                                     ! represented by the data points

   real,    intent (out) :: yc       ! Center y coordinate calculated here

   real,    intent (out) :: rsq      ! Squared residual from the linear least
                                     ! squares fit performed by SUPERFIT1
!  Procedures:

   external  fminrc                  ! 1-D reverse-communication minimizer

!  Local constants:

   real, parameter   :: one = 1.

!  Local constants:

   integer,   parameter :: lunout =  6    ! Suppress FMINRC iteration printout
   integer,   parameter :: nfmax  = 50    ! Limit on # function evaluations
   character, parameter :: caller * 9 = 'SUPERFIT2'

!  Local variables:

   integer :: istat, lunerr, numfun
   real    :: tol

!  Execution:

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

   numfun = nfmax      ! Limit; FMINRC typically takes about 6 iterations
   lunerr = abs (lunout)
   istat = 2           ! Initialize the minimization

10 continue

      call fminrc (pmin, pmax, p, rsq, tol, numfun, caller, lunout, istat)

      if (istat < -1) then

         write (lunerr, '(/, 2a)') caller, ': FMINRC fatal error'
         stop

      else if (istat < 0) then ! Iteration limit; may be usable

         write (lunerr, '(/, 2a)') caller, ': Iteration limit.'

      else if (istat > 0) then ! Evaluate the objective function at p

         call superfit1 (n, x, y, p, a, b, yc, rsq)
         go to 10

      else ! istat = 0 (success).

      end if

!  Ensure that everything matches p(best), not p(last):

   call superfit1 (n, x, y, p, a, b, yc, rsq)

   end subroutine superfit2
