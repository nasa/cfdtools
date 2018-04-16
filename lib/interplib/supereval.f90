!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine supereval (p, a, b, xc, yc, n, x, y)

!     SUPEREVAL evaluates the indicated superellipse at one or more abscissas.
!  If minor axis b is input positive, the upper portion of the curve is assumed,
!  else b < 0 and the lower portion is assumed (for all x(i)).
!
!                                       p             p
!                              (x - xc)   +  (y - yc)
!     The curve has the form   --------      --------   =   1
!                                   p             p
!                                 a             b
!
!     See SUPERFIT1 and SUPERFIT2 for fitting a [super]ellipse to a set of data
!  points.
!
!  07/29/2008  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:


   real,    intent (in)  :: p        ! Exponent of the superellipse

   real,    intent (in)  :: a        ! Major axis, with center x = 0 assumed

   real,    intent (in)  :: b        ! Minor axis; b > 0. means evaluate the
                                     ! upper portion of the curve, else b < 0.
                                     ! and the lower portion is evaluated

   real,    intent (in)  :: xc, yc   ! Center coordinates

   integer, intent (in)  :: n        ! Number of evaluation points

   real,    intent (in)  :: x(n)     ! Target abscissas on either the upper or
                                     ! lower half (but not both)

   real,    intent (out) :: y(n)     ! Evaluated ordinates, on upper half if
                                     ! b > 0. else on the lower half

!  Local constants:

   real, parameter :: one = 1., zero = 0.

!  Local variables:

   integer         :: i
   real            :: oneoverp, term

!  Execution:

   if (p == 2.) then

      do i = 1, n
         term = max (one - ((x(i) - xc) / a) ** 2, zero)
         y(i) = yc + b * sqrt (term)
      end do

   else

      oneoverp = one / p

      do i = 1, n
         term = max (one - (abs ((x(i) - xc) / a)) ** p, zero)
         y(i) = yc + b * term ** oneoverp
      end do

   end if

   end subroutine supereval
