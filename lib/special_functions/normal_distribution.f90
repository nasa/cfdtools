!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine normal_distribution (mu, sigma, n, x, f)
!
!  Evaluate the normal distribution function, f(x) = a exp (-.5((x - b)/c)**2),
!  at n abscissas, where b = mu, c = sigma, and a = 1/(sigma*sqrt (2pi)).  This
!  choice of a gives an area of 1 if the function is integrated on (-inf, +inf).
!
!  05/15/2014  D.A.Saunders  Companion for the Gaussian function, both now
!                            among the options of driving program GEN1D.
!
!  Author:  David Saunders, NASA Ames Research Center/ERC. Inc, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)  :: mu, sigma    ! Mean and standard deviation of f(x)
   integer, intent (in)  :: n            ! Number of evaluations to perform
   real,    intent (in)  :: x(n)         ! Target abscissa(s)
   real,    intent (out) :: f(n)         ! Requested evaluations

!  Local constants:

   real, parameter :: one = 1.0, twopi = 6.28318530717958647692528, zero = 0.0

!  Local variables:

   real :: a

!  Execution:

   a = one / sqrt (twopi)

   call gaussian (a, mu, sigma, zero, n, x, f)

   end subroutine normal_distribution
