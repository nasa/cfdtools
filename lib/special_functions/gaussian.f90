!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine gaussian (a, b, c, d, n, x, f)
!
!  Evaluate the Gaussian function  f(x) = a exp(-.5((x - b)/c)**2) + d at n
!  abscissas (n >= 1).  As Wikipedia explains, this is a symmetric bell-shaped
!  curve with height a above d, centered at x = b, with standard deviation c.
!  See companion routine normal_distribution, and driving program GEN1D.
!
!  05/15/2014  D.A.Saunders  Coded to construct plausible thruster pulses.
!
!  Author:  David Saunders, NASA Ames Research Center/ERC. Inc, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)  :: a, b, c, d   ! Parameters described above
   integer, intent (in)  :: n            ! Number of evaluations to perform
   real,    intent (in)  :: x(n)         ! Target abscissa(s)
   real,    intent (out) :: f(n)         ! Requested evaluations

!  Local constants:

   real, parameter :: half = 0.5, one = 1.0

!  Local variables:

   integer :: i
   real    :: constant

!  Execution:

   constant = one / c  !  [(x - b)/c]**2 should be preferable to squaring both

   do i = 1, n
      f(i) = a * exp (-half * ((x(i) - b)*constant)**2) + d
   end do

   end subroutine gaussian
