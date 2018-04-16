!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine trapezoid_recursion (n, t1, tprev, fprev, tnew, fnew,            &
                                   fmean, variance)
!
!  Incorporate one more data point in the mean & variance of a time-series-type
!  dataset via recursions based on the trapezoidal rule for quadrature with
!  respect to the time-like variable.  The variance is calculated as for the
!  mean, but on the squared deviation from the mean rather than the function.
!
!  The recursion is initialized here if input argument n is 1, just setting
!  fmean = fnew and variance = 0.  For n > 1, the mean and variance are updated
!  in place.  The value of n is otherwise not needed, but it may be useful in
!  the calling program, so n is used in favor of some other "initialize" flag.
!
!  The initial abscissa cannot be assumed to be zero, even if that is likely.
!
!  04/30/10  D. A. Saunders  Initial implementation, at Todd White's request.
!  05/03/10     "     "      The variance recursion is NOT simply that for the
!                            mean of the squared deviation from the mean,
!                            because that should be from the CURRENT mean,
!                            not the previous mean.  The algebra gets ugly,
!                            but it simplifies quite nicely in the end.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n         ! n = 1 => first data pt. in the series;
                                        ! n > 1 (or not 1) is not otherwise used
   real,    intent (in)    :: t1        ! The abscissa of data point 1
   real,    intent (in)    :: tprev, &  ! Data point from preceding call
                              fprev
   real,    intent (in)    :: tnew,  &  ! Current data point
                              fnew
   real,    intent (inout) :: fmean, &  ! Input with previous values;
                              variance  ! output with updated values

!  Local constants:

   real, parameter :: half = 0.5, zero = 0.

!  Local variables:

   real :: dmean, dtnew, dtprev, gnew, gprev, prevmean, term

!  Execution:

   if (n == 1) then  ! Initialize this recursion

      fmean     =  fnew
      variance  =  zero

   else

      dtprev    =  tprev - t1
      dtnew     =  tnew  - t1
      term      =  half * (tnew - tprev)
      prevmean  =  fmean

      fmean     = (fmean * dtprev    + term * (fprev + fnew)) / dtnew
      dmean     =  fmean - prevmean

      gprev     = (fprev - prevmean) ** 2           ! g = (f - previous fmean)^2
      gnew      = (fnew  - prevmean) ** 2

      variance  = (variance * dtprev + term * (gprev + gnew)) / dtnew  & ! Orig.
                   - dmean ** 2                                ! Correction term

   end if

   end subroutine trapezoid_recursion
