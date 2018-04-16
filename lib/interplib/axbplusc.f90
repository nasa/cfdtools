!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine AxBplusC (n, x, y, A, B, C, ssqmin, ier)
!
!  Purpose:
!
!     Calculate the nonlinear least squares fit of y = Ax^B + C to the given
!     (x,y) dataset representing growth or decay according to some power law.
!     A > 0 is assumed; B and C may be positive or negative.
!
!  Method:
!
!     Non-derivative minimizer FMINRC is provided with a likely search interval
!     for exponent B as follows:
!
!     If y(n) > y(1), B > 0 is assumed, and vice versa.
!
!     Slightly less than the minimum y is subtracted from all y to approximate
!     C = 0, and a linear least squares fit of log A + B log x ~ log (y - yshft)
!     provides an estimate for both A and B.  Then a 1-D minimization in the
!     original space is performed with linear least squares updates of A and C.
!
!     The C term may be suppressed by entering it as 0.
!
!     There is no requirement for the data abscissas to be uniformly spaced, but
!     they are assumed to increase monotonically and to be strictly positive.
!
!  History:
!
!     11/22/2010  D.A.Saunders  Adaptation of AeBxplusC and AxBeCx, q.v.
!
!  Author:
!
!     David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n          ! Number of data points
   real,    intent (in)    :: x(n), y(n) ! Data pt. coords.; x >= 0, monotonic
   real,    intent (out)   :: A, B       ! Coefs. of fitted y = Ax^B [+ C]; A>0
   real,    intent (inout) :: C          ! C = 0. on input means exclude it
   real,    intent (out)   :: ssqmin     ! Minimized sum of squared deviations
   integer, intent (inout) :: ier        ! Input 0 suppresses iteration printing
                                         ! Output 0 => no error detected;
                                         !       -1 => 1-D minimizer fatal error
                                         !       -2 => too few points or x <= 0
!  Constants:

   real,      parameter :: one    = 1.
   real,      parameter :: zero   = 0.
   real,      parameter :: scale1 = 5.   ! Applied to the estimate of B to give
   real,      parameter :: scale2 = 0.2  ! the search interval for optimizing it
   integer,   parameter :: nfmax  = 50   ! Limit on # function evaluations
   character, parameter :: caller * 8 = 'AxBplusC'

!  Variables:

   integer :: istat, lunerr, lunout, numfun
   logical :: excludeC
   real    :: Bhi, Blo, tol
   real    :: A1(n,3), coefs(n)

!  Procedures:

   external :: fminrc  ! Reverse-communication 1-D minimizer
   external :: hdesol  ! Linear least squares via orthogonal factorization

!  Execution:

!  Check that the inputs are reasonable:

   if (n < 3 .or. x(1) <= zero) then
      ier = -2
      go to 99
   end if

   excludeC = C == zero

   lunout = -6
   if (ier /= 0) lunout = -lunout

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

!  Calculate an estimate for B (and A and C) in log space:

   call estimateABC ()

   if (lunout > 0) write (lunout, '(a, 1p, 4e17.8)') &
      'Log-space estimates for A, B, C, ssq: ', A, B, C, ssqmin

   if (B > zero) then  ! Search interval for optimized B
      Bhi = scale1 * B
      Blo = scale2 * B
   else
      Blo = scale1 * B
      Bhi = scale2 * B
   end if

   lunerr = abs (lunout)
   numfun = nfmax      ! Limit; FMINRC typically takes 12-14 iterations
   istat  = 2          ! Initialize the minimization
   ier    = 0          ! Upon exit, normally

   do while (istat > 0)

      call fminrc (Blo, Bhi, B, ssqmin, tol, numfun, caller, lunout, istat)

      if (istat < -1) then

         write (lunerr, '(/, 2a)') caller, ': FMINRC fatal error'
         ier = -1

      else if (istat < 0) then ! Iteration limit; may be usable

         write (lunerr, '(/, 2a)') caller, ': Iteration limit.'
         istat = 0

      else  ! istat >= 0      ! Evaluate the objective function at B

         call updateAC ()     ! Linear least sqrs. soln. for A [, C] for this B;
                              ! the sum of squared deviations is a by-product
         if (istat == 0) exit ! Success
      end if

   end do

99 return

   contains

!     Internal procedures for subroutine AeBxplusC:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine estimateABC ()

!     Estimate A and B via linear least squares in log space after effectively
!     shifting the data to approximate C = 0 while avoiding log (0).
!     A can only be positive from this approach.  Distinguishing, say,
!     an A < 0, B < 0 case from A > 0, 0 < B < 1 doesn't seem possible,
!     so we don't try to detect a negative A case.

!     Local variables:

      integer :: i, i1, nsolve
      real    :: ymax, ymin, ymargin, yshift

!     Execution:

      ymax = maxval (y)
      ymin = minval (y)
      ymargin = (ymax - ymin) * tol
      yshift  = ymin - ymargin  ! Avoid getting too close to log (0)

!     Avoiding the closest to log (0) improves results:

      nsolve  = 0
      do i = 1, n
         if (y(i) == ymin) cycle
         nsolve = nsolve + 1
         A1(nsolve,1) = one
         A1(nsolve,2) = log (x(i))
         A1(nsolve,3) = log (y(i) - yshift)  ! RHS
      end do

      if (nsolve < 2) then
         ier = -2
         go to 99
      end if

      call hdesol (n, nsolve, 3, A1, coefs, ssqmin)

      A = exp (coefs(1))
      B = coefs(2)
      if (.not. excludeC) C = yshift

   99 return

      end subroutine estimateABC

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine updateAC ()  ! Linear least squares calculation for given B

!     Local variables:

      integer :: i, ncp1

!     Execution:

!     Set up the overdetermined system with the RHS as column nc + 1:

      if (excludeC) then
         ncp1 = 2  ! Only one coefficient to calculate
         do i = 1, n
            A1(i,1) = x(i) ** B
            A1(i,2) = y(i)
         end do
      else
         ncp1 = 3
         do i = 1, n
            A1(i,1) = x(i) ** B
            A1(i,2) = one
            A1(i,3) = y(i)
         end do
      end if

      call hdesol (n, n, ncp1, A1, coefs, ssqmin) ! Orthogonal factorzn. & soln.

      A = coefs(1)
      if (.not. excludeC) C = coefs(2)

      end subroutine updateAC

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine AxBplusC
